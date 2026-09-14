"""Adaptive spacecraft propagation with burn-separated dense solutions."""

import numpy as np
from scipy.integrate import solve_ivp

from .interplanetary import UnsafeTrajectoryError, make_ephemeris_derivative
from .results import MissionResult


def simulate_adaptive_mission(definition, candidate, ephemerides, bodies, *,
                              rtol=1e-10, position_atol_km=1e-4,
                              velocity_atol_km_s=1e-10, max_step_s=86400.0):
    """DOP853 integration in the provider's barycentric, fixed-axis frame.

    max_step_s limits event-search gaps, not saved-plot resolution. Event
    sign changes can still miss two crossings in one step; the evaluator
    additionally refines closest approaches on the dense trajectory.
    """
    if not np.isfinite([rtol, position_atol_km, velocity_atol_km_s, max_step_s]).all() or min(
            rtol, position_atol_km, velocity_atol_km_s, max_step_s) <= 0:
        raise ValueError("integration tolerances and maximum step must be positive and finite")
    if definition.departure_epoch != ephemerides.departure_epoch:
        raise ValueError("mission and ephemeris departure epochs must match")
    if definition.reference_frame != ephemerides.frame:
        raise ValueError("mission and ephemeris reference frames must match")
    if not definition.arrival_time_bounds_s[0] <= candidate.arrival_time_s <= definition.arrival_time_bounds_s[1]:
        raise ValueError("candidate arrival must lie within mission bounds")
    bodies = tuple(bodies)
    names = {body.name.strip().lower() for body in bodies}
    required = {name.strip().lower() for name in definition.flyby_body_names}
    required |= {"sun", definition.destination_body_name.strip().lower()}
    if not required <= names:
        raise ValueError("bodies must include Sun, destination and required flyby bodies")
    limits = {name.strip().lower(): value for name, value in definition.minimum_altitudes_km.items()}
    derivative = make_ephemeris_derivative(bodies, ephemerides, limits, check_clearance=False)
    initial = np.concatenate((definition.initial_spacecraft_state.position_km,
                              definition.initial_spacecraft_state.velocity_km_s))
    make_ephemeris_derivative(bodies, ephemerides, limits)(0.0, initial)

    events = []
    for body in bodies:
        def clearance(t, state, body=body):
            planet = ephemerides.get_state(body.name, t)
            return (np.linalg.norm(state[:3]-planet.position_km) - body.radius_km
                    - limits.get(body.name.strip().lower(), 0.0))
        clearance.terminal = True
        clearance.direction = -1
        events.append(clearance)

    times, states = [0.0], [initial.copy()]
    segments = []
    current = initial.copy()
    current_time = 0.0
    burns = sorted(candidate.manoeuvres, key=lambda burn: burn.time_s)
    for burn in [*burns, None]:
        stop = candidate.arrival_time_s if burn is None else burn.time_s
        if stop > current_time:
            solution = solve_ivp(
                derivative, (current_time, stop), current, method="DOP853",
                rtol=rtol, atol=[position_atol_km]*3+[velocity_atol_km_s]*3,
                max_step=max_step_s, dense_output=True, events=events,
            )
            if solution.status == 1:
                index = next(i for i, values in enumerate(solution.t_events) if len(values))
                raise UnsafeTrajectoryError(
                    f"clearance violation for {bodies[index].name} at {solution.t[-1]:.6f} s"
                )
            if not solution.success:
                raise RuntimeError(f"integration failed: {solution.message}")
            segments.append((current_time, stop, solution.sol))
            times.extend(solution.t[1:])
            states.extend(solution.y.T[1:].copy())
            current = solution.y[:, -1].copy()
            current_time = stop
        if burn is not None:
            current[3:] += burn.delta_velocity_km_s
            states[-1] = current.copy()

    starts = np.array([segment[0] for segment in segments])
    final = current.copy()

    def dense_state(t):
        if t == candidate.arrival_time_s:
            return final.copy()
        index = max(0, int(np.searchsorted(starts, t, side="right"))-1)
        return segments[index][2](t)

    return MissionResult(times, states, tuple(burns), state_layout="spacecraft",
                         reference_frame=ephemerides.frame,
                         departure_epoch=definition.departure_epoch, dense_state=dense_state)
