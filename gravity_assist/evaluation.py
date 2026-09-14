"""Mission measurements and explicit constraint residuals, without a search algorithm."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class EncounterMeasurement:
    body_name: str
    time_s: float
    distance_km: float
    altitude_km: float
    relative_speed_km_s: float


@dataclass(frozen=True)
class MissionEvaluation:
    total_delta_v_km_s: float
    arrival_distance_km: float
    arrival_relative_velocity_km_s: tuple[float, float, float]
    closest_approaches: tuple[EncounterMeasurement, ...]
    flybys: tuple[EncounterMeasurement, ...]
    # Every numeric residual uses <= 0 for satisfaction; units are in keys.
    constraint_residuals: dict[str, float]
    feasible: bool

    @property
    def constraint_violations(self):
        """Positive excesses only; each retains the unit in its key."""
        return {name: max(0.0, value) for name, value in self.constraint_residuals.items()}


@dataclass(frozen=True)
class CandidateAssessment:
    status: str
    total_delta_v_km_s: float
    evaluation: MissionEvaluation | None
    result: object | None
    failure_reason: str | None = None


def assess_candidate(definition, candidate, ephemerides, bodies, **integration_options):
    """One simulation and evaluation, ready for a future search to call.

    Clearance failures are explicitly unsuccessful. Configuration errors still
    raise, so a bad frame or missing body cannot masquerade as a bad trajectory.
    """
    from .adaptive import simulate_adaptive_mission
    from .interplanetary import UnsafeTrajectoryError

    bodies = tuple(bodies)
    try:
        result = simulate_adaptive_mission(definition, candidate, ephemerides,
                                           bodies, **integration_options)
    except UnsafeTrajectoryError as error:
        return CandidateAssessment("unsafe", candidate.total_delta_v_km_s,
                                   None, None, str(error))
    evaluation = evaluate_mission(definition, candidate, result, ephemerides, bodies)
    return CandidateAssessment("feasible" if evaluation.feasible else "infeasible",
                               candidate.total_delta_v_km_s, evaluation, result)


def closest_approach(result, ephemerides, body, start_s, end_s):
    """Refine extrema using roots of r_relative dot v_relative.

    Search each adaptive mesh interval in four subintervals, including burn
    boundaries and window endpoints. This is a numerical search, not a proof
    that no unresolved extrema exist; verify mesh/tolerance convergence.
    """
    if not result.times_s[0] <= start_s <= end_s <= result.times_s[-1]:
        raise ValueError("encounter window outside trajectory")

    def relative(t):
        spacecraft = result.spacecraft_state_at(t)
        planet = ephemerides.get_state(body.name, t)
        return spacecraft[:3]-planet.position_km, spacecraft[3:]-planet.velocity_km_s

    def radial(t):
        position, velocity = relative(t)
        return float(np.dot(position, velocity))

    mesh = np.unique(np.r_[start_s, result.times_s[
        (result.times_s > start_s) & (result.times_s < end_s)], end_s])
    candidates = list(mesh)
    for left, right in zip(mesh[:-1], mesh[1:]):
        grid = np.linspace(left, right, 5)
        values = [radial(t) for t in grid]
        for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if fa < 0 < fb:
                candidates.append(brentq(radial, a, b, xtol=1e-5))
            elif fa == 0:
                candidates.append(a)
        candidates.extend(grid)
    best_time = min(candidates, key=lambda t: np.linalg.norm(relative(t)[0]))
    position, velocity = relative(best_time)
    distance = float(np.linalg.norm(position))
    return EncounterMeasurement(body.name, float(best_time), distance,
                                distance-body.radius_km, float(np.linalg.norm(velocity)))


def evaluate_mission(definition, candidate, result, ephemerides, bodies):
    """Evaluate an actual propagated proposal; no weighted objective is imposed."""
    bodies = tuple(bodies)
    catalogue = {body.name.strip().lower(): body for body in bodies}
    if len(catalogue) != len(bodies):
        raise ValueError("duplicate body names")
    if (result.reference_frame != definition.reference_frame
            or ephemerides.frame != definition.reference_frame
            or result.departure_epoch != definition.departure_epoch
            or ephemerides.departure_epoch != definition.departure_epoch):
        raise ValueError("result, definition and ephemerides must share epoch and frame")
    if result.times_s[0] != 0 or result.times_s[-1] != candidate.arrival_time_s:
        raise ValueError("trajectory must cover departure through candidate arrival")
    if len(result.manoeuvres) != len(candidate.manoeuvres) or any(
        a.time_s != b.time_s or not np.array_equal(a.delta_velocity_km_s, b.delta_velocity_km_s)
        for a, b in zip(result.manoeuvres, candidate.manoeuvres)
    ):
        raise ValueError("trajectory manoeuvres do not match candidate")
    initial = np.r_[definition.initial_spacecraft_state.position_km,
                    definition.initial_spacecraft_state.velocity_km_s].copy()
    for burn in candidate.manoeuvres:
        if burn.time_s == 0:
            initial[3:] += burn.delta_velocity_km_s
    if not np.allclose(result.spacecraft_states[0], initial, rtol=0, atol=1e-9):
        raise ValueError("trajectory initial state does not match mission definition")
    requirements = definition.flyby_requirements
    if tuple(r.body_name.strip().lower() for r in requirements) != tuple(
        name.strip().lower() for name in definition.flyby_body_names
    ):
        raise ValueError("each named flyby needs a matching ordered FlybyRequirement")
    required_names = {definition.destination_body_name.strip().lower()}
    required_names |= {name.strip().lower() for name in definition.flyby_body_names}
    required_names |= {name.strip().lower() for name in definition.minimum_altitudes_km}
    if not required_names <= catalogue.keys():
        raise ValueError("missing body parameters for mission evaluation")

    arrival = candidate.arrival_time_s
    limits = {name.strip().lower(): altitude for name, altitude in definition.minimum_altitudes_km.items()}
    destination = catalogue[definition.destination_body_name.strip().lower()]
    if definition.destination_max_distance_km <= destination.radius_km + limits.get(destination.name.lower(), 0.):
        raise ValueError("destination target lies inside its clearance boundary")
    spacecraft = result.spacecraft_state_at(arrival)
    planet = ephemerides.get_state(destination.name, arrival)
    arrival_distance = float(np.linalg.norm(spacecraft[:3]-planet.position_km))
    relative_velocity = tuple(float(v) for v in spacecraft[3:]-planet.velocity_km_s)
    residuals = {
        "arrival_distance_km": arrival_distance-definition.destination_max_distance_km,
        "arrival_earliest_s": definition.arrival_time_bounds_s[0]-arrival,
        "arrival_latest_s": arrival-definition.arrival_time_bounds_s[1],
    }
    closest = tuple(closest_approach(result, ephemerides, body, 0., arrival) for body in bodies)
    for measurement in closest:
        name = measurement.body_name.lower()
        residuals[f"clearance_{name}_km"] = limits.get(name, 0.)-measurement.altitude_km
    flybys = []
    previous_time = 0.
    for index, requirement in enumerate(requirements):
        body = catalogue[requirement.body_name.strip().lower()]
        if requirement.max_distance_km <= body.radius_km + limits.get(body.name.lower(), 0.):
            raise ValueError("flyby target lies inside its clearance boundary")
        start, end = requirement.time_bounds_s
        stop = min(end, arrival)
        if start >= stop:
            raise ValueError("flyby window must overlap trajectory before arrival")
        measurement = closest_approach(result, ephemerides, body, start, stop)
        flybys.append(measurement)
        prefix = f"flyby_{index}_{body.name.lower()}"
        residuals[prefix+"_distance_km"] = measurement.distance_km-requirement.max_distance_km
        residuals[prefix+"_window_start_s"] = start+1e-6-measurement.time_s
        residuals[prefix+"_window_end_s"] = measurement.time_s+1e-6-stop
        # Require a strictly intermediate event; 1 microsecond is a numerical margin.
        residuals[prefix+"_order_s"] = previous_time+1e-6-measurement.time_s
        residuals[prefix+"_before_arrival_s"] = measurement.time_s+1e-6-arrival
        previous_time = measurement.time_s
    feasible = all(value <= 0 for value in residuals.values())
    # Touching a physical surface is a collision, even when minimum altitude is zero.
    feasible = feasible and all(item.altitude_km > 0 for item in closest)
    return MissionEvaluation(candidate.total_delta_v_km_s, arrival_distance,
                             relative_velocity, closest, tuple(flybys), residuals, feasible)
