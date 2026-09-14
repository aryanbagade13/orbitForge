"""Spacecraft propagation using prescribed planetary ephemerides.

This Newtonian model uses an inertial barycentric frame. Clearance is checked
at RK4 stages and output samples, not with continuous event detection. A
coarse timestep can miss a close encounter. JPL giant-planet barycentres are
only approximations to planet centres for surface-clearance purposes.
"""

from collections.abc import Mapping, Sequence

import numpy as np

from .ephemerides import EphemerisProvider
from .forces import total_gravitational_acceleration
from .mission import propagate_with_manoeuvres
from .mission_candidate import MissionCandidate
from .mission_definition import MissionDefinition
from .models import CelestialBody
from .results import MissionResult


class UnsafeTrajectoryError(ValueError):
    """An evaluated state touches a body or violates its clearance."""


def make_ephemeris_derivative(
    bodies: Sequence[CelestialBody],
    ephemerides: EphemerisProvider,
    minimum_altitudes_km: Mapping[str, float] | None = None,
    *, check_clearance: bool = True,
):
    """Build dr/dt=v, dv/dt=sum(gravity) for a six-component spacecraft."""
    bodies = tuple(bodies)
    names = [body.name.strip().lower() for body in bodies]
    if not bodies or len(set(names)) != len(names):
        raise ValueError("bodies must be nonempty with unique names")
    limits = {}
    for name, altitude in (minimum_altitudes_km or {}).items():
        key = name.strip().lower()
        if key not in names or key in limits:
            raise ValueError("clearance names must uniquely match supplied bodies")
        if not np.isfinite(altitude) or altitude < 0:
            raise ValueError("minimum altitudes must be non-negative and finite")
        limits[key] = float(altitude)

    def derivative(time_s, state):
        state = np.asarray(state, dtype=float)
        if state.shape != (6,) or not np.isfinite(state).all():
            raise ValueError("spacecraft state must be finite with shape (6,)")
        body_states = [ephemerides.get_state(name, time_s) for name in names]
        for body, name, body_state in zip(bodies, names, body_states):
            distance = np.linalg.norm(state[:3] - body_state.position_km)
            if check_clearance and (distance <= body.radius_km or distance < body.radius_km + limits.get(name, 0.0)):
                raise UnsafeTrajectoryError(
                    f"spacecraft violates clearance for {body.name} at {time_s} s"
                )
        acceleration = total_gravitational_acceleration(state[:3], bodies, body_states)
        return np.concatenate((state[3:6], acceleration))

    return derivative


def simulate_mission(
    definition: MissionDefinition,
    candidate: MissionCandidate,
    ephemerides: EphemerisProvider,
    bodies: Sequence[CelestialBody],
    dt_s: float,
) -> MissionResult:
    """Propagate a candidate, without claiming it reaches its required flybys.

    Destination miss distance and encounter sequence belong to the future
    evaluator. Burn vectors use the same fixed axes as the spacecraft state.
    Caller supplies body masses/radii consistently with the ephemeris model.
    """
    if definition.departure_epoch != ephemerides.departure_epoch:
        raise ValueError("mission and ephemeris departure epochs must match")
    if definition.reference_frame != ephemerides.frame:
        raise ValueError("mission and ephemeris reference frames must match")
    earliest, latest = definition.arrival_time_bounds_s
    if not earliest <= candidate.arrival_time_s <= latest:
        raise ValueError("candidate arrival must lie within mission bounds")
    bodies = tuple(bodies)
    names = {body.name.strip().lower() for body in bodies}
    required = {name.strip().lower() for name in definition.flyby_body_names}
    required.add(definition.destination_body_name.strip().lower())
    required.add("sun")
    if not required <= names:
        raise ValueError("bodies must include Sun, destination and required flyby bodies")
    derivative = make_ephemeris_derivative(
        bodies, ephemerides, definition.minimum_altitudes_km
    )
    initial = np.concatenate((
        definition.initial_spacecraft_state.position_km,
        definition.initial_spacecraft_state.velocity_km_s,
    ))
    derivative(0.0, initial)
    result = propagate_with_manoeuvres(
        initial, 0.0, candidate.arrival_time_s, dt_s,
        derivative, candidate.manoeuvres,
        state_layout="spacecraft", reference_frame=ephemerides.frame,
    )
    # Weighted RK4 endpoints differ from its internal stage states.
    for time_s, state in zip(result.times_s, result.states):
        derivative(time_s, state)
    return result
