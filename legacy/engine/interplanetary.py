"""Earlier fixed-step ephemeris propagation; current missions use adaptive.py."""

from collections.abc import Sequence

import numpy as np

from gravity_assist.ephemerides import EphemerisProvider
from gravity_assist.interplanetary import make_ephemeris_derivative
from gravity_assist.mission_candidate import MissionCandidate
from gravity_assist.mission_definition import MissionDefinition
from gravity_assist.models import CelestialBody
from gravity_assist.results import MissionResult
from .mission import propagate_with_manoeuvres


def simulate_mission(
    definition: MissionDefinition,
    candidate: MissionCandidate,
    ephemerides: EphemerisProvider,
    bodies: Sequence[CelestialBody],
    dt_s: float,
) -> MissionResult:
    """Propagate a candidate with fixed-step RK4 and sampled clearance checks."""
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
