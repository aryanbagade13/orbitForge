from datetime import datetime, timezone

import numpy as np
import pytest

from gravity_assist.constants import G
from legacy.engine.integrators import propagate_fixed_step
from gravity_assist.interplanetary import (
    UnsafeTrajectoryError, make_ephemeris_derivative,
)
from legacy.engine.interplanetary import simulate_mission
from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from legacy.engine.mission import propagate_with_manoeuvres
from gravity_assist.mission_candidate import MissionCandidate
from gravity_assist.mission_definition import MissionDefinition
from gravity_assist.models import CelestialBody, OrbitalState


class LinearEphemeris:
    departure_epoch = datetime(2000, 1, 1, tzinfo=timezone.utc)
    frame = "ICRS axes; Solar System barycentric origin"

    def __init__(self, velocity=None):
        self.velocity = np.zeros(3) if velocity is None else np.asarray(velocity)

    def get_state(self, body_name, time_s):
        return OrbitalState(self.velocity * time_s, self.velocity)


def unit_sun():
    return CelestialBody("Sun", 1 / G, 0.01)


def test_derivative_uses_moving_body_and_preserves_input():
    state = np.array([4., 0., 0., 0., 2., 0.])
    original = state.copy()
    derivative = make_ephemeris_derivative([unit_sun()], LinearEphemeris([1., 0., 0.]))
    np.testing.assert_allclose(derivative(2., state), [0., 2., 0., -0.25, 0., 0.])
    np.testing.assert_array_equal(state, original)


def test_multiple_bodies_accelerations_add():
    bodies = [unit_sun(), CelestialBody("Other", 2/G, 0.01)]
    derivative = make_ephemeris_derivative(bodies, LinearEphemeris())
    np.testing.assert_allclose(derivative(0., [2., 0., 0., 0., 0., 0.]),
                               [0., 0., 0., -0.75, 0., 0.])


def test_clearance_rejects_unsafe_state():
    derivative = make_ephemeris_derivative([unit_sun()], LinearEphemeris(), {"Sun": 1.})
    with pytest.raises(UnsafeTrajectoryError, match="Sun"):
        derivative(0., [1., 0., 0., 0., 0., 0.])


def test_six_state_burns_at_start_off_grid_and_arrival():
    def coast(t, state):
        return np.concatenate((state[3:], np.zeros(3)))

    burns = [ImpulsiveManoeuvre(t, [1., 0., 0.]) for t in [4., 0., 2.5]]
    initial = np.zeros(6)
    result = propagate_with_manoeuvres(initial, 0., 4., 1., coast, burns,
                                       state_layout="spacecraft")
    np.testing.assert_allclose(result.states[-1], [5.5, 0., 0., 3., 0., 0.])
    np.testing.assert_array_equal(initial, np.zeros(6))
    assert np.all(np.diff(result.times_s) > 0)
    assert [b.time_s for b in result.manoeuvres] == [0., 2.5, 4.]
    np.testing.assert_array_equal(result.spacecraft_states, result.states)


def test_circular_orbit_converges_when_step_halved():
    derivative = make_ephemeris_derivative([unit_sun()], LinearEphemeris())
    initial = np.array([1., 0., 0., 0., 1., 0.])
    errors = []
    for dt in [0.1, 0.05]:
        _, states = propagate_fixed_step(initial, 0., 2*np.pi, dt, derivative)
        errors.append(np.linalg.norm(states[-1] - initial))
    assert errors[1] < errors[0] / 10
    assert errors[1] < 1e-5


def definition(**overrides):
    args = dict(departure_epoch=LinearEphemeris.departure_epoch,
                initial_spacecraft_state=OrbitalState([1., 0., 0.], [0., 1., 0.]),
                flyby_body_names=(), destination_body_name="Sun",
                arrival_time_bounds_s=(0.1, 10.), destination_max_distance_km=2.,
                minimum_altitudes_km={})
    return MissionDefinition(**(args | overrides))


def test_mission_entry_matches_direct_propagation():
    candidate = MissionCandidate([], 1.)
    result = simulate_mission(definition(), candidate, LinearEphemeris(), [unit_sun()], .01)
    _, expected = propagate_fixed_step(
        [1., 0., 0., 0., 1., 0.], 0., 1., .01,
        make_ephemeris_derivative([unit_sun()], LinearEphemeris()),
    )
    np.testing.assert_allclose(result.states, expected)
    assert result.state_layout == "spacecraft"
    assert result.reference_frame == LinearEphemeris.frame
    assert result.times_s[-1] == 1.


@pytest.mark.parametrize("overrides,message", [
    ({"reference_frame": "heliocentric"}, "reference frames"),
    ({"departure_epoch": datetime(2001, 1, 1, tzinfo=timezone.utc)}, "epochs"),
    ({"arrival_time_bounds_s": (2., 3.)}, "arrival"),
    ({"destination_body_name": "Saturn"}, "bodies"),
    ({"minimum_altitudes_km": {"Unknown": 1.}}, "clearance names"),
])
def test_mission_rejects_incompatible_inputs(overrides, message):
    with pytest.raises(ValueError, match=message):
        simulate_mission(definition(**overrides), MissionCandidate([], 1.),
                         LinearEphemeris(), [unit_sun()], .01)
