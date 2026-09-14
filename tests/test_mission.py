import numpy as np
import pytest

from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from gravity_assist.mission import propagate_with_manoeuvres
from gravity_assist.results import MissionResult


def constant_velocity_derivative(time_s, state):
    derivative = np.zeros_like(state)
    derivative[0:3] = state[3:6]
    derivative[6:9] = state[9:12]
    return derivative


def test_manoeuvre_changes_motion_at_scheduled_time():
    initial_state = np.zeros(12)
    initial_state[9] = 1.0
    manoeuvre = ImpulsiveManoeuvre(
        time_s=2.5,
        delta_velocity_km_s=np.array([1.0, 0.0, 0.0]),
    )

    result = propagate_with_manoeuvres(
        initial_state=initial_state,
        start_time_s=0.0,
        end_time_s=4.0,
        dt_s=1.0,
        derivative_function=constant_velocity_derivative,
        manoeuvres=[manoeuvre],
    )

    assert isinstance(result, MissionResult)
    assert len(result.manoeuvres) == 1
    assert result.manoeuvres[0] is manoeuvre
    times_s = result.times_s
    states = result.states

    np.testing.assert_allclose(
        times_s,
        [0.0, 1.0, 2.0, 2.5, 3.5, 4.0],
    )
    np.testing.assert_allclose(states[3, 6:9], [2.5, 0.0, 0.0])
    np.testing.assert_allclose(states[3, 9:12], [2.0, 0.0, 0.0])
    np.testing.assert_allclose(states[-1, 6:9], [5.5, 0.0, 0.0])
    np.testing.assert_allclose(states[-1, 9:12], [2.0, 0.0, 0.0])
    np.testing.assert_allclose(initial_state[9:12], [1.0, 0.0, 0.0])


def test_manoeuvre_time_must_lie_within_mission():
    initial_state = np.zeros(12)
    manoeuvre = ImpulsiveManoeuvre(
        time_s=5.0,
        delta_velocity_km_s=np.zeros(3),
    )

    with pytest.raises(
        ValueError,
        match="manoeuvre time must lie within the mission",
    ):
        propagate_with_manoeuvres(
            initial_state=initial_state,
            start_time_s=0.0,
            end_time_s=4.0,
            dt_s=1.0,
            derivative_function=constant_velocity_derivative,
            manoeuvres=[manoeuvre],
        )
