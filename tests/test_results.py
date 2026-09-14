import numpy as np
import pytest

from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from gravity_assist.results import MissionResult


def test_mission_result_converts_inputs_to_expected_types():
    manoeuvre = ImpulsiveManoeuvre(
        time_s=5.0,
        delta_velocity_km_s=[0.0, 0.01, 0.0],
    )

    result = MissionResult(
        times_s=[0.0, 1.0],
        states=[[0.0] * 12, [1.0] * 12],
        manoeuvres=[manoeuvre],
    )

    assert isinstance(result.times_s, np.ndarray)
    assert isinstance(result.states, np.ndarray)
    assert result.times_s.dtype == float
    assert result.states.dtype == float
    assert result.manoeuvres == (manoeuvre,)


def test_mission_result_copies_input_arrays():
    original_times_s = np.array([0.0, 1.0])
    original_states = np.zeros((2, 12))

    result = MissionResult(
        times_s=original_times_s,
        states=original_states,
        manoeuvres=[],
    )

    original_times_s[0] = 99.0
    original_states[0, 0] = 99.0

    assert result.times_s[0] == 0.0
    assert result.states[0, 0] == 0.0


@pytest.mark.parametrize(
    "invalid_states",
    [
        np.zeros(12),
        np.zeros((2, 6)),
        np.zeros((2, 12, 1)),
    ],
)
def test_mission_result_rejects_states_with_wrong_shape(invalid_states):
    with pytest.raises(ValueError, match="states must have shape"):
        MissionResult(
            times_s=[0.0, 1.0],
            states=invalid_states,
            manoeuvres=[],
        )


def test_mission_result_rejects_non_one_dimensional_times():
    with pytest.raises(ValueError, match="times_s must be one-dimensional"):
        MissionResult(
            times_s=[[0.0], [1.0]],
            states=np.zeros((2, 12)),
            manoeuvres=[],
        )


def test_mission_result_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="matching lengths"):
        MissionResult(
            times_s=[0.0, 1.0],
            states=np.zeros((3, 12)),
            manoeuvres=[],
        )


@pytest.mark.parametrize(
    "invalid_times_s",
    [
        [0.0, 0.0],
        [1.0, 0.0],
    ],
)
def test_mission_result_rejects_times_that_are_not_increasing(
    invalid_times_s,
):
    with pytest.raises(ValueError, match="strictly increasing"):
        MissionResult(
            times_s=invalid_times_s,
            states=np.zeros((2, 12)),
            manoeuvres=[],
        )


@pytest.mark.parametrize("invalid_time", [np.nan, np.inf, -np.inf])
def test_mission_result_rejects_non_finite_times(invalid_time):
    with pytest.raises(ValueError, match="times_s.*finite"):
        MissionResult(
            times_s=[0.0, invalid_time],
            states=np.zeros((2, 12)),
            manoeuvres=[],
        )


@pytest.mark.parametrize("invalid_state", [np.nan, np.inf, -np.inf])
def test_mission_result_rejects_non_finite_states(invalid_state):
    states = np.zeros((2, 12))
    states[1, 4] = invalid_state

    with pytest.raises(ValueError, match="states.*finite"):
        MissionResult(
            times_s=[0.0, 1.0],
            states=states,
            manoeuvres=[],
        )
