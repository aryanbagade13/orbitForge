import numpy as np
import pytest

from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from gravity_assist.mission_candidate import MissionCandidate


def test_candidate_orders_burns_without_reordering_input():
    early = ImpulsiveManoeuvre(0.0, [0.01, 0.0, 0.0])
    late = ImpulsiveManoeuvre(10.0, [-0.01, 0.0, 0.0])
    burns = [late, early]
    candidate = MissionCandidate(burns, 10.0)

    assert isinstance(candidate.manoeuvres, tuple)
    assert candidate.manoeuvres[0] is early
    assert candidate.manoeuvres[1] is late
    assert burns[0] is late
    assert candidate.total_delta_v_km_s == pytest.approx(0.02)


def test_candidate_without_burns_has_zero_cost():
    candidate = MissionCandidate([], 10.0)
    assert candidate.manoeuvres == ()
    assert candidate.total_delta_v_km_s == 0.0


@pytest.mark.parametrize("arrival", [0.0, -1.0, np.nan, np.inf, -np.inf])
def test_candidate_rejects_invalid_arrival(arrival):
    with pytest.raises(ValueError, match="positive and finite"):
        MissionCandidate([], arrival)


@pytest.mark.parametrize("arrival", ["10", True, None])
def test_candidate_rejects_non_numeric_arrival(arrival):
    with pytest.raises(TypeError, match="real number"):
        MissionCandidate([], arrival)


def test_candidate_rejects_invalid_burn_object():
    with pytest.raises(TypeError, match="ImpulsiveManoeuvre"):
        MissionCandidate([object()], 10.0)


def test_candidate_rejects_burn_after_arrival():
    burn = ImpulsiveManoeuvre(11.0, [0.0, 0.01, 0.0])
    with pytest.raises(ValueError, match="after arrival"):
        MissionCandidate([burn], 10.0)
