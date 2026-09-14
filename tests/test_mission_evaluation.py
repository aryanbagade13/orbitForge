from datetime import datetime, timezone
from dataclasses import replace

import numpy as np
import pytest

from gravity_assist.models import CelestialBody, OrbitalState
from gravity_assist.mission_definition import MissionDefinition, FlybyRequirement
from gravity_assist.mission_candidate import MissionCandidate
from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from gravity_assist.results import MissionResult
from gravity_assist.evaluation import evaluate_mission, closest_approach
from gravity_assist.adaptive import simulate_adaptive_mission
from gravity_assist.departure import make_departure_state
from gravity_assist.ephemerides import TabulatedEphemerisProvider


class TestEphemeris:
    __test__ = False
    departure_epoch = datetime(2000, 1, 1, tzinfo=timezone.utc)
    frame = "test inertial frame"
    positions = {"sun": [-1000., 0., 0.], "jupiter": [2., 0., 0.], "saturn": [10., 0., 0.]}

    def get_state(self, name, time):
        return OrbitalState(self.positions[name.lower()], np.zeros(3))


def setup():
    provider = TestEphemeris()
    bodies = tuple(CelestialBody(name, 1., .1) for name in provider.positions)
    initial = OrbitalState([0., 1., 0.], [1., 0., 0.])
    definition = MissionDefinition(
        provider.departure_epoch, initial, ("Jupiter",), "Saturn", (9., 11.),
        2., {}, provider.frame, (FlybyRequirement("Jupiter", (1., 3.), 2.),)
    )
    candidate = MissionCandidate([], 10.)
    def straight(t):
        return np.array([t, 1., 0., 1., 0., 0.])
    result = MissionResult([0., 10.], [straight(0.), straight(10.)], (),
                           "spacecraft", provider.frame, provider.departure_epoch, straight)
    return provider, bodies, definition, candidate, result


def test_refinement_finds_encounter_between_two_samples():
    p, bodies, d, c, result = setup()
    measurement = closest_approach(result, p, bodies[1], 0., 10.)
    assert measurement.time_s == pytest.approx(2., abs=1e-5)
    assert measurement.distance_km == pytest.approx(1.)
    assert measurement.relative_speed_km_s == pytest.approx(1.)


def test_evaluator_reports_known_success_and_zero_cost():
    p, bodies, d, c, result = setup()
    evaluation = evaluate_mission(d, c, result, p, bodies)
    assert evaluation.feasible
    assert evaluation.arrival_distance_km == pytest.approx(1.)
    assert evaluation.total_delta_v_km_s == 0.


def test_evaluator_reports_missed_flyby():
    p, bodies, d, c, result = setup()
    d = replace(d, flyby_requirements=(FlybyRequirement("Jupiter", (1., 3.), .5),))
    evaluation = evaluate_mission(d, c, result, p, bodies)
    assert not evaluation.feasible
    assert evaluation.constraint_residuals["flyby_0_jupiter_distance_km"] == pytest.approx(.5)


def test_evaluator_rejects_ambiguous_named_flyby():
    p, bodies, d, c, result = setup()
    with pytest.raises(ValueError, match="matching ordered"):
        evaluate_mission(replace(d, flyby_requirements=()), c, result, p, bodies)


def test_evaluator_detects_between_sample_clearance_failure():
    p, bodies, d, c, result = setup()
    d = replace(d, minimum_altitudes_km={"Jupiter": 1.2})
    evaluation = evaluate_mission(d, c, result, p, bodies)
    assert not evaluation.feasible
    assert evaluation.constraint_residuals["clearance_jupiter_km"] == pytest.approx(.3)


def test_evaluator_checks_encounter_order():
    p, bodies, d, c, result = setup()
    d = replace(d, flyby_body_names=("Saturn", "Jupiter"), flyby_requirements=(
        FlybyRequirement("Saturn", (7., 9.), 2.),
        FlybyRequirement("Jupiter", (1., 3.), 2.),
    ))
    assert not evaluate_mission(d, c, result, p, bodies).feasible


def test_adaptive_burn_dense_states_and_final_burn():
    p, bodies, d, _, _ = setup()
    candidate = MissionCandidate([ImpulsiveManoeuvre(2.5, [1.,0.,0.]),
                                  ImpulsiveManoeuvre(10., [-1.,0.,0.])], 10.)
    result = simulate_adaptive_mission(d, candidate, p, bodies, max_step_s=1.)
    assert result.spacecraft_state_at(2.5)[3] == pytest.approx(2.)
    assert result.spacecraft_state_at(2.5-1e-5)[3] == pytest.approx(1.)
    assert result.spacecraft_state_at(10.)[0] == pytest.approx(17.5)
    assert result.spacecraft_state_at(10.)[3] == pytest.approx(1.)
    with pytest.raises(ValueError, match="outside"):
        result.spacecraft_state_at(11.)


def test_adaptive_stops_at_clearance_crossing():
    p, bodies, d, _, _ = setup()
    d = replace(d, minimum_altitudes_km={"Jupiter": 1.4})
    from gravity_assist.interplanetary import UnsafeTrajectoryError
    with pytest.raises(UnsafeTrajectoryError, match="jupiter"):
        simulate_adaptive_mission(d, MissionCandidate([],10.), p, bodies, max_step_s=.1)


def test_departure_offsets_are_in_ephemeris_frame():
    state = make_departure_state(TestEphemeris(), "Jupiter", [3,4,5], [1,2,3])
    np.testing.assert_array_equal(state.position_km, [5,4,5])
    np.testing.assert_array_equal(state.velocity_km_s, [1,2,3])


def test_ephemeris_table_agrees_with_cubic_motion():
    class Cubic(TestEphemeris):
        def get_state(self, name, time):
            return OrbitalState([time**3,0,0], [3*time**2,0,0])
    table = TabulatedEphemerisProvider(Cubic(), ["Sun"], 10., 2.)
    np.testing.assert_allclose(table.get_state("Sun",3.).position_km, [27,0,0])
    np.testing.assert_allclose(table.get_state("Sun",3.).velocity_km_s, [27,0,0])
    with pytest.raises(ValueError, match="coverage"):
        table.get_state("Sun",11.)


def test_table_rejects_excessive_interpolation_error():
    class Quartic(TestEphemeris):
        def get_state(self, name, time):
            return OrbitalState([time**4,0,0], [4*time**3,0,0])
    with pytest.raises(ValueError, match="too coarse"):
        TabulatedEphemerisProvider(Quartic(), ["Sun"], 10., 10.)


def test_assessment_reports_unsafe_instead_of_success():
    from gravity_assist.evaluation import assess_candidate
    p, bodies, d, _, _ = setup()
    d = replace(d, minimum_altitudes_km={"Jupiter": 1.4})
    assessment = assess_candidate(d, MissionCandidate([],10.), p, bodies, max_step_s=.1)
    assert assessment.status == "unsafe"
    assert assessment.evaluation is None
    assert "jupiter" in assessment.failure_reason


def test_assessment_reports_feasible_known_trajectory():
    from gravity_assist.evaluation import assess_candidate
    p, bodies, d, c, _ = setup()
    assessment = assess_candidate(d, c, p, bodies, max_step_s=1.)
    assert assessment.status == "feasible"
    assert all(v == 0. for v in assessment.evaluation.constraint_violations.values())


def test_adaptive_unit_orbit_matches_analytic_solution():
    from gravity_assist.constants import G
    class FixedSun(TestEphemeris):
        positions = {"sun": [0.,0.,0.]}
    p = FixedSun()
    d = MissionDefinition(p.departure_epoch, OrbitalState([1,0,0],[0,1,0]), (),
                          "Sun", (1.,10.), 2., {}, p.frame)
    result = simulate_adaptive_mission(d, MissionCandidate([],2*np.pi), p,
                                       [CelestialBody("Sun",1/G,.01)],
                                       rtol=1e-11, position_atol_km=1e-12,
                                       velocity_atol_km_s=1e-12, max_step_s=.1)
    np.testing.assert_allclose(result.states[-1], [1,0,0,0,1,0], atol=1e-9, rtol=0)


def test_adaptive_zero_time_and_simultaneous_burns():
    p, bodies, d, _, _ = setup()
    c = MissionCandidate([ImpulsiveManoeuvre(0., [1.,0.,0.]),
                          ImpulsiveManoeuvre(0., [-.5,0.,0.])], 10.)
    result = simulate_adaptive_mission(d,c,p,bodies,max_step_s=1.)
    assert result.spacecraft_state_at(0.)[3] == pytest.approx(1.5)
    assert result.states[-1,0] == pytest.approx(15.)
    assert np.all(np.diff(result.times_s)>0)


def test_evaluator_rejects_mismatched_candidate():
    p, bodies, d, _, result = setup()
    c = MissionCandidate([ImpulsiveManoeuvre(5.,[1.,0.,0.])],10.)
    with pytest.raises(ValueError,match="do not match"):
        evaluate_mission(d,c,result,p,bodies)


@pytest.mark.parametrize("bounds", [(-1.,2.), (2.,2.), (3.,2.), (0.,np.inf)])
def test_flyby_requirements_reject_invalid_time_windows(bounds):
    with pytest.raises(ValueError):
        FlybyRequirement("Jupiter",bounds,1e6)
