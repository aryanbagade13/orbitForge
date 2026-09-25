import numpy as np
import pytest

from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from gravity_assist.models import CelestialBody, OrbitalState
from legacy.engine.simulation import (
    apply_spacecraft_manoeuvre_to_system_state,
    pack_system_state,
    restricted_three_body_derivative,
)


def make_test_system(spacecraft_position_km):
    sun_body = CelestialBody(
        name="Sun",
        mass_kg=1.0,
        radius_km=10.0,
    )
    sun_state = OrbitalState(
        position_km=np.array([0.0, 0.0, 0.0]),
        velocity_km_s=np.array([0.0, 0.0, 0.0]),
    )
    planet_body = CelestialBody(
        name="Test planet",
        mass_kg=1.0,
        radius_km=5.0,
    )
    planet_state = OrbitalState(
        position_km=np.array([100.0, 0.0, 0.0]),
        velocity_km_s=np.array([0.0, 0.0, 0.0]),
    )
    spacecraft_state = OrbitalState(
        position_km=np.asarray(spacecraft_position_km, dtype=float),
        velocity_km_s=np.array([0.0, 0.0, 0.0]),
    )

    return (
        sun_body,
        sun_state,
        planet_body,
        pack_system_state(planet_state, spacecraft_state),
    )


def test_derivative_rejects_spacecraft_inside_sun():
    sun_body, sun_state, planet_body, system_state = make_test_system(
        [5.0, 0.0, 0.0]
    )

    with pytest.raises(ValueError, match="spacecraft intersects the Sun"):
        restricted_three_body_derivative(
            0.0,
            system_state,
            sun_body,
            sun_state,
            planet_body,
        )


def test_derivative_rejects_spacecraft_inside_planet():
    sun_body, sun_state, planet_body, system_state = make_test_system(
        [103.0, 0.0, 0.0]
    )

    with pytest.raises(
        ValueError,
        match="spacecraft intersects Test planet",
    ):
        restricted_three_body_derivative(
            0.0,
            system_state,
            sun_body,
            sun_state,
            planet_body,
        )


def test_spacecraft_manoeuvre_changes_only_spacecraft_velocity():
    original_system_state = np.array([
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0,
        10.0, 11.0, 12.0,
    ])
    manoeuvre = ImpulsiveManoeuvre(
        time_s=100.0,
        delta_velocity_km_s=np.array([0.1, -0.2, 0.3]),
    )

    updated_system_state = apply_spacecraft_manoeuvre_to_system_state(
        original_system_state,
        manoeuvre,
    )

    np.testing.assert_allclose(
        updated_system_state,
        [
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
            10.1, 10.8, 12.3,
        ],
    )
    np.testing.assert_allclose(
        original_system_state,
        [
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
            10.0, 11.0, 12.0,
        ],
    )
