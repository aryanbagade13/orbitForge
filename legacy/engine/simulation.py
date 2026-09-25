import numpy as np
from gravity_assist.collisions import position_intersects_body
from gravity_assist.manoeuvres import ImpulsiveManoeuvre, apply_impulsive_manoeuvre
from gravity_assist.models import OrbitalState, CelestialBody
from gravity_assist.forces import total_gravitational_acceleration

def pack_system_state(
    planet_state: OrbitalState,
    spacecraft_state: OrbitalState,
) -> np.ndarray:
    return np.concatenate((
        planet_state.position_km,
        planet_state.velocity_km_s,
        spacecraft_state.position_km,
        spacecraft_state.velocity_km_s,
    ))


def apply_spacecraft_manoeuvre_to_system_state(
    system_state: np.ndarray,
    manoeuvre: ImpulsiveManoeuvre,
) -> np.ndarray:
    updated_system_state = np.asarray(system_state, dtype=float).copy()

    if updated_system_state.shape != (12,):
        raise ValueError("system_state must have shape (12,)")
    if not np.isfinite(updated_system_state).all():
        raise ValueError("system_state must contain only finite values")

    spacecraft_state = OrbitalState(
        position_km=updated_system_state[6:9],
        velocity_km_s=updated_system_state[9:12],
    )
    manoeuvred_spacecraft_state = apply_impulsive_manoeuvre(
        spacecraft_state,
        manoeuvre.delta_velocity_km_s,
    )

    updated_system_state[6:9] = manoeuvred_spacecraft_state.position_km
    updated_system_state[9:12] = manoeuvred_spacecraft_state.velocity_km_s

    return updated_system_state


def restricted_three_body_derivative(time, system_state, sun_body, sun_state, planet_body):
    if system_state.shape != (12,):
        raise ValueError("system_state must have shape (12,)")

    planet_state = OrbitalState(
        position_km=system_state[0:3],
        velocity_km_s=system_state[3:6],
    )

    spacecraft_state = OrbitalState(
        position_km=system_state[6:9],
        velocity_km_s=system_state[9:12],
    )

    if position_intersects_body(
        spacecraft_state.position_km,
        sun_body,
        sun_state,
    ):
        raise ValueError("spacecraft intersects the Sun")

    if position_intersects_body(
        spacecraft_state.position_km,
        planet_body,
        planet_state,
    ):
        raise ValueError(f"spacecraft intersects {planet_body.name}")

    planet_acceleration = total_gravitational_acceleration(
        target_position_km=planet_state.position_km,
        bodies=[sun_body],
        body_states=[sun_state],
    )

    spacecraft_acceleration = total_gravitational_acceleration(
        target_position_km=spacecraft_state.position_km,
        bodies=[sun_body, planet_body],
        body_states=[sun_state, planet_state],
    )

    return np.concatenate((
        planet_state.velocity_km_s,
        planet_acceleration,
        spacecraft_state.velocity_km_s,
        spacecraft_acceleration,
    ))
