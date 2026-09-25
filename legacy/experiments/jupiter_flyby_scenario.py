import numpy as np

from legacy.engine.mission import propagate_with_manoeuvres
from gravity_assist.results import MissionResult
from gravity_assist.models import CelestialBody, OrbitalState
from legacy.engine.simulation import (
    pack_system_state,
    restricted_three_body_derivative,
)

EPOCH = "2000-01-01 12:00:00 TDB"
SECONDS_PER_DAY = 24 * 60 * 60

sun_body = CelestialBody(
    name="Sun",
    mass_kg=1.989e30,
    radius_km=695_700.0,
)

sun_state = OrbitalState(
    position_km=[0.0, 0.0, 0.0],
    velocity_km_s=[0.0, 0.0, 0.0],
)

jupiter_body = CelestialBody(
    name="Jupiter",
    mass_kg=1.8986e27,
    radius_km=69_911.0,
)

jupiter_state = OrbitalState(
    position_km=[
        5.985676246570644e8,
        4.396046799481729e8,
        -1.522686167298746e7,
    ],
    velocity_km_s=[
        -7.909860292172008,
        11.15621752636729,
        0.1308656815823666,
    ],
)

spacecraft_position_relative_to_jupiter_km = np.array([
    -10_000_000.0,
    1_500_000.0,
    300_000.0,
])

hyperbolic_excess_speed_km_s = 5.0
initial_distance_from_jupiter_km = np.linalg.norm(
    spacecraft_position_relative_to_jupiter_km
)
initial_speed_relative_to_jupiter_km_s = np.sqrt(
    hyperbolic_excess_speed_km_s**2
    + 2 * jupiter_body.mu / initial_distance_from_jupiter_km
)
incoming_direction = np.array([1.0, 0.0, 0.0])
spacecraft_velocity_relative_to_jupiter_km_s = (
    initial_speed_relative_to_jupiter_km_s * incoming_direction
)

spacecraft_state = OrbitalState(
    position_km=(
        jupiter_state.position_km
        + spacecraft_position_relative_to_jupiter_km
    ),
    velocity_km_s=(
        jupiter_state.velocity_km_s
        + spacecraft_velocity_relative_to_jupiter_km_s
    ),
)

system_state = pack_system_state(jupiter_state, spacecraft_state)

start_time_s = 0.0
end_time_s = 50 * SECONDS_PER_DAY
dt_s = 60.0


def derivative_for_this_flyby(time_s, state):
    return restricted_three_body_derivative(
        time_s,
        state,
        sun_body,
        sun_state,
        jupiter_body,
    )


def propagate_jupiter_flyby(manoeuvres=()) -> MissionResult:
    return propagate_with_manoeuvres(
        initial_state=system_state,
        start_time_s=start_time_s,
        end_time_s=end_time_s,
        dt_s=dt_s,
        derivative_function=derivative_for_this_flyby,
        manoeuvres=manoeuvres,
    )


def calculate_flyby_comparison_metrics(times_s, system_states):
    planet_positions_km = system_states[:, 0:3]
    planet_velocities_km_s = system_states[:, 3:6]
    spacecraft_positions_km = system_states[:, 6:9]
    spacecraft_velocities_km_s = system_states[:, 9:12]

    relative_positions_km = spacecraft_positions_km - planet_positions_km
    relative_velocities_km_s = (
        spacecraft_velocities_km_s - planet_velocities_km_s
    )
    distances_km = np.linalg.norm(relative_positions_km, axis=1)
    closest_index = int(np.argmin(distances_km))

    if closest_index == 0 or closest_index == len(distances_km) - 1:
        raise ValueError("flyby closest approach must occur inside the mission")

    comparison_distance_km = 10_000_000.0
    incoming_index = int(np.argmin(
        np.abs(distances_km[:closest_index] - comparison_distance_km)
    ))
    outgoing_index = closest_index + int(np.argmin(
        np.abs(distances_km[closest_index:] - comparison_distance_km)
    ))

    incoming_relative_velocity_km_s = relative_velocities_km_s[incoming_index]
    outgoing_relative_velocity_km_s = relative_velocities_km_s[outgoing_index]
    incoming_relative_speed_km_s = np.linalg.norm(
        incoming_relative_velocity_km_s
    )
    outgoing_relative_speed_km_s = np.linalg.norm(
        outgoing_relative_velocity_km_s
    )
    turning_angle_cosine = np.dot(
        incoming_relative_velocity_km_s,
        outgoing_relative_velocity_km_s,
    ) / (
        incoming_relative_speed_km_s * outgoing_relative_speed_km_s
    )
    measured_turning_angle_degrees = np.degrees(
        np.arccos(np.clip(turning_angle_cosine, -1.0, 1.0))
    )

    incoming_heliocentric_speed_km_s = np.linalg.norm(
        spacecraft_velocities_km_s[incoming_index]
    )
    outgoing_heliocentric_speed_km_s = np.linalg.norm(
        spacecraft_velocities_km_s[outgoing_index]
    )

    return {
        "closest_approach_time_days": (
            times_s[closest_index] / SECONDS_PER_DAY
        ),
        "closest_approach_altitude_km": (
            distances_km[closest_index] - jupiter_body.radius_km
        ),
        "turning_angle_degrees": measured_turning_angle_degrees,
        "outgoing_relative_speed_km_s": outgoing_relative_speed_km_s,
        "heliocentric_speed_change_km_s": (
            outgoing_heliocentric_speed_km_s
            - incoming_heliocentric_speed_km_s
        ),
    }
