import matplotlib.pyplot as plt
import numpy as np

from gravity_assist.manoeuvres import ImpulsiveManoeuvre
from legacy.engine.mission import propagate_with_manoeuvres
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


system_state = pack_system_state(
    jupiter_state,
    spacecraft_state,
)


def derivative_for_this_flyby(time_s, state):
    return restricted_three_body_derivative(
        time_s,
        state,
        sun_body,
        sun_state,
        jupiter_body,
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
        np.abs(
            distances_km[:closest_index] - comparison_distance_km
        )
    ))
    outgoing_index = closest_index + int(np.argmin(
        np.abs(
            distances_km[closest_index:] - comparison_distance_km
        )
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


start_time_s = 0.0
end_time_s = 50 * SECONDS_PER_DAY
dt_s = 60.0

result = propagate_with_manoeuvres(
    initial_state=system_state,
    start_time_s=start_time_s,
    end_time_s=end_time_s,
    dt_s=dt_s,
    derivative_function=derivative_for_this_flyby,
    manoeuvres=[],
)

times = result.times_s
states = result.states
baseline_metrics = calculate_flyby_comparison_metrics(times, states)

trial_manoeuvre = ImpulsiveManoeuvre(
    time_s=5 * SECONDS_PER_DAY,
    delta_velocity_km_s=np.array([0.0, 0.01, 0.0]),
)
manoeuvred_result = propagate_with_manoeuvres(
    initial_state=system_state,
    start_time_s=start_time_s,
    end_time_s=end_time_s,
    dt_s=dt_s,
    derivative_function=derivative_for_this_flyby,
    manoeuvres=[trial_manoeuvre],
)
manoeuvred_times = manoeuvred_result.times_s
manoeuvred_states = manoeuvred_result.states
manoeuvred_metrics = calculate_flyby_comparison_metrics(
    manoeuvred_times,
    manoeuvred_states,
)

jupiter_positions_km = states[:, 0:3]
spacecraft_positions_km = states[:, 6:9]
jupiter_velocities_km_s = states[:, 3:6]
spacecraft_velocities_km_s = states[:, 9:12]

spacecraft_velocities_relative_to_jupiter_km_s = (
    spacecraft_velocities_km_s - jupiter_velocities_km_s
)

spacecraft_positions_relative_to_jupiter_km = (
    spacecraft_positions_km - jupiter_positions_km
)
spacecraft_distances_from_jupiter_km = np.linalg.norm(
    spacecraft_positions_relative_to_jupiter_km,
    axis=1,
)

closest_approach_index = np.argmin(
    spacecraft_distances_from_jupiter_km
)

closest_approach_distance_km = (
    spacecraft_distances_from_jupiter_km[closest_approach_index]
)

closest_approach_altitude_km = (
    closest_approach_distance_km - jupiter_body.radius_km
)

closest_approach_time_days = (
    times[closest_approach_index] / SECONDS_PER_DAY
)

# Compare the incoming and outgoing motion at the same Jupiter distance.
comparison_distance_km = 10_000_000.0

incoming_index = np.argmin(
    np.abs(
        spacecraft_distances_from_jupiter_km[:closest_approach_index]
        - comparison_distance_km
    )
)

outgoing_index = closest_approach_index + np.argmin(
    np.abs(
        spacecraft_distances_from_jupiter_km[closest_approach_index:]
        - comparison_distance_km
    )
)

incoming_velocity_relative_to_jupiter_km_s = (
    spacecraft_velocities_relative_to_jupiter_km_s[incoming_index]
)
outgoing_velocity_relative_to_jupiter_km_s = (
    spacecraft_velocities_relative_to_jupiter_km_s[outgoing_index]
)

incoming_speed_relative_to_jupiter_km_s = np.linalg.norm(
    incoming_velocity_relative_to_jupiter_km_s
)
outgoing_speed_relative_to_jupiter_km_s = np.linalg.norm(
    outgoing_velocity_relative_to_jupiter_km_s
)

turning_angle_cosine = np.dot(
    incoming_velocity_relative_to_jupiter_km_s,
    outgoing_velocity_relative_to_jupiter_km_s,
) / (
    incoming_speed_relative_to_jupiter_km_s
    * outgoing_speed_relative_to_jupiter_km_s
)
turning_angle_degrees = np.degrees(
    np.arccos(np.clip(turning_angle_cosine, -1.0, 1.0))
)

theoretical_eccentricity = 1 + closest_approach_distance_km*hyperbolic_excess_speed_km_s**2 / jupiter_body.mu

theoretical_turning_angle_degrees = np.degrees(2*np.arcsin(1/theoretical_eccentricity))

incoming_heliocentric_speed_km_s = np.linalg.norm(
    spacecraft_velocities_km_s[incoming_index]
)
outgoing_heliocentric_speed_km_s = np.linalg.norm(
    spacecraft_velocities_km_s[outgoing_index]
)
heliocentric_speed_change_km_s = (
    outgoing_heliocentric_speed_km_s
    - incoming_heliocentric_speed_km_s
)

# Construct a correctly scaled spherical surface for Jupiter.
longitude = np.linspace(0, 2 * np.pi, 80)
latitude = np.linspace(0, np.pi, 40)

# Use spherical coordinates.
jupiter_surface_x_km = (
    jupiter_body.radius_km
    * np.outer(np.cos(longitude), np.sin(latitude))
)

jupiter_surface_y_km = (
    jupiter_body.radius_km
    * np.outer(np.sin(longitude), np.sin(latitude))
)

jupiter_surface_z_km = (
    jupiter_body.radius_km
    * np.outer(np.ones_like(longitude), np.cos(latitude))
)

# Select only the nearby section of the trajectory for the close-up.
close_up_limit_km = 500_000.0
close_up_mask = (
    spacecraft_distances_from_jupiter_km <= close_up_limit_km
)
close_up_positions_km = (
    spacecraft_positions_relative_to_jupiter_km[close_up_mask]
)

closest_position_relative_to_jupiter_km = (
    spacecraft_positions_relative_to_jupiter_km[
        closest_approach_index
    ]
)

figure = plt.figure(figsize=(9, 8))
axes = figure.add_subplot(projection="3d")

axes.plot_surface(
    jupiter_surface_x_km,
    jupiter_surface_y_km,
    jupiter_surface_z_km,
    color="orange",
    alpha=0.8,
    linewidth=0,
)

axes.plot(
    close_up_positions_km[:, 0],
    close_up_positions_km[:, 1],
    close_up_positions_km[:, 2],
    color="blue",
    label="Spacecraft trajectory",
)

axes.scatter(
    closest_position_relative_to_jupiter_km[0],
    closest_position_relative_to_jupiter_km[1],
    closest_position_relative_to_jupiter_km[2],
    color="red",
    s=50,
    label="Closest approach",
)

axes.set_xlim(-close_up_limit_km, close_up_limit_km)
axes.set_ylim(-close_up_limit_km, close_up_limit_km)
axes.set_zlim(-close_up_limit_km, close_up_limit_km)
axes.set_box_aspect((1, 1, 1))
axes.set_xlabel("x relative to Jupiter (km)")
axes.set_ylabel("y relative to Jupiter (km)")
axes.set_zlabel("z relative to Jupiter (km)")
axes.set_title("Close-up of Jupiter gravity-assist trajectory")
axes.legend()
figure.tight_layout()
plt.show()

print(f"Epoch: {EPOCH}")
print(
    "Chosen hyperbolic-excess speed: "
    f"{hyperbolic_excess_speed_km_s:.3f} km/s"
)
print(
    "Initial Jupiter-relative speed: "
    f"{initial_speed_relative_to_jupiter_km_s:.3f} km/s"
)
print(f"Closest-approach time: {closest_approach_time_days:.3f} days")
print(f"Distance from Jupiter's centre: {closest_approach_distance_km:,.1f} km")
print(f"Altitude above Jupiter's surface: {closest_approach_altitude_km:,.1f} km")
print(
    "Incoming comparison distance: "
    f"{spacecraft_distances_from_jupiter_km[incoming_index]:,.1f} km"
)
print(
    "Outgoing comparison distance: "
    f"{spacecraft_distances_from_jupiter_km[outgoing_index]:,.1f} km"
)
print(
    "Incoming Jupiter-relative speed: "
    f"{incoming_speed_relative_to_jupiter_km_s:.6f} km/s"
)
print(
    "Outgoing Jupiter-relative speed: "
    f"{outgoing_speed_relative_to_jupiter_km_s:.6f} km/s"
)
print(f"Turning angle: {turning_angle_degrees:.3f} degrees")
print(
    "Heliocentric speed change: "
    f"{heliocentric_speed_change_km_s:+.6f} km/s"
)
print(f"Stored times shape: {times.shape}")
print(f"Stored states shape: {states.shape}")
print(
    "Jupiter-relative velocity shape:",
    spacecraft_velocities_relative_to_jupiter_km_s.shape,
)

print(f"Absolute difference between theoretical turning angle and measured is: {abs(theoretical_turning_angle_degrees - turning_angle_degrees): .3f} degrees")

print("\nTimed-manoeuvre comparison")
print(
    "Burn: day "
    f"{trial_manoeuvre.time_s / SECONDS_PER_DAY:.1f}, "
    f"delta-v = {trial_manoeuvre.magnitude_km_s * 1000:.1f} m/s "
    "in the +y direction"
)
print(
    "Closest-approach altitude: "
    f"{baseline_metrics['closest_approach_altitude_km']:,.1f} km "
    "without burn, "
    f"{manoeuvred_metrics['closest_approach_altitude_km']:,.1f} km "
    "with burn"
)
print(
    "Change in closest-approach altitude: "
    f"{manoeuvred_metrics['closest_approach_altitude_km'] - baseline_metrics['closest_approach_altitude_km']:+,.1f} km"
)
print(
    "Turning angle: "
    f"{baseline_metrics['turning_angle_degrees']:.3f} degrees "
    "without burn, "
    f"{manoeuvred_metrics['turning_angle_degrees']:.3f} degrees "
    "with burn"
)
print(
    "Outgoing Jupiter-relative speed: "
    f"{baseline_metrics['outgoing_relative_speed_km_s']:.6f} km/s "
    "without burn, "
    f"{manoeuvred_metrics['outgoing_relative_speed_km_s']:.6f} km/s "
    "with burn"
)
print(
    "Heliocentric speed change: "
    f"{baseline_metrics['heliocentric_speed_change_km_s']:+.6f} km/s "
    "without burn, "
    f"{manoeuvred_metrics['heliocentric_speed_change_km_s']:+.6f} km/s "
    "with burn"
)
