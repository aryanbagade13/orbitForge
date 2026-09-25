import numpy as np

def heliocentric_state(provider, body_name, time_s):
    body = provider.get_state(body_name, time_s)
    sun = provider.get_state("Sun", time_s)

    position_km = body.position_km - sun.position_km
    velocity_km_s = body.velocity_km_s - sun.velocity_km_s

    return position_km, velocity_km_s


def transfer_endpoints(provider, departure_body, arrival_body, flight_time_s):
    """Return Sun-relative planet states at departure and arrival.

    Positions are in km; velocities are in km/s.
    """
    if flight_time_s <= 0:
        raise ValueError("flight_time_s must be positive")

    departure_position, departure_velocity = heliocentric_state(
        provider, departure_body, 0.0
    )
    arrival_position, arrival_velocity = heliocentric_state(
        provider, arrival_body, flight_time_s
    )

    return (
        departure_position,
        departure_velocity,
        arrival_position,
        arrival_velocity,
    )

def relative_speed(spacecraft_velocity_km_s: np.ndarray, planet_velocity_km_s: np.ndarray):
    v_rel = planet_velocity_km_s - spacecraft_velocity_km_s
    magnitude_v_rel = np.linalg.norm(v_rel)
    return magnitude_v_rel


def hohmann_transfer_time(mu, departure_radius_km, arrival_radius_km):
    """Returns the half-ellipse transfer time between circular, coplanar orbits.

    mu is the central body's gravitational parameter in km^3/s^2.
    Both orbital radii are in km; the returned duration is in seconds.
    """
    if mu <= 0 or departure_radius_km <= 0 or arrival_radius_km <= 0:
        raise ValueError("mu and both orbital radii must be positive")

    semi_major_axis = (departure_radius_km + arrival_radius_km) / 2
    return np.pi * np.sqrt(semi_major_axis**3 / mu)
