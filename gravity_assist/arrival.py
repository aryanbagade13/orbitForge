import numpy as np

def periapsis_speed(mu, periapsis_km, apoapsis_km):
    if mu <= 0 or periapsis_km <=0 or apoapsis_km <= 0:
        raise ValueError("Mu, periapsis_km and apoapsis_km must be positive")
    if apoapsis_km < periapsis_km:
        raise ValueError("apoapsis_km must be larger than or equal to periapsis_km")

    semi_major_axis = (periapsis_km+apoapsis_km)/2
    velocity_periapsis = np.sqrt(mu*(2/periapsis_km - 1/semi_major_axis))
    return velocity_periapsis

def arrival_periapsis_speed(mu, periapsis_km, v_infinity_km_s):
    if v_infinity_km_s < 0:
        raise ValueError("v_infinity_km_s must be positive or 0")
    if mu <= 0 or periapsis_km <= 0:
        raise ValueError("Mu, periapsis_km must be positive")

    v_incoming_km_s = np.sqrt(v_infinity_km_s**2 + 2*mu/periapsis_km)
    return v_incoming_km_s

def capture_delta_v(mu, periapsis_km, apoapsis_km, v_infinity_km_s):
    delta_v = arrival_periapsis_speed(mu, periapsis_km, v_infinity_km_s) - periapsis_speed(mu, periapsis_km, apoapsis_km)
    return delta_v

def arrival_delta_v(arrival_mode: str, mu, periapsis_km, apoapsis_km, v_infinity_km_s):
    if arrival_mode == "flyby":
        return 0.0
    elif arrival_mode == "orbit_insertion":
        return capture_delta_v(mu, periapsis_km, apoapsis_km, v_infinity_km_s)
    else:
        raise ValueError("Arrival mode must be either flyby or orbit_insertion")