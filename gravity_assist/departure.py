"""Construct post-departure states in the ephemeris's reference frame."""

from .models import OrbitalState


def make_departure_state(ephemerides, body_name, relative_position_km,
                         relative_velocity_km_s):
    """Offsets use the same fixed axes as the ephemeris, not a rotating frame.

    Relative velocity is the actual velocity at this position, not v-infinity.
    This does not model launch, parking-orbit escape, or their propellant cost.
    """
    offset = OrbitalState(relative_position_km, relative_velocity_km_s)
    planet = ephemerides.get_state(body_name, 0.0)
    return OrbitalState(planet.position_km + offset.position_km,
                        planet.velocity_km_s + offset.velocity_km_s)
