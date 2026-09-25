"""Shared ephemeris-driven forces and clearance checks for the current engine.

Use adaptive.simulate_adaptive_mission or evaluation.assess_candidate to run
a mission. The old fixed-step wrapper lives in legacy.engine.interplanetary.
"""

from collections.abc import Mapping, Sequence

import numpy as np

from .ephemerides import EphemerisProvider
from .forces import total_gravitational_acceleration
from .models import CelestialBody


class UnsafeTrajectoryError(ValueError):
    """An evaluated state touches a body or violates its clearance."""


def make_ephemeris_derivative(
    bodies: Sequence[CelestialBody],
    ephemerides: EphemerisProvider,
    minimum_altitudes_km: Mapping[str, float] | None = None,
    *, check_clearance: bool = True,
):
    """Build dr/dt=v, dv/dt=sum(gravity) for a six-component spacecraft."""
    bodies = tuple(bodies)
    names = [body.name.strip().lower() for body in bodies]
    if not bodies or len(set(names)) != len(names):
        raise ValueError("bodies must be nonempty with unique names")
    limits = {}
    for name, altitude in (minimum_altitudes_km or {}).items():
        key = name.strip().lower()
        if key not in names or key in limits:
            raise ValueError("clearance names must uniquely match supplied bodies")
        if not np.isfinite(altitude) or altitude < 0:
            raise ValueError("minimum altitudes must be non-negative and finite")
        limits[key] = float(altitude)

    def derivative(time_s, state):
        state = np.asarray(state, dtype=float)
        if state.shape != (6,) or not np.isfinite(state).all():
            raise ValueError("spacecraft state must be finite with shape (6,)")
        body_states = [ephemerides.get_state(name, time_s) for name in names]
        for body, name, body_state in zip(bodies, names, body_states):
            distance = np.linalg.norm(state[:3] - body_state.position_km)
            if check_clearance and (distance <= body.radius_km or distance < body.radius_km + limits.get(name, 0.0)):
                raise UnsafeTrajectoryError(
                    f"spacecraft violates clearance for {body.name} at {time_s} s"
                )
        acceleration = total_gravitational_acceleration(state[:3], bodies, body_states)
        return np.concatenate((state[3:6], acceleration))

    return derivative
