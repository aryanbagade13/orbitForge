from dataclasses import dataclass, field
from numbers import Real
from datetime import datetime
from .models import OrbitalState
import numpy as np


@dataclass(frozen=True)
class FlybyRequirement:
    """An intermediate encounter, within a mission-relative time window."""

    body_name: str
    time_bounds_s: tuple[float, float]
    max_distance_km: float

    def __post_init__(self):
        if not isinstance(self.body_name, str) or not self.body_name.strip():
            raise ValueError("flyby body_name must be nonblank")
        if len(self.time_bounds_s) != 2:
            raise ValueError("flyby time bounds must have two values")
        lo, hi = self.time_bounds_s
        if not np.isfinite([lo, hi]).all() or not 0 <= lo < hi:
            raise ValueError("flyby bounds must satisfy 0 <= start < end")
        if not np.isfinite(self.max_distance_km) or self.max_distance_km <= 0:
            raise ValueError("flyby max_distance_km must be positive and finite")
        object.__setattr__(self, "time_bounds_s", (float(lo), float(hi)))

@dataclass(frozen=True)
class MissionDefinition:
    departure_epoch: datetime
    initial_spacecraft_state: OrbitalState
    flyby_body_names: tuple[str, ...]
    destination_body_name: str
    arrival_time_bounds_s: tuple[float, float]
    destination_max_distance_km: float
    minimum_altitudes_km: dict[str, float]
    reference_frame: str = "ICRS axes; Solar System barycentric origin"
    flyby_requirements: tuple[FlybyRequirement, ...] = field(default_factory=tuple)

    def __post_init__(self):
        requirements = tuple(self.flyby_requirements)
        if not all(isinstance(item, FlybyRequirement) for item in requirements):
            raise TypeError("flyby_requirements must contain FlybyRequirement objects")
        object.__setattr__(self, "flyby_requirements", requirements)
        if not isinstance(self.reference_frame, str) or not self.reference_frame.strip():
            raise ValueError("reference_frame must be a nonblank string")
        if len(self.arrival_time_bounds_s) != 2:
            raise ValueError(
                "arrival_time_bounds_s must contain two values"
            )

        earliest_s, latest_s = self.arrival_time_bounds_s

        if not np.isfinite([earliest_s, latest_s]).all():
            raise ValueError("arrival time bounds must be finite")

        if not 0 < earliest_s <= latest_s:
            raise ValueError(
                "arrival time bounds must satisfy 0 < earliest <= latest"
            )

        object.__setattr__(
            self,
            "arrival_time_bounds_s",
            (float(earliest_s), float(latest_s)),
        )

        if not isinstance(self.departure_epoch, datetime):
            raise TypeError("departure_epoch must be a datetime")

        if (
                self.departure_epoch.tzinfo is None
                or self.departure_epoch.utcoffset() is None
        ):
            raise ValueError("departure_epoch must be timezone-aware")

        if not isinstance(self.initial_spacecraft_state, OrbitalState):
            raise TypeError(
                "initial_spacecraft_state must be an OrbitalState"
            )

        if not isinstance(self.destination_body_name, str):
            raise TypeError("destination_body_name must be a string")

        if not self.destination_body_name.strip():
            raise ValueError("destination_body_name must not be blank")

        if isinstance(self.flyby_body_names, str):
            raise TypeError(
                "flyby_body_names must be a collection of names"
            )

        flyby_names = tuple(self.flyby_body_names)
        for name in flyby_names:
            if not isinstance(name, str):
                raise TypeError("each flyby body name must be a string")
            if not name.strip():
                raise ValueError("flyby body names must not be blank")

        object.__setattr__(self, "flyby_body_names", flyby_names)

        if (
            not np.isfinite(self.destination_max_distance_km)
            or self.destination_max_distance_km <= 0
        ):
            raise ValueError(
                "destination_max_distance_km must be positive and finite"
            )

        if not isinstance(self.minimum_altitudes_km, dict):
            raise TypeError("minimum_altitudes_km must be a dictionary")

        for body_name, altitude_km in self.minimum_altitudes_km.items():
            if not isinstance(body_name, str):
                raise TypeError("minimum-altitude body names must be strings")
            if not body_name.strip():
                raise ValueError("minimum-altitude body names must not be blank")
            if isinstance(altitude_km, (bool, np.bool_)) or not isinstance(
                altitude_km, Real
            ):
                raise TypeError("minimum altitudes must be real numbers")
            if not np.isfinite(altitude_km) or altitude_km < 0:
                raise ValueError(
                    "minimum altitudes must be non-negative and finite"
                )
