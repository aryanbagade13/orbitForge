"""A proposed burn schedule and arrival time for mission evaluation."""

from dataclasses import dataclass
from numbers import Real

import numpy as np

from .manoeuvres import ImpulsiveManoeuvre


@dataclass(frozen=True)
class MissionCandidate:
    """Times are seconds since departure; delta-v is measured in km/s.

    Frozen fields prevent reassignment, but contained manoeuvres remain mutable.
    """

    manoeuvres: tuple[ImpulsiveManoeuvre, ...]
    arrival_time_s: float

    def __post_init__(self):
        if isinstance(self.arrival_time_s, (bool, np.bool_)) or not isinstance(
            self.arrival_time_s, Real
        ):
            raise TypeError("arrival_time_s must be a real number")
        if not np.isfinite(self.arrival_time_s) or self.arrival_time_s <= 0:
            raise ValueError("arrival_time_s must be positive and finite")

        manoeuvres = tuple(self.manoeuvres)
        for manoeuvre in manoeuvres:
            if not isinstance(manoeuvre, ImpulsiveManoeuvre):
                raise TypeError("each manoeuvre must be an ImpulsiveManoeuvre")
            if not np.isfinite(manoeuvre.time_s) or manoeuvre.time_s < 0:
                raise ValueError("manoeuvre times must be non-negative and finite")
            if manoeuvre.time_s > self.arrival_time_s:
                raise ValueError("manoeuvre time must not be after arrival")

        object.__setattr__(self, "arrival_time_s", float(self.arrival_time_s))
        object.__setattr__(
            self,
            "manoeuvres",
            tuple(sorted(manoeuvres, key=lambda manoeuvre: manoeuvre.time_s)),
        )

    @property
    def total_delta_v_km_s(self) -> float:
        """Sum burn magnitudes, including burns in opposing directions."""
        return sum(
            (manoeuvre.magnitude_km_s for manoeuvre in self.manoeuvres),
            0.0,
        )
