from dataclasses import dataclass, field
from collections.abc import Callable
from datetime import datetime

import numpy as np

from gravity_assist.manoeuvres import ImpulsiveManoeuvre


@dataclass
class MissionResult:
    times_s: np.ndarray
    states: np.ndarray
    manoeuvres: tuple[ImpulsiveManoeuvre, ...]
    state_layout: str = "planet_spacecraft"
    reference_frame: str | None = None
    departure_epoch: datetime | None = None
    dense_state: Callable | None = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        self.times_s = np.asarray(self.times_s, dtype=float).copy()
        self.states = np.asarray(self.states, dtype=float).copy()
        self.manoeuvres = tuple(self.manoeuvres)

        if self.times_s.ndim != 1:
            raise ValueError("times_s must be one-dimensional")

        widths = {"planet_spacecraft": 12, "spacecraft": 6}
        if self.state_layout not in widths:
            raise ValueError("unknown state_layout")
        width = widths[self.state_layout]
        if self.states.ndim != 2 or self.states.shape[1] != width:
            raise ValueError(f"states must have shape (n, {width})")

        if len(self.times_s) != len(self.states):
            raise ValueError(
                "times_s and states must have matching lengths"
            )

        if not np.isfinite(self.times_s).all():
            raise ValueError("times_s must contain only finite values")

        if not np.isfinite(self.states).all():
            raise ValueError("states must contain only finite values")

        if np.any(np.diff(self.times_s) <= 0):
            raise ValueError("times_s must be strictly increasing")

    @property
    def spacecraft_states(self) -> np.ndarray:
        """Spacecraft [x, y, z, vx, vy, vz], independent of stored layout."""
        return self.states if self.state_layout == "spacecraft" else self.states[:, 6:12]

    def spacecraft_state_at(self, time_s: float) -> np.ndarray:
        """Continuous six-component state, with right-hand velocity at burns."""
        if not np.isfinite(time_s) or not self.times_s[0] <= time_s <= self.times_s[-1]:
            raise ValueError("requested time lies outside the trajectory")
        if self.dense_state is None:
            raise ValueError("result has no continuous trajectory; use adaptive propagation")
        return np.asarray(self.dense_state(time_s), dtype=float).copy()
