"""Mission propagation with scheduled impulsive manoeuvres."""

from collections.abc import Sequence

from .results import MissionResult

import numpy as np

from .integrators import DerivativeFunction, StateVector, propagate_fixed_step
from .manoeuvres import ImpulsiveManoeuvre
from .simulation import apply_spacecraft_manoeuvre_to_system_state


def propagate_with_manoeuvres(
    initial_state: StateVector,
    start_time_s: float,
    end_time_s: float,
    dt_s: float,
    derivative_function: DerivativeFunction,
    manoeuvres: Sequence[ImpulsiveManoeuvre],
    *,
    state_layout: str = "planet_spacecraft",
    reference_frame: str | None = None,
) -> MissionResult:
    """Propagate a declared state layout; burn-time samples are post-burn."""
    current_state = np.asarray(initial_state, dtype=float).copy()

    widths = {"planet_spacecraft": 12, "spacecraft": 6}
    if state_layout not in widths:
        raise ValueError("unknown state_layout")
    if current_state.shape != (widths[state_layout],):
        raise ValueError(f"initial_state must have shape ({widths[state_layout]},)")
    if not np.isfinite(current_state).all():
        raise ValueError("initial_state must contain only finite values")
    if not np.isfinite(start_time_s) or not np.isfinite(end_time_s):
        raise ValueError("mission times must be finite")
    if end_time_s < start_time_s:
        raise ValueError("end_time_s must not be before start_time_s")
    if not np.isfinite(dt_s) or dt_s <= 0:
        raise ValueError("dt_s must be positive and finite")

    ordered_manoeuvres = sorted(
        manoeuvres,
        key=lambda manoeuvre: manoeuvre.time_s,
    )
    for manoeuvre in ordered_manoeuvres:
        if not start_time_s <= manoeuvre.time_s <= end_time_s:
            raise ValueError("manoeuvre time must lie within the mission")

    current_time_s = start_time_s
    all_times_s = [start_time_s]
    all_states = [current_state.copy()]

    for manoeuvre in ordered_manoeuvres:
        segment_times_s, segment_states = propagate_fixed_step(
            initial_state=current_state,
            start_time_s=current_time_s,
            end_time_s=manoeuvre.time_s,
            dt_s=dt_s,
            derivative_function=derivative_function,
        )

        all_times_s.extend(segment_times_s[1:])
        all_states.extend(segment_states[1:])

        if state_layout == "planet_spacecraft":
            current_state = apply_spacecraft_manoeuvre_to_system_state(
                segment_states[-1], manoeuvre
            )
        else:
            current_state = segment_states[-1].copy()
            current_state[3:6] += manoeuvre.delta_velocity_km_s
        current_time_s = manoeuvre.time_s

        # Store the post-burn state at the manoeuvre time.
        all_states[-1] = current_state.copy()

    final_times_s, final_states = propagate_fixed_step(
        initial_state=current_state,
        start_time_s=current_time_s,
        end_time_s=end_time_s,
        dt_s=dt_s,
        derivative_function=derivative_function,
    )
    all_times_s.extend(final_times_s[1:])
    all_states.extend(final_states[1:])

    return MissionResult(
        times_s=all_times_s,
        states=all_states,
        manoeuvres=ordered_manoeuvres,
        state_layout=state_layout,
        reference_frame=reference_frame,
    )
