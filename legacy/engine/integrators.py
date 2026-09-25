"""Archived fixed-step RK4 implementation, retained for learning and checks."""

from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray


StateVector = NDArray[np.float64]
DerivativeFunction = Callable[[float, StateVector], StateVector]


def _evaluate_derivative(
    derivative_function: DerivativeFunction,
    time_s: float,
    state: StateVector,
) -> StateVector:
    derivative = np.asarray(derivative_function(time_s, state), dtype=float)

    if derivative.shape != state.shape:
        raise ValueError("derivative must have the same shape as state")
    if not np.isfinite(derivative).all():
        raise ValueError("derivative must contain only finite values")

    return derivative


def rk4_step(
    state: StateVector,
    time_s: float,
    dt_s: float,
    derivative_function: DerivativeFunction,
) -> StateVector:
    """Advance a state vector by one fourth-order Runge-Kutta step."""
    state = np.asarray(state, dtype=float)

    if state.ndim != 1:
        raise ValueError("state must be a one-dimensional array")
    if not np.isfinite(state).all():
        raise ValueError("state must contain only finite values")
    if not np.isfinite(time_s):
        raise ValueError("time_s must be finite")
    if not np.isfinite(dt_s) or dt_s <= 0:
        raise ValueError("dt_s must be positive and finite")

    k1 = _evaluate_derivative(derivative_function, time_s, state)
    k2 = _evaluate_derivative(
        derivative_function,
        time_s + dt_s / 2,
        state + dt_s * k1 / 2,
    )
    k3 = _evaluate_derivative(
        derivative_function,
        time_s + dt_s / 2,
        state + dt_s * k2 / 2,
    )
    k4 = _evaluate_derivative(
        derivative_function,
        time_s + dt_s,
        state + dt_s * k3,
    )

    return state + dt_s * (k1 + 2 * k2 + 2 * k3 + k4) / 6

def propagate_fixed_step(initial_state, start_time_s, end_time_s, dt_s, derivative_function: DerivativeFunction):
    current_state = np.asarray(initial_state, dtype=float).copy()
    current_time = start_time_s
    times = [start_time_s]
    states = [current_state.copy()]

    while current_time < end_time_s:
        step_s = min(dt_s, end_time_s - current_time)
        current_state = rk4_step(current_state, current_time, step_s, derivative_function)

        current_time += step_s

        times.append(current_time)
        states.append(current_state.copy())

    return np.asarray(times), np.asarray(states)
