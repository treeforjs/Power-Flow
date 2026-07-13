"""Trace loading and drive fitting."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Trace:
    time_s: np.ndarray
    values: dict[str, np.ndarray]

    def sample(self, name: str, t: float) -> float:
        return float(np.interp(t, self.time_s, self.values[name]))


def load_trace_csv(path: str | Path) -> Trace:
    table = np.genfromtxt(path, delimiter=",", names=True)
    if table.dtype.names is None or "time_s" not in table.dtype.names:
        raise ValueError("trace CSV must include a time_s column")
    time = np.asarray(table["time_s"], dtype=float)
    values = {name: np.asarray(table[name], dtype=float) for name in table.dtype.names if name != "time_s"}
    return Trace(time_s=time, values=values)


def fit_effective_rl(time_s: np.ndarray, voltage_v: np.ndarray, current_a: np.ndarray) -> dict[str, float]:
    """Least-squares fit V = R I + L dI/dt."""
    time_s = np.asarray(time_s, dtype=float)
    voltage_v = np.asarray(voltage_v, dtype=float)
    current_a = np.asarray(current_a, dtype=float)
    didt = np.gradient(current_a, time_s, edge_order=1)
    a = np.column_stack([current_a, didt])
    coeff, *_ = np.linalg.lstsq(a, voltage_v, rcond=None)
    pred = a @ coeff
    rms = float(np.sqrt(np.mean((pred - voltage_v) ** 2)))
    return {"resistance_ohm": float(coeff[0]), "inductance_h": float(coeff[1]), "rms_voltage_error_v": rms}


def synthesize_current_from_rl(
    time_s: np.ndarray,
    voltage_v: np.ndarray,
    resistance_ohm: float,
    inductance_h: float,
) -> np.ndarray:
    """Forward Euler current reconstruction for an R-L series element."""
    current = np.zeros_like(time_s, dtype=float)
    inductance_h = max(inductance_h, 1e-18)
    for k in range(1, len(time_s)):
        dt = time_s[k] - time_s[k - 1]
        di = (voltage_v[k - 1] - resistance_ohm * current[k - 1]) / inductance_h
        current[k] = current[k - 1] + dt * di
    return current
