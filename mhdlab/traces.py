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


@dataclass
class DriveProfile:
    time_s: np.ndarray
    current_a: np.ndarray
    voltage_v: np.ndarray
    metadata: dict[str, float | str]

    def sample_current(self, t: float) -> float:
        return float(np.interp(t, self.time_s, self.current_a))

    def sample_voltage(self, t: float) -> float:
        return float(np.interp(t, self.time_s, self.voltage_v))


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


def parametric_current(
    time_s: np.ndarray,
    peak_current_a: float,
    rise_time_s: float,
    waveform: str = "half_sine",
) -> np.ndarray:
    """Build a current waveform when measured load current is unavailable."""
    time_s = np.asarray(time_s, dtype=float)
    peak_current_a = float(peak_current_a)
    rise_time_s = float(rise_time_s)
    if rise_time_s <= 0.0:
        raise ValueError("rise_time_s must be positive")

    name = waveform.lower().replace("-", "_")
    tau = np.maximum(time_s / rise_time_s, 0.0)
    if name in {"half_sine", "quarter_period_sine", "quarter_sine"}:
        current = peak_current_a * np.sin(0.5 * np.pi * tau)
        return np.where((tau >= 0.0) & (tau <= 2.0), np.maximum(current, 0.0), 0.0)
    if name in {"quarter_sine_hold", "sine_rise_hold"}:
        return peak_current_a * np.sin(0.5 * np.pi * np.minimum(tau, 1.0))
    if name in {"linear_ramp_hold", "linear"}:
        return peak_current_a * np.minimum(tau, 1.0)
    raise ValueError(f"unsupported parametric current waveform: {waveform}")


def drive_from_trace(
    trace: Trace,
    voltage_column: str,
    current_column: str,
    mode: str,
) -> DriveProfile:
    rl_fit = fit_effective_rl(trace.time_s, trace.values[voltage_column], trace.values[current_column])
    current = trace.values[current_column]
    if mode == "fit_rl":
        current = synthesize_current_from_rl(
            trace.time_s,
            trace.values[voltage_column],
            rl_fit["resistance_ohm"],
            rl_fit["inductance_h"],
        )
    metadata: dict[str, float | str] = {"mode": mode, **rl_fit}
    return DriveProfile(
        time_s=trace.time_s,
        current_a=np.asarray(current, dtype=float),
        voltage_v=np.asarray(trace.values[voltage_column], dtype=float),
        metadata=metadata,
    )


def drive_from_parametric_current(time_s: np.ndarray, config: dict) -> DriveProfile:
    peak_current_a = float(config["peak_current_a"])
    rise_time_s = float(config["rise_time_s"])
    waveform = str(config.get("waveform", "half_sine"))
    current = parametric_current(time_s, peak_current_a, rise_time_s, waveform)
    voltage = np.full_like(time_s, float(config.get("voltage_v", 0.0)), dtype=float)
    return DriveProfile(
        time_s=np.asarray(time_s, dtype=float),
        current_a=current,
        voltage_v=voltage,
        metadata={
            "mode": "parametric_current",
            "waveform": waveform,
            "peak_current_a": peak_current_a,
            "rise_time_s": rise_time_s,
            "voltage_v": float(config.get("voltage_v", 0.0)),
            "resistance_ohm": 0.0,
            "inductance_h": 0.0,
            "rms_voltage_error_v": 0.0,
        },
    )
