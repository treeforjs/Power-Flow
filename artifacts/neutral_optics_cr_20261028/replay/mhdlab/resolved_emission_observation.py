"""Observe externally supplied level populations through explicit ray weights.

This module evolves no population, thermal or MHD state. Its only optical model
is optically thin emission with authored Gaussian line profiles and a rectangular
gate over the piecewise-linear interpolant of supplied population samples.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Real

import numpy as np
from scipy.special import ndtr


_H = 6.62607015e-34
_C = 299792458.0


def _label(value, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _number(value, name: str, *, positive: bool = False, signed: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not np.isfinite(value) or (value <= 0 if positive else not signed and value < 0):
        bound = " and positive" if positive else "" if signed else " and nonnegative"
        raise ValueError(f"{name} must be finite{bound}")
    return value


def _array(value, name: str, ndim: int) -> np.ndarray:
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.ndim != ndim or not raw.size:
        raise ValueError(f"{name} must be a nonempty {ndim}-dimensional numeric array")
    array = np.array(raw, dtype=float, copy=True)
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class ResolvedTransition:
    transition_id: str
    upper_state_id: str
    vacuum_wavelength_nm: float
    einstein_a_per_s: float

    def __post_init__(self):
        _label(self.transition_id, "transition_id")
        _label(self.upper_state_id, "upper_state_id")
        _number(self.vacuum_wavelength_nm, "vacuum_wavelength_nm", positive=True)
        _number(self.einstein_a_per_s, "einstein_a_per_s", positive=True)


def project_resolved_emission(
    time_s,
    upper_state_density_m3: dict[str, np.ndarray],
    path_weights_m,
    transitions: list[ResolvedTransition],
    bin_edges_nm,
    gate_start_s: float,
    gate_end_s: float,
    fwhm_nm,
    *,
    time_registration: str,
    population_provenance: str,
    geometry_provenance: str,
    status: str,
    optical_model: str,
) -> dict:
    """Return gate-averaged bin energy radiance from mapped upper-level states.

    ``path_weights_m`` contains simulation ray-cell lengths multiplied by
    dimensionless normalized radiance-collection weights, applied once. These
    are not detector etendue or a fitted depth. Geometry is fixed during the
    gate; the caller owns collection-weight normalization. Each
    population array is ``(time, cell)`` in m^-3. Upper states may feed multiple
    uniquely identified transitions; no band-total population is substituted.
    The caller registers the source clock and supplies all population physics.
    Registered times and gate endpoints may be signed; no origin shift occurs.
    """
    for name, value in (("time_registration", time_registration),
                        ("population_provenance", population_provenance),
                        ("geometry_provenance", geometry_provenance),
                        ("status", status), ("optical_model", optical_model)):
        _label(value, name)
    if status not in {"manufactured", "owner_supplied", "source_evaluated"}:
        raise ValueError("status must be manufactured, owner_supplied or source_evaluated")
    if optical_model != "optically_thin":
        raise ValueError("only explicit optical_model='optically_thin' is supported")
    times = _array(time_s, "time_s", 1)
    if times.size < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("time_s needs at least two finite, strictly increasing samples")
    start = _number(gate_start_s, "gate_start_s", signed=True)
    stop = _number(gate_end_s, "gate_end_s", signed=True)
    if stop <= start or start < times[0] or stop > times[-1]:
        raise ValueError("a positive-width gate must be fully contained in time_s; extrapolation is unsupported")
    weights = _array(path_weights_m, "path_weights_m", 1)
    if np.any(weights < 0):
        raise ValueError("path_weights_m must be nonnegative")
    edges = _array(bin_edges_nm, "bin_edges_nm", 1)
    if edges.size < 2 or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
        raise ValueError("bin_edges_nm must be positive and strictly increasing")
    if not isinstance(transitions, list) or not transitions or any(
        not isinstance(item, ResolvedTransition) for item in transitions
    ):
        raise ValueError("transitions must be a nonempty list of ResolvedTransition values")
    identifiers = [item.transition_id for item in transitions]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("transition IDs must be unique")
    if not isinstance(upper_state_density_m3, dict) or not upper_state_density_m3:
        raise ValueError("explicit keyed upper-state populations are required")
    populations = {}
    for state_id, values in upper_state_density_m3.items():
        _label(state_id, "upper-state population ID")
        array = _array(values, f"population {state_id}", 2)
        if array.shape != (times.size, weights.size) or np.any(array < 0):
            raise ValueError(f"population {state_id} must be nonnegative with shape (time, cell)")
        populations[state_id] = array
    missing = {item.upper_state_id for item in transitions} - populations.keys()
    if missing:
        raise ValueError(f"missing resolved upper-state populations: {sorted(missing)}")
    widths_raw = np.asarray(fwhm_nm)
    if widths_raw.ndim == 0:
        widths = np.full(len(transitions), _number(widths_raw.item(), "fwhm_nm", positive=True))
    else:
        widths = _array(fwhm_nm, "fwhm_nm", 1)
        if widths.shape != (len(transitions),) or np.any(widths <= 0):
            raise ValueError("fwhm_nm requires one positive width per transition")

    knots = np.concatenate(([start], times[(times > start) & (times < stop)], [stop]))
    indices = np.clip(np.searchsorted(times, knots, side="right") - 1, 0, times.size - 2)
    fractions = (knots - times[indices]) / (times[indices + 1] - times[indices])
    intervals = np.diff(knots) / (stop - start)
    means = {}
    bin_radiance = np.zeros(edges.size - 1)
    lines = []
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            for state_id, values in populations.items():
                sampled = ((1 - fractions[:, None]) * values[indices]
                           + fractions[:, None] * values[indices + 1])
                means[state_id] = np.sum(
                    (0.5 * sampled[:-1] + 0.5 * sampled[1:]) * intervals[:, None], axis=0,
                )
            for item, width in zip(transitions, widths):
                column = _number(float(np.dot(means[item.upper_state_id], weights)), "upper-state column")
                radiance = _number(float((_H * _C / (4 * np.pi * item.vacuum_wavelength_nm * 1e-9))
                                         * item.einstein_a_per_s * column), "integrated line radiance")
                sigma = width / np.sqrt(8 * np.log(2))
                z_left = (edges[:-1] - item.vacuum_wavelength_nm) / sigma
                z_right = (edges[1:] - item.vacuum_wavelength_nm) / sigma
                area = np.where(z_left > 0, ndtr(-z_left) - ndtr(-z_right),
                                ndtr(z_right) - ndtr(z_left))
                bin_radiance += radiance * area
                lines.append({"transition_id": item.transition_id, "upper_state_id": item.upper_state_id,
                              "vacuum_wavelength_nm": float(item.vacuum_wavelength_nm),
                              "einstein_a_per_s": float(item.einstein_a_per_s),
                              "fwhm_nm": float(width), "gate_mean_upper_column_m2": column,
                              "integrated_radiance_w_m2_sr": radiance})
            bin_mean = bin_radiance / np.diff(edges)
    except FloatingPointError as exc:
        raise ValueError("resolved emission projection exceeds the finite numerical domain") from exc
    if not np.isfinite(bin_mean).all() or not np.isfinite(bin_radiance).all():
        raise ValueError("resolved emission bins exceed the finite numerical domain")
    return {
        "schema": "mhdlab-resolved-emission-observation-v1", "status": status,
        "optical_model": optical_model, "time_registration": time_registration,
        "population_provenance": population_provenance, "geometry_provenance": geometry_provenance,
        "gate": {"start_s": start, "end_s": stop, "averaging": "rectangular_piecewise_linear_population_history"},
        "time_s": times.tolist(), "path_weights_m": weights.tolist(),
        "path_weight_convention": "ray-cell lengths multiplied by dimensionless normalized radiance-collection weights once; not detector etendue or fitted scalar depth; geometry fixed during gate and normalization owned by caller",
        "gate_mean_upper_state_density_m3": {key: value.tolist() for key, value in means.items()},
        "lines": lines, "bin_edges_nm": edges.tolist(),
        "bin_integrated_radiance_w_m2_sr": bin_radiance.tolist(),
        "bin_mean_radiance_w_m2_sr_nm": bin_mean.tolist(),
        "profile": "area_normalized_gaussian; finite-window losses are not renormalized",
        "population_closure_supplied": False, "density_temperature_fit_performed": False,
        "flash_simulation_performed": False,
    }
