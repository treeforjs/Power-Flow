"""Coverage regression against externally supplied, already binned predictions.

This module synthesizes no spectrum and evaluates no plasma or chemistry model.
One monolayer is an explicit areal-density definition supplied with provenance;
it is never inferred from a desorption reference coverage. A declared known
background is added without calibration scaling. The covariance is used in
full, without diagonal replacement, regularization, or estimated noise floors.

Coverage is interpolated linearly only inside the supplied finite grid. Every
segment's endpoints, analytic profile stationary point and active scale-zero
transition are considered. Calculations use float64 Cholesky whitening and
long-double scalar reductions; unresolved or nonfinite arithmetic is rejected.
The fitted grid and numerical diagnostics establish no confidence interval or
physical validation of the supplied predictions, calibration or layer definition.
"""

from __future__ import annotations

from dataclasses import KW_ONLY, dataclass, field
import math
from numbers import Real

import numpy as np
from scipy.linalg import solve_triangular


SCHEMA = "mhdlab-coverage-regression-v1"
STATUSES = frozenset({"manufactured", "observed", "reconstructed", "sourced", "source_evaluated", "owner_supplied"})
LD = np.longdouble
EPS = float(np.finfo(float).eps)
ILL_CONDITIONED_THRESHOLD = 1.e12
SEARCH_METHOD = "all_piecewise_linear_segments_analytic_profile_stationary_points_endpoints_and_scale_zero_transitions_v1"


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _status(value):
    if not isinstance(value, str) or value not in STATUSES:
        raise ValueError(f"status must be one of {sorted(STATUSES)}")
    return value


def _real(value, name, *, nonnegative=False, positive=False):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar, not a boolean or masked value")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise ValueError(f"{name} must be representable as a finite real scalar") from error
    if not math.isfinite(result) or (nonnegative and result < 0) or (positive and result <= 0):
        raise ValueError(f"{name} is outside its finite {'positive' if positive else 'nonnegative' if nonnegative else 'real'} domain")
    return result


def _array(value, name, ndim, *, nonnegative=False):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        objects = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real array") from error
    if objects.ndim != ndim or objects.size == 0:
        raise ValueError(f"{name} must be a nonempty {ndim}-dimensional array")
    for item in objects.flat:
        _real(item, name, nonnegative=nonnegative)
    result = np.asarray(objects, dtype=np.float64)
    return np.frombuffer(np.ascontiguousarray(result).tobytes(), dtype=np.float64).reshape(result.shape)


def _edges(value, bins):
    edges = _array(value, "bin_edges", 1)
    if len(edges) != bins + 1 or np.any(edges[1:] <= edges[:-1]):
        raise ValueError("bin_edges must have nbin+1 finite strictly increasing edges")
    return edges


def _float(value, name):
    if not np.isfinite(value) or abs(value) > np.finfo(float).max:
        raise ValueError(f"{name} exceeds supported finite float64 output range")
    result = float(value)
    if value != 0 and result == 0:
        raise ValueError(f"{name} is positive in magnitude but below supported float64 output range")
    return result


def _immutable(value):
    value = np.asarray(value, dtype=np.float64)
    return np.frombuffer(np.ascontiguousarray(value).tobytes(), dtype=np.float64).reshape(value.shape)


@dataclass(frozen=True, slots=True, eq=False)
class CoveragePredictionGrid:
    coverage_ml: np.ndarray
    predictions: np.ndarray
    _: KW_ONLY
    monolayer_areal_density_m2: float
    monolayer_definition: str
    monolayer_provenance: str
    signal_units: str
    bin_edges: np.ndarray
    bin_axis_units: str
    bin_support: str
    status: str
    provenance: str

    def __post_init__(self):
        coverage = _array(self.coverage_ml, "coverage_ml", 1, nonnegative=True)
        predictions = _array(self.predictions, "predictions", 2, nonnegative=True)
        if len(coverage) < 2 or np.any(coverage[1:] <= coverage[:-1]):
            raise ValueError("coverage_ml must contain at least two finite strictly increasing nonnegative knots")
        if predictions.shape[0] != len(coverage):
            raise ValueError("predictions must have shape (ncoverage, nbin)")
        density = _real(self.monolayer_areal_density_m2, "monolayer_areal_density_m2", positive=True)
        _float(LD(coverage[-1])*LD(density), "maximum areal inventory")
        object.__setattr__(self, "coverage_ml", coverage)
        object.__setattr__(self, "predictions", predictions)
        object.__setattr__(self, "monolayer_areal_density_m2", density)
        object.__setattr__(self, "bin_edges", _edges(self.bin_edges, predictions.shape[1]))
        for name in ("monolayer_definition", "monolayer_provenance", "signal_units", "bin_axis_units", "bin_support", "provenance"):
            _text(getattr(self, name), name)
        _status(self.status)


@dataclass(frozen=True, slots=True, eq=False)
class BinnedObservation:
    values: np.ndarray
    covariance: np.ndarray
    known_background: np.ndarray
    _: KW_ONLY
    signal_units: str
    bin_edges: np.ndarray
    bin_axis_units: str
    bin_support: str
    status: str
    provenance: str
    _cholesky: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        # Signed calibrated/dark-subtracted bins are retained rather than censored.
        values = _array(self.values, "values", 1)
        background = _array(self.known_background, "known_background", 1)
        covariance = _array(self.covariance, "covariance", 2)
        if background.shape != values.shape or covariance.shape != (len(values), len(values)):
            raise ValueError("background and covariance must match the observation bins")
        if not np.array_equal(covariance, covariance.T):
            raise ValueError("covariance must be exactly symmetric; no implicit symmetrization is applied")
        try:
            factor = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as error:
            raise ValueError("covariance must be symmetric positive definite, without regularization") from error
        if not np.all(np.isfinite(factor)) or np.any(np.diag(factor) <= 0):
            raise ValueError("covariance Cholesky factor is not finite positive definite")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "known_background", background)
        object.__setattr__(self, "covariance", covariance)
        object.__setattr__(self, "_cholesky", _immutable(factor))
        object.__setattr__(self, "bin_edges", _edges(self.bin_edges, len(values)))
        for name in ("signal_units", "bin_axis_units", "bin_support", "provenance"):
            _text(getattr(self, name), name)
        _status(self.status)


@dataclass(frozen=True, slots=True, eq=False)
class CalibrationScale:
    mode: str
    _: KW_ONLY
    fixed_scale: float | None = None
    prior_mean: float | None = None
    prior_sigma: float | None = None
    provenance: str

    def __post_init__(self):
        if not isinstance(self.mode, str) or self.mode not in {"fixed", "gaussian_prior", "free"}:
            raise ValueError("calibration mode must be fixed, gaussian_prior or free")
        _text(self.provenance, "calibration provenance")
        if self.mode == "fixed":
            if self.fixed_scale is None or self.prior_mean is not None or self.prior_sigma is not None:
                raise ValueError("fixed calibration requires only an explicit fixed_scale")
            object.__setattr__(self, "fixed_scale", _real(self.fixed_scale, "fixed_scale", nonnegative=True))
        elif self.mode == "gaussian_prior":
            if self.fixed_scale is not None or self.prior_mean is None or self.prior_sigma is None:
                raise ValueError("gaussian_prior calibration requires only explicit prior_mean and prior_sigma")
            object.__setattr__(self, "prior_mean", _real(self.prior_mean, "prior_mean", nonnegative=True))
            object.__setattr__(self, "prior_sigma", _real(self.prior_sigma, "prior_sigma", positive=True))
        elif any(value is not None for value in (self.fixed_scale, self.prior_mean, self.prior_sigma)):
            raise ValueError("free calibration cannot supply fixed or prior parameters")


def _whiten(factor, values):
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        result = solve_triangular(factor, values, lower=True, check_finite=False)
    if not np.all(np.isfinite(result)):
        raise ValueError("Cholesky whitening exceeds the supported finite float64 range")
    result_ld, factor_ld, values_ld = (np.asarray(value, dtype=LD) for value in (result, factor, values))
    reconstructed = factor_ld @ result_ld
    backward_scale = np.abs(factor_ld) @ np.abs(result_ld)+np.abs(values_ld)
    if np.any(np.abs(reconstructed-values_ld) > LD(256*EPS)*backward_scale):
        raise ValueError("Cholesky whitening loses a resolved input component, including possible underflow")
    return result_ld


def _stationary_fraction(a, b, y, precision, prior_mean):
    """Rotate the segment direction, retaining weak perpendicular information.

    Subtracting products of Gram entries loses a small but measured shape
    direction beside a large amplitude. A Householder reflection instead gives
    model=(A+t*B, w), target=(Y, v). The profile stationary condition is
    Y*||w||^2-(w.v)*(A+t*B)=0, with no difference of large Gram products.
    """
    if precision > 0:
        root_precision = np.sqrt(precision)
        a = np.concatenate((a, np.array([root_precision], dtype=LD)))
        b = np.concatenate((b, np.array([0], dtype=LD)))
        y = np.concatenate((y, np.array([root_precision*prior_mean], dtype=LD)))
    norm = np.sqrt(np.dot(b, b))
    if norm == 0:
        return None
    pivot = int(np.argmax(np.abs(b)))
    signed_norm = norm if b[pivot] >= 0 else -norm
    reflector = b.copy()
    reflector[pivot] += signed_norm
    reflector_norm = np.dot(reflector, reflector)
    rotated_a = a-reflector*(2*np.dot(reflector, a)/reflector_norm)
    rotated_y = y-reflector*(2*np.dot(reflector, y)/reflector_norm)
    perpendicular = np.arange(len(a)) != pivot
    w, v = rotated_a[perpendicular], rotated_y[perpendicular]
    projection = np.dot(w, v)
    if projection == 0:
        return None
    parallel_derivative = -signed_norm
    return (rotated_y[pivot]*np.dot(w, w)-projection*rotated_a[pivot])/(projection*parallel_derivative)


def fit_surface_coverage(grid: CoveragePredictionGrid, observation: BinnedObservation,
                         calibration: CalibrationScale) -> dict:
    """Return a finite-grid conditional fit and explicitly qualified diagnostics.

    The Gaussian prior is a quadratic penalty on linear calibration scale,
    constrained to scale>=0. It is not a log-scale or truncated-normal density
    normalization. The reported chi-square omits parameter-independent constants.
    """
    if not isinstance(grid, CoveragePredictionGrid) or not isinstance(observation, BinnedObservation) or not isinstance(calibration, CalibrationScale):
        raise ValueError("fit requires validated CoveragePredictionGrid, BinnedObservation and CalibrationScale records")
    if (grid.signal_units != observation.signal_units or grid.bin_axis_units != observation.bin_axis_units
            or grid.bin_support != observation.bin_support or not np.array_equal(grid.bin_edges, observation.bin_edges)):
        raise ValueError("prediction and observation signal units, bin edges, axis units and declared support must match exactly")
    if grid.predictions.shape[1] != len(observation.values):
        raise ValueError("prediction and observation bin counts must match")
    difference = np.asarray(observation.values, dtype=LD)-np.asarray(observation.known_background, dtype=LD)
    if np.any(np.abs(difference) > np.finfo(float).max):
        raise ValueError("observation minus known background exceeds finite whitening range")
    y = _whiten(observation._cholesky, np.asarray(difference, dtype=float))
    models = _whiten(observation._cholesky, grid.predictions.T).T
    knots = np.asarray(grid.coverage_ml, dtype=LD)
    precision = LD(0) if calibration.mode != "gaussian_prior" else LD(1)/LD(calibration.prior_sigma)**2
    prior_mean = LD(0) if calibration.mode != "gaussian_prior" else LD(calibration.prior_mean)

    def evaluate(interval, fraction, reason):
        # Retain the same representable fraction in the candidate record and
        # final replay, so its reported objective is exactly reproducible.
        t = LD(_float(LD(fraction), "interval fraction"))
        if not np.isfinite(t) or t < 0 or t > 1:
            raise ValueError("candidate lies outside a supplied prediction segment")
        model = (1-t)*models[interval]+t*models[interval+1]
        physical_model = (1-t)*np.asarray(grid.predictions[interval], dtype=LD)+t*np.asarray(grid.predictions[interval+1], dtype=LD)
        norm = np.dot(model, model)
        numerator = np.dot(model, y)+precision*prior_mean
        denominator = norm+precision
        unidentified = calibration.mode == "free" and denominator == 0
        if calibration.mode == "fixed":
            scale = LD(calibration.fixed_scale)
        elif denominator == 0:
            scale = LD(0)
        else:
            scale = max(LD(0), numerator/denominator)
        residual = y-scale*model
        data_cost = np.dot(residual, residual)
        prior_cost = precision*(scale-prior_mean)**2
        coverage = (1-t)*knots[interval]+t*knots[interval+1]
        represented_coverage = _float(coverage, "coverage")
        if 0 < t < 1 and represented_coverage in (grid.coverage_ml[interval], grid.coverage_ml[interval+1]):
            raise ValueError("interior coverage rounds onto a grid endpoint; the coverage coordinate is unresolved")
        result = {
            "coverage_ml": represented_coverage,
            "areal_inventory_m2": _float(coverage*LD(grid.monolayer_areal_density_m2), "areal inventory"),
            "calibration_scale": _float(scale, "calibration scale"),
            "chi2": _float(data_cost+prior_cost, "chi-square"),
            "data_chi2": _float(data_cost, "data chi-square"),
            "prior_chi2": _float(prior_cost, "prior chi-square"),
            "zero_scale_boundary": bool(scale == 0),
            "zero_model": bool(norm == 0),
            "scale_unidentified": bool(unidentified),
            "interval_index": int(interval), "interval_fraction": _float(t, "interval fraction"),
            "candidate_reason": reason,
            "_residual_ld": residual, "_scale_ld": scale,
        }
        return result, physical_model, residual

    candidates = []
    knot_records = []
    for index in range(len(knots)):
        interval, fraction = (index, LD(0)) if index < len(knots)-1 else (index-1, LD(1))
        record, _, _ = evaluate(interval, fraction, "supplied_grid_knot")
        candidates.append(record)
        knot_records.append(record)
    for interval in range(len(knots)-1):
        a, b = models[interval], models[interval+1]-models[interval]
        fractions = []
        if calibration.mode == "fixed":
            scale = LD(calibration.fixed_scale)
            denominator = scale*np.dot(b, b)
            if denominator != 0:
                fractions.append((np.dot(b, y-scale*a)/denominator, "fixed_scale_quadratic_stationary"))
        else:
            b0, b1 = np.dot(a, y)+precision*prior_mean, np.dot(b, y)
            stationary = _stationary_fraction(a, b, y, precision, prior_mean)
            if stationary is not None:
                fractions.append((stationary, "profile_rational_quadratic_stationary"))
            if b1 != 0:
                fractions.append((-b0/b1, "scale_zero_transition"))
        for fraction, reason in fractions:
            if np.isfinite(fraction) and 0 < fraction < 1:
                record, _, _ = evaluate(interval, fraction, reason)
                candidates.append(record)
    # A deterministic lowest-coverage tie choice is not a uniqueness assertion.
    candidates.sort(key=lambda item: (item["coverage_ml"], item["candidate_reason"]))
    anchor_residual, anchor_scale = candidates[0]["_residual_ld"], candidates[0]["_scale_ld"]
    for record in candidates:
        residual, scale = record["_residual_ld"], record["_scale_ld"]
        record["_contrast_ld"] = np.dot(residual-anchor_residual, residual+anchor_residual)+precision*(scale-anchor_scale)*(scale+anchor_scale-2*prior_mean)
    best_record = min(candidates, key=lambda item: (item["_contrast_ld"], item["coverage_ml"]))
    best_contrast = best_record["_contrast_ld"]
    best, physical, residual = evaluate(best_record["interval_index"], best_record["interval_fraction"], best_record["candidate_reason"])
    best["coverage_boundary"] = ("lower" if best["coverage_ml"] == grid.coverage_ml[0] else "upper" if best["coverage_ml"] == grid.coverage_ml[-1] else "interior")
    # Use the same profiled scale as the objective; rounding that scalar before
    # multiplying can introduce a resolved error for very precise observations.
    signal = best["_scale_ld"]*physical
    prediction_with_background = signal+np.asarray(observation.known_background, dtype=LD)
    best["prediction_signal"] = [_float(value, "prediction signal") for value in signal]
    best["prediction_with_background"] = [_float(value, "prediction plus background") for value in prediction_with_background]
    best["whitened_residual"] = [_float(value, "whitened residual") for value in residual]
    costs = np.asarray([record["_contrast_ld"] for record in candidates], dtype=LD)
    cost_span = np.max(costs)-np.min(costs)
    flat_tolerance = LD(256*EPS)*max(LD(1), np.max(np.abs(costs)))
    flat = bool(cost_span <= flat_tolerance)
    norms = np.sqrt(np.sum(models*models, axis=1, dtype=LD))
    positive = norms > 0
    shape_difference = LD(0)
    shape_rank = 0
    if np.any(positive):
        directions = models[positive]/norms[positive, None]
        reference = directions[int(np.argmax(norms[positive]))]
        shape_difference = np.max(np.sqrt(np.sum((directions-reference)**2, axis=1, dtype=LD)))
        singular = np.linalg.svd(np.asarray(directions, dtype=float), compute_uv=False)
        shape_rank = int(np.count_nonzero(singular > max(directions.shape)*EPS*singular[0]))
    proportional = bool(shape_difference <= LD(256*EPS))
    degeneracy = bool(calibration.mode == "free" and proportional)
    try:
        condition = float(np.linalg.cond(observation.covariance))
    except np.linalg.LinAlgError:
        condition = math.inf
    condition_output = condition if math.isfinite(condition) else None
    grid_profile = {name: [record[name] for record in knot_records] for name in (
        "coverage_ml", "areal_inventory_m2", "calibration_scale", "chi2", "data_chi2", "prior_chi2", "zero_scale_boundary")}
    grid_profile["delta_chi2"] = [_float(max(LD(0), record["_contrast_ld"]-best_contrast), "profile delta chi-square") for record in knot_records]
    rounded_ties = any(record["chi2"] == best["chi2"] and record["_contrast_ld"] != best_contrast for record in candidates)
    for record in candidates:
        record["delta_chi2"] = _float(max(LD(0), record["_contrast_ld"]-best_contrast), "candidate delta chi-square")
        for key in ("_residual_ld", "_scale_ld", "_contrast_ld"):
            del record[key]
    del best["_residual_ld"], best["_scale_ld"]
    return {
        "schema": SCHEMA,
        "best_fit": best, "grid_profile": grid_profile, "search_candidates": candidates,
        "monolayer": {"areal_density_m2": grid.monolayer_areal_density_m2, "definition": grid.monolayer_definition, "provenance": grid.monolayer_provenance},
        "calibration": {"mode": calibration.mode, "fixed_scale": calibration.fixed_scale, "prior_mean": calibration.prior_mean, "prior_sigma": calibration.prior_sigma, "provenance": calibration.provenance},
        "bin_support": {"signal_units": grid.signal_units, "bin_edges": grid.bin_edges.tolist(), "bin_axis_units": grid.bin_axis_units, "declaration": grid.bin_support, "known_background": observation.known_background.tolist()},
        "inputs": {"prediction_status": grid.status, "prediction_provenance": grid.provenance, "observation_status": observation.status, "observation_provenance": observation.provenance},
        "diagnostics": {
            "flat_profile": flat, "flat_profile_tolerance": _float(flat_tolerance, "profile flatness tolerance"), "profile_span_chi2": _float(cost_span, "profile span"),
            "coverage_amplitude_degeneracy": degeneracy,
            "proportional_prediction_shapes": proportional, "maximum_normalized_shape_difference": _float(shape_difference, "prediction shape difference"),
            "prediction_shape_rank": shape_rank, "zero_prediction_knots": np.flatnonzero(~positive).tolist(),
            "zero_model_at_best": best["zero_model"], "scale_unidentified_at_best": best["scale_unidentified"],
            "coverage_unconstrained": bool(flat or degeneracy or best["zero_scale_boundary"]),
            "calibration_mode": calibration.mode, "covariance_condition_number": condition_output,
            "ill_conditioned_covariance": bool(condition >= ILL_CONDITIONED_THRESHOLD), "ill_conditioned_covariance_threshold": ILL_CONDITIONED_THRESHOLD,
            "search_method": SEARCH_METHOD, "segments_searched": len(knots)-1, "candidates_evaluated": len(candidates),
            "reported_chi2_rounding_ties": bool(rounded_ties),
            "objective_comparison": "long-double centered residual and prior contrasts; reported total chi-square may round equal while delta chi-square remains resolved",
            "parameter_output_rounding": "scalar parameters and arrays are float64 summaries of the internal profiled fit; use supplied prediction arrays rather than multiplying a rounded calibration scalar when its rounding is resolved by the covariance",
            "coverage_domain_ml": [float(knots[0]), float(knots[-1])], "extrapolation_performed": False,
        },
        "limitations": [
            "Externally supplied prediction grid only; no spectral, chemistry, plasma or gap solver evaluated.",
            "Piecewise-linear coverage interpolation on the supplied finite grid; no extrapolation or prediction-grid refinement bound.",
            "Full supplied symmetric positive-definite observation covariance is fixed; no covariance, background or noise model fitted.",
            "Calibration scale multiplies predictions only; explicit known background remains unscaled.",
            "Free scale and proportional prediction shapes confound coverage with amplitude; a selected numerical minimum is not an absolute monolayer measurement.",
            "Flatness and rank are floating-point diagnostics, not physical identifiability proofs; zero-scale/zero-model boundaries may leave coverage unconstrained.",
            "No automatic confidence interval, likelihood threshold, degrees-of-freedom claim or physical model validation.",
        ],
    }
