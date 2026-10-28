"""Explicit binned radiance/count conversion with measurement covariance.

The supplied photon response excludes exposure, etendue and wavelength-bin
integration. No throughput, molecular population, spectral transport or actual
instrument calibration is inferred. Values are bin means in the declared units.
"""

from __future__ import annotations

from dataclasses import KW_ONLY, dataclass
from fractions import Fraction
import math
from numbers import Real
from types import MappingProxyType
from collections.abc import Mapping

import numpy as np


PLANCK_CONSTANT_J_S = 6.62607015e-34
SPEED_OF_LIGHT_M_S = 299792458.0
SCHEMA = "mhdlab-radiometric-response-v1"
RADIANCE_UNITS = MappingProxyType({
    "photon_spectral_radiance": "photons s^-1 m^-2 sr^-1 nm^-1",
    "energy_spectral_radiance": "W m^-2 sr^-1 nm^-1",
})
COVARIANCE_PSD_RELATIVE_TOLERANCE = 1e-12
_H_EXACT = Fraction("6.62607015e-34")
_C_EXACT = Fraction(299792458)
_NM_TO_M_EXACT = Fraction(1, 10**9)
_INTEGRATION_IDENTITIES = frozenset({
    "exposure_time", "effective_exposure_time", "effective_exposure_s", "gate_duration_s",
    "etendue", "etendue_m2_sr", "bin_width_nm", "wavelength_bin_width",
})


def _real(value, name, *, positive=False):
    if (np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_))
            or not isinstance(value, Real)):
        raise ValueError(f"{name} must be a finite real number, without bool, mask, string or complex values")
    try:
        result = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} exceeds finite real numeric support") from error
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and strictly positive" if positive else ""))
    return result


def _array(value, name, ndim):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        raw = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real array") from error
    if raw.ndim != ndim or raw.size == 0:
        raise ValueError(f"{name} must be a nonempty {ndim}-dimensional array")
    return np.array([_real(item, name) for item in raw.flat], dtype=float).reshape(raw.shape)


def _immutable(value):
    array = np.ascontiguousarray(value, dtype=float)
    return np.frombuffer(array.tobytes(), dtype=float).reshape(array.shape)


def _label(value, name):
    if (not isinstance(value, str) or not value.strip() or value != value.strip()
            or any(ord(character) < 32 for character in value)):
        raise ValueError(f"{name} must be a nonempty label without surrounding whitespace or control characters")
    return value


def _kind(value):
    if not isinstance(value, str) or value not in RADIANCE_UNITS:
        raise ValueError("radiance_kind must explicitly be photon_spectral_radiance or energy_spectral_radiance; ambiguous ordinate units are unsupported")
    return value


def _rounded(exact, name):
    try:
        result = float(exact)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} overflows finite numeric support") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} overflows finite numeric support")
    if exact and result == 0:
        raise ValueError(f"positive or signed nonzero {name} underflows numeric support")
    return result


def _freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(item) for item in value)
    return value


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    return value


def _covariance(value, n, name):
    covariance = _array(value, name, 2)
    if covariance.shape != (n, n):
        raise ValueError(f"{name} must have shape ({n},{n}); a variance vector is not a full covariance")
    if not np.array_equal(covariance, covariance.T):
        raise ValueError(f"{name} must be exactly symmetric; no hidden symmetrization is applied")
    diagonal = np.diag(covariance)
    if np.any(diagonal < 0):
        raise ValueError(f"{name} cannot contain a negative variance")
    positive = diagonal > 0
    if np.any(covariance[~positive] != 0):
        raise ValueError(f"{name} zero-variance rows and columns must have zero covariance")
    minimum, maximum = 0., 0.
    if np.any(positive):
        selected = covariance[np.ix_(positive, positive)]
        deviation = np.sqrt(diagonal[positive])
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                correlation = selected / deviation[:, None] / deviation[None, :]
            # A diagonal is exactly one in correlation coordinates. This changes
            # only the validation matrix, never the supplied covariance.
            np.fill_diagonal(correlation, 1.)
            eigenvalues = np.linalg.eigvalsh(correlation)
        except (FloatingPointError, np.linalg.LinAlgError) as error:
            raise ValueError(f"{name} has unsupported or non-PSD correlation scaling") from error
        if not np.all(np.isfinite(eigenvalues)):
            raise ValueError(f"{name} has unsupported finite-range PSD qualification")
        minimum, maximum = float(eigenvalues[0]), float(eigenvalues[-1])
        if minimum < -COVARIANCE_PSD_RELATIVE_TOLERANCE * max(1., maximum):
            raise ValueError(f"{name} is not positive semidefinite within the declared correlation-eigenvalue tolerance")
    qualification = {
        "symmetry": "exact input symmetry required; no symmetrization",
        "psd_coordinates": "positive-variance correlation submatrix; exact zero-variance rows checked separately",
        "psd_relative_eigenvalue_tolerance": COVARIANCE_PSD_RELATIVE_TOLERANCE,
        "minimum_correlation_eigenvalue": minimum,
        "maximum_correlation_eigenvalue": maximum,
        "singular_covariance_supported": True,
        "jitter_or_projection_applied": False,
        "qualification": "numerical PSD qualification, not a positive-definite likelihood guarantee",
    }
    return covariance, qualification


@dataclass(frozen=True, slots=True, eq=False, init=False)
class CalibratedRadiance:
    """Validated output of ``RadiometricResponse.calibrate_counts``.

    Direct construction is deliberately unavailable: conversion and covariance
    qualification must run before this result record is created.
    """

    bin_edges_nm: np.ndarray
    radiance: np.ndarray
    covariance: np.ndarray
    radiance_kind: str
    radiance_units: str
    _metadata: Mapping

    def __init__(self, *args, **kwargs):
        raise TypeError("use RadiometricResponse.calibrate_counts to construct validated radiance results")

    @classmethod
    def _from_validated(cls, bin_edges_nm, radiance, covariance, radiance_kind, radiance_units, metadata):
        result = object.__new__(cls)
        for name, value in (("bin_edges_nm", bin_edges_nm), ("radiance", radiance), ("covariance", covariance)):
            object.__setattr__(result, name, _immutable(value))
        object.__setattr__(result, "radiance_kind", radiance_kind)
        object.__setattr__(result, "radiance_units", radiance_units)
        object.__setattr__(result, "_metadata", _freeze(metadata))
        return result

    @property
    def metadata(self):
        """Detached JSON-compatible conversion and qualification metadata."""
        return _plain(self._metadata)


@dataclass(frozen=True, slots=True, eq=False)
class RadiometricResponse:
    """Supplied counts per collected photon, with explicit geometry and time.

    The response and selected radiance are constant within each bin. Effective
    exposure is one separately supplied integration time, including the caller's
    temporal acceptance convention. Gain may exceed one count per photon.
    """

    bin_edges_nm: np.ndarray
    etendue_m2_sr: float
    effective_exposure_s: float
    counts_per_collected_photon: np.ndarray
    _: KW_ONLY
    reference_plane: str
    included_factors: tuple[str, ...]
    status: str
    provenance: str
    measurement_conditions_id: str
    response_includes_exposure: bool
    response_includes_etendue: bool
    response_includes_bin_width: bool

    def __post_init__(self):
        edges = _array(self.bin_edges_nm, "bin_edges_nm", 1)
        if len(edges) < 2 or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
            raise ValueError("bin_edges_nm must be positive, strictly increasing and define at least one bin")
        response = _array(self.counts_per_collected_photon, "counts_per_collected_photon", 1)
        if response.shape != (len(edges) - 1,) or np.any(response <= 0):
            raise ValueError("counts_per_collected_photon must have one strictly positive value per bin")
        object.__setattr__(self, "bin_edges_nm", _immutable(edges))
        object.__setattr__(self, "counts_per_collected_photon", _immutable(response))
        for name in ("etendue_m2_sr", "effective_exposure_s"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        for name in ("reference_plane", "status", "provenance", "measurement_conditions_id"):
            object.__setattr__(self, name, _label(getattr(self, name), name))
        for name in ("response_includes_exposure", "response_includes_etendue", "response_includes_bin_width"):
            value = getattr(self, name)
            if not isinstance(value, (bool, np.bool_)) or bool(value):
                raise ValueError(f"{name} must explicitly be False; integration is performed separately")
            object.__setattr__(self, name, False)
        if not isinstance(self.included_factors, tuple) or not self.included_factors:
            raise ValueError("included_factors must explicitly be a nonempty tuple of supplied response factor labels")
        factors = tuple(_label(item, "included_factors") for item in self.included_factors)
        if len(factors) != len(set(factors)):
            raise ValueError("included_factors cannot repeat the same factor")
        if any(item.casefold() in _INTEGRATION_IDENTITIES for item in factors):
            raise ValueError("included_factors declares geometry/time/bin integration already included in photon response")
        object.__setattr__(self, "included_factors", factors)
        # Validate numerical support for both declared kinds at construction.
        for kind in RADIANCE_UNITS:
            self.conversion_factors(kind)

    @property
    def metadata(self):
        return {
            "schema": SCHEMA,
            "bin_edges_nm": self.bin_edges_nm.tolist(),
            "etendue_m2_sr": self.etendue_m2_sr,
            "effective_exposure_s": self.effective_exposure_s,
            "counts_per_collected_photon": self.counts_per_collected_photon.tolist(),
            "reference_plane": self.reference_plane,
            "included_factors": list(self.included_factors),
            "status": self.status,
            "provenance": self.provenance,
            "measurement_conditions_id": self.measurement_conditions_id,
            "response_includes_exposure": False,
            "response_includes_etendue": False,
            "response_includes_bin_width": False,
            "assumptions": [
                "selected photon or energy spectral radiance is a constant bin mean, per nm",
                "counts per collected photon is constant within each bin at the declared reference plane",
                "etendue and effective exposure are constant across these bins and supplied separately",
                "energy-radiance conversion integrates wavelength/(h*c) across each bin exactly under these assumptions",
                "signed values describe a linear background-subtracted signal, not a nonnegative count likelihood",
                "condition-ID equality is a caller-supplied association, not independent recovery of instrument conditions",
            ],
            "radiance_units": dict(RADIANCE_UNITS),
            "planck_constant_j_s": PLANCK_CONSTANT_J_S,
            "speed_of_light_m_s": SPEED_OF_LIGHT_M_S,
            "constants_reference": "exact SI defining values; factors rounded to float64 once",
            "automatic_throughput_or_gate_multiplication": False,
            "system_calibration_uncertainty_included": False,
        }

    def conversion_factors(self, radiance_kind):
        """Counts per unit of the selected bin-mean spectral radiance."""
        kind = _kind(radiance_kind)
        common = Fraction(self.etendue_m2_sr) * Fraction(self.effective_exposure_s)
        values = []
        for lower, upper, response in zip(self.bin_edges_nm[:-1], self.bin_edges_nm[1:],
                                           self.counts_per_collected_photon):
            lo, hi = Fraction(float(lower)), Fraction(float(upper))
            factor = (hi - lo) * common * Fraction(float(response))
            if kind == "energy_spectral_radiance":
                factor *= (lo + hi) / 2 * _NM_TO_M_EXACT / (_H_EXACT * _C_EXACT)
            values.append(_rounded(factor, "radiometric conversion factor"))
        return _immutable(values)

    def _conditions(self, measurement_conditions_id):
        if _label(measurement_conditions_id, "measurement_conditions_id") != self.measurement_conditions_id:
            raise ValueError("measurement conditions do not match the supplied response")

    def _bins(self, value, name):
        result = _array(value, name, 1)
        if result.shape != self.counts_per_collected_photon.shape:
            raise ValueError(f"{name} must have one value per response bin; no wavelength extrapolation or interpolation")
        return result

    def counts_from_radiance(self, radiance, *, radiance_kind, measurement_conditions_id):
        """Convert a signed bin-mean radiance signal to integrated ADU counts."""
        self._conditions(measurement_conditions_id)
        factors = self.conversion_factors(radiance_kind)
        values = self._bins(radiance, "radiance")
        return _immutable([_rounded(Fraction(float(value)) * Fraction(float(factor)), "counts")
                           for value, factor in zip(values, factors)])

    def calibrate_counts(self, counts, covariance_adu2, *, radiance_kind, measurement_conditions_id):
        """Apply D=diag(1/factor) to counts and D C D to full covariance.

        The input covariance is measurement covariance in ADU squared. Singular
        PSD matrices are supported under the reported numerical qualification;
        no calibration nuisance uncertainty, Poisson term, jitter or projection
        is added. A downstream likelihood may require a stricter SPD contract.
        """
        self._conditions(measurement_conditions_id)
        kind = _kind(radiance_kind)
        factors = self.conversion_factors(kind)
        values = self._bins(counts, "counts")
        covariance, qualification = _covariance(covariance_adu2, len(values), "covariance_adu2")
        exact_factors = [Fraction(float(factor)) for factor in factors]
        radiance = [_rounded(Fraction(float(value)) / factor, "calibrated radiance")
                    for value, factor in zip(values, exact_factors)]
        transformed = np.empty(covariance.shape)
        for i in range(len(values)):
            for j in range(i, len(values)):
                value = _rounded(Fraction(float(covariance[i, j])) / (exact_factors[i] * exact_factors[j]),
                                 "calibrated covariance component")
                transformed[i, j] = transformed[j, i] = value
        _, output_qualification = _covariance(transformed, len(values), "calibrated covariance")
        metadata = self.metadata
        metadata.update({
            "radiance_kind": kind,
            "radiance_units": RADIANCE_UNITS[kind],
            "covariance_units": f"({RADIANCE_UNITS[kind]})^2",
            "counts_units": "ADU",
            "input_covariance_units": "ADU^2",
            "counts_per_radiance_unit": factors.tolist(),
            "input_covariance_qualification": qualification,
            "output_covariance_qualification": output_qualification,
            "uncertainty_scope": "supplied measurement covariance only; calibration systematics require a separate nuisance model",
        })
        return CalibratedRadiance._from_validated(self.bin_edges_nm, radiance, transformed, kind,
                                                  RADIANCE_UNITS[kind], metadata)
