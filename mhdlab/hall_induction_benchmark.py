"""Periodic Hall-induction verification with prescribed stationary ions.

This is an exact constant-coefficient induction-subsystem reference, not a
plasma closure, material model or full Hall-MHD solution. Scalar resistivity,
permeability and electron density are supplied, with no electron pressure or
inertia. The periodic helical perturbation rotates without Hall dissipation.
Scalar resistivity damps it and deposits exactly its lost magnetic energy.

Cell-averaged fields are analytic; Joule power is averaged independently of J.
The uniform axial background energy stays separate from perturbation energy.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real

import numpy as np

from .constants import QE
from .thermal_release import _deplete


def _real(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number without bool/string coercion")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result) or (positive and result <= 0) or (nonnegative and result < 0):
        raise ValueError(f"{name} is outside its finite physical domain")
    return result


def _scaled(numerators, denominators=(), *, attenuation=0.0):
    """Combine finite factors before restoring their common binary scale."""
    if any(value == 0 for value in numerators):
        return 0.0
    mantissa, exponent = 1.0, 0
    for value in numerators:
        fraction, power = math.frexp(value)
        mantissa *= fraction
        exponent += power
    for value in denominators:
        fraction, power = math.frexp(value)
        mantissa /= fraction
        exponent -= power
    if attenuation == math.inf:
        return 0.0
    if attenuation > 700:
        # Combine the decay with the full binary scale before exponentiation;
        # exp(-attenuation) alone may underflow while this result survives.
        logarithm = math.fsum((math.log(abs(mantissa)), exponent * math.log(2.0), -attenuation))
        try:
            result = math.copysign(math.exp(logarithm), mantissa)
        except OverflowError as exc:
            raise ValueError("Hall reference result exceeds finite float64 range") from exc
        return result
    if attenuation:
        fraction, power = math.frexp(math.exp(-attenuation))
        mantissa *= fraction
        exponent += power
    try:
        result = math.ldexp(mantissa, exponent)
    except OverflowError as exc:
        raise ValueError("Hall reference result exceeds finite float64 range") from exc
    if not math.isfinite(result):
        raise ValueError("Hall reference result exceeds finite float64 range")
    return result


def _array(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must not be masked")
    raw = np.asarray(value, dtype=object)
    if raw.ndim != 1 or raw.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional array")
    return np.array([_real(item, name) for item in raw], dtype=float)


def _immutable(value):
    array = np.ascontiguousarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("Hall reference fields must remain finite")
    return np.frombuffer(array.tobytes(), dtype=float).reshape(array.shape)


@dataclass(frozen=True, slots=True)
class PeriodicHallInduction:
    resistivity_ohm_m: float
    permeability_h_m: float
    electron_density_m3: float
    background_field_t: float
    perturbation_field_t: float
    wavelength_m: float
    helicity: int
    phase_rad: float

    def __post_init__(self):
        for name in ("resistivity_ohm_m", "permeability_h_m", "electron_density_m3",
                     "background_field_t", "perturbation_field_t", "wavelength_m", "phase_rad"):
            object.__setattr__(self, name, _real(getattr(self, name), name,
                positive=name in {"permeability_h_m", "electron_density_m3", "wavelength_m"},
                nonnegative=name in {"resistivity_ohm_m", "perturbation_field_t"}))
        if isinstance(self.helicity, (bool, np.bool_)) or not isinstance(self.helicity, (int, np.integer)) or self.helicity not in {-1, 1}:
            raise ValueError("helicity must be explicitly integer -1 or +1")
        # Reject an unrepresentable reference before generating partial output.
        for value in (self.wavenumber_m_inv, self.decay_rate_s, self.angular_frequency_s):
            if not math.isfinite(value):
                raise ValueError("derived Hall coefficients must be finite")

    @property
    def wavenumber_m_inv(self):
        return _scaled((float(self.helicity), math.tau), (self.wavelength_m,))

    @property
    def decay_rate_s(self):
        k = self.wavenumber_m_inv
        return _scaled((self.resistivity_ohm_m, k, k), (self.permeability_h_m,))

    @property
    def angular_frequency_s(self):
        k = self.wavenumber_m_inv
        return _scaled((self.background_field_t, k, k),
                       (self.permeability_h_m, self.electron_density_m3, QE))

    def _amplitude(self, time_s):
        time = _real(time_s, "time_s", nonnegative=True)
        attenuation = self._attenuation(time)
        return time, _scaled((self.perturbation_field_t,), attenuation=attenuation)

    def _attenuation(self, time):
        k = self.wavenumber_m_inv
        try:
            return _scaled((self.resistivity_ohm_m, k, k, time), (self.permeability_h_m,))
        except ValueError:
            # A positive exponent beyond the finite range proves complete
            # representable damping; it need not itself be stored finitely.
            return math.inf

    def _fields(self, z, time, amplitude, averaging):
        k, mu, ne = self.wavenumber_m_inv, self.permeability_h_m, self.electron_density_m3
        attenuation = self._attenuation(time)
        initial = self.perturbation_field_t
        current_amplitude = _scaled((-k, initial), (mu,), attenuation=attenuation)
        resistive_amplitude = _scaled((-self.resistivity_ohm_m, k, initial), (mu,), attenuation=attenuation)
        hall_scale = _scaled((self.background_field_t, k, initial), (mu, ne, QE), attenuation=attenuation)
        electron_amplitude = _scaled((k, initial), (mu, ne, QE), attenuation=attenuation)
        if all(value == 0 for value in (amplitude, current_amplitude, resistive_amplitude, hall_scale, electron_amplitude)):
            cosine, sine = np.zeros_like(z), np.zeros_like(z)
        else:
            with np.errstate(over="ignore", invalid="ignore"):
                temporal_phase = _scaled((self.background_field_t, k, k, time), (mu, ne, QE))
                phase = k * z + temporal_phase + self.phase_rad
            if not np.all(np.isfinite(phase)):
                raise ValueError("requested phase exceeds finite float64 range")
            cosine, sine = np.cos(phase) * averaging, np.sin(phase) * averaging
        zero = np.zeros_like(z)
        transverse = np.column_stack((cosine, sine, zero))
        magnetic = amplitude * transverse
        magnetic[:, 2] = self.background_field_t
        current = current_amplitude * transverse
        resistive = resistive_amplitude * transverse
        hall = hall_scale * np.column_stack((-sine, cosine, zero))
        electron = electron_amplitude * transverse
        power = _scaled((self.resistivity_ohm_m, k, k, initial, initial), (mu, mu), attenuation=2 * attenuation)
        return {name: _immutable(value) for name, value in {
            "magnetic_field_t": magnetic,
            "current_density_a_m2": current,
            "resistive_electric_field_v_m": resistive,
            "hall_electric_field_v_m": hall,
            "electric_field_v_m": resistive + hall,
            "electron_velocity_m_s": electron,
            "joule_power_w_m3": np.full_like(z, power),
            "hall_power_w_m3": zero,
        }.items()}

    def fields(self, z_m, time_s):
        """Point fields at periodic coordinates and one nonnegative time."""
        z = _array(z_m, "z_m")
        time, amplitude = self._amplitude(time_s)
        return self._fields(z, time, amplitude, 1.0)

    def cell_averages(self, z_edges_m, time_s):
        """Exact local cell means; mean Joule power is not rho*|mean J|²."""
        edges = _array(z_edges_m, "z_edges_m")
        with np.errstate(over="ignore", invalid="ignore"):
            widths = np.diff(edges)
        if len(edges) < 2 or not np.all(np.isfinite(widths)) or np.any(widths <= 0):
            raise ValueError("z_edges_m must be finite strictly increasing cell edges")
        centers = edges[:-1] / 2 + edges[1:] / 2
        time, amplitude = self._amplitude(time_s)
        # np.sinc(x) = sin(pi*x)/(pi*x), with its regular limit at zero.
        with np.errstate(over="ignore", invalid="ignore"):
            argument = widths / self.wavelength_m
        if not np.all(np.isfinite(argument)):
            raise ValueError("cell phase extent exceeds finite float64 range")
        return self._fields(centers, time, amplitude, np.sinc(argument))

    def energy_budget(self, time_s):
        """Periodic-volume energy densities; Hall power and net flux are zero.

        Ions are prescribed stationary, so the induction subsystem has no
        fluid-work reservoir. The constant background energy is separate.
        """
        time = _real(time_s, "time_s", nonnegative=True)
        initial = _scaled((0.5, self.perturbation_field_t, self.perturbation_field_t), (self.permeability_h_m,))
        if time == 0 or self.resistivity_ohm_m == 0:
            remaining, heat = initial, 0.0
        else:
            k = self.wavenumber_m_inv
            log_hazard = math.fsum((math.log(2.0), math.log(self.resistivity_ohm_m),
                                   2 * math.log(abs(k)), math.log(time), -math.log(self.permeability_h_m)))
            remaining, heat = _deplete(initial, log_hazard)
        return {
            "initial_perturbation_energy_j_m3": initial,
            "magnetic_perturbation_energy_j_m3": remaining,
            "joule_heat_j_m3": heat,
            "background_magnetic_energy_j_m3": _scaled((0.5, self.background_field_t, self.background_field_t), (self.permeability_h_m,)),
            "hall_work_j_m3": 0.0,
            "balance_residual_j_m3": math.fsum((remaining, -initial, heat)),
        }
