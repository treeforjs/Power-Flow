"""Temperature-limited electron supply; no field, sheath, or space-charge closure.

All surface parameters are explicit. The ideal free-electron Richardson constant
is a reference, not a calibrated effective constant for a particular foil.
"""
from __future__ import annotations
import math
from numbers import Real
import numpy as np
from .constants import KB, QE

ME_REFERENCE_KG = 9.1093837015e-31
PLANCK_J_S = 6.62607015e-34
IDEAL_RICHARDSON_A_M2_K2 = 4 * math.pi * ME_REFERENCE_KG * QE * KB**2 / PLANCK_J_S**3


def log_thermionic_supply(temperature_k, *, work_function_ev, richardson_a_m2_k2):
    """Return ln(J/[A m^-2]) for the zero-field Richardson-Dushman supply."""
    raw = np.asarray(temperature_k)
    if raw.dtype.kind not in 'iuf':
        raise ValueError('surface temperatures must be real numeric values')
    temperature = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(temperature)) or np.any(temperature <= 0):
        raise ValueError('surface temperatures must be finite and positive')
    for name, value in [('work_function_ev', work_function_ev), ('richardson_a_m2_k2', richardson_a_m2_k2)]:
        if isinstance(value, (bool,np.bool_)) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be explicitly finite and positive')
    return math.log(richardson_a_m2_k2) + 2*np.log(temperature) - work_function_ev*QE/(KB*temperature)


def thermionic_supply(temperature_k, *, work_function_ev, richardson_a_m2_k2):
    """Return the unrestricted supply in A/m², with no collection-current cap.

    Numerical underflow can yield zero; the companion logarithmic function
    retains its magnitude. Nonfinite large-current outputs are refused.
    """
    logarithm = log_thermionic_supply(temperature_k, work_function_ev=work_function_ev,
                                    richardson_a_m2_k2=richardson_a_m2_k2)
    with np.errstate(under='ignore',over='ignore'):
        current = np.exp(logarithm)
    if not np.all(np.isfinite(current)):
        raise ValueError('current exceeds the finite output range; use logarithmic supply')
    return current
