import numpy as np

from mhdlab.traces import fit_effective_rl


def test_fit_effective_rl_recovers_synthetic_values():
    t = np.linspace(0.0, 1.0e-6, 200)
    current = 1.0e6 * (1.0 - np.exp(-t / 2.0e-7))
    r_true = 0.08
    l_true = 12.0e-9
    voltage = r_true * current + l_true * np.gradient(current, t)
    fit = fit_effective_rl(t, voltage, current)
    assert np.isclose(fit["resistance_ohm"], r_true, rtol=1e-3)
    assert np.isclose(fit["inductance_h"], l_true, rtol=1e-3)
