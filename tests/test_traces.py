import numpy as np

from mhdlab.traces import drive_from_parametric_current, fit_effective_rl, parametric_current


def test_fit_effective_rl_recovers_synthetic_values():
    t = np.linspace(0.0, 1.0e-6, 200)
    current = 1.0e6 * (1.0 - np.exp(-t / 2.0e-7))
    r_true = 0.08
    l_true = 12.0e-9
    voltage = r_true * current + l_true * np.gradient(current, t)
    fit = fit_effective_rl(t, voltage, current)
    assert np.isclose(fit["resistance_ohm"], r_true, rtol=1e-3)
    assert np.isclose(fit["inductance_h"], l_true, rtol=1e-3)


def test_parametric_half_sine_current_uses_rise_time_as_quarter_period():
    rise = 125.0e-9
    t = np.asarray([0.0, rise, 2.0 * rise, 3.0 * rise])
    current = parametric_current(t, peak_current_a=850.0e3, rise_time_s=rise, waveform="half_sine")
    assert np.isclose(current[0], 0.0)
    assert np.isclose(current[1], 850.0e3)
    assert np.isclose(current[2], 0.0, atol=1.0e-9)
    assert np.isclose(current[3], 0.0)


def test_parametric_drive_profile_defaults_to_zero_voltage():
    t = np.linspace(0.0, 250.0e-9, 11)
    drive = drive_from_parametric_current(
        t,
        {
            "peak_current_a": 850.0e3,
            "rise_time_s": 125.0e-9,
            "waveform": "half_sine",
        },
    )
    assert np.isclose(drive.sample_current(125.0e-9), 850.0e3)
    assert np.isclose(drive.sample_voltage(125.0e-9), 0.0)
    assert drive.metadata["mode"] == "parametric_current"
