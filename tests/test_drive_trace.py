"""Analytic interpolation, timing, and piecewise-linear integral checks."""

import hashlib

import numpy as np
import pytest

from mhdlab.drive_trace import DriveTrace, load_current_csv


def trace(**overrides):
    values = dict(time_s=[-1e-9, 1e-9, 4e-9], current_a=[-2.0, 2.0, 5.0],
                  source_kind="synthetic", shot_id="analytic-check", source_reference="analytic linear trace",
                  diagnostic_time_offset_s=0.0)
    values.update(overrides)
    return DriveTrace(**values)


def loader(path, **overrides):
    options = dict(time_column="delay_ns", current_column="current_kA", time_scale_to_s=1e-9,
                   current_scale_to_a=1e3, diagnostic_time_offset_s=-2e-9,
                   source_kind="reconstructed", shot_id="figure-check", source_reference="published figure test fixture",
                   figure_reference="Figure X, current curve")
    options.update(overrides)
    return load_current_csv(path, **options)


def test_explicit_csv_units_time_offset_and_figure_provenance(tmp_path):
    path = tmp_path / "trace.csv"
    content = b"delay_ns,current_kA,other\n1,-2,ignored\n3,2,ignored\n6,5,ignored\n"
    path.write_bytes(content)
    loaded = loader(path)
    np.testing.assert_allclose(loaded.time_s, [-1e-9, 1e-9, 4e-9], rtol=1e-15, atol=0.0)
    np.testing.assert_array_equal(loaded.current_a, [-2000.0, 2000.0, 5000.0])
    assert loaded.source_kind == "reconstructed"
    assert loaded.figure_reference == "Figure X, current curve"
    assert loaded.diagnostic_time_offset_s == -2e-9
    assert loaded.csv_import_metadata["sha256"] == hashlib.sha256(content).hexdigest()
    assert loaded.csv_import_metadata["current_scale_to_a"] == 1e3


def test_interpolation_preserves_sign_and_finite_support():
    sample = trace()
    np.testing.assert_allclose(sample.current_at([-1e-9, 0.0, 2e-9, 4e-9]), [-2.0, 0.0, 3.0, 5.0],
                               rtol=2e-15, atol=5.0 * np.finfo(float).eps)
    for query in (-1.01e-9, 4.01e-9, np.nan, True, 0j, [], [[0.0]]):
        with pytest.raises(ValueError):
            sample.current_at(query)


def test_piecewise_slopes_require_explicit_knot_policy():
    sample = trace()
    np.testing.assert_allclose(sample.slope_a_s([0.0, 2e-9], knot_policy="raise"), [2e9, 1e9], rtol=2e-15, atol=0.0)
    assert sample.slope_a_s(1e-9, knot_policy="left") == pytest.approx(2e9, rel=2e-15, abs=0.0)
    assert sample.slope_a_s(1e-9, knot_policy="right") == pytest.approx(1e9, rel=2e-15, abs=0.0)
    np.testing.assert_allclose(sample.slope_a_s([-1e-9, 4e-9], knot_policy="raise"), [2e9, 1e9], rtol=2e-15, atol=0.0)
    with pytest.raises(ValueError, match="interior knot"):
        sample.slope_a_s(1e-9, knot_policy="raise")
    with pytest.raises(ValueError, match="knot_policy"):
        sample.slope_a_s(0.0, knot_policy="average")


def test_current_and_squared_current_integrals_match_analytic_polynomials():
    sample = trace(time_s=[0.0, 2.0, 5.0], current_a=[0.0, 4.0, 10.0])
    # I=2t over both unequal segments: int I dt=t^2; int I^2 dt=4t^3/3.
    assert sample.integrate_current_a_s(0.5, 4.0) == pytest.approx(16.0 - .25, rel=2e-15, abs=0.0)
    assert sample.integrate_current_squared_a2_s(0.5, 4.0) == pytest.approx(4.0 * (64.0 - .125) / 3.0, rel=2e-15, abs=0.0)
    total = sample.integrate_current_squared_a2_s(0.0, 5.0)
    pieces = sample.integrate_current_squared_a2_s(0.0, 2.0) + sample.integrate_current_squared_a2_s(2.0, 5.0)
    assert total == pytest.approx(pieces, rel=2e-15, abs=0.0)


def test_signed_zero_crossing_has_zero_current_integral_and_positive_square_integral():
    sample = trace(time_s=[0.0, 1e-9], current_a=[-2e5, 2e5])
    assert sample.integrate_current_a_s(0.0, 1e-9) == 0.0
    assert sample.integrate_current_squared_a2_s(0.0, 1e-9) == pytest.approx(4e10 * 1e-9 / 3.0, rel=2e-15, abs=0.0)
    assert sample.integrate_current_a_s(.5e-9, .5e-9) == 0.0
    assert sample.integrate_current_squared_a2_s(.5e-9, .5e-9) == 0.0


def test_opposite_large_currents_interpolate_without_overflowing_their_difference():
    sample = trace(time_s=[0.0, 1.0], current_a=[-1e308, 1e308])
    assert sample.current_at(0.5) == 0.0
    assert sample.current_at(0.25) == pytest.approx(-5e307, rel=2e-15, abs=0.0)
    assert sample.integrate_current_a_s(0.0, 1.0) == 0.0


@pytest.mark.parametrize("start,end", [(-2e-9, 0.0), (0.0, 5e-9), (2e-9, 0.0), (np.nan, 1e-9), (False, 1e-9)])
def test_integrals_reject_unsupported_or_reversed_intervals(start, end):
    for operation in (trace().integrate_current_a_s, trace().integrate_current_squared_a2_s):
        with pytest.raises(ValueError):
            operation(start, end)


@pytest.mark.parametrize("overrides", [
    {"time_s": [0.0]}, {"time_s": [0.0, 0.0]}, {"time_s": [1.0, 0.0]},
    {"time_s": [0.0, np.nan]}, {"time_s": [[0.0, 1.0]]}, {"current_a": [1.0, 2.0]},
    {"current_a": [True, 1.0, 2.0]}, {"current_a": [1.0, 2j, 3.0]},
    {"current_a": [1.0, np.inf, 3.0]}, {"source_kind": "digitized"}, {"shot_id": " "},
    {"source_reference": ""}, {"source_kind": "reconstructed", "figure_reference": ""},
    {"diagnostic_time_offset_s": True}, {"diagnostic_time_offset_s": np.inf},
])
def test_invalid_trace_contract_is_rejected(overrides):
    with pytest.raises(ValueError):
        trace(**overrides)


def test_arrays_and_csv_provenance_are_copied_and_read_only():
    time, current = np.array([0.0, 1.0]), np.array([1.0, 2.0])
    metadata = {"column": "I"}
    sample = trace(time_s=time, current_a=current, csv_import_metadata=metadata)
    current[0] = 50.0
    time[0] = -5.0
    metadata["column"] = "changed"
    np.testing.assert_array_equal(sample.time_s, [0.0, 1.0])
    np.testing.assert_array_equal(sample.current_a, [1.0, 2.0])
    assert sample.csv_import_metadata["column"] == "I"
    with pytest.raises(ValueError):
        sample.current_a[0] = 100.0
    with pytest.raises(TypeError):
        sample.csv_import_metadata["column"] = "I2"


@pytest.mark.parametrize("content", [
    "time,current\n0,1\n1,2\n", "delay_ns,current_kA\n0,1\n0,2\n",
    "delay_ns,current_kA\n0,1\n", "delay_ns,current_kA\n0,nan\n1,2\n",
    "delay_ns,current_kA\n0,\n1,2\n", "delay_ns,current_kA\n0,1,extra\n1,2\n",
    "delay_ns,current_kA,current_kA\n0,1,2\n1,3,4\n",
])
def test_malformed_csv_is_rejected(tmp_path, content):
    path = tmp_path / "bad.csv"
    path.write_text(content)
    with pytest.raises(ValueError):
        loader(path)


@pytest.mark.parametrize("override", [{"time_scale_to_s": 0.0}, {"current_scale_to_a": -1.0},
                                     {"time_scale_to_s": True}, {"diagnostic_time_offset_s": np.nan}])
def test_loader_requires_explicit_valid_conversion_and_alignment(tmp_path, override):
    path = tmp_path / "trace.csv"
    path.write_text("delay_ns,current_kA\n0,1\n1,2\n")
    with pytest.raises(ValueError):
        loader(path, **override)
