"""Finite-volume heating comparisons and provenance-preserving SI exchange."""

from dataclasses import FrozenInstanceError
import csv
import hashlib
import json
import math

import numpy as np
import pytest

from mhdlab.heating_benchmarks import ConstantConductivityHalfSpace
from mhdlab.heating_comparison import (
    HeatingProfile,
    compare_heating_profile,
    load_heating_profile_csv,
    write_heating_profile_csv,
)


FIELDS = ("magnetic_field_t", "current_density_a_m2", "joule_power_w_m3")
CSV_COLUMNS = ("depth_left_m", "depth_right_m", "bx_t", "jz_a_m2", "joule_w_m3")


def benchmark(**overrides):
    parameters = dict(conductivity_s_m=3e6, permeability_h_m=4e-7 * math.pi,
                      volumetric_heat_capacity_j_m3_k=2.8e6, boundary_field_t=5.)
    parameters.update(overrides)
    return ConstantConductivityHalfSpace(**parameters)


def profile(**overrides):
    parameters = dict(time_s=2e-6, depth_edges_m=[0., 1e-4, 4e-4],
                      magnetic_field_t=[-2., 3.], current_density_a_m2=[-4., 5.],
                      joule_power_w_m3=[6., 7.], source_kind="synthetic",
                      run_id="manufactured-profile", source_reference="manufactured SI fixture")
    parameters.update(overrides)
    return HeatingProfile(**parameters)


def analytic_profile(model=None, **overrides):
    model = benchmark() if model is None else model
    time = 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    edges = np.array([0., .07, .4, 1.1, 3.5, 8.]) * scale
    cells = model.cell_averages(edges, time)
    parameters = dict(time_s=time, depth_edges_m=edges,
                      magnetic_field_t=cells.magnetic_field_t,
                      current_density_a_m2=cells.current_density_a_m2,
                      joule_power_w_m3=cells.joule_power_w_m3)
    parameters.update(overrides)
    return profile(**parameters)


def independent_norms(values, edges):
    """Straight polynomial moments suffice for the ordinary-sized fixtures."""
    widths = np.diff(edges)
    span = edges[-1] - edges[0]
    return {"l1": float(np.dot(widths, np.abs(values)) / span),
            "l2": math.sqrt(float(np.dot(widths, np.square(values)) / span)),
            "linf": float(np.max(np.abs(values)))}


def test_exact_cell_averages_have_zero_error_and_json_safe_metadata():
    model = benchmark(boundary_field_t=-5.)
    sample = analytic_profile(model)
    result = compare_heating_profile(sample, model)
    assert set(result["errors"]) == set(FIELDS)
    for name in FIELDS:
        record = result["errors"][name]
        assert record["absolute"] == dict(l1=0., l2=0., linf=0.)
        assert record["relative"] == dict(l1=0., l2=0., linf=0.)
        expected = independent_norms(getattr(sample, name), sample.depth_edges_m)
        assert record["reference_norms"] == pytest.approx(expected, rel=3e-13, abs=0.)
    encoded = json.dumps(result, allow_nan=False)
    assert "temperature" not in result["errors"]
    assert sample.source_reference in encoded
    assumptions = " ".join(result["reference"]["assumptions"])
    for condition in ("half-space", "constant", "conductivity", "permeability",
                      "initial magnetic field zero", "boundary field", "motion", "phase change"):
        assert condition in assumptions
    assert result["reference"]["conductivity_s_m"] == model.conductivity_s_m
    assert result["reference"]["boundary_field_t"] == model.boundary_field_t
    assert result["grid"]["sampling"] == "cell_average"


def test_comparison_uses_mean_of_joule_power_not_square_of_mean_current():
    model = benchmark()
    exact = analytic_profile(model)
    wrong_q = exact.current_density_a_m2**2 / model.conductivity_s_m
    assert wrong_q[-1] < .5 * exact.joule_power_w_m3[-1]
    sample = analytic_profile(model, joule_power_w_m3=wrong_q)
    result = compare_heating_profile(sample, model)
    assert result["errors"]["magnetic_field_t"]["absolute"]["linf"] == 0.
    assert result["errors"]["current_density_a_m2"]["absolute"]["linf"] == 0.
    expected = independent_norms(wrong_q - exact.joule_power_w_m3, exact.depth_edges_m)
    assert result["errors"]["joule_power_w_m3"]["absolute"] == pytest.approx(expected, rel=3e-13, abs=0.)


def test_constant_perturbations_have_constant_absolute_norms_and_correct_relatives():
    model = benchmark()
    exact = analytic_profile(model)
    perturbations = dict(magnetic_field_t=-.25, current_density_a_m2=3.5,
                         joule_power_w_m3=13.)
    sample = analytic_profile(model, **{name: getattr(exact, name) + delta
                                       for name, delta in perturbations.items()})
    result = compare_heating_profile(sample, model)
    for name, delta in perturbations.items():
        expected_reference = independent_norms(getattr(exact, name), exact.depth_edges_m)
        for norm in ("l1", "l2", "linf"):
            assert result["errors"][name]["absolute"][norm] == pytest.approx(abs(delta), rel=1e-10, abs=0.)
            assert result["errors"][name]["relative"][norm] == pytest.approx(
                abs(delta) / expected_reference[norm], rel=1e-10, abs=0.)


def test_nonuniform_cells_are_weighted_by_physical_width():
    model = benchmark(boundary_field_t=0.)
    edges = np.array([0., 1., 4., 10.])
    sample = profile(depth_edges_m=edges, magnetic_field_t=[-2., 5., -1.],
                     current_density_a_m2=[3., -4., 8.], joule_power_w_m3=[0., 7., 2.])
    result = compare_heating_profile(sample, model)
    for name in FIELDS:
        expected = independent_norms(getattr(sample, name), edges)
        assert result["errors"][name]["absolute"] == pytest.approx(expected, rel=2e-14, abs=0.)
        assert result["errors"][name]["reference_norms"] == dict(l1=0., l2=0., linf=0.)
        assert result["errors"][name]["relative"] == dict(l1=None, l2=None, linf=None)
    assert json.loads(json.dumps(result, allow_nan=False))["errors"][FIELDS[0]]["relative"]["l2"] is None


def test_zero_profiles_do_not_turn_zero_reference_relative_error_into_zero():
    sample = profile(magnetic_field_t=[0., 0.], current_density_a_m2=[0., 0.],
                     joule_power_w_m3=[0., 0.])
    result = compare_heating_profile(sample, benchmark(boundary_field_t=0.))
    for name in FIELDS:
        assert result["errors"][name]["absolute"] == dict(l1=0., l2=0., linf=0.)
        assert result["errors"][name]["relative"] == dict(l1=None, l2=None, linf=None)


def test_norms_keep_large_values_without_overflowing_squared_moments():
    sample = profile(depth_edges_m=[0., 1., 2.], magnetic_field_t=[1e250, -1e250],
                     current_density_a_m2=[-1e250, 1e250], joule_power_w_m3=[1e250, 1e250])
    result = compare_heating_profile(sample, benchmark(boundary_field_t=0.))
    for name in FIELDS:
        for value in result["errors"][name]["absolute"].values():
            assert value == pytest.approx(1e250, rel=3e-13, abs=0.)
    json.dumps(result, allow_nan=False)


def test_norms_keep_tiny_width_contributions_that_restore_finite_moments():
    sample = profile(depth_edges_m=[0., 1e-300, 1e300], magnetic_field_t=[1e300, 0.],
                     current_density_a_m2=[1e300, 0.], joule_power_w_m3=[1e300, 0.])
    result = compare_heating_profile(sample, benchmark(boundary_field_t=0.))
    for name in FIELDS:
        assert result["errors"][name]["absolute"] == pytest.approx(
            dict(l1=1e-300, l2=1., linf=1e300), rel=6e-13, abs=0.)


def test_constant_maximum_finite_profile_does_not_overflow_width_normalized_norms():
    maximum = np.finfo(float).max
    sample = profile(depth_edges_m=[0., 1e-300, 1e300],
                     magnetic_field_t=[maximum, maximum], current_density_a_m2=[0., 0.],
                     joule_power_w_m3=[0., 0.])
    result = compare_heating_profile(sample, benchmark(boundary_field_t=0.))
    for value in result["errors"]["magnetic_field_t"]["absolute"].values():
        assert value == pytest.approx(maximum, rel=5e-13, abs=0.)


def test_underflowed_reference_moments_keep_nonzero_reference_and_relative_ratios():
    model = benchmark(conductivity_s_m=1., permeability_h_m=1.,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_t=1e-300)
    sample = profile(time_s=1., depth_edges_m=[0., 1e-300, 1e300],
                     magnetic_field_t=[0., 0.], current_density_a_m2=[0., 0.],
                     joule_power_w_m3=[0., 0.])
    result = compare_heating_profile(sample, model)
    for name in ("magnetic_field_t", "current_density_a_m2"):
        record = result["errors"][name]
        assert record["reference_is_zero"] is False
        assert record["reference_norms"]["l1"] == 0.
        assert record["reference_norms"]["l2"] == 0.
        assert record["reference_norms"]["linf"] > 0.
        assert record["relative"] == dict(l1=1., l2=1., linf=1.)


def test_unrepresentable_opposite_signed_error_is_rejected_instead_of_infinity():
    model = benchmark(conductivity_s_m=1e-308, permeability_h_m=1e308,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_t=-1e308)
    sample = profile(time_s=1., depth_edges_m=[0., 1e-10, 2e-10],
                     magnetic_field_t=[1e308, 1e308], current_density_a_m2=[0., 0.],
                     joule_power_w_m3=[0., 0.])
    with pytest.raises(ValueError):
        compare_heating_profile(sample, model)


def test_unrepresentable_relative_ratio_is_rejected_instead_of_json_infinity():
    model = benchmark(conductivity_s_m=1., permeability_h_m=1.,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_t=1e-300)
    sample = profile(time_s=1., depth_edges_m=[0., 1., 2.],
                     magnetic_field_t=[1e300, 1e300], current_density_a_m2=[0., 0.],
                     joule_power_w_m3=[0., 0.])
    with pytest.raises(ValueError):
        compare_heating_profile(sample, model)


def test_arrays_are_copied_into_an_immutable_profile():
    edges = np.array([0., 1., 2.])
    magnetic = np.array([-2., 3.])
    sample = profile(depth_edges_m=edges, magnetic_field_t=magnetic)
    edges[1], magnetic[0] = .2, 10.
    np.testing.assert_array_equal(sample.depth_edges_m, [0., 1., 2.])
    np.testing.assert_array_equal(sample.magnetic_field_t, [-2., 3.])
    for name in ("depth_edges_m",) + FIELDS:
        with pytest.raises(ValueError):
            getattr(sample, name).setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        sample.time_s = 10.


@pytest.mark.parametrize("overrides", [
    {"time_s": 0.}, {"time_s": -1.}, {"time_s": np.nan}, {"time_s": np.inf},
    {"time_s": True}, {"time_s": 2e-6 + 0j}, {"time_s": "2e-6"},
    {"depth_edges_m": []}, {"depth_edges_m": [0.]}, {"depth_edges_m": [0., 0., 1.]},
    {"depth_edges_m": [0., 2., 1.]}, {"depth_edges_m": [-1., 0., 1.]},
    {"depth_edges_m": [0., np.nan, 1.]}, {"depth_edges_m": [0., 1., np.inf]},
    {"depth_edges_m": [[0., 1., 2.]]}, {"depth_edges_m": [False, 1., 2.]},
    {"depth_edges_m": [0., 1.+0j, 2.]}, {"depth_edges_m": ["0", "1", "2"]},
    {"magnetic_field_t": [1.]}, {"current_density_a_m2": [[1., 2.]]},
    {"joule_power_w_m3": [1., -1.]}, {"magnetic_field_t": [1., np.inf]},
    {"current_density_a_m2": [1., np.nan]}, {"joule_power_w_m3": [1., False]},
    {"magnetic_field_t": [0j, 2.]}, {"current_density_a_m2": ["1", "2"]},
    {"source_kind": "measured"}, {"source_kind": ""}, {"run_id": " "},
    {"source_reference": ""}, {"run_id": None}, {"source_reference": True},
    {"source_kind": "flash", "flash_version": ""}, {"flash_version": 4.8},
    {"depth_edges_m": np.ma.array([0., 1., 2.], mask=[False, True, False])},
    {"joule_power_w_m3": np.ma.array([1., 2.], mask=[False, True])},
])
def test_profile_rejects_ambiguous_geometry_invalid_values_and_missing_provenance(overrides):
    with pytest.raises(ValueError):
        profile(**overrides)


def test_flash_provenance_is_required_and_preserved():
    sample = profile(source_kind="flash", flash_version="FLASH 4.8 custom application")
    assert sample.source_kind == "flash"
    assert sample.flash_version == "FLASH 4.8 custom application"


def test_canonical_csv_round_trip_is_lossless_and_requires_metadata(tmp_path):
    sample = profile(source_kind="flash", flash_version="FLASH 4.8",
                     magnetic_field_t=[-np.finfo(float).max, np.nextafter(1., 2.)],
                     current_density_a_m2=[np.nextafter(0., 1.), -np.nextafter(2., 3.)],
                     joule_power_w_m3=[np.nextafter(0., 1.), np.finfo(float).max])
    csv_path, metadata_path = write_heating_profile_csv(tmp_path / "lineout.csv", sample)
    assert csv_path == tmp_path / "lineout.csv"
    assert metadata_path == tmp_path / "lineout.csv.json"
    with csv_path.open(newline="") as handle:
        assert tuple(next(csv.reader(handle))) == CSV_COLUMNS
    loaded = load_heating_profile_csv(csv_path)
    for name in ("depth_edges_m",) + FIELDS:
        np.testing.assert_array_equal(getattr(loaded, name), getattr(sample, name))
    for name in ("time_s", "source_kind", "run_id", "source_reference", "flash_version"):
        assert getattr(loaded, name) == getattr(sample, name)
    metadata_path.unlink()
    with pytest.raises((ValueError, FileNotFoundError)):
        load_heating_profile_csv(csv_path)


def test_metadata_copy_cannot_change_the_profile_or_future_serialization():
    sample = profile()
    metadata = sample.metadata
    metadata["units"]["bx_t"] = "G"
    metadata["source_kind"] = "flash"
    assert sample.metadata["units"]["bx_t"] == "T"
    assert sample.source_kind == "synthetic"


def test_csv_and_metadata_cannot_overwrite_each_other(tmp_path):
    path = tmp_path / "lineout.csv"
    path.write_text("existing-content")
    with pytest.raises(ValueError):
        write_heating_profile_csv(path, profile(), metadata_path=path)
    assert path.read_text() == "existing-content"


def test_explicit_metadata_path_and_csv_content_hash(tmp_path):
    csv_path = tmp_path / "lineout.csv"
    metadata_path = tmp_path / "provenance.json"
    written = write_heating_profile_csv(csv_path, profile(), metadata_path=metadata_path)
    assert written == (csv_path, metadata_path)
    load_heating_profile_csv(csv_path, metadata_path=metadata_path)
    # A numerically plausible edit still needs the original content digest.
    content = csv_path.read_text().replace("-2", "-3", 1)
    csv_path.write_text(content)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path, metadata_path=metadata_path)


def read_metadata(path):
    return json.loads(path.read_text())


def write_metadata(path, metadata):
    path.write_text(json.dumps(metadata, allow_nan=False))


def refresh_digest(csv_path, metadata_path):
    metadata = read_metadata(metadata_path)
    metadata["csv_sha256"] = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    write_metadata(metadata_path, metadata)


def test_metadata_records_fixed_units_sampling_time_and_provenance(tmp_path):
    sample = profile(source_kind="flash", flash_version="FLASH 4.8")
    csv_path, metadata_path = write_heating_profile_csv(tmp_path / "lineout.csv", sample)
    metadata = read_metadata(metadata_path)
    assert metadata["schema"] == "mhdlab-heating-profile-v1"
    assert metadata["units"] == dict(time_s="s", depth_left_m="m", depth_right_m="m",
                                      bx_t="T", jz_a_m2="A m^-2", joule_w_m3="W m^-3")
    assert metadata["cell_sampling"] == "cell_average"
    assert metadata["coordinate_system"] == "y=conductor_inward_depth;Bx=tangential;Jz=-dBx/dy/mu"
    assert metadata["csv_sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    for name in ("time_s", "source_kind", "run_id", "source_reference", "flash_version"):
        assert metadata[name] == getattr(sample, name)


@pytest.mark.parametrize("key,value", [
    ("schema", "unversioned"), ("time_s", 0.), ("time_s", True),
    ("source_kind", "measured"), ("run_id", ""), ("source_reference", " "),
    ("flash_version", ""), ("coordinate_system", "y=vacuum_gap_depth"),
    ("cell_sampling", "cell_center"), ("csv_sha256", "0" * 64),
])
def test_loader_rejects_forged_metadata_and_missing_required_keys(tmp_path, key, value):
    csv_path, metadata_path = write_heating_profile_csv(
        tmp_path / "lineout.csv", profile(source_kind="flash", flash_version="FLASH 4.8"))
    metadata = read_metadata(metadata_path)
    metadata[key] = value
    write_metadata(metadata_path, metadata)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)
    del metadata[key]
    write_metadata(metadata_path, metadata)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)


@pytest.mark.parametrize("column,unit", [
    ("depth_left_m", "mm"), ("depth_right_m", "cm"), ("bx_t", "G"),
    ("jz_a_m2", "A cm^-2"), ("joule_w_m3", "W cm^-3"), ("time_s", "ns"),
])
def test_loader_rejects_unit_reinterpretation(tmp_path, column, unit):
    csv_path, metadata_path = write_heating_profile_csv(tmp_path / "lineout.csv", profile())
    metadata = read_metadata(metadata_path)
    metadata["units"][column] = unit
    write_metadata(metadata_path, metadata)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)


@pytest.mark.parametrize("rows", [
    # A gap and an overlap must not turn into an apparently contiguous edge array.
    [["0", "1", "1", "1", "1"], ["1.1", "2", "1", "1", "1"]],
    [["0", "1", "1", "1", "1"], ["0.9", "2", "1", "1", "1"]],
    [["1", "2", "1", "1", "1"], ["0", "1", "1", "1", "1"]],
    [["0", "1", "1", "1", "1"], ["1", "1", "1", "1", "1"]],
    [["0", "1", "NaN", "1", "1"], ["1", "2", "1", "1", "1"]],
    [["0", "1", "1", "inf", "1"], ["1", "2", "1", "1", "1"]],
    [["0", "1", "1", "1", "-1"], ["1", "2", "1", "1", "1"]],
    [["0", "1", "1", "1"], ["1", "2", "1", "1", "1"]],
    [["0", "1", "1", "1", "1", "0"], ["1", "2", "1", "1", "1"]],
    [],
])
def test_loader_revalidates_cells_even_if_digest_is_recomputed(tmp_path, rows):
    csv_path, metadata_path = write_heating_profile_csv(tmp_path / "lineout.csv", profile())
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        writer.writerows(rows)
    refresh_digest(csv_path, metadata_path)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)


@pytest.mark.parametrize("header", [
    ["depth_left_m", "depth_right_m", "bx_t", "jz_a_m2", "joule_w_m3", "temperature_k"],
    ["depth_left_m", "depth_right_m", "bx_gauss", "jz_a_m2", "joule_w_m3"],
    ["depth_right_m", "depth_left_m", "bx_t", "jz_a_m2", "joule_w_m3"],
    ["depth_left_m", "depth_right_m", "bx_t", "bx_t", "joule_w_m3"],
])
def test_loader_requires_exact_canonical_header(tmp_path, header):
    csv_path, metadata_path = write_heating_profile_csv(tmp_path / "lineout.csv", profile())
    rows = list(csv.reader(csv_path.read_text().splitlines()))
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows[1:])
    refresh_digest(csv_path, metadata_path)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)


def test_json_duplicate_keys_and_nonfinite_constants_are_rejected(tmp_path):
    csv_path, metadata_path = write_heating_profile_csv(tmp_path / "lineout.csv", profile())
    original = metadata_path.read_text()
    duplicate = original.rstrip()[:-1] + ',"time_s":2e-6}'
    assert json.loads(duplicate)["time_s"] == 2e-6
    metadata_path.write_text(duplicate)
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)
    metadata = json.loads(original)
    metadata["time_s"] = float("nan")
    metadata_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError):
        load_heating_profile_csv(csv_path)
