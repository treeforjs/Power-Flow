"""Analytical overlap and independently indexed conservation checks."""

import numpy as np
import pytest

from mhdlab.surface_sources import (
    conservative_remap,
    matched_uniform_moments,
    matched_uniform_source,
    source_contrast,
)


def test_uneven_bin_remap_matches_analytic_overlap_integrals():
    source_edges = [0.0, 0.1, 0.4, 1.0]
    target_edges = [0.0, 0.2, 0.8, 1.0]
    values = [10.0, 20.0, 50.0]
    result = conservative_remap(values, source_edges, target_edges)
    np.testing.assert_allclose(result, [15.0, 40.0, 50.0], rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(np.sum(result * np.diff(target_edges)), 37.0, rtol=2e-15, atol=0.0)


def test_remap_preserves_each_time_reservoir_component_with_spatial_middle_axis():
    values = np.arange(1.0, 25.0).reshape(2, 3, 4)
    source_edges = np.array([-0.002, -0.0018, 0.0008, 0.003])
    target_edges = np.array([-0.002, -0.001, 0.0, 0.0005, 0.003])
    result = conservative_remap(values, source_edges, target_edges, axis=1)
    assert result.shape == (2, 4, 4)
    before = np.sum(values * np.diff(source_edges)[None, :, None], axis=1)
    after = np.sum(result * np.diff(target_edges)[None, :, None], axis=1)
    np.testing.assert_allclose(after, before, rtol=3e-15, atol=0.0)


def test_single_bin_constant_limit_and_identity_do_not_alias_inputs():
    values = np.array([[7.0], [0.0]])
    result = conservative_remap(values, [-1.0, 2.0], [-1.0, -0.9, 0.5, 2.0])
    np.testing.assert_allclose(result, [[7.0] * 3, [0.0] * 3], rtol=1e-15, atol=0.0)
    identity = conservative_remap(values, [-1.0, 2.0], [-1.0, 2.0])
    identity[0, 0] = 100.0
    assert values[0, 0] == 7.0


def test_matched_uniform_number_preserves_every_separate_integral():
    edges = np.array([0.0, 0.1, 0.4, 1.0])
    source = np.array([[[0.0, 10.0, 20.0], [5.0, 1.0, 2.0]],
                       [[0.0, 0.0, 0.0], [3.0, 3.0, 3.0]]])
    result = matched_uniform_source(source, edges)
    expected = np.sum(source * np.diff(edges), axis=-1)
    np.testing.assert_allclose(result, np.repeat(expected[..., None], 3, axis=-1), rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(np.sum(result * np.diff(edges), axis=-1), expected, rtol=2e-15, atol=0.0)
    result[0, 0, 0] = 1000.0
    assert source[0, 0, 0] == 0.0


def test_moment_constrained_null_preserves_number_energy_and_explicit_signed_momentum():
    edges = np.array([0.0, 0.1, 0.4, 1.0])
    moments = {"number": np.array([10.0, 20.0, 50.0]),
               "energy": np.array([40.0, 10.0, 80.0]),
               "momentum_x": np.array([-1.0, 2.0, -3.0])}
    result = matched_uniform_moments(moments, edges, allow_signed={"momentum_x"})
    for name, values in moments.items():
        expected_integral = np.dot(values, np.diff(edges))
        np.testing.assert_allclose(result[name], expected_integral, rtol=2e-15, atol=0.0)
        np.testing.assert_allclose(np.dot(result[name], np.diff(edges)), expected_integral, rtol=2e-15, atol=0.0)
    assert result["momentum_x"][0] < 0.0
    assert not np.isclose(result["energy"][0] / result["number"][0], moments["energy"][0] / moments["number"][0])


def test_signed_remap_requires_an_explicit_opt_in():
    values = [-1.0, 1.0]
    with pytest.raises(ValueError, match="nonnegative"):
        conservative_remap(values, [0, 1, 2], [0, 2])
    result = conservative_remap(values, [0, 1, 2], [0, 0.5, 2], allow_signed=True)
    np.testing.assert_allclose(result, [-1.0, 1.0 / 3.0], rtol=2e-15, atol=0.0)


def test_contrast_uses_partial_bin_overlaps():
    edges = [0.0, 0.05, 0.3, 0.55, 0.9, 1.0]
    result = source_contrast([100.0, 200.0, 300.0, 400.0, 500.0], edges)
    np.testing.assert_allclose(result["left_edge_mean"], 150.0, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(result["right_edge_mean"], 500.0, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(result["center_mean"], 325.0, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(result["combined_edge_mean"], 325.0, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(result["edge_to_center_ratio"], 1.0, rtol=2e-15, atol=0.0)


def test_constant_source_and_zero_center_have_explicit_ratio_masks():
    edges = [0.0, 0.1, 0.4, 0.6, 0.9, 1.0]
    values = np.array([[7.0] * 5, [0.0] * 5, [2.0, 0.0, 0.0, 0.0, 4.0]])
    result = source_contrast(values, edges)
    ratio = result["edge_to_center_ratio"]
    np.testing.assert_array_equal(np.ma.getmaskarray(ratio), [False, True, True])
    assert ratio[0] == 1.0
    assert np.ma.is_masked(ratio[1])
    assert np.ma.is_masked(ratio[2])
    assert np.all(np.isfinite(ratio.data))
    np.testing.assert_allclose(result["edge_minus_center"], [0.0, 0.0, 3.0], atol=2e-15, rtol=2e-15)


def test_contrast_preserves_nonspatial_dimensions_and_accepts_touching_regions():
    result = source_contrast(np.full((2, 1, 3), 4.0), [-0.5, 0.5], axis=1,
                             edge_fraction=0.25, center_fraction=0.5)
    assert result["center_mean"].shape == (2, 3)
    np.testing.assert_array_equal(result["center_mean"], 4.0)
    np.testing.assert_array_equal(result["edge_minus_center"], 0.0)


@pytest.mark.parametrize("target", [[-0.1, 1.0], [0.0, 0.9], [0.0, np.nextafter(1.0, 0.0)]])
def test_domain_mismatch_never_trims_or_extrapolates(target):
    with pytest.raises(ValueError, match="domain endpoints"):
        conservative_remap([1.0], [0.0, 1.0], target)


@pytest.mark.parametrize("values", [
    [True, 1.0], np.array([True, False]), [1.0, 2.0j], np.array([1.0 + 0j, 2.0 + 0j]),
    ["1", "2"], [np.nan, 1.0], [np.inf, 1.0], [-1.0, 1.0], np.array([1, 2], dtype=object),
])
def test_nonphysical_source_values_are_rejected(values):
    for operation in (matched_uniform_source, source_contrast):
        with pytest.raises(ValueError):
            operation(values, [0.0, 0.5, 1.0])
    with pytest.raises(ValueError):
        conservative_remap(values, [0.0, 0.5, 1.0], [0.0, 1.0])


@pytest.mark.parametrize("edges", [[0.0], [0.0, 0.0], [1.0, 0.0], [0.0, np.inf],
                                    [False, 1.0], [[0.0, 1.0]], [0.0, 1.0 + 0j]])
def test_invalid_edges_are_rejected(edges):
    with pytest.raises(ValueError):
        matched_uniform_source([1.0], edges)


@pytest.mark.parametrize("axis", [True, 0.5, 2, -3])
def test_invalid_spatial_axes_are_rejected(axis):
    with pytest.raises(ValueError):
        matched_uniform_source([[1.0, 2.0]], [0.0, 0.5, 1.0], axis=axis)


def test_shape_and_scalar_sources_are_rejected():
    with pytest.raises(ValueError, match="bin count"):
        matched_uniform_source([1.0, 2.0], [0.0, 1.0])
    with pytest.raises(ValueError, match="existing spatial"):
        matched_uniform_source(1.0, [0.0, 1.0])


@pytest.mark.parametrize("edge,center", [(0.0, .2), (-.1, .2), (.1, 0.0), (.3, .5),
                                         (True, .2), (.1, np.nan), ([.1], .2), (.1, .2j)])
def test_invalid_or_overlapping_contrast_regions_are_rejected(edge, center):
    with pytest.raises(ValueError):
        source_contrast([1.0], [0.0, 1.0], edge_fraction=edge, center_fraction=center)


def test_unknown_signed_moments_and_unsigned_negative_energy_are_rejected():
    with pytest.raises(ValueError, match="unknown moment"):
        matched_uniform_moments({"number": [1.0]}, [0, 1], allow_signed={"momentum"})
    with pytest.raises(ValueError, match="collection"):
        matched_uniform_moments({"number": [1.0]}, [0, 1], allow_signed="number")
    with pytest.raises(ValueError, match="nonnegative"):
        matched_uniform_moments({"number": [1.0], "energy": [-1.0]}, [0, 1])
    with pytest.raises(ValueError, match="nonempty mapping"):
        matched_uniform_moments({}, [0, 1])
