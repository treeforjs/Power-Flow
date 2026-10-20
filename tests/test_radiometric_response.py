"""Independent unit, bin-integration and covariance checks for supplied response."""

from dataclasses import FrozenInstanceError
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
import pytest
from scipy.integrate import quad

from mhdlab.radiometric_response import RadiometricResponse


def response(**overrides):
    parameters = dict(
        bin_edges_nm=np.array([100., 101., 104.]),
        etendue_m2_sr=2e-8,
        effective_exposure_s=3e-9,
        counts_per_collected_photon=np.array([.2, 5.]),
        reference_plane="manufactured spectrometer entrance plane",
        included_factors=("mirror_transmission", "detector_quantum_efficiency", "camera_gain_at_matching_gate_condition"),
        status="manufactured",
        provenance="manufactured response; no instrument calibration is inferred",
        measurement_conditions_id="manufactured-shot-gate-and-gain-setting-1",
        response_includes_exposure=False,
        response_includes_etendue=False,
        response_includes_bin_width=False,
    )
    parameters.update(overrides)
    return RadiometricResponse(**parameters)


def options(kind="photon_spectral_radiance"):
    return dict(radiance_kind=kind, measurement_conditions_id="manufactured-shot-gate-and-gain-setting-1")


def test_unequal_bin_photon_counts_include_width_etendue_exposure_and_supplied_gain_once():
    calibration = response()
    # Direct dimensional arithmetic, with gain >1 in the wider second bin.
    photons = np.array([2e20, 4e19])
    actual = calibration.counts_from_radiance(photons, **options())
    np.testing.assert_allclose(actual, [2400., 36000.], rtol=3e-15, atol=0.)
    factors = calibration.conversion_factors("photon_spectral_radiance")
    np.testing.assert_allclose(factors, np.array([1., 3.]) * 2e-8 * 3e-9 * [.2, 5.], rtol=3e-15)
    # Neither mirror nor gain labels trigger any additional multiplication.
    renamed = response(included_factors=("supplied_total_photon_response",))
    np.testing.assert_array_equal(renamed.counts_from_radiance(photons, **options()), actual)


def test_energy_radiance_integrates_lambda_over_hc_using_exact_si_values():
    calibration = response()
    energy = np.array([.001, .002])
    actual = calibration.counts_from_radiance(energy, **options("energy_spectral_radiance"))
    expected = []
    with localcontext() as context:
        context.prec = 65
        h, c = Decimal("6.62607015e-34"), Decimal(299792458)
        for lower, upper, gain, radiance in zip([100, 101], [101, 104], ["0.2", "5"], ["0.001", "0.002"]):
            integral = (Decimal(upper)**2 - Decimal(lower)**2) / 2 * Decimal("1e-9")
            expected.append(float(Decimal(radiance) * Decimal("2e-8") * Decimal("3e-9")
                                  * Decimal(gain) * integral / (h * c)))
    np.testing.assert_allclose(actual, expected, rtol=4e-15, atol=0.)
    numerical = [quad(lambda wavelength: value * wavelength * 1e-9 / (6.62607015e-34 * 299792458.)
                      * 2e-8 * 3e-9 * gain, lower, upper, epsabs=0., epsrel=1e-12)[0]
                 for lower, upper, gain, value in zip([100, 101], [101, 104], [.2, 5.], energy)]
    np.testing.assert_allclose(actual, numerical, rtol=4e-15)
    photon_bin_means = energy * np.array([100.5, 102.5]) * 1e-9 / (6.62607015e-34 * 299792458.)
    np.testing.assert_allclose(actual, calibration.counts_from_radiance(photon_bin_means, **options()), rtol=4e-15)


@pytest.mark.parametrize("kind", ["photon_spectral_radiance", "energy_spectral_radiance"])
def test_signed_counts_and_correlated_full_covariance_propagate_without_diagonalization(kind):
    calibration = response()
    counts = np.array([-12., 30.])
    covariance = np.array([[4., -1.5], [-1.5, 9.]])
    result = calibration.calibrate_counts(counts, covariance, **options(kind))
    factors = calibration.conversion_factors(kind)
    inverse = np.diag(1 / factors)
    np.testing.assert_allclose(result.radiance, inverse @ counts, rtol=4e-15)
    np.testing.assert_allclose(result.covariance, inverse @ covariance @ inverse, rtol=4e-15)
    assert result.radiance[0] < 0 and result.covariance[0, 1] < 0
    before_correlation = covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])
    after_correlation = result.covariance[0, 1] / np.sqrt(result.covariance[0, 0] * result.covariance[1, 1])
    assert after_correlation == pytest.approx(before_correlation, rel=3e-15)
    np.testing.assert_allclose(calibration.counts_from_radiance(result.radiance, **options(kind)), counts, rtol=4e-15)
    np.testing.assert_allclose(np.diag(factors) @ result.covariance @ np.diag(factors), covariance, rtol=4e-15)
    assert result.metadata["system_calibration_uncertainty_included"] is False
    assert result.metadata["input_covariance_units"] == "ADU^2"
    assert result.metadata["radiance_kind"] == kind


def test_energy_inverse_and_photon_inverse_have_consistent_bin_mean_units():
    calibration = response()
    counts, covariance = np.array([12., -30.]), np.array([[4., 1.5], [1.5, 9.]])
    photons = calibration.calibrate_counts(counts, covariance, **options())
    energy = calibration.calibrate_counts(counts, covariance, **options("energy_spectral_radiance"))
    photon_energy_conversion = np.array([100.5, 102.5]) * 1e-9 / (6.62607015e-34 * 299792458.)
    np.testing.assert_allclose(photons.radiance, energy.radiance * photon_energy_conversion, rtol=4e-15)
    np.testing.assert_allclose(photons.covariance, energy.covariance * np.outer(photon_energy_conversion, photon_energy_conversion), rtol=4e-15)
    assert photons.radiance_units == "photons s^-1 m^-2 sr^-1 nm^-1"
    assert energy.radiance_units == "W m^-2 sr^-1 nm^-1"


def test_zero_signal_and_singular_covariance_are_supported_without_poisson_or_jitter():
    calibration = response()
    zero = calibration.calibrate_counts([0., 0.], np.zeros((2, 2)), **options())
    np.testing.assert_array_equal(zero.radiance, 0.)
    np.testing.assert_array_equal(zero.covariance, 0.)
    singular = np.array([[4., -2.], [-2., 1.]])
    actual = calibration.calibrate_counts([-100., 100.], singular, **options())
    factors = calibration.conversion_factors("photon_spectral_radiance")
    np.testing.assert_allclose(actual.covariance, singular / np.outer(factors, factors), rtol=3e-15)
    qualification = actual.metadata["input_covariance_qualification"]
    assert qualification["singular_covariance_supported"] is True
    assert qualification["jitter_or_projection_applied"] is False
    assert qualification["psd_relative_eigenvalue_tolerance"] == 1e-12


def test_declared_psd_roundoff_tolerance_is_not_an_unreported_projection():
    calibration = response()
    covariance = np.array([[1., 1. + 5e-13], [1. + 5e-13, 1.]])
    actual = calibration.calibrate_counts([0., 0.], covariance, **options())
    factors = calibration.conversion_factors("photon_spectral_radiance")
    np.testing.assert_allclose(actual.covariance * np.outer(factors, factors), covariance, rtol=3e-15)
    assert actual.metadata["input_covariance_qualification"]["minimum_correlation_eigenvalue"] < 0
    assert actual.metadata["input_covariance_qualification"]["jitter_or_projection_applied"] is False


def test_widely_separated_variances_preserve_correlation_and_unit_scaling():
    calibration = response()
    covariance = np.array([[1e300, .5], [.5, 1e-300]])
    # Use energy factors to keep the converted covariance within finite range.
    actual = calibration.calibrate_counts([0., 0.], covariance, **options("energy_spectral_radiance"))
    factors = calibration.conversion_factors("energy_spectral_radiance")
    np.testing.assert_allclose(actual.covariance, covariance / np.outer(factors, factors), rtol=4e-15)
    assert actual.covariance[0, 1] > 0


def test_response_condition_and_integration_metadata_prevent_duplicate_time_geometry_or_bins():
    calibration = response()
    for name in ("response_includes_exposure", "response_includes_etendue", "response_includes_bin_width"):
        with pytest.raises(ValueError, match="integration"):
            response(**{name: True})
        with pytest.raises(ValueError):
            response(**{name: 0})
    for factor in ("effective_exposure_s", "etendue_m2_sr", "wavelength_bin_width"):
        with pytest.raises(ValueError, match="integration"):
            response(included_factors=(factor,))
    # Matching-gate gain is a response condition, not a second time integral.
    assert "camera_gain_at_matching_gate_condition" in calibration.included_factors
    for method, args in ((calibration.counts_from_radiance, ([1., 2.],)),
                         (calibration.calibrate_counts, ([1., 2.], np.eye(2)))):
        with pytest.raises(ValueError, match="conditions"):
            method(*args, radiance_kind="photon_spectral_radiance", measurement_conditions_id="another-gate-setting")
    assert calibration.metadata["automatic_throughput_or_gate_multiplication"] is False


@pytest.mark.parametrize("kind", ["µW/(sr nm mm)", "μW/(sr nm mm)", "uW/(sr nm mm)",
                                 "µW/(sr nm mm²)", "W/(m² sr nm)", "photon", "energy", "", True])
def test_ambiguous_thesis_ordinate_and_unknown_radiance_kinds_are_not_area_units(kind):
    calibration = response()
    with pytest.raises(ValueError, match="radiance_kind"):
        calibration.counts_from_radiance([1., 2.], **options(kind))
    with pytest.raises(ValueError, match="radiance_kind"):
        calibration.calibrate_counts([1., 2.], np.eye(2), **options(kind))


@pytest.mark.parametrize("overrides", [
    {"bin_edges_nm": []}, {"bin_edges_nm": [100.]}, {"bin_edges_nm": [100., 100., 104.]},
    {"bin_edges_nm": [101., 100., 104.]}, {"bin_edges_nm": [0., 101., 104.]},
    {"bin_edges_nm": [100., np.inf, 104.]}, {"bin_edges_nm": [100., True, 104.]},
    {"bin_edges_nm": np.ma.array([100., 101., 104.], mask=False)},
    {"counts_per_collected_photon": [.2]}, {"counts_per_collected_photon": [0., 5.]},
    {"counts_per_collected_photon": [-.2, 5.]}, {"counts_per_collected_photon": [.2, np.nan]},
    {"counts_per_collected_photon": [.2, True]}, {"counts_per_collected_photon": [.2, 5. + 0j]},
    {"etendue_m2_sr": 0.}, {"etendue_m2_sr": True}, {"etendue_m2_sr": -1.},
    {"effective_exposure_s": 0.}, {"effective_exposure_s": np.inf},
    {"effective_exposure_s": np.ma.array(1., mask=False)}, {"effective_exposure_s": "1"},
    {"reference_plane": ""}, {"provenance": ""}, {"status": ""},
    {"measurement_conditions_id": ""}, {"included_factors": []},
    {"included_factors": ("detector_gain", "detector_gain")},
])
def test_invalid_response_inputs_fail_explicitly(overrides):
    with pytest.raises(ValueError):
        response(**overrides)


@pytest.mark.parametrize("values", [[1.], [[1., 2.]], [1., True], [1., "2"], [1., 2. + 0j],
                                   [1., np.nan], [1., np.inf], [], np.ma.array([1., 2.], mask=False)])
def test_signal_vectors_reject_missing_wrong_shape_mask_bool_complex_or_nonfinite(values):
    calibration = response()
    with pytest.raises(ValueError):
        calibration.counts_from_radiance(values, **options())
    with pytest.raises(ValueError):
        calibration.calibrate_counts(values, np.eye(2), **options())


@pytest.mark.parametrize("covariance", [
    [1., 1.], np.eye(3), [[1., .1], [.2, 1.]], [[-1., 0.], [0., 1.]],
    [[0., .1], [.1, 1.]], [[1., 2.], [2., 1.]], [[1., True], [True, 1.]],
    [[1., 0j], [0j, 1.]], [[1., np.nan], [np.nan, 1.]],
    np.ma.array(np.eye(2), mask=False),
])
def test_covariance_shape_symmetry_psd_and_real_guards(covariance):
    with pytest.raises(ValueError):
        response().calibrate_counts([1., -1.], covariance, **options())


@pytest.mark.parametrize("overrides", [
    {"etendue_m2_sr": 1e-200, "effective_exposure_s": 1e-200},
    {"etendue_m2_sr": 1e200, "effective_exposure_s": 1e200},
])
def test_unrepresentable_nonzero_conversion_factors_fail_without_zero_response(overrides):
    with pytest.raises(ValueError, match="numeric support"):
        response(**overrides)


def test_final_unit_conversion_avoids_intermediate_outer_factor_overflow():
    calibration = response(bin_edges_nm=[100., 101.], etendue_m2_sr=1e-300,
                           effective_exposure_s=1., counts_per_collected_photon=[1.])
    actual = calibration.calibrate_counts([0.], [[1e-300]], **options())
    factor = calibration.conversion_factors("photon_spectral_radiance")[0]
    exact = Fraction(1e-300) / Fraction(float(factor))**2
    assert actual.covariance[0, 0] == float(exact)
    assert np.isfinite(actual.covariance[0, 0])
    np.testing.assert_allclose(calibration.counts_from_radiance([1e300], **options()), [1.], rtol=3e-15)


def test_nonzero_covariance_underflow_is_rejected_instead_of_discarding_information():
    calibration = response(bin_edges_nm=[100., 101.], etendue_m2_sr=1e200,
                           effective_exposure_s=1., counts_per_collected_photon=[1.])
    with pytest.raises(ValueError, match="underflows"):
        calibration.calibrate_counts([0.], [[1e-300]], **options())


def test_arrays_and_nested_metadata_are_immutable_or_detached():
    edges, gain = np.array([100., 101., 104.]), np.array([.2, 5.])
    calibration = response(bin_edges_nm=edges, counts_per_collected_photon=gain)
    edges[0], gain[0] = 1., 100.
    assert calibration.bin_edges_nm[0] == 100.
    assert calibration.counts_per_collected_photon[0] == .2
    counts = np.array([-2., 3.])
    covariance = np.array([[1., -.1], [-.1, 2.]])
    result = calibration.calibrate_counts(counts, covariance, **options())
    counts[:] = 0.
    covariance[:] = 0.
    assert result.radiance[0] < 0 and result.covariance[0, 1] < 0
    for array in (calibration.bin_edges_nm, calibration.counts_per_collected_photon,
                  calibration.conversion_factors("photon_spectral_radiance"),
                  calibration.counts_from_radiance([1., 2.], **options()),
                  result.bin_edges_nm, result.radiance, result.covariance):
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.flat[0] = 0.
        with pytest.raises(ValueError):
            array.flags.writeable = True
    metadata = result.metadata
    metadata["input_covariance_qualification"]["jitter_or_projection_applied"] = True
    metadata["bin_edges_nm"][0] = 0.
    assert result.metadata["input_covariance_qualification"]["jitter_or_projection_applied"] is False
    assert result.metadata["bin_edges_nm"][0] == 100.
    separate = calibration.metadata
    separate["included_factors"].clear()
    assert calibration.metadata["included_factors"]
    with pytest.raises((FrozenInstanceError, AttributeError)):
        calibration.etendue_m2_sr = 1.


def test_result_cannot_bypass_calibration_and_covariance_validation_by_direct_construction():
    result = response().calibrate_counts([0., 0.], np.eye(2), **options())
    with pytest.raises(TypeError, match="calibrate_counts"):
        type(result)(bin_edges_nm=[True, False], radiance=[True], covariance=[[True]],
                     radiance_kind="ambiguous", radiance_units="ambiguous", _metadata={})
