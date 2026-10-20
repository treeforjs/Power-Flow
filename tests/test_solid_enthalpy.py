"""Independent checks of the qualified Mills austenitic-solid input contract.

Numerical references below use the printed 2004 equations, not implementation
helpers.  The selected branch is an evaluated mass-specific Cp/H relation; the
fixtures do not establish the dissertation foil's composition or phase.
"""

from dataclasses import FrozenInstanceError
from decimal import Decimal, localcontext
import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from mhdlab.solid_enthalpy import load_mills_2004_solid_enthalpy


MANIFEST = Path(__file__).resolve().parents[1] / "data/materials/mills_2004_source_manifest.json"
PDF_SHA256 = "bf120dcbdf20d5f9b00494efd5e05643793edc0838deee8ba2c52988f9074836"


def _parameters(**overrides):
    parameters = dict(
        temperature_min_k=300.0,
        temperature_max_k=1700.0,
        reference_temperature_k=800.0,
        reference_convention="printed_eq11",
        expected_material_family="austenitic stainless steel",
        expected_phase="solid",
        material_qualification="Selected evaluated 304/316-based austenitic correlation; specimen chemistry remains unmeasured.",
        phase_qualification="Caller restricts this candidate to solid austenite over the explicitly selected interval.",
    )
    parameters.update(overrides)
    return parameters


def model(**overrides):
    return load_mills_2004_solid_enthalpy(MANIFEST, **_parameters(**overrides))


def _decimal_temperature(temperature):
    # Keep the exact binary input, including a one-ULP displacement from Tref.
    return Decimal.from_float(float(temperature))


def _printed_cp(temperature):
    with localcontext() as context:
        context.prec = 70
        t = _decimal_temperature(temperature)
        return Decimal(472) + Decimal("0.136") * t - Decimal("2.82e6") / t**2


def _printed_h(temperature):
    with localcontext() as context:
        context.prec = 70
        t = _decimal_temperature(temperature)
        return (Decimal(472) * t + Decimal("0.068") * t**2
                + Decimal("2.82e6") / t - Decimal(156000))


def _delta_h(temperature, reference):
    with localcontext() as context:
        context.prec = 70
        return _printed_h(temperature) - _printed_h(reference)


def test_printed_cp_enthalpy_and_reference_difference_against_decimal():
    relation = model()
    temperatures = np.array([[300., 350., 500., 800.], [1000., 1250., 1600., 1700.]])
    values = relation.evaluate(temperatures)
    expected_cp = np.array([float(_printed_cp(t)) for t in temperatures.flat]).reshape(temperatures.shape)
    expected_h = np.array([float(_printed_h(t)) for t in temperatures.flat]).reshape(temperatures.shape)
    expected_delta = np.array([float(_delta_h(t, 800.)) for t in temperatures.flat]).reshape(temperatures.shape)
    np.testing.assert_allclose(values.cp_j_kg_k, expected_cp, rtol=4e-15, atol=0.)
    np.testing.assert_allclose(values.printed_enthalpy_j_kg, expected_h, rtol=4e-15, atol=0.)
    np.testing.assert_allclose(values.enthalpy_j_kg, expected_h, rtol=4e-15, atol=0.)
    np.testing.assert_allclose(values.delta_enthalpy_from_reference_j_kg, expected_delta, rtol=4e-15, atol=0.)
    assert np.all(values.cp_j_kg_k > 0)


@pytest.mark.parametrize("convention", ["printed_eq11", "rezeroed_eq11"])
def test_enthalpy_increment_is_independent_integral_of_source_cp(convention):
    relation = model(reference_temperature_k=375., reference_convention=convention)
    for temperature in (300., 375., 550., 1000., 1700.):
        # Independent adaptive integration of printed Cp verifies the sign and
        # inverse-square contribution without calling the implementation Cp.
        expected = quad(lambda t: 472. + .136 * t - 2.82e6 / t**2,
                        375., temperature, epsabs=1e-8, epsrel=3e-13)[0]
        assert relation.delta_enthalpy_from_reference_j_kg(temperature) == pytest.approx(expected, rel=8e-14, abs=1e-9)
    for temperature in (400., 900., 1600.):
        step = .02
        derivative = (relation.enthalpy_j_kg(temperature + step)
                      - relation.enthalpy_j_kg(temperature - step)) / (2 * step)
        assert derivative == pytest.approx(float(_printed_cp(temperature)), rel=3e-9)


def test_rounded_source_offset_is_preserved_until_explicit_rezeroing():
    printed = model(reference_temperature_k=300., reference_convention="printed_eq11")
    rezeroed = model(reference_temperature_k=300., reference_convention="rezeroed_eq11")
    assert printed.printed_enthalpy_j_kg(300.) == pytest.approx(1120., rel=2e-14)
    assert printed.enthalpy_j_kg(300.) == pytest.approx(1120., rel=2e-14)
    assert rezeroed.enthalpy_j_kg(300.) == 0.
    assert rezeroed.printed_enthalpy_j_kg(300.) == printed.printed_enthalpy_j_kg(300.)
    assert rezeroed.printed_enthalpy_at_reference_j_kg == printed.printed_enthalpy_at_reference_j_kg
    # The rounded printed source is not secretly corrected to zero at 298 K or
    # to Table 3's separately rounded 298.2 K row.
    manifest = json.loads(MANIFEST.read_text())
    assert float(manifest["reference_sanity"]["printed_eq11_at_298_j_kg"]) == pytest.approx(float(_printed_h(298.)), rel=2e-15)
    assert float(manifest["reference_sanity"]["printed_eq11_at_298_2_j_kg"]) == pytest.approx(float(_printed_h(298.2)), rel=2e-13)
    assert float(_printed_h(298.)) > 150.
    assert float(_printed_h(298.2)) > 250.
    temperatures = np.array([300., 800., 1700.])
    np.testing.assert_allclose(printed.enthalpy_j_kg(temperatures) - rezeroed.enthalpy_j_kg(temperatures),
                               1120., rtol=0., atol=2e-10)
    np.testing.assert_array_equal(printed.cp_j_kg_k(temperatures), rezeroed.cp_j_kg_k(temperatures))


@pytest.mark.parametrize("reference", [300., 800., 1699.])
def test_factored_tiny_reference_deltas_retain_both_signs(reference):
    relation = model(reference_temperature_k=reference, reference_convention="rezeroed_eq11")
    temperatures = [np.nextafter(reference, np.inf)]
    if reference > 300.:
        temperatures.append(np.nextafter(reference, -np.inf))
    for temperature in temperatures:
        expected = float(_delta_h(temperature, reference))
        actual = relation.delta_enthalpy_from_reference_j_kg(temperature)
        assert actual != 0.
        assert actual == pytest.approx(expected, rel=5e-15, abs=0.)
        assert relation.enthalpy_j_kg(temperature) == actual
        assert np.sign(actual) == np.sign(temperature - reference)
    assert relation.delta_enthalpy_from_reference_j_kg(reference) == 0.


@pytest.mark.parametrize("convention", ["printed_eq11", "rezeroed_eq11"])
def test_inverse_is_bounded_monotone_and_matches_independent_source_values(convention):
    relation = model(reference_convention=convention)
    temperatures = np.linspace(300., 1700., 153).reshape(9, 17)
    offset = _printed_h(800.) if convention == "rezeroed_eq11" else Decimal(0)
    enthalpies = np.array([float(_printed_h(t) - offset) for t in temperatures.flat]).reshape(temperatures.shape)
    actual = relation.temperature_from_enthalpy_k(enthalpies)
    np.testing.assert_allclose(actual, temperatures, rtol=0., atol=4e-10)
    assert np.all(np.diff(actual.ravel()) > 0)
    assert np.all(actual >= 300.) and np.all(actual <= 1700.)
    lower, upper = relation.supported_enthalpy_range_j_kg
    assert relation.temperature_from_enthalpy_k(lower) == 300.
    assert relation.temperature_from_enthalpy_k(upper) == 1700.
    for outside in (np.nextafter(lower, -np.inf), np.nextafter(upper, np.inf)):
        with pytest.raises(ValueError):
            relation.temperature_from_enthalpy_k(outside)


def test_rezeroed_inverse_supports_negative_enthalpy_and_resolves_tiny_offsets():
    relation = model(reference_convention="rezeroed_eq11")
    assert relation.enthalpy_j_kg(300.) < 0
    assert relation.temperature_from_enthalpy_k(0.) == 800.
    for temperature in (np.nextafter(800., -np.inf), np.nextafter(800., np.inf)):
        enthalpy = float(_delta_h(temperature, 800.))
        assert relation.temperature_from_enthalpy_k(enthalpy) == pytest.approx(temperature, rel=0., abs=np.spacing(800.))


def test_arrays_are_detached_immutable_and_scalars_remain_scalars():
    relation = model()
    temperatures = np.array([300., 800., 1700.])
    values = relation.evaluate(temperatures)
    temperatures[0] = 900.
    for name in ("temperature_k", "cp_j_kg_k", "printed_enthalpy_j_kg",
                 "delta_enthalpy_from_reference_j_kg", "enthalpy_j_kg"):
        array = getattr(values, name)
        assert isinstance(array, np.ndarray) and array.shape == (3,)
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array[0] = 0.
        with pytest.raises(ValueError):
            array.flags.writeable = True
    assert values.temperature_k[0] == 300.
    with pytest.raises((FrozenInstanceError, AttributeError)):
        values.enthalpy_j_kg = np.ones(3)
    for operation in (relation.cp_j_kg_k, relation.printed_enthalpy_j_kg,
                      relation.delta_enthalpy_from_reference_j_kg, relation.enthalpy_j_kg):
        assert isinstance(operation(500.), float)
    assert isinstance(relation.temperature_from_enthalpy_k(relation.enthalpy_j_kg(500.)), float)


def test_nondefault_narrower_domain_and_reference_are_not_replaced():
    relation = model(temperature_min_k=450., temperature_max_k=1200.,
                     reference_temperature_k=1000., reference_convention="rezeroed_eq11")
    assert relation.supported_temperature_range_k == (450., 1200.)
    assert relation.enthalpy_j_kg(1000.) == 0.
    for temperature in (449.999, 1200.001):
        with pytest.raises(ValueError):
            relation.evaluate(temperature)


@pytest.mark.parametrize("name", list(_parameters()))
def test_all_physical_selection_keywords_are_required(name):
    parameters = _parameters()
    del parameters[name]
    with pytest.raises(TypeError):
        load_mills_2004_solid_enthalpy(MANIFEST, **parameters)
    assert inspect.signature(load_mills_2004_solid_enthalpy).parameters[name].default is inspect.Parameter.empty


@pytest.mark.parametrize("overrides", [
    {"temperature_min_k": 299.999}, {"temperature_max_k": 1700.001},
    {"temperature_min_k": 1700.}, {"temperature_min_k": 800., "temperature_max_k": 700.},
    {"reference_temperature_k": 298.}, {"reference_temperature_k": 1700.001},
    {"temperature_min_k": True}, {"temperature_max_k": np.inf},
    {"reference_temperature_k": np.nan}, {"reference_temperature_k": 800. + 0j},
    {"temperature_min_k": np.ma.array(300., mask=False)},
    {"temperature_max_k": "1700"}, {"reference_temperature_k": np.bool_(True)},
    {"reference_convention": "zero_at_298"}, {"reference_convention": ""},
    {"expected_material_family": "ferritic stainless steel"},
    {"expected_material_family": "304"}, {"expected_phase": "liquid"},
    {"expected_phase": "solid-liquid mixture"}, {"material_qualification": ""},
    {"material_qualification": "   "}, {"phase_qualification": ""},
])
def test_unsupported_source_phase_reference_and_missing_qualification_fail(overrides):
    with pytest.raises((ValueError, TypeError)):
        model(**overrides)


INVALID_VALUES = [True, np.bool_(False), np.nan, np.inf, -np.inf, 500. + 0j,
                  [500., True], [500., "600"], np.array([500. + 0j]),
                  np.ma.array([500.], mask=[False]), [], np.array([])]


@pytest.mark.parametrize("value", INVALID_VALUES)
def test_temperature_inputs_reject_bool_mask_complex_nonfinite_or_empty(value):
    relation = model()
    for operation in (relation.cp_j_kg_k, relation.printed_enthalpy_j_kg,
                      relation.delta_enthalpy_from_reference_j_kg,
                      relation.enthalpy_j_kg, relation.evaluate):
        with pytest.raises((ValueError, TypeError)):
            operation(value)


@pytest.mark.parametrize("value", INVALID_VALUES)
def test_inverse_rejects_nonreal_or_missing_inputs(value):
    with pytest.raises((ValueError, TypeError)):
        model().temperature_from_enthalpy_k(value)


def test_temperature_endpoints_are_supported_without_extrapolation():
    relation = model()
    relation.evaluate(np.array([300., 1700.]))
    for outside in (np.nextafter(300., -np.inf), np.nextafter(1700., np.inf), 1800., 1900.):
        with pytest.raises(ValueError):
            relation.evaluate(outside)


def test_source_manifest_and_pdf_identity_are_pinned_and_tampering_fails(tmp_path):
    relation = model()
    assert PDF_SHA256 in repr(relation.source_metadata)
    assert relation.source_manifest_sha256 == hashlib.sha256(MANIFEST.read_bytes()).hexdigest()
    assert relation.metadata["source_manifest_sha256"] == relation.source_manifest_sha256
    assert "10.2355/isijinternational.44.1661" in repr(relation.source_metadata)
    copied = tmp_path / "manifest.json"
    copied.write_bytes(MANIFEST.read_bytes())
    identical = load_mills_2004_solid_enthalpy(copied, **_parameters())
    assert identical.enthalpy_j_kg(800.) == relation.enthalpy_j_kg(800.)
    payload = json.loads(copied.read_text())
    payload["source"]["pdf_sha256"] = "0" * 64
    copied.write_text(json.dumps(payload))
    with pytest.raises(ValueError):
        load_mills_2004_solid_enthalpy(copied, **_parameters())


def test_metadata_is_deeply_immutable_or_detached_and_cannot_change_physics():
    relation = model()
    before = repr(relation.source_metadata)
    metadata = relation.metadata
    assert isinstance(metadata, dict)
    original_metadata = relation.metadata
    metadata["source_metadata"]["source"]["doi"] = "forged"
    metadata["source_metadata"]["source"]["authors"][0] = "forged"
    metadata["supported_temperature_range_k"][0] = 0.
    assert relation.metadata == original_metadata
    metadata.clear()
    assert relation.metadata
    assert repr(relation.source_metadata) == before
    with pytest.raises(TypeError):
        relation.source_metadata["forged"] = "value"
    # Check a nested container rather than merely the top mapping's proxy.
    def nested_containers(value):
        if hasattr(value, "items"):
            for item in value.values():
                if hasattr(item, "items") or isinstance(item, (list, tuple)):
                    yield item
                    yield from nested_containers(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                if hasattr(item, "items") or isinstance(item, (list, tuple)):
                    yield item
                    yield from nested_containers(item)
    containers = list(nested_containers(relation.source_metadata))
    assert containers
    for container in containers:
        if hasattr(container, "items"):
            with pytest.raises(TypeError):
                container["forged"] = "value"
        elif len(container):
            with pytest.raises(TypeError):
                container[0] = "forged"
    with pytest.raises((FrozenInstanceError, AttributeError)):
        relation.reference_convention = "rezeroed_eq11"
    assert repr(relation.source_metadata) == before
    assert relation.enthalpy_j_kg(300.) == pytest.approx(1120., rel=2e-14)


def test_contract_exposes_no_cv_density_internal_energy_or_liquid_branch():
    relation = model()
    for unsupported in ("cv_j_kg_k", "volumetric_heat_capacity_j_m3_k", "density_kg_m3",
                        "internal_energy_j_kg", "liquid_enthalpy_j_kg", "latent_heat_j_kg"):
        assert not hasattr(relation, unsupported)
    assert "mass" in repr(relation.metadata).lower()
    assert "enthalpy" in repr(relation.metadata).lower()
    assert relation.metadata["cp_is_cv"] is False
    assert relation.metadata["enthalpy_is_internal_energy"] is False
    assert relation.metadata["actual_dissertation_foil_validated"] is False
    assert relation.metadata["density_or_fusion_supplied"] is False
    assert relation.metadata["adopted_as_solver_default"] is False
    assert relation.parameter_status == "evaluated_source_correlation"
    assert relation.material_qualification == _parameters()["material_qualification"]
    assert relation.phase_qualification == _parameters()["phase_qualification"]
    with pytest.raises(TypeError):
        type(relation)(**_parameters())
