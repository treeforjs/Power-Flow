"""Independent parent/state/pool inventory and coverage-sweep checks."""

from dataclasses import FrozenInstanceError, replace
import json
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mhdlab.constants import KB, QE
from mhdlab.coverage_inventory import (
    InventoryPoolAllocation,
    MonolayerDefinition,
    ParentSurfaceInventory,
    build_coverage_source_sweep,
    prepare_coverage_inventory,
)
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_release import TemkinReservoir


def reservoir(name="pool-a", state_id="M:v0", **changes):
    arguments = dict(name=name, state_id=state_id, attempt_frequency_s=1e7,
                     activation_energy_ev=.08, desorption_energy_ev=.05,
                     internal_energy_ev=.02, coverage_coefficient=.2,
                     reference_coverage_m2=2e18, max_coverage_m2=2e18,
                     kinetic_provenance="Manufactured fixed independent-pool parameters",
                     energy_provenance="Manufactured energy inputs independent of the layer convention")
    arguments.update(changes)
    return TemkinReservoir(**arguments)


def monolayer(density=1e17, **changes):
    arguments = dict(areal_density_m2=density,
                     definition="Explicit manufactured parent-particle layer convention",
                     provenance="Manufactured layer density, separate from each kinetic reference coverage",
                     status="manufactured")
    arguments.update(changes)
    return MonolayerDefinition(**arguments)


def allocation(pool=None, fraction=1., **changes):
    arguments = dict(reservoir=reservoir() if pool is None else pool, fraction=fraction,
                     provenance="Explicit manufactured state/pool population fraction", status="manufactured")
    arguments.update(changes)
    return InventoryPoolAllocation(**arguments)


def parent(*, allocations=None, patch_multiplier=1., layer=None, **changes):
    arguments = dict(parent_id="M", monolayer=monolayer() if layer is None else layer,
                     allocations=[allocation()] if allocations is None else allocations,
                     patch_multiplier=patch_multiplier,
                     profile_provenance="Explicit manufactured patch multipliers, without normalization",
                     status="manufactured", provenance="Manufactured one-parent inventory partition")
    arguments.update(changes)
    return ParentSurfaceInventory(**arguments)


def surface(*, x_edges=(-2., -1., 2.), z_edges=(0., 2., 5.), temperatures=None):
    base = FlashSurfaceHistory([0., 1e-7], [0., 1.], [[300.], [300.]], "synthetic",
                               "coverage-inventory-test", "Manufactured prescribed temperature patches",
                               "Explicit manufactured stationary boundary")
    attributes = base.attributes
    attributes["inventory_test_marker"] = {"meaning": "retained upstream temperature provenance"}
    if temperatures is None:
        temperatures = np.full((2, len(z_edges)-1, len(x_edges)-1), 600.)
    return PlanarSurfaceHistory([0., 1e-7], x_edges, z_edges, temperatures, attributes)


def heterogeneous_parent(*, constant_hazard=False):
    first = np.array([[.1, .2], [.3, .4]])
    second = np.array([[.6, .2], [.1, .2]])
    third = np.array([[.3, .6], [.6, .4]])
    fractions = np.array([first, second, third])
    overrides = dict(activation_energy_ev=0.) if constant_hazard else {}
    pools = [reservoir("pool-a", "M:v0", **overrides),
             reservoir("pool-b", "M:v0", desorption_energy_ev=.1, **overrides),
             reservoir("pool-c", "M:v1", internal_energy_ev=.07, **overrides)]
    multipliers = np.array([[1., 2.], [4., .5]])
    sample = parent(allocations=[allocation(pool, fraction) for pool, fraction in zip(pools, fractions)],
                    patch_multiplier=multipliers)
    return sample, pools, fractions, multipliers


def test_inventory_is_partitioned_once_across_states_and_pools_using_actual_patch_areas():
    sample, pools, fractions, multipliers = heterogeneous_parent()
    coverage = np.array([0., .5, 2.])
    plan = prepare_coverage_inventory(surface(), sample, coverage)
    expected_parent = coverage[:, None, None] * 1e17 * multipliers
    expected_pools = expected_parent[:, None] * fractions[None]
    np.testing.assert_allclose(plan.initial_parent_coverage_m2, expected_parent, rtol=3e-15, atol=0.)
    np.testing.assert_allclose(plan.initial_pool_coverage_m2, expected_pools, rtol=3e-15, atol=0.)
    np.testing.assert_allclose(np.sum(plan.initial_pool_coverage_m2, axis=1), expected_parent,
                               rtol=3e-15, atol=0.)
    # Area=[[2,6],[3,9]]. The unnormalized multiplier has area integral30.5,
    # and pool-specific integrals8,5.7,16.8. Total face area is20.
    np.testing.assert_allclose(plan.parent_number, coverage * 1e17 * 30.5, rtol=3e-15, atol=0.)
    np.testing.assert_allclose(plan.pool_number, coverage[:, None] * 1e17 * np.array([8., 5.7, 16.8]),
                               rtol=4e-15, atol=0.)
    np.testing.assert_allclose(np.sum(plan.pool_number, axis=1), plan.parent_number, rtol=4e-15, atol=0.)
    np.testing.assert_array_equal(plan.patch_multiplier, multipliers)
    np.testing.assert_array_equal(plan.allocation_fractions, fractions)
    assert pools[0].state_id == pools[1].state_id
    assert pools[0].name != pools[1].name
    assert plan.parent_number[-1] != coverage[-1] * 1e17 * 20.
    json.dumps(plan.manifest(), allow_nan=False)


def test_monolayer_conversion_never_substitutes_the_kinetic_reference_coverage():
    pool = reservoir(reference_coverage_m2=7e18, max_coverage_m2=2e18)
    sample = parent(allocations=[allocation(pool)], layer=monolayer(3e17))
    plan = prepare_coverage_inventory(surface(), sample, [1., 2.])
    np.testing.assert_array_equal(plan.initial_pool_coverage_m2[:, 0, 0, 0], [3e17, 6e17])
    assert pool.reference_coverage_m2 == 7e18
    assert pool.max_coverage_m2 == 2e18
    other = parent(allocations=[allocation(pool)], layer=monolayer(6e17))
    changed = prepare_coverage_inventory(surface(), other, [1., 2.])
    np.testing.assert_array_equal(changed.initial_pool_coverage_m2, 2 * plan.initial_pool_coverage_m2)
    assert pool.reference_coverage_m2 == 7e18


def test_source_sweep_preserves_local_mixture_and_independent_analytic_depletion():
    sample, pools, fractions, multipliers = heterogeneous_parent(constant_hazard=True)
    temperatures = np.array([[[300., 600.], [900., 1200.]]] * 2)
    result = build_coverage_source_sweep(surface(temperatures=temperatures), sample, [0., 1., 2.],
                                         parameter_provenance="Manufactured analytic constant-hazard verification",
                                         subdivisions_per_interval=2)
    areas = np.array([[2., 6.], [3., 9.]])
    emitted_fraction = -math.expm1(-1.)
    for coverage, history in zip([0., 1., 2.], result.histories):
        initial = coverage * 1e17 * multipliers[None] * fractions
        emitted = initial * emitted_fraction
        np.testing.assert_allclose(history.coverage_m2[0], initial, rtol=3e-15, atol=0.)
        np.testing.assert_allclose(np.sum(history.released_m2, axis=0), emitted, rtol=3e-14, atol=0.)
        np.testing.assert_allclose(history.coverage_m2[-1], initial * math.exp(-1.), rtol=3e-14, atol=0.)
        np.testing.assert_allclose(np.sum(history.extensive_totals(), axis=0),
                                   np.sum(emitted * areas, axis=(-2, -1)), rtol=3e-14, atol=0.)
        expected_energy = np.sum(emitted * 2 * KB * temperatures[0] * areas, axis=(-2, -1))
        np.testing.assert_allclose(np.sum(history.extensive_totals("translational_energy_j_m2"), axis=0),
                                   expected_energy, rtol=3e-14, atol=0.)
        assert all(actual is supplied for actual, supplied in zip(history.reservoirs, pools))
        assert [item.reference_coverage_m2 for item in history.reservoirs] == [2e18] * 3
        assert history.surface_attributes["inventory_test_marker"]["meaning"] == "retained upstream temperature provenance"
    np.testing.assert_array_equal(result.histories[2].released_m2, 2 * result.histories[1].released_m2)
    json.dumps(result.manifest(), allow_nan=False)


def test_nonlinear_sweep_uses_one_fixed_reference_coverage_against_an_independent_ode():
    pool = reservoir(activation_energy_ev=.1, coverage_coefficient=.8,
                     reference_coverage_m2=1e18, max_coverage_m2=1e18)
    sample = parent(layer=monolayer(1e18), allocations=[allocation(pool)])
    result = build_coverage_source_sweep(surface(), sample, [.25, .75],
                                         parameter_provenance="Manufactured fixed-reference ODE comparison")
    beta = .1 * QE / (KB * 600.)
    emitted_fractions = []
    for initial_fraction, history in zip([.25, .75], result.histories):
        def rhs(_unit_time, state):
            return -state * np.exp(-beta * (1. - .8 * state))

        numerical = solve_ivp(rhs, (0., 1.), [initial_fraction], method="DOP853",
                              rtol=1e-12, atol=1e-14)
        assert numerical.success
        expected_remaining = numerical.y[0, -1] * 1e18
        np.testing.assert_allclose(history.coverage_m2[-1, 0], expected_remaining,
                                   rtol=4e-11, atol=0.)
        assert history.reservoirs[0] is pool
        emitted_fractions.append(float(np.sum(history.released_m2[:, 0, 0, 0])) /
                                 (initial_fraction * 1e18))
    assert emitted_fractions[1] > 1.3 * emitted_fractions[0]


def test_coverage_effect_preflight_keeps_temperature_and_population_on_the_same_patch():
    pool = reservoir(activation_energy_ev=1., coverage_coefficient=1.,
                     reference_coverage_m2=1e18, max_coverage_m2=1e18)
    temperatures = np.array([[[1., 1000.], [1000., 1000.]]] * 2)
    sample = parent(layer=monolayer(1e18), allocations=[allocation(pool)],
                    patch_multiplier=[[0., 1.], [0., 0.]])
    plan = prepare_coverage_inventory(surface(temperatures=temperatures), sample, [1.])
    assert plan.max_dimensionless_coverage_effect[0] == pytest.approx(QE / (KB * 1000.),
                                                                    rel=3e-14, abs=0.)


def test_accepted_fraction_roundoff_is_retained_without_renormalization():
    second_fraction = .5 + 2.**-47
    pools = [reservoir("a", reference_coverage_m2=128., max_coverage_m2=128.),
             reservoir("b", "M:v1", reference_coverage_m2=128., max_coverage_m2=128.)]
    sample = parent(layer=monolayer(16.), allocations=[allocation(pools[0], .5), allocation(pools[1], second_fraction)])
    plan = prepare_coverage_inventory(surface(x_edges=(0., 1.), z_edges=(0., 1.)), sample, [1.])
    assert plan.allocation_fractions[1, 0, 0] == second_fraction
    assert plan.initial_pool_coverage_m2[0, 1, 0, 0] == 16 * second_fraction
    assert float(np.sum(plan.initial_pool_coverage_m2[0, :, 0, 0])) > 16.
    assert np.any(plan.allocation_balance_residual_m2 != 0)


def test_a_zero_parent_inventory_is_a_valid_zero_source_not_a_population_prior():
    pool = reservoir()
    result = build_coverage_source_sweep(surface(), parent(allocations=[allocation(pool)]), [0.],
                                         parameter_provenance="Manufactured zero inventory check")
    assert len(result.histories) == 1
    assert np.all(result.plan.initial_parent_coverage_m2 == 0)
    assert np.all(result.plan.parent_number == 0)
    assert np.all(result.histories[0].coverage_m2 == 0)
    assert np.all(result.histories[0].released_m2 == 0)


def test_nonunit_patch_profile_is_explicitly_not_area_normalized():
    plan = prepare_coverage_inventory(surface(), parent(patch_multiplier=3.), [1.])
    assert np.all(plan.initial_parent_coverage_m2 == 3e17)
    assert plan.parent_number[0] == 60e17


def test_zero_patch_and_zero_pool_fractions_preserve_structural_zeros():
    fractions = np.array([[1., 0.], [0., 1.]])
    pools = [reservoir("a"), reservoir("b", "M:v1")]
    sample = parent(allocations=[allocation(pools[0], fractions), allocation(pools[1], 1.-fractions)],
                    patch_multiplier=[[1., 0.], [2., 0.]])
    plan = prepare_coverage_inventory(surface(), sample, [1.])
    np.testing.assert_array_equal(plan.initial_pool_coverage_m2[0],
                                  [[[1e17, 0.], [0., 0.]], [[0., 0.], [2e17, 0.]]])
    assert np.all(plan.initial_parent_coverage_m2[0, :, 1] == 0)


def test_duplicate_pool_names_reject_even_if_states_and_fractions_differ():
    with pytest.raises(ValueError):
        parent(allocations=[allocation(reservoir("same", "M:v0"), .5),
                            allocation(reservoir("same", "M:v1"), .5)])


def test_multiple_parent_records_are_not_silently_pooled():
    with pytest.raises(ValueError):
        prepare_coverage_inventory(surface(), [parent(), parent(parent_id="Q")], [1.])


def test_evidence_labels_and_exact_state_ids_are_retained_separately():
    layer = monolayer(status="sourced", definition="Explicit supplied parent-count convention")
    sample = parent(layer=layer, status="owner_supplied",
                    allocations=[allocation(reservoir("a", "M:v0"), .25, status="unvalidated"),
                                 allocation(reservoir("b", "M:v1"), .75, status="reconstructed")])
    manifest = prepare_coverage_inventory(surface(), sample, [1.]).manifest()
    assert manifest["parent_status"] == "owner_supplied"
    assert manifest["monolayer"]["status"] == "sourced"
    assert manifest["monolayer"]["definition"] == layer.definition
    assert [item["state_id"] for item in manifest["pools"]] == ["M:v0", "M:v1"]
    assert [item["allocation_status"] for item in manifest["pools"]] == ["unvalidated", "reconstructed"]
    assert [item["reference_coverage_m2"] for item in manifest["pools"]] == [2e18, 2e18]


@pytest.mark.parametrize("changes", [{"definition": ""}, {"provenance": ""},
                                      {"status": "unknown"}, {"status": True}])
def test_monolayer_definition_and_evidence_status_must_be_explicit(changes):
    with pytest.raises(ValueError):
        monolayer(**changes)


@pytest.mark.parametrize("fractions", [[.4, .4], [.6, .6]])
def test_parent_partition_must_sum_to_one_in_every_patch(fractions):
    with pytest.raises(ValueError):
        sample = parent(allocations=[allocation(reservoir("a"), fractions[0]),
                                    allocation(reservoir("b", "M:v1"), fractions[1])])
        prepare_coverage_inventory(surface(), sample, [1.])


def test_fraction_partition_validation_is_local_not_only_area_averaged():
    with pytest.raises(ValueError):
        sample = parent(allocations=[allocation(reservoir("a"), [[.2, .8], [.2, .8]]),
                                    allocation(reservoir("b", "M:v1"), .5)])
        prepare_coverage_inventory(surface(), sample, [1.])


@pytest.mark.parametrize("bad", [True, "1", -1., np.nan, np.inf,
                                 np.ma.array([[1., 1.], [1., 1.]], mask=False),
                                 np.ones((2, 2), dtype=complex), np.ones((1, 2))])
def test_patch_multiplier_rejects_coercion_masks_nonfinite_negative_and_wrong_geometry(bad):
    with pytest.raises(ValueError):
        prepare_coverage_inventory(surface(), parent(patch_multiplier=bad), [1.])


@pytest.mark.parametrize("bad", [True, "1", -.1, 1.1, np.nan, np.inf,
                                 np.ma.array([[1., 1.], [1., 1.]], mask=False),
                                 np.ones((2, 2), dtype=complex), np.ones((1, 2))])
def test_allocation_fraction_rejects_invalid_values_or_geometry(bad):
    with pytest.raises(ValueError):
        prepare_coverage_inventory(surface(), parent(allocations=[allocation(fraction=bad)]), [1.])


@pytest.mark.parametrize("bad", [True, "1e17", 0., -1., np.nan, np.inf])
def test_monolayer_density_is_explicit_positive_finite_and_not_coerced(bad):
    with pytest.raises(ValueError):
        monolayer(bad)


@pytest.mark.parametrize("bad", [[-1., 0.], [1., 1.], [2., 1.], [], [False, True], ["0", "1"],
                                 [0., np.inf], [0., np.nan], np.ma.array([0., 1.], mask=False)])
def test_coverage_sweep_knots_require_sorted_finite_nonnegative_real_values(bad):
    with pytest.raises(ValueError):
        prepare_coverage_inventory(surface(), parent(), bad)


@pytest.mark.parametrize("bad", [True, 0, 1.5, "2"])
def test_time_subdivision_is_an_explicit_positive_integer(bad):
    with pytest.raises(ValueError):
        prepare_coverage_inventory(surface(), parent(), [1.], subdivisions_per_interval=bad)


def test_inventory_above_a_pool_domain_is_rejected_before_any_source_integration(monkeypatch):
    import mhdlab.coverage_inventory as implementation

    def forbidden(*args, **kwargs):
        raise AssertionError("source evolution started before validating every coverage knot")

    monkeypatch.setattr(implementation, "evolve_temkin_planar_release", forbidden)
    with pytest.raises(ValueError):
        build_coverage_source_sweep(surface(), parent(layer=monolayer(1e18)), [0., 1., 3.],
                                     parameter_provenance="Manufactured preflight check")


def test_unsupported_nonlinear_coverage_effect_is_rejected_during_preparation():
    pool = reservoir(activation_energy_ev=1000., coverage_coefficient=.01,
                     reference_coverage_m2=1e18, max_coverage_m2=1e18)
    cold = surface(temperatures=np.ones((2, 2, 2)))
    with pytest.raises(ValueError):
        prepare_coverage_inventory(cold, parent(layer=monolayer(1e18), allocations=[allocation(pool)]), [1.])


@pytest.mark.parametrize("density,coverage,multiplier", [(1e300, [1e10], 1.),
                                                        (1e-300, [1e-100], 1.),
                                                        (1., [1.], 1e-300)])
def test_unrepresentable_positive_inventory_or_patch_counts_are_rejected(density, coverage, multiplier):
    pool = reservoir(activation_energy_ev=0., reference_coverage_m2=1e308, max_coverage_m2=1e308)
    tiny = multiplier == 1e-300
    geometry = surface(x_edges=(0., 1e-100), z_edges=(0., 1e-100)) if tiny else surface()
    with pytest.raises(ValueError):
        prepare_coverage_inventory(geometry, parent(layer=monolayer(density), patch_multiplier=multiplier,
                                                    allocations=[allocation(pool)]), coverage)


@pytest.mark.parametrize("density,coverage,multiplier,expected", [
    (1e-300, 1e-100, 1e300, 1e-100), (1e300, 1e100, 1e-300, 1e100),
])
def test_finite_final_inventory_is_preserved_across_intermediate_float64_range(density, coverage, multiplier, expected):
    pool = reservoir(activation_energy_ev=0., reference_coverage_m2=1e308, max_coverage_m2=1e308)
    plan = prepare_coverage_inventory(surface(), parent(layer=monolayer(density), patch_multiplier=multiplier,
                                                        allocations=[allocation(pool)]), [coverage])
    np.testing.assert_allclose(plan.initial_parent_coverage_m2, expected, rtol=4e-15, atol=0.)
    assert plan.parent_number[0] == pytest.approx(20. * expected, rel=4e-15, abs=0.)


def test_large_density_on_a_tiny_patch_retains_its_finite_extensive_population():
    pool = reservoir(activation_energy_ev=0., reference_coverage_m2=1e308, max_coverage_m2=1e308)
    geometry = surface(x_edges=(0., 1e-200, 1e200), z_edges=(0., 1.))
    sample = parent(layer=monolayer(1.), allocations=[allocation(pool)], patch_multiplier=[[1e300, 0.]])
    plan = prepare_coverage_inventory(geometry, sample, [1.])
    assert plan.parent_number[0] == pytest.approx(1e100, rel=4e-15, abs=0.)
    assert plan.pool_number[0, 0] == pytest.approx(1e100, rel=4e-15, abs=0.)


@pytest.mark.parametrize("field", ["parent_number", "pool_number", "initial_pool_coverage_m2"])
def test_derived_plan_inventory_cannot_be_replaced_with_an_inconsistent_total(field):
    plan = prepare_coverage_inventory(surface(), parent(), [1.])
    with pytest.raises(ValueError):
        replace(plan, **{field: 2 * getattr(plan, field)})


def test_input_arrays_and_metadata_are_detached_and_plan_arrays_are_immutable():
    multiplier = np.array([[1., 2.], [3., 4.]])
    fraction = np.array([[.2, .4], [.6, .8]])
    allocations = [allocation(reservoir("a"), fraction), allocation(reservoir("b", "M:v1"), 1.-fraction)]
    sample = parent(allocations=allocations, patch_multiplier=multiplier)
    knots = np.array([0., 1.])
    plan = prepare_coverage_inventory(surface(), sample, knots)
    multiplier[:] = 99.
    fraction[:] = 99.
    knots[:] = 99.
    allocations.clear()
    np.testing.assert_array_equal(plan.patch_multiplier, [[1., 2.], [3., 4.]])
    np.testing.assert_array_equal(plan.allocation_fractions[0], [[.2, .4], [.6, .8]])
    np.testing.assert_array_equal(plan.coverage_ml, [0., 1.])
    for values in (plan.coverage_ml, plan.patch_multiplier, plan.allocation_fractions,
                   plan.initial_parent_coverage_m2, plan.initial_pool_coverage_m2,
                   plan.parent_number, plan.pool_number, plan.allocation_balance_residual_m2):
        with pytest.raises(ValueError):
            values.setflags(write=True)
    with pytest.raises((FrozenInstanceError, AttributeError)):
        sample.parent_id = "changed"
    manifest = plan.manifest()
    manifest.clear()
    assert plan.manifest()
