"""The fitted monolayer coordinate owns one physical H2O source inventory."""
import numpy as np
import pytest

from mhdlab.coverage_inventory import InventoryPoolAllocation, MonolayerDefinition
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.neutral_cr_water_source import bind_initial_h2o_monolayers, build_common_water_source_sweep
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_release import TemkinReservoir


def inputs():
    base=FlashSurfaceHistory([0.,1e-7],[0.,1.],[[600.],[600.]],'synthetic',
        'water-inventory-control','Authored constant temperature','Authored stationary face')
    surface=PlanarSurfaceHistory([0.,1e-7],[0.,.001,.003],[0.,.01],np.full((2,1,2),600.),base.attributes)
    layer=MonolayerDefinition(1e17,'Authored water-layer convention','Independent manufactured areal density','manufactured')
    pool=TemkinReservoir(name='water_parent',state_id='H2O:authored',attempt_frequency_s=1e7,
        activation_energy_ev=.08,desorption_energy_ev=.05,internal_energy_ev=.02,
        coverage_coefficient=.2,reference_coverage_m2=2e18,max_coverage_m2=2e18,
        kinetic_provenance='Authored coverage-dependent release',energy_provenance='Authored energy channels')
    allocation=InventoryPoolAllocation(pool,1.,'One complete water parent pool','manufactured')
    return surface, dict(monolayer=layer,allocations=[allocation],state_species={'H2O:authored':'H2O'},
        provenance='Manufactured common-water parameter binding',status='manufactured')


def test_monolayer_candidates_set_actual_shared_initial_inventory_once():
    surface,kw=inputs()
    plan=bind_initial_h2o_monolayers(surface,[0.,1.,2.],**kw)
    np.testing.assert_allclose(plan.parent_number,[0.,3e12,6e12],rtol=2e-15)
    np.testing.assert_array_equal(plan.initial_pool_coverage_m2.sum(axis=1),plan.initial_parent_coverage_m2)
    assert plan.parent.parent_id=='H2O'
    assert plan.parent.monolayer.areal_density_m2==1e17
    assert plan.parent.allocations[0].reservoir.reference_coverage_m2==2e18
    scalar=bind_initial_h2o_monolayers(surface,.5,**kw)
    np.testing.assert_allclose(scalar.parent_number,[1.5e12],rtol=2e-15)


def test_non_linear_release_is_solved_again_instead_of_amplitude_scaling():
    surface,kw=inputs()
    sweep=build_common_water_source_sweep(surface,[0.,1.,2.],parameter_provenance='Authored nonlinear source cases',**kw)
    n=[float(np.sum(h.released_m2*h.patch_areas_m2)) for h in sweep.histories]
    assert n[0]==0
    assert n[2]>2*n[1]
    assert not np.isclose(n[2],2*n[1],rtol=1e-3,atol=0)
    assert all(h.reservoirs[0].reference_coverage_m2==2e18 for h in sweep.histories)


@pytest.mark.parametrize('value',[-1.,float('nan'),float('inf'),True,[1.,1.],[2.,1.]])
def test_invalid_free_parameter_candidates_reject(value):
    surface,kw=inputs()
    with pytest.raises(ValueError):
        bind_initial_h2o_monolayers(surface,value,**kw)


def test_fragment_allocation_and_missing_monolayer_definition_reject():
    surface,kw=inputs()
    kw['state_species']={'H2O:authored':'OH'}
    with pytest.raises(ValueError,match='fragments come from chemistry'):
        bind_initial_h2o_monolayers(surface,1.,**kw)
    surface,kw=inputs();kw['monolayer']=None
    with pytest.raises(ValueError,match='MonolayerDefinition'):
        bind_initial_h2o_monolayers(surface,1.,**kw)
