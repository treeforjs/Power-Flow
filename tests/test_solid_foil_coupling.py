import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from mhdlab.thermal_expansion import load_ss304_material
from mhdlab.solid_foil_mechanics import structured_hex_mesh,transient_foil
from mhdlab.solid_foil_coupling import thermoelastic_exchange,coupled_cell_heat_step,von_mises_Pa,project_exchange_to_thermal_xy
from mhdlab.solid_foil_caloric import SealedNominalCaloric


class CaloricFixture:
    def __init__(self,material):
        self.receipt={"reference_density_kg_m3":material.density_at_stress_free_temperature(300.),
                      "basis":"manufactured constant-capacity stress-free caloric model"}
    def energy_J_m3(self,t):return 3.8e6*(np.asarray(t)-300.)
    def capacity_J_m3_K(self,t):return np.asarray(t)*0+3.8e6
    def temperature_K(self,e):return np.asarray(e)/3.8e6+300.


def fixtures():
    mat=load_ss304_material()
    mesh=structured_hex_mesh([-.0025,0,.0025],[.0025,.0026,.0027],[-.0125,0,.0125])
    caloric=CaloricFixture(mat)
    args=dict(material=mat,stress_free_temperature_K=300.,poisson_ratio=.29,boundary="free",
              caloric=caloric,state_revision=0,mesh_revision="mesh0",material_sha256="a"*64,
              source_sha256="b"*64,transaction_id="heat0",proof_stress_screen_Pa=lambda t:np.full_like(t,190e6))
    return mesh,args


def test_nominal_corrector_conserves_combined_internal_kinetic_heat_energy():
    mesh,args=fixtures(); t0=np.full(len(mesh.cells),300.); heat=mesh.reference_volumes_m3*3.8e6*10
    trial=coupled_cell_heat_step(mesh,t0,heat,dt_s=1e-9,**args)
    packet=trial["exchange"]
    assert packet["corrector_iterations"]>1
    assert np.max(abs(np.array(packet["end_temperature_K"])-310.))>1e-3
    assert abs(packet["coupled_energy_residual_J"])<1e-10*heat.sum()
    assert np.max(np.abs(packet["cell_energy_residual_J"]))<1e-10*heat.sum()
    np.testing.assert_allclose(np.array(packet["cell_heat_to_thermal_J"]),-np.array(packet["cell_thermoelastic_debit_J"]))
    assert packet["physical_EOS_ready"] is False and packet["accepted_state_mutated"] is False
    np.testing.assert_array_equal(t0,300.)
    transfer=project_exchange_to_thermal_xy(mesh,packet,target_x_edges_m=np.linspace(-.0025,.0025,5),
        target_y_edges_m=np.linspace(.0025,.0027,5),expected_revision=0,expected_mesh_revision='mesh0')
    assert transfer['cell_heat_J'].shape==(4,4)
    assert transfer['target_total_J']==pytest.approx(sum(packet['cell_heat_to_thermal_J']),rel=5e-13)
    assert transfer['reference_axial_extent_m']==.025
    with pytest.raises(ValueError,match='stale'):
        project_exchange_to_thermal_xy(mesh,packet,target_x_edges_m=mesh.x_edges_m,
            target_y_edges_m=mesh.y_edges_m,expected_revision=1,expected_mesh_revision='mesh0')


def test_zero_heat_no_motion_no_exchange():
    mesh,args=fixtures()
    trial=coupled_cell_heat_step(mesh,np.full(len(mesh.cells),300.),np.zeros(len(mesh.cells)),dt_s=1e-9,**args)
    np.testing.assert_array_equal(trial["mechanics"]["displacement_m"],0.)
    np.testing.assert_array_equal(trial["exchange"]["cell_thermoelastic_debit_J"],0.)


def test_coefficient_crossing_rejects_whole_coupled_trial():
    mesh,args=fixtures(); t0=np.full(len(mesh.cells),300.)
    with pytest.raises(ValueError,match="knots"):
        coupled_cell_heat_step(mesh,t0,mesh.reference_volumes_m3*3.8e6*80,dt_s=1e-9,
                              **dict(args,proof_stress_screen_Pa=lambda t:np.full_like(t,1e12)))
    np.testing.assert_array_equal(t0,300.)


def test_yield_screen_and_density_mismatch_are_not_ignored():
    mesh,args=fixtures(); t0=np.full(len(mesh.cells),300.)
    with pytest.raises(ValueError,match="proof-stress"):
        coupled_cell_heat_step(mesh,t0,mesh.reference_volumes_m3*3.8e6*10,dt_s=1e-9,
                              **dict(args,proof_stress_screen_Pa=lambda t:np.full_like(t,1.)))
    args["caloric"].receipt["reference_density_kg_m3"]=8000.
    with pytest.raises(ValueError,match="density"):
        coupled_cell_heat_step(mesh,t0,mesh.reference_volumes_m3*3.8e6,dt_s=1e-9,**args)


def test_von_mises_respects_physical_shear_and_hydrostatic_stress():
    assert von_mises_Pa([100.,100.,100.,0.,0.,0.])==0.
    assert von_mises_Pa([0.,0.,0.,100.,0.,0.])==pytest.approx(np.sqrt(3)*100.)


def test_thermal_source_seal_mismatch_rejected_before_loading_library(tmp_path):
    seal=tmp_path/'seal.json';seal.write_text('{}')
    with pytest.raises(ValueError,match="seal SHA256"):
        SealedNominalCaloric(seal,expected_source_seal_sha256='0'*64,reference_density_kg_m3=7900.,
                            qualification_mode="evaluated_solid_1500")


def test_exchange_uses_sum_of_subcycle_stress_work():
    mesh,args=fixtures(); mat=args['material']; times=np.linspace(0.,2e-8,11)
    temperatures=np.repeat(np.linspace(300.,310.,11)[:,None],len(mesh.cells),axis=1)
    result=transient_foil(mesh,times,temperatures,material=mat,stress_free_temperature_K=300.,
                          poisson_ratio=.29,boundary="free")
    keys={k:v for k,v in args.items() if k not in ['poisson_ratio','boundary']}
    packet=thermoelastic_exchange(mesh,result,**keys)
    assert abs(packet['kinetic_energy_change_J']+sum(packet['cell_stress_work_J']))<1e-15
    assert abs(packet['mechanical_plus_debit_residual_J'])<1e-14
    plain=dict(packet);checksum=plain.pop('payload_sha256')
    assert checksum==hashlib.sha256(json.dumps(plain,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def test_mismatched_material_reference_and_canceling_local_work_rejected():
    mesh,args=fixtures();times=[0.,1e-9]
    result=transient_foil(mesh,times,np.repeat([[300.],[310.]],len(mesh.cells),axis=1),
        material=args['material'],stress_free_temperature_K=300.,poisson_ratio=.29,boundary='free')
    keys={k:v for k,v in args.items() if k not in ['poisson_ratio','boundary']}
    with pytest.raises(ValueError,match='reference binding'):
        thermoelastic_exchange(mesh,result,**dict(keys,stress_free_temperature_K=301.))
    bad=dict(result);bad['gauss_stress_work_increment_J']=result['gauss_stress_work_increment_J'].copy()
    bad['gauss_stress_work_increment_J'][0,0,0]+=1e-8
    bad['gauss_stress_work_increment_J'][0,1,0]-=1e-8
    with pytest.raises(ValueError,match='local stress-work'):
        thermoelastic_exchange(mesh,bad,**keys)
    wrong=structured_hex_mesh([-.003,0,.003],[.0025,.0026,.0027],[-.0125,0,.0125])
    with pytest.raises(ValueError,match='identity mismatch'):
        thermoelastic_exchange(wrong,result,**keys)
