#!/usr/bin/env python3
"""Actual sealed source/water + native conduction + mechanics component cases."""
from pathlib import Path
import argparse,hashlib,json,os,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from mhdlab.thermal_expansion import load_ss304_material
from mhdlab.solid_foil_mechanics import structured_hex_mesh
from mhdlab.solid_foil_caloric_v2 import SealedNominalCaloricV2,SOURCE_SEAL_SHA256
from mhdlab.solid_foil_native_heat_v4 import NativeSourceConductiveHeat,electron_action
from mhdlab.solid_foil_native_coupling_v4 import NativeSourceFoilCorrector,digest,file_digest

MATERIAL_SEAL=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029/source-manifest.json')
SCENARIO={'schema':'sealed-source-native-mechanics-component-v4','geometry_m':[2e-6,1e-6,2e-6],'cells_xyz':[4,2,2],
 'thermal_state':'mean +0.2*cos(pi*x/width), uniform y,z; open smooth coefficient segment',
 'initial_mechanics':'affine thermal strain at mean T; signed0.2m/s x dilation at faces',
 'surface_partitions':'LowY first half oxide-water; second half clean4.3eV electron component hypothesis',
 'source_weight':1000.,'source_weight_is_predicted_event_count':False,
 'water':'sealed ideal-oxide surrogate;0.33 monolayer; explicit range selection',
 'source_surface':'stationary fixed reference surface approximation; no moving-wall momentum/work coupling',
 'field':'explicit prescribed one-sided total E/B for component verification; no field solve',
 'Joule_authored_heat_J':0.,'production_run':False}
SOURCE_SHA=digest(SCENARIO)


def make_case(executable_path,*,temperature_K=700.,velocity_sign=1.,water_range=0,boundary='free',water_enabled=True):
    material=load_ss304_material();rho=material.density_at_stress_free_temperature(300.)
    caloric=SealedNominalCaloricV2(MATERIAL_SEAL,expected_source_seal_sha256=SOURCE_SEAL_SHA256,
        reference_density_kg_m3=rho,qualification_mode='nominal_1800')
    x=np.linspace(0.,2e-6,5);y=np.linspace(0.,1e-6,3);z=np.linspace(-1e-6,1e-6,3)
    mesh=structured_hex_mesh(x,y,z);t=np.broadcast_to(temperature_K+.2*np.cos(np.pi*x/2e-6),(3,5)).copy()
    patches=[{'id':'oxide-water','gas_state_id':'H2O(g,0)','parent_id':'adsorbed-H2O','scenario_sha256':SOURCE_SHA,
        'face':0,'begin':0.,'end':.5,'water_enabled':water_enabled,'water_range':water_range,'initial_coverage':.33,
        'registry_mass_kg':2.9915076246e-26},
        {'id':'clean-electron','gas_state_id':'H2O(g,0)','parent_id':'unused-water','scenario_sha256':SOURCE_SHA,
        'face':0,'begin':.5,'end':1.,'water_enabled':False,'water_range':water_range,'initial_coverage':0.,
        'registry_mass_kg':2.9915076246e-26}]
    native=NativeSourceConductiveHeat(executable_path,expected_executable_sha256=file_digest(executable_path),
        nx=4,ny=2,x0_m=0.,y0_m=0.,width_m=2e-6,thickness_m=1e-6,height_m=2e-6,
        initial_temperature_K=t,reference_density_kg_m3=rho,mesh_revision='native-component-4x2-v4',
        caloric_sha256=SOURCE_SEAL_SHA256,source_sha256=SOURCE_SHA,electrode_id='component-solid',
        patches=patches,electron_state_id='e-',electron_parent_id='electrode-electrons')
    proof=json.loads((ROOT/'data/materials/ss304-v1.json').read_text())['proof_stress_screening']
    screen=lambda t:np.interp(t,proof['temperature_K'],proof['value_Pa'])
    velocity=np.zeros_like(mesh.nodes_m);velocity[:,0]=velocity_sign*.2*(mesh.nodes_m[:,0]-1e-6)/1e-6
    displacement=material.strain(temperature_K,stress_free_temperature_K=300.)*(mesh.nodes_m-[1e-6,.5e-6,0.])
    if boundary!='free':velocity[:]=0;displacement[:]=0
    try:
        return NativeSourceFoilCorrector(native,mesh,caloric=caloric,material=material,stress_free_temperature_K=300.,
            poisson_ratio=.29,boundary=boundary,mechanical_mesh_revision='solid-component-4x2x2-v4',
            material_sha256=file_digest(ROOT/'data/materials/ss304-v1.json'),source_sha256=SOURCE_SHA,
            proof_stress_screen_Pa=screen,proof_stress_screen_sha256=digest(proof),
            initial_displacement_m=displacement,initial_velocity_m_s=velocity)
    except BaseException:native.close();raise


def inputs(case,dt_s,*,begin=False,field_V_m=1e9,B_T=0.,layer_id=1,weight=1000.):
    state=case.native.state();revision=state['revision'];source=case.native.source_state()
    stage={'common_state':revision+1,'field_state':revision+1,'source_mesh_revision':1,'source_material_revision':1,
        'field_time_s':state['time_s'],'E_x_V_m':0.,'E_y_V_m':field_V_m,'E_z_V_m':0.,
        'B_x_T':B_T,'B_y_T':0.,'B_z_T':0.,'material_phase_admitted':True,'planar_barrier_admitted':True,
        'electron_supply_admitted':True,'field_basis':'prescribed total one-sided component fields; no field solve'}
    stage['field_sha256']=digest(stage)
    actions=[]
    common=dict(patch_index=1,patch_id=2,origin_m=[1.5e-6,0.,0.],work_function_eV=4.3,local_linear_vacuum_extent_m=25e-6)
    if begin:
        actions.append(electron_action(kind=1,layer_id=layer_id,weight=weight,duration_s=0.,**common))
        actions.append(electron_action(kind=2,layer_id=layer_id,weight=weight,duration_s=dt_s,**common))
    for l in source['layers']:
        if int(l['status'])==0:
            actions.append(electron_action(kind=2,layer_id=int(l['layer_id']),weight=l['weight'],duration_s=dt_s,**common))
    draws=[dict(patch_index=i,normal=.4,radius=.6,angle=.3,along=.5,z=.5,time=.5)
           for i,p in enumerate(case.native.binding['patches']) if p['water_enabled']]
    return dict(source_stage=stage,source_actions=actions,water_draws=draws)


def run(executable_path,out):
    out.mkdir(parents=True,exist_ok=False)
    if sorted(os.sched_getaffinity(0))!=[7] or os.getpriority(os.PRIO_PROCESS,0)<15:
        raise RuntimeError('CPU7 shared lock, nice15 and frozen runtime required')
    summary={'schema':'sealed-source-native-mechanics-v4-validation','scenario':SCENARIO,'source_sha256':SOURCE_SHA,
        'native_executable_sha256':file_digest(executable_path),'cases':{},'Joule_authored_heat_J':0.,
        'gas_reference_J':0.,'mechanics_max_K':773.15,'moving_source_surface_coupled':False,'production_run_performed':False,
        'physical_EOS_ready':False,'GPU_used':False,'resources':{'CPU_affinity':sorted(os.sched_getaffinity(0)),
        'nice':os.getpriority(os.PRIO_PROCESS,0),'thread_jobs':1}}
    for sign in (1.,-1.):
        key='expansion' if sign>0 else 'compression';c=make_case(executable_path,velocity_sign=sign);history=[]
        before=c.state_digest();first_inputs=inputs(c,1e-15,begin=True)
        discarded=c.prepare(1e-15,transaction_id=key+'-rollback',**first_inputs);c.discard()
        assert c.state_digest()==before
        first=c.prepare(1e-15,transaction_id=key+'-step1',**first_inputs);c.commit(first);history.append(first.receipt)
        c.save(out/(key+'-checkpoint'));r=make_case(executable_path,velocity_sign=sign);r.load(out/(key+'-checkpoint'))
        assert c.state_digest()==r.state_digest()
        checkpoint_layers=c.native.source_state()['layers']
        for step,dt in enumerate((1e-12,1e-10,1e-10),2):
            source=inputs(c,dt);a=c.prepare(dt,transaction_id=f'{key}-step{step}',**source)
            b=r.prepare(dt,transaction_id=f'{key}-step{step}',**source)
            assert a.receipt==b.receipt
            c.commit(a);r.commit(b);assert c.state_digest()==r.state_digest();history.append(a.receipt)
        (out/(key+'-history.json')).write_text(json.dumps(history,indent=2,allow_nan=False)+'\n')
        summary['cases'][key]={'rollback_preserved_complete_state':True,'restart_bitwise_equal':True,
            'checkpoint_layer_statuses':[int(l['status']) for l in checkpoint_layers],
            'max_joint_energy_residual_J':max(abs(h['global_combined_energy_residual_J']) for h in history),
            'max_native_nodal_energy_residual_J':max(h['native_budget']['nodal_residual_max_abs_J'] for h in history),
            'max_proof_screen_ratio':max(h['mechanical_exchange']['maximum_nominal_proof_screen_ratio'] for h in history),
            'corrector_iterations':[h['corrector_iterations'] for h in history],
            'mechanical_exchange_J':[h['native_budget']['mechanical_exchange_J'] for h in history],
            'electron_wall_J':[h['native_budget']['source_wall_J'] for h in history],
            'water_adsorbed_change_J':[h['native_budget']['adsorbed_reference_change_J'] for h in history],
            'source_field_work_J':[h['native_budget']['source_field_work_J'] for h in history],
            'final_temperature_range_K':[float(c.native.state()['temperature_K'].min()),float(c.native.state()['temperature_K'].max())],
            'final_state_digest':c.state_digest(),'final_source_state':c.native.source_state()}
        c.native.close();r.native.close()
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='cases'},indent=2))
    for key,c in summary['cases'].items():print(key,json.dumps({k:v for k,v in c.items() if k!='final_source_state'}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native-executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.native_executable,a.output)
