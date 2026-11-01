#!/usr/bin/env python3
"""Bounded native-conduction / mechanical-exchange demonstration, no EM heat."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from mhdlab.thermal_expansion import load_ss304_material
from mhdlab.solid_foil_mechanics import structured_hex_mesh
from mhdlab.solid_foil_caloric_v2 import SealedNominalCaloricV2,SOURCE_SEAL_SHA256
from mhdlab.solid_foil_native_heat import NativeConductiveHeat
from mhdlab.solid_foil_native_coupling import NativeFoilCorrector,file_digest

MATERIAL_SEAL=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029/source-manifest.json')
SCENARIO={'schema':'manufactured-native-thermoelastic-v3','geometry_m':[2e-6,1e-6,2e-6],
    'cells_xyz':[4,2,2],'temperature':'299.5+0.6*cos(pi*x/width), uniform y,z',
    'initial_velocity':'signed 0.2m/s linear dilation at x faces, centered',
    'initial_displacement':'affine thermal strain at mean initial T, centered',
    'boundary':'free','dt_s':1e-10,'external_heat_J':0.,'reference_temperature_K':300.,
    'interpretation':'microscale numerical binding benchmark; not the 5mm foil or production water breakdown'}
SOURCE_SHA=hashlib.sha256(json.dumps(SCENARIO,sort_keys=True).encode()).hexdigest()


def make_case(executable_path,*,velocity_sign=1.,initial_shift_K=0.,boundary='free'):
    material=load_ss304_material();rho=material.density_at_stress_free_temperature(300.)
    caloric=SealedNominalCaloricV2(MATERIAL_SEAL,expected_source_seal_sha256=SOURCE_SEAL_SHA256,
        reference_density_kg_m3=rho,qualification_mode='nominal_1800')
    x=np.linspace(0.,2e-6,5);y=np.linspace(0.,1e-6,3);z=np.linspace(-1e-6,1e-6,3)
    mesh=structured_hex_mesh(x,y,z)
    t=np.broadcast_to(299.5+initial_shift_K+.6*np.cos(np.pi*x/2e-6),(3,5)).copy()
    native=NativeConductiveHeat(executable_path,expected_executable_sha256=file_digest(executable_path),
        nx=4,ny=2,x0_m=0.,y0_m=0.,width_m=2e-6,thickness_m=1e-6,height_m=2e-6,
        initial_temperature_K=t,reference_density_kg_m3=rho,mesh_revision='native-manufactured-4x2-v3',
        caloric_sha256=SOURCE_SEAL_SHA256,source_sha256=SOURCE_SHA,electrode_id='manufactured-solid')
    record=json.loads((ROOT/'data/materials/ss304-v1.json').read_text())['proof_stress_screening']
    screen=lambda t:np.interp(t,record['temperature_K'],record['value_Pa'])
    velocity=np.zeros_like(mesh.nodes_m);velocity[:,0]=velocity_sign*.2*(mesh.nodes_m[:,0]-1e-6)/1e-6
    displacement=material.strain(299.5+initial_shift_K,stress_free_temperature_K=300.)*(mesh.nodes_m-[1e-6,.5e-6,0.])
    if boundary!='free':velocity[:]=0;displacement[:]=0
    return NativeFoilCorrector(native,mesh,caloric=caloric,material=material,stress_free_temperature_K=300.,
        poisson_ratio=.29,boundary=boundary,mechanical_mesh_revision='solid-manufactured-4x2x2-v3',
        material_sha256=file_digest(ROOT/'data/materials/ss304-v1.json'),source_sha256=SOURCE_SHA,
        proof_stress_screen_Pa=screen,proof_stress_screen_sha256=hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest(),
        initial_displacement_m=displacement,initial_velocity_m_s=velocity)


def run(executable_path,out):
    out.mkdir(parents=True,exist_ok=False)
    if sorted(os.sched_getaffinity(0))!=[7] or os.getpriority(os.PRIO_PROCESS,0)<15:
        raise RuntimeError('demo must use the CPU7 lock/nice15 frozen launcher')
    histories={};checks={}
    for sign in (1.,-1.):
        key='expansion' if sign>0 else 'compression';case=make_case(executable_path,velocity_sign=sign)
        before=case.state_digest();initial=case.native.state()['temperature_K'].copy()
        trial=case.prepare(1e-10,transaction_id=key+'-discard')
        assert case.state_digest()==before
        discarded=trial.receipt
        case.discard();assert case.state_digest()==before
        trial=case.prepare(1e-10,transaction_id=key+'-step1');case.commit(trial)
        histories[key]=[trial.receipt]
        case.save(out/(key+'-checkpoint'))
        resumed=make_case(executable_path,velocity_sign=sign);resumed.load(out/(key+'-checkpoint'))
        assert resumed.state_digest()==case.state_digest()
        for step in (2,3,4):
            a=case.prepare(1e-10,transaction_id=f'{key}-step{step}')
            b=resumed.prepare(1e-10,transaction_id=f'{key}-step{step}')
            assert a.receipt==b.receipt
            case.commit(a);resumed.commit(b)
            assert case.state_digest()==resumed.state_digest()
            histories[key].append(a.receipt)
        checks[key]={'rollback_preserved_state':True,'restart_bitwise_equal':True,
            'initial_minmax_K':[float(initial.min()),float(initial.max())],
            'final_minmax_K':[float(case.native.state()['temperature_K'].min()),float(case.native.state()['temperature_K'].max())],
            'max_global_energy_residual_J':max(abs(r['global_internal_plus_kinetic_residual_J']) for r in histories[key]),
            'max_native_nodal_energy_residual_J':max(r['native_budget']['nodal_residual_max_abs_J'] for r in histories[key]),
            'peak_native_conduction_node_J':max(float(np.max(np.abs(r['native_budget']['nodal_ledger']['conduction_out_J']))) for r in histories[key]),
            'signed_mechanical_transfer_by_step_J':[r['native_budget']['mechanical_exchange_J'] for r in histories[key]],
            'max_proof_screen_ratio':max(r['mechanical_exchange']['maximum_nominal_proof_screen_ratio'] for r in histories[key]),
            'corrector_iterations':[r['corrector_iterations'] for r in histories[key]],
            'final_state_digest':case.state_digest()}
        (out/(key+'-history.json')).write_text(json.dumps(histories[key],indent=2,allow_nan=False)+'\n')
        case.native.close();resumed.native.close()
    summary={'schema':'native-thermoelastic-v3-validation','scenario':SCENARIO,'source_sha256':SOURCE_SHA,
        'native_executable_sha256':file_digest(executable_path),'caloric_source_seal_sha256':SOURCE_SEAL_SHA256,
        'cases':checks,'external_heat_J':0.,'joule_J':0.,'authored_heat_J':0.,
        'native_conduction_solved':True,'joint_rollback_restart_verified':True,
        'mechanics_max_K':773.15,'reference_mesh_fixed':True,'physical_EOS_ready':False,
        'production_run_performed':False,'GPU_used':False,
        'limitations':['manufactured microscale benchmark, not prescribed source foil history',
            'fixed reference geometry; no deformed Maxwell or native moving mesh',
            'conditional Cp-as-Cv caloric energy; no physical EOS',
            'same XY grid; cell caloric averaging reduces subcell thermal variation',
            'piecewise thermoelastic coefficient-knot crossings rejected',
            'nominal proof-stress screen, no plasticity or measured fixture qualification'],
        'resources':{'affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0)}}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native-executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();run(args.native_executable,args.output)
