#!/usr/bin/env python3
"""Bounded 3D support scenarios and conservative conditional work coupling."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mhdlab.thermal_expansion import load_ss304_material
from mhdlab.thermal_expansion_source import load_matched_history
from mhdlab.solid_foil_caloric import SealedNominalCaloric
from mhdlab.solid_foil_mechanics import structured_hex_mesh,transient_foil
from mhdlab.solid_foil_geometry import deform_hex_mesh,midplane_gap,parallel_plate_shape_sensitivity,remap_extensive_xy
from mhdlab.solid_foil_coupling import thermoelastic_exchange,coupled_cell_heat_step,von_mises_Pa,project_exchange_to_thermal_xy

SOURCE='/home/treefor/projects/Power-Flow/runs/quick_heating_256x32_domain40mm_height25mm_20261029'
SEAL='/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-20261029/source-manifest.json'
SEAL_SHA='8fc064a3341519c8f6254f2f0e971bff6638a6d6160d69176fa00c5d1b1a17fc'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,value):p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def run(out):
    out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
                         'font.size':10,'axes.grid':True,'grid.alpha':.22})
    mat=load_ss304_material(); tref=300.; nu=.29
    record=json.loads((ROOT/'data/materials/ss304-v1.json').read_text())
    proof=record['proof_stress_screening']
    screen=lambda t:np.interp(t,proof['temperature_K'],proof['value_Pa'])
    rho=mat.density_at_stress_free_temperature(tref)
    caloric=SealedNominalCaloric(SEAL,expected_source_seal_sha256=SEAL_SHA,reference_density_kg_m3=rho,
                                qualification_mode='evaluated_solid_1500')
    # Real C-ABI forward/inverse checks bind the thermal-owner model, no refit.
    checks=np.array([300.,301.,330.,368.])
    np.testing.assert_allclose(caloric.temperature_K(caloric.energy_J_m3(checks)),checks,rtol=0,atol=2e-13)
    h=load_matched_history(SOURCE,lower_volume='mirror_with_surface_validation')
    face_max=np.max(np.stack([v.max(axis=1) for v in h.face_temperature_k.values()]),axis=0)
    first_bad=np.flatnonzero(face_max>=368.)
    end_limit=h.face_time_s[first_bad[0]] if len(first_bad) else np.inf
    ids=np.flatnonzero((h.time_s<end_limit)&(h.upper_temperature_k.max(axis=(1,2))<368.))
    end=int(ids[-1]); source_times=h.time_s[:end+1]; source_T=h.upper_temperature_k[:end+1]
    if len(source_times)<2:raise ValueError('no smooth source prefix')
    source_energy=caloric.energy_J_m3(source_T)
    source_xyvolume=np.diff(h.depth_edges_m)[:,None]*np.diff(h.x_edges_m)[None,:]*.025
    source_extensive=source_energy*source_xyvolume
    summaries={};profiles={};exchanges={};mapping_receipts={};rejected={}
    base_args=dict(material=mat,stress_free_temperature_K=tref,poisson_ratio=nu)
    identities=dict(material_sha256=sha(ROOT/'data/materials/ss304-v1.json'),source_sha256=h.manifest['source_run_sha256'],
                    state_revision=0,caloric=caloric,proof_stress_screen_Pa=screen)

    def solve_case(label,nx,ny,nz,dt,boundary,edge_resolution=0):
        print('RUN '+label,flush=True)
        x=np.linspace(-.0025,.0025,nx+1);depth=np.linspace(0.,.0002,ny+1)
        if edge_resolution:
            ix=([0,1,2,4,8,16,32,64,128,192,224,240,248,252,254,255,256] if edge_resolution==1 else
                [0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,160,192,208,224,232,240,244,248,250,252,253,254,255,256])
            iy=([0,1,2,4,8,16,24,28,30,31,32] if edge_resolution==1 else
                [0,1,2,3,4,6,8,12,16,20,24,26,28,29,30,31,32])
            x=h.x_edges_m[ix];depth=h.depth_edges_m[iy];nx=len(x)-1;ny=len(depth)-1
        y=.0025+depth; z=np.linspace(-.0125,.0125,nz+1)
        mesh=structured_hex_mesh(x,y,z)
        mapped=[]; errors=[]
        target_xyvolume=np.diff(depth)[:,None]*np.diff(x)[None,:]*.025
        for values in source_extensive:
            transfer=remap_extensive_xy(h.x_edges_m,h.depth_edges_m,x,depth,values)
            mapped.append(transfer['values']/target_xyvolume)
            errors.append(transfer['relative_conservation_error'])
        mapped=np.array(mapped)
        count=int(np.ceil((source_times[-1]-source_times[0])/dt))
        times=np.linspace(source_times[0],source_times[-1],count+1)
        # Linear interpolation of the declared caloric energy, then inversion.
        energy=np.stack([np.interp(times,source_times,mapped[:,j,i]) for j in range(ny) for i in range(nx)],axis=1)
        temp2=caloric.temperature_K(energy).reshape(len(times),ny,nx)
        temp=np.repeat(temp2[:,None,:,:],nz,axis=1)
        result=transient_foil(mesh,times,temp,boundary=boundary,**base_args)
        try:
            exchange=thermoelastic_exchange(mesh,result,material=mat,stress_free_temperature_K=tref,
                mesh_revision=label,transaction_id=label+'-work',**identities)
        except ValueError as error:
            if 'proof-stress screening exceeded' not in str(error):raise
            ratios=von_mises_Pa(result['gauss_stress_Pa'])/screen(temp.reshape(len(times),-1))[...,None]
            first=int(np.flatnonzero(ratios.max(axis=(1,2))>1)[0])
            rejected[label]={'reason':str(error),'mesh_cells_xyz':[nx,ny,nz],
                'first_rejected_time_s':float(times[first]),'maximum_nominal_proof_screen_ratio':float(ratios.max()),
                'max_von_mises_Pa':float(von_mises_Pa(result['gauss_stress_Pa']).max()),
                'geometry_handoff_performed':False,'thermal_exchange_emitted':False,
                'interpretation':'elastic trial fails nominal proof screen; plastic response unsupported'}
            dump(out/(label+'-rejected.json'),rejected[label]);print('REJECT '+label+': '+str(error),flush=True)
            return
        heat_projection=project_exchange_to_thermal_xy(mesh,exchange,target_x_edges_m=h.x_edges_m,
            target_y_edges_m=h.upper_y_edges_m,expected_revision=0,expected_mesh_revision=label)
        heat_projection['cell_heat_J']=heat_projection['cell_heat_J'].tolist()
        dump(out/(label+'-native-thermal-load.json'),heat_projection)
        upper=deform_hex_mesh(mesh.nodes_m,mesh.cells,result['displacement_m'][-1],
            reference_density_kg_m3=rho,cell_temperature_k=temp[-1].reshape(-1),
            admissible_temperature_range_k=[293.15,773.15])
        lower_mesh=structured_hex_mesh(x,-y[::-1],z)
        lower_u=result['displacement_m'][-1].reshape(nz+1,ny+1,nx+1,3)[:,::-1,:,:]*[1,-1,1]
        lower=deform_hex_mesh(lower_mesh.nodes_m,lower_mesh.cells,lower_u.reshape(-1,3),
            reference_density_kg_m3=rho,cell_temperature_k=temp[-1,:,::-1,:].reshape(-1),
            admissible_temperature_range_k=[293.15,773.15])
        np.savez_compressed(out/(label+'-geometry.npz'),reference_nodes_m=mesh.nodes_m,cells=mesh.cells,
            displacement_m=result['displacement_m'][-1],temperature_K=temp[-1],lower_reference_nodes_m=lower_mesh.nodes_m,
            lower_displacement_m=lower_u.reshape(-1,3))
        gap=midplane_gap(lower,upper,z_m=0.)
        sensitivity=parallel_plate_shape_sensitivity(1.,gap['gap_m'])
        vm=von_mises_Pa(result['gauss_stress_Pa'])
        mass0=float(result['reference_mass_kg'])
        summaries[label]={'mesh_cells_xyz':[nx,ny,nz],'max_dt_s':float(np.diff(times).max()),'boundary':boundary,
            'max_displacement_m':float(np.max(np.linalg.norm(result['displacement_m'],axis=-1))),
            'max_von_mises_Pa':float(vm.max()),'max_nominal_proof_ratio':exchange['maximum_nominal_proof_screen_ratio'],
            'max_gap_narrowing_m':float(np.max(.005-gap['gap_m'])),
            'minimum_gap_m':float(np.min(gap['gap_m'])),
            'mass_relative_error':float(np.max(abs(result['current_mass_kg']/mass0-1))),
            'max_mapping_energy_relative_error':max(errors),
            'mechanical_work_residual_J':exchange['mechanical_work_residual_J'],
            'thermoelastic_debit_J':sum(exchange['cell_thermoelastic_debit_J']),
            'algorithmic_free_energy_debit_J':sum(exchange['cell_recoverable_strain_energy_change_J'])-sum(exchange['cell_stress_work_J']),
            'constitutive_internal_change_J':sum(exchange['cell_elastic_internal_energy_change_J']),
            'fixed_voltage_parallel_plate_max_relative_E_change':float(np.max(.005/np.asarray(gap['gap_m'])-1)),
            'field_solve_performed':False}
        profiles[label]={k:np.asarray(v).tolist() if isinstance(v,np.ndarray) else v for k,v in gap.items()}
        exchanges[label]=exchange
        dump(out/(label+'-exchange.json'),exchange)
        mapping_receipts[label]={'source_total_model_energy_J':source_extensive.sum(axis=(1,2)).tolist(),
            'maximum_relative_error':max(errors),'method':'extensive reference XY cell energy overlap then scalar caloric inversion; total height once',
            'source_energy_interpretation':'chosen caloric model evaluated on prescribed native T; not original constant-property heating energy ledger'}

    for boundary in ('free','short_plane_clamped','both_end_clamped'):
        solve_case(boundary,8,4,4,2.5e-10,boundary)
    solve_case('free_time_refined',8,4,4,1.25e-10,'free')
    solve_case('free_mesh_refined',16,8,8,2.5e-10,'free')
    solve_case('clamped_mesh_refined',16,8,8,2.5e-10,'both_end_clamped')
    solve_case('free_edge_resolved',16,10,4,2.5e-10,'free',edge_resolution=1)
    solve_case('free_edge_refined',28,16,4,2.5e-10,'free',edge_resolution=2)
    solve_case('clamped_edge_refined',28,16,4,2.5e-10,'both_end_clamped',edge_resolution=2)
    left=max(p['x_m'][0] for p in profiles.values());right=min(p['x_m'][-1] for p in profiles.values())
    xcommon=np.unique(np.r_[left,right,*[p['x_m'] for p in profiles.values()]])
    xcommon=xcommon[(xcommon>=left)&(xcommon<=right)]
    def compare(label):
        values=np.interp(xcommon,profiles[label]['x_m'],profiles[label]['gap_m'],left=np.nan,right=np.nan)
        if not np.isfinite(values).all():raise ValueError('profile comparison cannot extrapolate')
        return values
    base=compare('free')
    refined=compare('free_mesh_refined')
    comparison={
        'support_gap_span_Linf_m':float(np.ptp(np.stack([compare(b) for b in ('free','short_plane_clamped','both_end_clamped')]),axis=0).max()),
        'time_refinement_gap_Linf_m':float(np.max(abs(compare('free_time_refined')-base))),
        'mesh_refinement_gap_Linf_m':float(np.max(abs(refined-base))),
        'refined_free_to_clamped_gap_Linf_m':float(np.max(abs(compare('clamped_mesh_refined')-refined))),
        'edge_resolution_refinement_gap_Linf_m':float(np.max(abs(compare('free_edge_refined')-compare('free_edge_resolved')))),
        'edge_refined_free_to_clamped_gap_Linf_m':(float(np.max(abs(compare('clamped_edge_refined')-compare('free_edge_refined')))) if 'clamped_edge_refined' in profiles else None),
        'physical_convergence_established':False}
    dump(out/'support-summaries.json',summaries);dump(out/'midplane-gap-profiles.json',profiles)
    dump(out/'mapping-receipts.json',mapping_receipts)
    for name,packet in exchanges.items():dump(out/(name+'-exchange.json'),packet)

    # Actual two-way temperature correction, independent of prescribed native T.
    mesh=structured_hex_mesh([-.0025,0,.0025],[.0025,.0026,.0027],[-.0125,0,.0125])
    t0=np.full(len(mesh.cells),300.)
    heat=mesh.reference_volumes_m3*(caloric.energy_J_m3(310.)-caloric.energy_J_m3(300.))
    corrected={}
    for bc in ('free','both_end_clamped'):
        case=coupled_cell_heat_step(mesh,t0,heat,dt_s=1e-9,boundary=bc,mesh_revision='corrector-'+bc,
            transaction_id='corrector-'+bc,**base_args,**identities)
        corrected[bc]=case['exchange'];dump(out/('coupled-corrector-'+bc+'.json'),case['exchange'])
    wave=np.sqrt(float(mat.modulus(300.))*(1-nu)/(rho*(1+nu)*(1-2*nu)))
    summary={'schema':'solid-foil-support-and-work-v2','source_last_time_s':float(source_times[-1]),
        'raw_source_max_K':float(source_T.max()),'surface_prefix_cap_K':368.,'stress_free_temperature_K':tref,
        'poisson_ratio':nu,'poisson_basis':'nominal304 strip datum; held constant as explicit scenario',
        'caloric':caloric.receipt,'source_manifest':h.manifest,'supports':summaries,'comparison':comparison,'rejected_scenarios':rejected,
        'fastest_nominal_longitudinal_speed_m_s':float(wave),'short_to_midplane_wave_time_s':.0125/float(wave),
        'support_interpretation':'scenario bracket, not proven displacement bounds or measured fixture; no structural short thickness invented',
        'lower_body_interpretation':'explicit reflection of upper mechanical solution under symmetric brackets; native lower bulk unavailable, surfaces checked',
        'corrector':{k:{q:v for q,v in p.items() if q in ['end_temperature_K','corrector_iterations','coupled_energy_residual_J','maximum_nominal_proof_screen_ratio']} for k,p in corrected.items()},
        'physical_EOS_ready':False,'native_heat_feedback_integrated':False,'physical_field_handoff_ready':False,
        'limitations':['full-integration Q1 thin-foil bending not converged','no fixture/short stiffness or joint law','nominal proof screen, no plasticity','thermodynamic mode rejects coefficient knots','conditional Cp-as-Cv reference caloric model','parallel-plate fixed-voltage sensitivity, no deformed Maxwell solve'],
        'resources':{'affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0),'GPU_used':False}}
    dump(out/'summary.json',summary)
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    for name in ('free','short_plane_clamped','both_end_clamped','free_mesh_refined','free_edge_resolved','free_edge_refined'):
        xp=np.asarray(profiles[name]['x_m']);gp=np.asarray(profiles[name]['gap_m'])
        axes[0].plot(xp*1e3,(.005-gp)*1e9,label=name.replace('_',' '))
        axes[1].plot(xp*1e3,(.005/gp-1)*1e6,label=name.replace('_',' '))
    names=('free','short_plane_clamped','both_end_clamped')
    axes[2].bar(np.arange(3)-.18,[summaries[n]['thermoelastic_debit_J']*1e6 for n in names],.36,label='internal-energy thermal debit')
    axes[2].bar(np.arange(3)+.18,[summaries[n]['algorithmic_free_energy_debit_J']*1e6 for n in names],.36,label='strain-energy reservoir debit')
    axes[2].set_xticks(range(3),['free','short clamped','both clamped'],rotation=15)
    axes[0].set(xlabel='Width coordinate at z=0 (mm)',ylabel='Gap narrowing (nm)')
    axes[1].set(xlabel='Width coordinate at z=0 (mm)',ylabel='Fixed-V parallel-plate E change (ppm)')
    axes[2].set(ylabel='Signed energy transfer (µJ)')
    axes[0].legend(fontsize=6);axes[2].legend(fontsize=6)
    fig.suptitle(f'Conditional SS304 3-D trials at {source_times[-1]*1e9:.3f} ns • no Maxwell solve or physical fixture claim')
    fig.tight_layout();fig.savefig(out/'support-gap-field-work.png',dpi=170,facecolor='white',bbox_inches='tight');plt.close(fig)
    dump(out/'files-sha256.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps({'comparison':comparison,'corrector':summary['corrector'],'output':str(out)},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
