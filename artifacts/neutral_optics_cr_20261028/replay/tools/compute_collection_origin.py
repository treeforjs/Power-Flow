"""Attribute a saved optical unit-response coefficient to emitting material faces."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.meta_path[:]=[f for f in sys.meta_path if type(f).__module__!='_editable_skbc_mhdlab']
from mhdlab.neutral_vacuum_map import observe_vacuum_packets
from compute_target_collection import gate_mean, write, overlap

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--transport-dir',type=Path,required=True)
    p.add_argument('--collection-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();out=a.output_dir
    if out.exists(): p.error('use a fresh output directory')
    source=a.transport_dir;collection=a.collection_dir
    meta=json.loads((source/'metadata.json').read_text())
    report=json.loads((collection/'collection_report.json').read_text())
    case=json.loads((collection/'case.json').read_text())
    with np.load(source/'neutral_maps.npz') as f: maps={k:f[k] for k in f.files}
    with np.load(collection/'collection_kernel.npz') as f: kernel={k:f[k] for k in f.files}
    # Keep exactly the original knots bracketing both gates: the resulting
    # interpolant agrees with the original saved-map gate convention.
    t=maps['time_s'];indices=[]
    for lo,hi in case['gates_s']:
        start=max(0,np.searchsorted(t,lo,side='right')-1)
        stop=min(len(t)-1,np.searchsorted(t,hi,side='left'))
        indices.extend(range(start,stop+1))
    indices=np.unique(indices);times=t[indices]
    packet_fields=['birth_time_s','position_m','velocity_m_s','particle_weight','species_index','mass_kg']
    with np.load(source/'source_packets.npz') as f: packets={k:f[k] for k in packet_fields+['face_index']}
    density_sum=np.zeros_like(maps['number_density_m3'][indices]);face_profiles=[]
    for i,record in enumerate(meta['source_records']):
        selected=packets['face_index']==i
        print('Observing '+record['face']+' origin',flush=True)
        arrays,_=observe_vacuum_packets(**{key:packets[key][selected] for key in packet_fields},time_s=times,
            x_edges_m=maps['x_edges_m'],y_edges_m=maps['y_edges_m'],foil_bounds_m=maps['foil_bounds_m'],
            invariant_current_length_m=meta['invariant_current_length_m'],initial_surface_particles=record['initial_particles'])
        density=arrays['number_density_m3'];density_sum+=density
        profiles=[]
        for lo,hi in case['gates_s']:
            profiles.append(np.einsum('ryx,syx->sr',kernel['kernel_m3'],gate_mean(times,density,lo,hi))*(hi-lo))
        face_profiles.append(np.array(profiles))
    # Differences from histogram summation order are measured, not hidden.
    original=maps['number_density_m3'][indices]
    volumes=np.diff(maps['y_edges_m'])[:,None]*np.diff(maps['x_edges_m'])[None,:]*meta['invariant_current_length_m']
    difference=float(np.max(np.sum(np.abs(density_sum-original)*volumes,axis=(-2,-1))/np.maximum(1,np.sum(original*volumes,axis=(-2,-1)))))
    if difference>1e-11: raise ValueError('face-origin density sum disagrees with saved map')
    profiles=np.array(face_profiles) # face, gate, species, row
    out.mkdir(parents=True)
    np.savez_compressed(out/'face_origin_profiles.npz',profiles=profiles,face_labels=np.array(meta['emitting_faces']),
        species=maps['species'],gates_s=np.array(case['gates_s']),detector_bin_edges_m=kernel['detector_bin_edges_m'])
    rows=[]
    for g,gate in enumerate(case['gates_s']):
        items=[]
        for s,species in enumerate(maps['species']):
            coeff=profiles[:,g,s].sum(axis=-1);total=coeff.sum()
            items.append({'species':str(species),'unit_rate_total_coefficient':float(total),
                'face_coefficients':dict(zip(meta['emitting_faces'],map(float,coeff))),
                'side_origin_fraction':None if total==0 else float(coeff[2:].sum()/total),
                'left_right_origin_profile_overlap':overlap(profiles[2,g,s],profiles[3,g,s]),
                'profile':profiles[:,g,s].sum(axis=0).tolist()})
        rows.append({'gate_s':gate,'species':items})
    result={'schema':'mhdlab-optical-source-origin-v1','status':'conditional_unit_photon_rate_per_neutral',
        'physical_detector_counts_predicted':False,'source_origin_distinct_from_current_edge_position':True,
        'foil_shadow_and_finite_pupil_included':True,'reference_scenario':case['reference_scenario'],
        'whole_density_reconstruction_maximum_relative_L1_error':difference,'gates':rows,
        'inputs_sha256':{str(path):sha(path) for path in [source/'source_packets.npz',source/'metadata.json',collection/'collection_kernel.npz']}}
    plot(out,profiles,kernel['detector_bin_edges_m'],meta['emitting_faces'],maps['species'])
    result['outputs_sha256']={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    write(out/'source_origin_report.json',result)
    print(json.dumps({'report':str(out/'source_origin_report.json'),'density_error':difference}))

def plot(out,profiles,edges,faces,species):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained');y=.5*(edges[1:]+edges[:-1])*1e3
    for s,(ax,name) in enumerate(zip(axes.flat,species)):
        for f,label in enumerate(faces):
            row=profiles[f,-1,s];ax.plot(y,np.where(row>0,row,np.nan),label=label)
        ax.set(title=str(name),xlabel='Conditional slit-image spatial coordinate [mm]',ylabel='Gate photons for rate 1 s$^{-1}$ per neutral',yscale='log')
        peak=float(profiles[:,-1,s].max())
        if peak>0: ax.set_ylim(peak*1e-8,peak*2)
        ax.yaxis.label.set_size(9);ax.legend();ax.grid(alpha=.2)
    fig.suptitle('57–60 ns collection by original emitting face\nConditional geometry and equal per-neutral rate; display spans eight decades')
    fig.savefig(out/'source_face_collection.png',dpi=150);plt.close(fig)

if __name__=='__main__': main()
