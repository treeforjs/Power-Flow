"""Conditional finite-pupil target collection and edge-overlap coefficients."""
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
from mhdlab.target_optics import RelayGeometry, build_collection_kernel

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value): path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def gate_mean(times,values,lo,hi):
    if not times[0]<=lo<hi<=times[-1]: raise ValueError('gate outside saved population support')
    knots=np.r_[lo,times[(times>lo)&(times<hi)],hi]
    indices=np.clip(np.searchsorted(times,knots,side='right')-1,0,len(times)-2)
    f=(knots-times[indices])/(times[indices+1]-times[indices])
    samples=(1-f[:,None,None,None])*values[indices]+f[:,None,None,None]*values[indices+1]
    return np.einsum('t,tsij->sij',np.diff(knots)/(hi-lo),.5*(samples[:-1]+samples[1:]))

def overlap(a,b):
    sa,sb=a.sum(),b.sum()
    return None if sa==0 or sb==0 else float(np.minimum(a/sa,b/sb).sum())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,default=ROOT/'runs/neutral_desorption_current_17133_perimeter_v2/neutral_maps.npz')
    p.add_argument('--metadata',type=Path,default=ROOT/'runs/neutral_desorption_current_17133_perimeter_v2/metadata.json')
    p.add_argument('--case',type=Path,default=ROOT/'configs/target_collection_scenarios.json')
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();out=args.output_dir
    if out.exists(): p.error('use a fresh output directory')
    case=json.loads(args.case.read_bytes());metadata=json.loads(args.metadata.read_bytes())
    if case['schema']!='mhdlab-target-collection-scenarios-v1' or case['physical_absolute_detector_prediction'] is not False:
        p.error('explicit conditional geometry required')
    with np.load(args.maps) as data: maps={k:data[k] for k in data.files}
    t=maps['time_s'];n=maps['number_density_m3'];x,y=maps['x_edges_m'],maps['y_edges_m']
    foil=maps['foil_bounds_m'];length=metadata['invariant_current_length_m']
    if not np.isfinite(n).all() or (n<0).any(): p.error('invalid densities')
    ydet=np.array(case['detector_bin_edges_m']);xc=.5*(x[:-1]+x[1:])
    regions={'far_edge':xc<=foil[0]+case['edge_band_width_m'],
             'near_edge':xc>=foil[1]-case['edge_band_width_m'],
             'center':np.abs(xc)<=case['central_band_half_width_m']}
    regions['other_interior']=~np.logical_or.reduce(list(regions.values()))
    means=[gate_mean(t,n,*gate) for gate in case['gates_s']]
    out.mkdir(parents=True);write(out/'case.json',case)
    results=[];saved={};reference=None
    for spec in case['scenarios']:
        geometry=RelayGeometry(**spec['geometry']);print('Tracing '+spec['id'],flush=True)
        k,diag=build_collection_kernel(x,y,ydet,foil,[-length/2,length/2],geometry,**case['quadrature'])
        ku,du=build_collection_kernel(x,y,ydet,foil,[-length/2,length/2],geometry,**case['quadrature'],opaque_foil=False)
        entries=[]
        for gate,mean in zip(case['gates_s'],means):
            duration=gate[1]-gate[0]
            profiles={name:np.einsum('ryx,syx,x->sr',k,mean,mask.astype(float))*duration for name,mask in regions.items()}
            unocc={name:np.einsum('ryx,syx,x->sr',ku,mean,mask.astype(float))*duration for name,mask in regions.items()}
            all_profile=sum(profiles.values());all_unocc=sum(unocc.values())
            records=[]
            for s,name in enumerate(maps['species']):
                total=float(all_profile[s].sum());unblocked=float(all_unocc[s].sum())
                records.append({'species':str(name),'photons_per_gate_for_unit_per_neutral_photon_rate_s1':total,
                    'unit_emission_rate_is_physical_prediction':False,
                    'profile':all_profile[s].tolist(),
                    'region_coefficients':{r:float(v[s].sum()) for r,v in profiles.items()},
                    'region_profiles':{r:v[s].tolist() for r,v in profiles.items()},
                    'edge_fraction':None if total==0 else float((profiles['near_edge'][s].sum()+profiles['far_edge'][s].sum())/total),
                    'near_far_profile_overlap':overlap(profiles['near_edge'][s],profiles['far_edge'][s]),
                    'combined_edge_center_profile_overlap':overlap(profiles['near_edge'][s]+profiles['far_edge'][s],profiles['center'][s]),
                    'foil_visible_fraction':None if unblocked==0 else total/unblocked,
                    'region_visible_fraction':{r:None if unocc[r][s].sum()==0 else float(v[s].sum()/unocc[r][s].sum()) for r,v in profiles.items()}})
            entries.append({'gate_s':gate,'species':records})
        results.append({'id':spec['id'],'assumptions':spec['assumption'],'geometry':spec['geometry'],
                        'kernel_diagnostics':diag,'unoccluded_diagnostics':du,'gates':entries})
        if spec['id']==case['reference_scenario']:
            reference=results[-1]
            saved={'kernel_m3':k,'unoccluded_kernel_m3':ku,'detector_bin_edges_m':ydet,'x_edges_m':x,'y_edges_m':y,'source_z_interval_m':np.array([-length/2,length/2])}
            np.savez_compressed(out/'collection_kernel.npz',**saved)
    if reference is None: raise ValueError('missing reference scenario')
    radius=.0254/2;distance=case['scenarios'][0]['geometry']['target_to_first_m']
    omega=2*np.pi*(1-distance/np.hypot(distance,radius));fraction=omega/(4*np.pi)
    report={'schema':'mhdlab-target-collection-report-v1','status':'conditional_geometry_and_density_weighted_unit_emission',
        'physical_detector_counts_predicted':False,'monolayer_fit_performed':False,
        'first_pupil_on_axis_point_source':{'solid_angle_sr':float(omega),'isotropic_fraction':float(fraction),
            'photons_in_3ns_per_1W_monochromatic_121_6nm_before_losses_or_slit':float(3e-9*121.6e-9/(6.62607015e-34*299792458)*fraction)},
        'response_equation':'N_intensifier(lambda,row)=sum_s gamma_s(lambda)[photons/(neutral s)]*C_s(lambda,row)[neutral s]*T_owned(lambda); gamma and T are not yet supplied',
        'state_resolved_equation':'N_intensifier_ul(row)=A_ul*integral_gate dt*sum_cells n_u(t,cell)*K(row,cell)*T_owned(lambda_ul)',
        'factor_ownership':'K contains isotropic1/(4pi),dxdy and accepted z integration and pupil/slit/speccone geometry once. No extra volume, depth, etendue or4pi. Spectral profiles/bin fractions and measured wavelength losses remain separate.',
        'source_limit':'Historical independent assumed neutral pools used only for density-weighted unit-rate coefficients; not a common-water-monolayer CR source or inferred emissivity.',
        'edge_definition':'Geometric current-position regions: outer0.5mm of each width side plus exterior; center|x|<=0.5mm; remaining interior separate. These are not packet origin labels.',
        'time_interpolation':'Exact gate integration of declared piecewise-linear saved-density interpolant; no population extrapolation.',
        'detector_plane_status':case['detector_coordinates'],
        'unmodeled_response':case['spectral_response'],'reference_scenario':case['reference_scenario'],
        'scenarios':results,'input_sha256':{str(path):sha(path) for path in [args.maps,args.metadata,args.case,ROOT/case['geometry_evidence']]},
        'owned_source_sha256':{name:sha(ROOT/name) for name in ['mhdlab/target_optics.py','tools/compute_target_collection.py']}}
    make_plots(out,report,saved,maps['species'])
    report['outputs_sha256']={path.name:sha(path) for path in out.iterdir() if path.is_file()}
    write(out/'collection_report.json',report)
    print(json.dumps({'report':str(out/'collection_report.json'),'first_pupil_fraction':fraction}))

def make_plots(out,report,saved,species):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    ref=next(r for r in report['scenarios'] if r['id']==report['reference_scenario'])
    ym=.5*(saved['detector_bin_edges_m'][:-1]+saved['detector_bin_edges_m'][1:])*1e3
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for ax,row in zip(axes.flat,ref['gates'][-1]['species']):
        for name,color in [('near_edge','tab:red'),('far_edge','tab:orange'),('center','tab:blue'),('other_interior','0.5')]:
            values=np.array(row['region_profiles'][name]);ax.plot(ym,np.where(values>0,values,np.nan),label=name.replace('_',' '),color=color)
        ax.set(title=row['species'],xlabel='Conditional slit-image spatial coordinate [mm]',ylabel='Gate photons for rate 1 s$^{-1}$ per neutral',yscale='log')
        ax.yaxis.label.set_size(9)
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Conditional 57–60 ns collection: edge and center overlap\nUser distances + focus-compatible relay scenario; no CR emissivity or throughput fit')
    fig.savefig(out/'edge_collection_overlap.png',dpi=150);plt.close(fig)
    x,y=saved['x_edges_m'],saved['y_edges_m'];volume=np.diff(y)[:,None]*np.diff(x)[None,:]*np.diff(saved['source_z_interval_m'])[0]
    fig,axes=plt.subplots(2,1,figsize=(11,5.8),layout='constrained')
    for ax,key,title in zip(axes,['unoccluded_kernel_m3','kernel_m3'],['Without foil shadow','Opaque foil shadow included']):
        probability=saved[key].sum(axis=0)/volume
        display=np.ma.masked_where(probability<=0,np.log10(np.where(probability>0,probability,1)))
        image=ax.pcolormesh(x*1e3,y*1e3,display,vmin=-8,vmax=-3,cmap='inferno',shading='flat')
        ax.add_patch(Rectangle((-2.5,2.5),5,.2,color='.5'))
        ax.set(xlabel='Foil-width optical depth x [mm]',ylabel='Gap normal y [mm]',title=title,aspect='equal')
        fig.colorbar(image,ax=ax,label='log10 accepted photon fraction')
    fig.suptitle('Finite pupils, slit, spectrometer f/4.5 and depth defocus\nConditional reference geometry; wavelength losses excluded')
    fig.savefig(out/'collection_acceptance.png',dpi=150);plt.close(fig)

if __name__=='__main__': main()
