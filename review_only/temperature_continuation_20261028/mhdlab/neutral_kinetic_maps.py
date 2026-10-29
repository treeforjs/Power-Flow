"""Stream local 3V neutral kinetic moments from saved ballistic packets.

The scalar temperature is the trace of the drift-subtracted velocity covariance
times m/(3 kB), not launch/wall temperature or a claim of local equilibrium.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np

from .constants import KB
from .neutral_vacuum_map import SPECIES, packet_loss_times

SCHEMA = 'mhdlab-neutral-kinetic-temperature-v1'


def _real_array(value, name, rank):
    if np.ma.isMaskedArray(value):
        raise ValueError(f'{name} cannot be masked')
    raw=np.asarray(value)
    if raw.ndim!=rank or raw.dtype.kind not in 'iuf':
        raise ValueError(f'{name} requires a real rank-{rank} array')
    array=np.asarray(raw,dtype=float)
    if not np.isfinite(array).all(): raise ValueError(f'{name} must be finite')
    return array


def _threshold(value,name,*,positive=False):
    if (isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,float,np.integer,np.floating))
            or not math.isfinite(value) or value<0 or positive and value==0):
        raise ValueError(f'{name} must be finite and {"positive" if positive else "nonnegative"}')
    return float(value)


def _grid(x_edges_m,y_edges_m,foil_bounds_m):
    x=_real_array(x_edges_m,'x_edges_m',1)
    y=_real_array(y_edges_m,'y_edges_m',1)
    foil=_real_array(foil_bounds_m,'foil_bounds_m',1)
    if any(len(a)<2 or np.any(np.diff(a)<=0) for a in (x,y)):
        raise ValueError('grid edges must strictly increase')
    if (foil.shape!=(4,) or not x[0]<foil[0]<foil[1]<x[-1]
            or not y[0]<foil[2]<foil[3]<y[-1]):
        raise ValueError('grid must strictly contain the foil rectangle')
    for boundary,edges in ((foil[0],x),(foil[1],x),(foil[2],y),(foil[3],y)):
        if not np.any(edges==boundary): raise ValueError('grid must align exactly with foil boundaries')
    solid=((.5*(y[:-1]+y[1:]))[:,None]>foil[2]) & ((.5*(y[:-1]+y[1:]))[:,None]<foil[3])
    solid=solid & ((.5*(x[:-1]+x[1:]))[None,:]>foil[0]) & ((.5*(x[:-1]+x[1:]))[None,:]<foil[1])
    return x,y,foil,solid


def _packets(position_m,velocity_m_s,particle_weight,species_index,mass_kg):
    p=_real_array(position_m,'position_m',2)
    v=_real_array(velocity_m_s,'velocity_m_s',2)
    w=_real_array(particle_weight,'particle_weight',1)
    mass=_real_array(mass_kg,'mass_kg',1)
    species=np.asarray(species_index)
    n=len(w)
    if (p.shape!=(n,3) or v.shape!=p.shape or mass.shape!=(n,)
            or species.shape!=(n,) or species.dtype.kind not in 'iu'
            or np.any((species<0)|(species>=4)) or np.any(w<=0) or np.any(mass<=0)):
        raise ValueError('packet arrays need consistent sizes, canonical species and positive weights/masses')
    masses=np.full(4,np.nan)
    for s in range(4):
        selected=mass[species==s]
        if selected.size:
            if not np.all(selected==selected[0]):
                raise ValueError('each species requires one stable particle mass')
            masses[s]=selected[0]
    return p,v,w,species,mass,masses


def coarsen_aligned_edges(edges,boundaries,factor):
    """Coarsen each vacuum/solid segment separately, retaining every boundary."""
    if isinstance(factor,(bool,np.bool_)) or not isinstance(factor,(int,np.integer)) or factor<1:
        raise ValueError('coarsen_factor must be a positive integer')
    edges=_real_array(edges,'edges',1)
    if len(edges)<2 or np.any(np.diff(edges)<=0): raise ValueError('edges must strictly increase')
    anchors=[0,len(edges)-1]
    for value in boundaries:
        found=np.flatnonzero(edges==value)
        if len(found)!=1: raise ValueError('protected boundary must be an exact grid edge')
        anchors.append(int(found[0]))
    anchors=sorted(set(anchors));indices=[]
    for lo,hi in zip(anchors,anchors[1:]):
        indices.extend(range(lo,hi,int(factor)))
    indices.append(anchors[-1])
    return edges[indices]


def binned_kinetic_moments(*,position_m,velocity_m_s,particle_weight,species_index,
        mass_kg,x_edges_m,y_edges_m,foil_bounds_m,invariant_current_length_m,
        minimum_effective_sample_size=50.,minimum_expected_cell_particles=1.):
    """Compute one frame's weighted, locally centered 3V population moments.

    Velocities are centered in two passes around a local anchor. Per-bin weight
    normalization avoids squaring raw tiny macro-particle weights for Neff.
    No Bessel correction or equilibrium closure is applied.
    """
    x,y,foil,solid=_grid(x_edges_m,y_edges_m,foil_bounds_m)
    length=_threshold(invariant_current_length_m,'invariant_current_length_m',positive=True)
    minimum_neff=_threshold(minimum_effective_sample_size,'minimum_effective_sample_size',positive=True)
    minimum_particles=_threshold(minimum_expected_cell_particles,'minimum_expected_cell_particles')
    p,v,w,species,mass,masses=_packets(position_m,velocity_m_s,particle_weight,species_index,mass_kg)
    ny,nx=len(y)-1,len(x)-1
    shape=(4,ny,nx);nbins=4*ny*nx
    points=p[:,:2].copy()
    # Match observe_vacuum_packets exactly: adjust histogram ownership only.
    for axis in (0,1):
        other=1-axis
        face_span=(points[:,other]>=foil[2*other]) & (points[:,other]<=foil[2*other+1])
        for endpoint,direction in ((foil[2*axis],-np.inf),(foil[2*axis+1],np.inf)):
            at_face=face_span & (points[:,axis]==endpoint)
            points[at_face,axis]=np.nextafter(endpoint,direction)
    ix=np.searchsorted(x,points[:,0],side='right')-1
    iy=np.searchsorted(y,points[:,1],side='right')-1
    ix[points[:,0]==x[-1]]=nx-1;iy[points[:,1]==y[-1]]=ny-1
    if np.any((ix<0)|(ix>=nx)|(iy<0)|(iy>=ny)):
        raise ValueError('live packet lies outside the observation grid')
    if np.any(solid[iy,ix]): raise ValueError('a live packet entered a solid cell')
    bins=(species*ny+iy)*nx+ix
    maximum=np.zeros(nbins)
    np.maximum.at(maximum,bins,w)
    a=w/maximum[bins]
    sum_a=np.bincount(bins,weights=a,minlength=nbins)
    sum_a2=np.bincount(bins,weights=a*a,minlength=nbins)
    populated=sum_a>0
    neff=np.zeros(nbins)
    np.divide(sum_a*sum_a,sum_a2,out=neff,where=sum_a2>0)
    number=maximum*sum_a
    if not np.isfinite(number).all(): raise ValueError('cell number exceeds finite range')
    first=np.full(nbins,len(w),dtype=np.int64)
    np.minimum.at(first,bins,np.arange(len(w),dtype=np.int64))
    anchor=np.zeros((nbins,3))
    anchor[populated]=v[first[populated]]
    relative=v-anchor[bins]
    mean_relative=np.zeros((nbins,3))
    for axis in range(3):
        np.divide(np.bincount(bins,weights=a*relative[:,axis],minlength=nbins),sum_a,
                  out=mean_relative[:,axis],where=populated)
    centered=relative-mean_relative[bins]
    variance=np.zeros((nbins,3))
    for axis in range(3):
        np.divide(np.bincount(bins,weights=a*centered[:,axis]**2,minlength=nbins),sum_a,
                  out=variance[:,axis],where=populated)
    mean=anchor+mean_relative
    mean[~populated]=np.nan
    components=variance.reshape(*shape,3)*masses[:,None,None,None]/KB
    supported=(neff.reshape(shape)>=minimum_neff) & (number.reshape(shape)>=minimum_particles) & ~solid[None]
    components[~supported]=np.nan
    temperature=np.mean(components,axis=-1)
    density=number.reshape(shape)/(np.diff(y)[:,None]*np.diff(x)[None,:]*length)
    if not np.isfinite(density).all(): raise ValueError('number density exceeds finite range')
    return {'translational_temperature_k':temperature,'temperature_components_k':components,
        'effective_sample_size':neff.reshape(shape),'number_density_m3':density,
        'mean_velocity_m_s':mean.reshape(*shape,3),'expected_cell_particles':number.reshape(shape)}


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def compute_neutral_temperature(packets_path,maps_path,metadata_path,output_dir,*,
        coarsen_factor=1,minimum_effective_sample_size=50.,minimum_expected_cell_particles=1.):
    """Re-observe saved trajectories and write one HDF frame at a time.

    Packet arrays are resident, but no full time-by-species map cube is loaded
    or allocated. Saved map header coordinates/times define the observation.
    """
    output=Path(output_dir)
    if output.exists(): raise ValueError('use a fresh output directory')
    metadata=json.loads(Path(metadata_path).read_bytes())
    for key in ('time_registration','model_description'):
        if not isinstance(metadata.get(key),str) or not metadata[key].strip():
            raise ValueError(f'input metadata requires {key}')
    length=_threshold(metadata.get('invariant_current_length_m'),'invariant_current_length_m',positive=True)
    minimum_neff=_threshold(minimum_effective_sample_size,'minimum_effective_sample_size',positive=True)
    minimum_particles=_threshold(minimum_expected_cell_particles,'minimum_expected_cell_particles')
    with np.load(maps_path,allow_pickle=False) as maps:
        times=_real_array(maps['time_s'],'time_s',1)
        original_x,original_y,foil,_=_grid(maps['x_edges_m'],maps['y_edges_m'],maps['foil_bounds_m'])
        if tuple(maps['species'].tolist())!=SPECIES: raise ValueError('maps species must be H,H2,H2O,OH')
    if len(times)<2 or np.any(np.diff(times)<=0): raise ValueError('times must strictly increase')
    x=coarsen_aligned_edges(original_x,foil[:2],coarsen_factor)
    y=coarsen_aligned_edges(original_y,foil[2:],coarsen_factor)
    with np.load(packets_path,allow_pickle=False) as packets:
        birth=_real_array(packets['birth_time_s'],'birth_time_s',1)
        p,v,w,species,mass,masses=_packets(*(packets[key] for key in
            ('position_m','velocity_m_s','particle_weight','species_index','mass_kg')))
    if birth.shape!=w.shape: raise ValueError('birth times must match packets')
    exits,hits=packet_loss_times(p,v,[x[0],x[-1],y[0],y[-1]],foil)
    loss=np.minimum(exits,hits)
    definition={'scalar_temperature':'m/(3*kB) times the trace of local weighted centered3V velocity covariance',
        'component_temperature':'m/kB times each local weighted centered velocity variance',
        'bulk_drift_removed':True,'velocity_components':['x','y','z'],
        'variance_convention':'weighted population moment; no Bessel correction or equilibrium assumption',
        'weight_normalization':'per cell and species a_i=w_i/max(w); Neff=(sum a)^2/sum(a^2)',
        'mean_algorithm':'local velocity anchor, first pass weighted centered mean, second pass centered variance',
        'support':'NaN if Neff below threshold, expected physical particle count below threshold, empty or solid',
        'zero_kelvin':'valid zero variance remains0K when support thresholds pass',
        'trajectory':'same saved births and packet_loss_times absorption/exit rules as neutral_vacuum_map',
        'birth_bin_ownership':'gas-side nextafter only; physical trajectories unchanged',
        'minimum_effective_sample_size':minimum_neff,'minimum_expected_cell_particles':minimum_particles}
    output.mkdir(parents=True)
    hdf_path=output/'neutral_temperature.h5'
    shape=(len(times),4,len(y)-1,len(x)-1)
    chunks=(1,1,min(64,shape[2]),min(128,shape[3]))
    supported_counts=[]
    with h5py.File(hdf_path,'w') as hdf:
        hdf.attrs['schema']=SCHEMA;hdf.attrs['status']='provisional'
        hdf.attrs['definition_metadata_json']=json.dumps(definition,allow_nan=False)
        for key,value in (('time_s',times),('x_edges_m',x),('y_edges_m',y),('foil_bounds_m',foil)):
            hdf.create_dataset(key,data=value)
        hdf.create_dataset('species',data=np.asarray(SPECIES,dtype=object),dtype=h5py.string_dtype('utf-8'))
        datasets={}
        for key in ('translational_temperature_k','effective_sample_size','number_density_m3'):
            datasets[key]=hdf.create_dataset(key,shape=shape,dtype='f8',chunks=chunks,compression='gzip',compression_opts=1)
        for key in ('mean_velocity_m_s','temperature_components_k'):
            datasets[key]=hdf.create_dataset(key,shape=(*shape,3),dtype='f8',chunks=(*chunks,3),compression='gzip',compression_opts=1)
        for frame,stamp in enumerate(times):
            age=stamp-birth
            live=(age>=0)&(age<loss)
            points=p[live]+age[live,None]*v[live]
            moments=binned_kinetic_moments(position_m=points,velocity_m_s=v[live],particle_weight=w[live],
                species_index=species[live],mass_kg=mass[live],x_edges_m=x,y_edges_m=y,foil_bounds_m=foil,
                invariant_current_length_m=length,minimum_effective_sample_size=minimum_neff,
                minimum_expected_cell_particles=minimum_particles)
            for key,dataset in datasets.items(): dataset[frame]=moments[key]
            supported_counts.append(np.sum(np.isfinite(moments['translational_temperature_k']),axis=(1,2)).tolist())
            if frame%10==0: print(f'Kinetic temperature frame{frame+1}/{len(times)}',flush=True)
    assumptions=list(metadata.get('assumptions',[]))
    assumptions += ['Drift-subtracted kinetic second moments of the saved collisionless packets; no local thermal equilibrium is inferred.',
        'No launch-temperature substitution, temperature floor, velocity boost, missing-data interpolation or bulk-drift heating.',
        f'Temperature unsupported cells are NaN; required Neff>={minimum_neff:g} and expected particles/cell>={minimum_particles:g}.']
    result={'schema':SCHEMA,'status':'provisional','input_status':metadata.get('status'),
        'case_title':'Neutral3V kinetic temperature | local bulk drift removed',
        'time_registration':metadata['time_registration'],'model_description':'Framewise local centered3V kinetic moments from saved neutral trajectories',
        'assumptions':assumptions,'invariant_current_length_m':length,'normal_axis_magnification':1.,
        'minimum_effective_sample_size':minimum_neff,'minimum_expected_cell_particles':minimum_particles,
        'definition':definition,'coarsen_factor':int(coarsen_factor),'shape':list(shape),
        'original_grid_shape':[len(original_y)-1,len(original_x)-1],
        'species_mass_kg':{name:float(value) if np.isfinite(value) else None for name,value in zip(SPECIES,masses)},
        'supported_cell_counts':supported_counts,'input_sha256':{str(path):_sha(path) for path in (packets_path,maps_path,metadata_path)},
        'hdf_file':hdf_path.name,'hdf_sha256':_sha(hdf_path),'producer_source_sha256':_sha(__file__),
        'new_heating_or_transport_simulation_performed':False}
    (output/'metadata.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result
