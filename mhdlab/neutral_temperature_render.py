"""Stream saved neutral kinetic-temperature maps into a verified GIF.

The input writer owns velocity moments, statistical support and physical scope.
This renderer preserves finite temperatures, unsupported NaNs and timestamps;
it supplies neither a temperature closure nor a missing estimate.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform

import h5py
import numpy as np

from .neutral_transport_render import SPECIES, _finite_array


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def _text(value,name):
    if not isinstance(value,str) or not value.strip():
        raise ValueError(f'{name} requires explicit text')
    return value


def _json(raw):
    def pairs(items):
        result={}
        for key,value in items:
            if key in result: raise ValueError(f'duplicate JSON key: {key}')
            result[key]=value
        return result
    def invalid(value): raise ValueError(f'nonfinite JSON constant: {value}')
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=invalid)


def _positive(value,name,*,allow_zero=False):
    if (isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value)
            or value < 0 or (not allow_zero and value == 0)):
        raise ValueError(f'{name} requires a finite {"nonnegative" if allow_zero else "positive"} value')
    return float(value)


def _plain(value):
    if isinstance(value,np.ndarray): return _plain(value.tolist())
    if isinstance(value,(list,tuple)): return [_plain(item) for item in value]
    if isinstance(value,np.generic): return _plain(value.item())
    if isinstance(value,bytes): return value.decode('utf-8')
    if isinstance(value,(str,int,bool)) or value is None: return value
    if isinstance(value,float) and math.isfinite(value): return value
    raise ValueError('HDF attributes must contain finite JSON-compatible values')


def _metadata(path):
    metadata=_json(Path(path).read_bytes())
    if metadata.get('status') not in {'provisional','manufactured','owner_supplied'}:
        raise ValueError('metadata requires provisional, manufactured or owner_supplied status')
    for name in ('model_description','time_registration'): _text(metadata.get(name),name)
    hypotheses=metadata.get('assumptions',metadata.get('source_hypotheses'))
    if not isinstance(hypotheses,list) or not hypotheses:
        raise ValueError('metadata requires assumptions or source_hypotheses')
    for value in hypotheses: _text(value,'source hypothesis')
    length=_positive(metadata.get('invariant_current_length_m'),'invariant_current_length_m')
    neff=_positive(metadata.get('minimum_effective_sample_size',
        metadata.get('min_effective_sample_size',metadata.get('min_Neff',50.))),
                   'min_effective_sample_size')
    minimum=_positive(metadata.get('minimum_expected_cell_particles',
        metadata.get('min_expected_cell_particles',1.)),'minimum_expected_cell_particles',allow_zero=True)
    if metadata.get('normal_axis_magnification',1.) != 1:
        raise ValueError('temperature views require equal physical x/y scales')
    return metadata,length,neff,minimum


def inspect_temperature_maps(hdf_path,metadata_path):
    """Validate framewise and return small headers/statistics, never a 4D cube."""
    metadata,length,min_neff,min_particles=_metadata(metadata_path)
    required=('time_s','x_edges_m','y_edges_m','species','foil_bounds_m',
              'translational_temperature_k','effective_sample_size','number_density_m3')
    with h5py.File(hdf_path,'r') as source:
        missing=set(required)-set(source)
        if missing: raise ValueError(f'missing kinetic-temperature datasets: {sorted(missing)}')
        headers={name:_finite_array(source[name][...],name,1)
                 for name in ('time_s','x_edges_m','y_edges_m','foil_bounds_m')}
        times,x,y,foil=(headers[name] for name in ('time_s','x_edges_m','y_edges_m','foil_bounds_m'))
        for name,values in (('time_s',times),('x_edges_m',x),('y_edges_m',y)):
            if len(values)<2 or np.any(np.diff(values)<=0):
                raise ValueError(f'{name} must strictly increase')
        if (foil.shape!=(4,) or foil[0]>=foil[1] or foil[2]>=foil[3]
                or foil[0]<x[0] or foil[1]>x[-1] or foil[2]<y[0] or foil[3]>y[-1]):
            raise ValueError('foil rectangle must lie inside the full saved map')
        for boundary,edges in ((foil[0],x),(foil[1],x),(foil[2],y),(foil[3],y)):
            if not np.any(edges==boundary): raise ValueError('grid must align with foil boundaries')
        if h5py.check_string_dtype(source['species'].dtype) is None:
            raise ValueError('species must be the ordered UTF8 string vector')
        species=source['species'].asstr()[...]
        if species.shape!=(4,) or tuple(species)!=SPECIES:
            raise ValueError('species must be H,H2,H2O,OH in that order')
        shape=(len(times),4,len(y)-1,len(x)-1)
        for name in required[-3:]:
            if source[name].shape!=shape or source[name].dtype.kind not in 'iuf':
                raise ValueError(f'{name} must have shape(time,4,normal,width) and real numeric values')
        attrs={key:_plain(value) for key,value in source.attrs.items()}
        if attrs.get('status')!=metadata['status']:
            raise ValueError('HDF and sidecar status must agree')
        if 'definition_metadata_json' not in attrs:
            raise ValueError('HDF requires definition_metadata_json')
        definition=_json(attrs['definition_metadata_json'])
        if not isinstance(definition,dict) or not definition:
            raise ValueError('kinetic temperature definition must be a nonempty JSON object')
        volumes=np.diff(y)[:,None]*np.diff(x)[None,:]*length
        solid=((.5*(y[:-1]+y[1:]))[:,None]>foil[2]) & ((.5*(y[:-1]+y[1:]))[:,None]<foil[3])
        solid=solid & ((.5*(x[:-1]+x[1:]))[None,:]>foil[0]) & ((.5*(x[:-1]+x[1:]))[None,:]<foil[1])
        counts=np.zeros((len(times),4),dtype=np.int64)
        ranges=[{'species':name,'valid_min_k':None,'valid_max_k':None} for name in SPECIES]
        for frame in range(len(times)):
            temperature=np.asarray(source['translational_temperature_k'][frame],dtype=float)
            effective=np.asarray(source['effective_sample_size'][frame],dtype=float)
            density=np.asarray(source['number_density_m3'][frame],dtype=float)
            valid=np.isfinite(temperature)
            if np.isinf(temperature).any() or np.any(temperature[valid]<0):
                raise ValueError('temperature permits unsupported NaN, but no infinity or negative finite values')
            if np.isinf(effective).any() or np.any(effective[np.isfinite(effective)]<0):
                raise ValueError('effective sample size must be nonnegative or unsupported NaN')
            if not np.isfinite(density).all() or np.any(density<0):
                raise ValueError('number density must be finite and nonnegative')
            supported=(np.isfinite(effective) & (effective>=min_neff)
                       & (density*volumes[None]>=min_particles*(1-64*np.finfo(float).eps)) & ~solid[None])
            if np.any(valid & ~supported):
                raise ValueError('finite temperature violates declared statistical/particle/solid support')
            counts[frame]=np.sum(valid,axis=(1,2))
            for species_index,record in enumerate(ranges):
                values=temperature[species_index][valid[species_index]]
                if values.size:
                    low,high=float(values.min()),float(values.max())
                    record['valid_min_k']=low if record['valid_min_k'] is None else min(record['valid_min_k'],low)
                    record['valid_max_k']=high if record['valid_max_k'] is None else max(record['valid_max_k'],high)
    return headers,metadata,{'shape':list(shape),'hdf_attributes':attrs,'definition':definition,
        'valid_cell_counts':counts.tolist(),'species_temperature_ranges':ranges,
        'minimum_effective_sample_size':min_neff,'minimum_expected_cell_particles':min_particles,
        'invariant_current_length_m':length}


def render_temperature(hdf_path,metadata_path,output_dir,*,fps=12.,dpi=110,
                       view_bounds_m=None,temperature_max_k=None):
    """Render supplied estimates, keeping unsupported cells masked and times exact."""
    from PIL import Image,__version__ as pillow_version
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator

    hdf_path,metadata_path,output=(Path(value).resolve() for value in (hdf_path,metadata_path,output_dir))
    if output.exists(): raise ValueError('use a new render output directory')
    fps=_positive(fps,'fps')
    if fps>100: raise ValueError('fps must be <=100')
    if type(dpi) is not int or not 60<=dpi<=240: raise ValueError('dpi must be an integer from60to240')
    if temperature_max_k is not None: temperature_max_k=_positive(temperature_max_k,'temperature_max_k')
    hashes={'hdf':_sha(hdf_path),'metadata':_sha(metadata_path)}
    headers,metadata,summary=inspect_temperature_maps(hdf_path,metadata_path)
    times,x,y,foil=(headers[name] for name in ('time_s','x_edges_m','y_edges_m','foil_bounds_m'))
    full=np.array([x[0],x[-1],y[0],y[-1]])
    view=full.copy() if view_bounds_m is None else _finite_array(view_bounds_m,'view_bounds_m',1)
    if (view.shape!=(4,) or view[0]>=view[1] or view[2]>=view[3]
            or view[0]<full[0] or view[1]>full[1] or view[2]<full[2] or view[3]>full[3]):
        raise ValueError('view_bounds_m must be ordered inside the full saved map')
    cropped=not np.array_equal(view,full)
    highs=[record['valid_max_k'] for record in summary['species_temperature_ranges'] if record['valid_max_k'] is not None]
    data_max=max(highs) if highs else None
    upper=temperature_max_k if temperature_max_k is not None else (data_max if data_max and data_max>0 else 1.)
    cmap=plt.get_cmap('inferno').with_extremes(bad='white')
    palette=Image.new('P',(1,1))
    colors=[tuple(round(255*channel) for channel in cmap(i/223)[:3]) for i in range(224)]
    colors.extend([(round(255*i/31),)*3 for i in range(32)])
    palette.putpalette([channel for color in colors for channel in color])
    fig,axes=plt.subplots(2,2,figsize=(10.4,8.6),layout='constrained')
    fig.set_dpi(dpi)
    title=metadata.get('temperature_case_title','Neutral H / H2 / H2O / OH')
    suffix=' | edge view' if cropped else ''
    fig.suptitle(f'{title}{suffix}\nKinetic translational T: bulk drift removed — {metadata["status"].upper()}\nCollisionless; not equilibrium or internal state T',fontsize=11.5)
    artists=[];annotations=[]
    norm=Normalize(0.,upper)
    for species_index,ax in enumerate(axes.flat):
        artist=ax.pcolormesh(x*1000,y*1000,np.ma.masked_all((len(y)-1,len(x)-1)),
            shading='flat',cmap=cmap,norm=norm,rasterized=True)
        artists.append(artist)
        ax.add_patch(Rectangle((foil[0]*1000,foil[2]*1000),(foil[1]-foil[0])*1000,
            (foil[3]-foil[2])*1000,facecolor='.35',edgecolor='black',linewidth=.8,zorder=5))
        ax.set(xlabel='Foil width x [mm]',ylabel='Thickness / vacuum normal y [mm]',
            xlim=view[:2]*1000,ylim=view[2:]*1000,title=SPECIES[species_index])
        ax.set_aspect(1.,adjustable='box')
        ax.xaxis.set_major_locator(MaxNLocator(5));ax.yaxis.set_major_locator(MaxNLocator(5))
        annotations.append(ax.text(.5,.95,'',transform=ax.transAxes,ha='center',va='top',fontsize=8,
            bbox={'facecolor':'white','edgecolor':'.8','alpha':.9},zorder=8))
        if data_max is not None:
            colorbar=fig.colorbar(artist,ax=ax,fraction=.046,pad=.035,
                                 extend='max' if data_max>upper else 'neither')
            colorbar.set_label('Kinetic translational T [K]',fontsize=9)
    maximum=float(np.max(np.abs(times)))
    unit,multiplier=('ns',1e9) if maximum<1e-6 else ('µs',1e6) if maximum<1e-3 else ('ms',1e3) if maximum<1 else ('s',1.)
    footer=fig.supxlabel('',fontsize=9)
    frames=[];over_counts=np.zeros((len(times),4),dtype=np.int64)
    with h5py.File(hdf_path,'r') as source:
        for frame,stamp in enumerate(times):
            values=np.asarray(source['translational_temperature_k'][frame],dtype=float)
            for species_index,artist in enumerate(artists):
                valid=np.isfinite(values[species_index])
                artist.set_array(np.ma.array(values[species_index],mask=~valid))
                annotations[species_index].set_text('' if valid.any() else 'No supported estimates at this time')
                annotations[species_index].set_visible(not valid.any())
                over_counts[frame,species_index]=int(np.count_nonzero(valid & (values[species_index]>upper)))
            footer.set_text(f't = {stamp*multiplier:.5g} {unit}    |    frame {frame+1}/{len(times)}\n'
                'White: unsupported estimate   •   Gray: solid foil   •   Equal x/y scale\n'
                f'Effective N ≥ {summary["minimum_effective_sample_size"]:g}; expected particles/cell ≥ {summary["minimum_expected_cell_particles"]:g}; no time interpolation')
            fig.canvas.draw()
            rgb=np.asarray(fig.canvas.buffer_rgba())[...,:3].copy()
            frames.append(Image.fromarray(rgb).quantize(palette=palette,dither=Image.Dither.NONE))
    plt.close(fig)
    output.mkdir(parents=True)
    gif=output/'neutral_temperature.gif'
    dwell=max(10,int(round(100/fps))*10)
    frames[0].save(gif,save_all=True,append_images=frames[1:],duration=dwell,loop=0,disposal=2,optimize=False)
    targets={'preview_first.png':0,'preview_middle.png':len(times)//2,'preview_last.png':len(times)-1}
    durations=[]
    with Image.open(gif) as decoded:
        frame_count=decoded.n_frames
        if frame_count!=len(times): raise ValueError('encoded GIF failed to retain every source time')
        dimensions=list(decoded.size)
        for frame in range(frame_count):
            decoded.seek(frame);durations.append(decoded.info.get('duration',0))
            for name,index in targets.items():
                if frame==index: decoded.convert('RGB').save(output/name)
    if hashes!={'hdf':_sha(hdf_path),'metadata':_sha(metadata_path)}:
        raise ValueError('source kinetic maps changed during rendering')
    report={'schema':'mhdlab-neutral-temperature-gif-v1','status':metadata['status'],
        'input_hdf':str(hdf_path),'input_metadata':str(metadata_path),'input_sha256':hashes,
        'source_metadata':metadata,'source_summary':summary,
        'frame_count':frame_count,'dimensions_pixels':dimensions,'physical_time_s':times.tolist(),
        'requested_fps':fps,'frame_duration_ms':durations,'total_animation_duration_ms':sum(durations),
        'color_scale':{'kind':'shared_fixed_linear_kelvin','limits_k':[0.,upper],
            'maximum_supported_temperature_k':data_max,'requested_maximum_k':temperature_max_k,
            'all_supported_values_zero':data_max==0.,'bins_above_color_range':over_counts.tolist()},
        'temperature_values_changed':False,'unsupported_nans_replaced':False,'temperature_floor_applied':False,
        'time_interpolation':False,'spatial_interpolation':False,'geometry_aspect':'equal_millimetre_units',
        'normal_axis_magnification':1.,'full_map_bounds_m':full.tolist(),'display_view_bounds_m':view.tolist(),
        'view_crops_axes_only':cropped,'data_access':'HDF frames streamed; no 4D temperature cube loaded',
        'definition_label':'kinetic translational T, species bulk drift removed; collisionless; not equilibrium/internal state T',
        'physical_model_qualified_by_renderer':False,'palette':'fixed224inferno+32grays,no dithering',
        'runtime':{'python':platform.python_version(),'numpy':np.__version__,'h5py':h5py.__version__,
                   'matplotlib':matplotlib.__version__,'pillow':pillow_version},
        'renderer_source_sha256':_sha(__file__),
        'output_sha256':{p.name:_sha(p) for p in sorted(output.iterdir()) if p.is_file()}}
    (output/'render_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report
