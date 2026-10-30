from pathlib import Path
import bisect
import csv
import datetime
import fcntl
import hashlib
import json
import math
import os
import shutil
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from PIL import Image

ROOT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile')
AUDIT=ROOT/'particle-motion-audit'
OUT=AUDIT/'trajectory-diagnostic-10ps'
RUN=AUDIT/'tracking-run-10ps'
LOG=RUN/'selected-particle-forces.csv'
COMPARISON=RUN/'full-comparison.json'
SCHEMA=ROOT/'cpu-independent-review/checkpoint-validation/abi-schema.json'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
with open(ROOT/'testing-core.lock','a') as lease:
    try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:pass
    else:raise RuntimeError('Testing-core lease is required')
OUT.mkdir(exist_ok=True);(OUT/'frames').mkdir(exist_ok=True)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save_json(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
input_paths=[LOG,COMPARISON,SCHEMA,AUDIT/'selected-ids.txt',AUDIT/'selected-persistent-ids.csv',RUN/'progress.json']
before_hashes={str(path):sha(path) for path in input_paths}
comparison=json.loads(COMPARISON.read_text())
assert comparison['status']=='passed' and comparison['comparison_complete'] and not comparison['failures']
assert comparison['paths']['baseline']==str(AUDIT/'baseline-track-10ps')
schema=json.loads(SCHEMA.read_text())
assert (schema['endian'],schema['options_size'],schema['diagnostics_size'],schema['particle_size'])==('little',168,200,80)
ids=sorted(int(value) for value in (AUDIT/'selected-ids.txt').read_text().split());assert len(ids)==len(set(ids))==32
index={identifier:i for i,identifier in enumerate(ids)}
def checkpoint(directory):
    meta=directory/'restart.meta';fields=meta.read_text().splitlines()[0].split();path=directory/fields[0]
    data=path.read_bytes();assert struct.unpack_from('<Q',data)[0]==0x5241494b50494332
    options={}
    for item in schema['options']:
        fmt={'f64':'d','i32':'i','u64':'Q','bool':'B'}[item['kind']]
        options[item['name']]=struct.unpack_from('<'+fmt,data,8+item['offset'])[0]
    elapsed=struct.unpack_from('<d',data,176)[0];next_id,count=struct.unpack_from('<QQ',data,376)
    assert len(data)==392+80*count
    selected={}
    for row in struct.iter_unpack('<Q9d',data[392:]):
        identifier=row[0]
        if identifier in index:selected[identifier]={'position':np.array(row[1:4]),'u':np.array(row[4:7]),'q':row[7],'m':row[8],'weight':row[9]}
    assert set(selected)==set(ids)
    input_paths.extend([meta,path]);before_hashes.update({str(meta):sha(meta),str(path):sha(path)})
    return {'options':options,'elapsed':elapsed,'step':int(fields[1]),'selected':selected,'checkpoint':str(path),'next_id':next_id}
seed=checkpoint(ROOT/'seed-current');endpoint=checkpoint(RUN)
options=seed['options'];assert options['finite_foils']==1 and options['finite_axial_extent']==1
assert options==endpoint['options']
nsteps=endpoint['step']-seed['step'];assert nsteps>0
nblocks=1+2*nsteps
positions=np.full((nblocks,32,3),np.nan);velocities=np.full_like(positions,np.nan)
times=np.full(nblocks,np.nan);dts=np.full(nblocks,np.nan)
electric=np.full((nsteps,32,3),np.nan);magnetic=np.full_like(electric,np.nan);u_before=np.full_like(electric,np.nan)
seen=np.zeros((nblocks,32),dtype=np.uint8)
expected_header=['event','time_relative_s','dt_s','id','present','x_m','y_m','z_m','Ex_V_m','Ey_V_m','Ez_V_m','Bx_T','By_T','Bz_T','ux_before_m_s','uy_before_m_s','uz_before_m_s','ux_after_m_s','uy_after_m_s','uz_after_m_s']
count=0
with LOG.open(newline='') as source:
    reader=csv.reader(source);assert next(reader)==expected_header
    for row_number,row in enumerate(reader):
        assert len(row)==20
        block=row_number//32;assert block<nblocks
        event='initial' if block==0 else ('force_midpoint' if block%2 else 'step_end')
        assert row[0]==event and row[4]=='1','Unexpected event, missing survivor or event ordering'
        identifier=int(row[3]);assert identifier in index;column=index[identifier]
        assert seen[block,column]==0;seen[block,column]=1
        time=float(row[1]);dt=float(row[2]);assert math.isfinite(time) and math.isfinite(dt) and dt>0
        if math.isnan(times[block]):times[block]=time;dts[block]=dt
        else:assert times[block]==time and dts[block]==dt
        coordinates=np.array([float(value) for value in row[5:8]])
        after=np.array([float(value) for value in row[17:20]])
        assert np.isfinite(coordinates).all() and np.isfinite(after).all()
        positions[block,column]=coordinates;velocities[block,column]=after
        if event=='force_midpoint':
            values=np.array([float(value) for value in row[8:17]]);assert np.isfinite(values).all()
            step=(block-1)//2;electric[step,column]=values[:3];magnetic[step,column]=values[3:6];u_before[step,column]=values[6:]
        else:assert row[8:17]==['']*9,'State rows must not invent force samples'
        count+=1
assert count==nblocks*32 and np.all(seen==1)
assert np.isfinite(positions).all() and np.isfinite(velocities).all()
assert times[0]==seed['elapsed'] and times[-1]==endpoint['elapsed']
for identifier,column in index.items():
    assert np.array_equal(positions[0,column],seed['selected'][identifier]['position'])
    assert np.array_equal(velocities[0,column],seed['selected'][identifier]['u'])
    assert np.array_equal(positions[-1,column],endpoint['selected'][identifier]['position'])
    assert np.array_equal(velocities[-1,column],endpoint['selected'][identifier]['u'])
    assert all(seed['selected'][identifier][key]==endpoint['selected'][identifier][key] for key in ['q','m','weight'])
previous_positions=positions[0:-2:2];mid_positions=positions[1::2];end_positions=positions[2::2]
previous_u=velocities[0:-2:2];force_u=velocities[1::2];end_u=velocities[2::2]
step_dt=dts[2::2];previous_times=times[0:-2:2]
assert np.array_equal(dts[1::2],step_dt) and dts[0]==step_dt[0]
assert np.array_equal(previous_u,u_before) and np.array_equal(force_u,end_u)
time_error_mid=float(np.max(np.abs(times[1::2]-(previous_times+.5*step_dt))))
time_error_end=float(np.max(np.abs(times[2::2]-(previous_times+step_dt))))
assert max(time_error_mid,time_error_end)<=4*np.spacing(float(np.max(times)))
C=299792458.
def norm2(value):return value[...,0]*value[...,0]+value[...,1]*value[...,1]+value[...,2]*value[...,2]
def gamma(value):return np.sqrt(1+norm2(value)/(C*C))
def drift(position,u):
    half=.5*step_dt[:,None];g=gamma(u);result=np.empty_like(position)
    result[...,:2]=position[...,:2]+(u[...,:2]/g[...,None])*half[...,None]
    result[...,2]=position[...,2]+half*u[...,2]/g
    return result
first_error=float(np.max(np.abs(mid_positions-drift(previous_positions,u_before))))
second_error=float(np.max(np.abs(end_positions-drift(mid_positions,force_u))))
roundoff_limit=64*np.finfo(float).eps*max(options['width_m'],options['gap_m'],options['depth_m'])
assert max(first_error,second_error)<=roundoff_limit
assert np.all((positions[...,0]>=0)&(positions[...,0]<=options['width_m']))
assert np.all((positions[...,1]>=0)&(positions[...,1]<=options['gap_m']))
assert np.all((positions[...,2]>=options['short_z_m'])&(positions[...,2]<=options['feed_z_m']))
# An open rectangle intersection would contradict straight, surviving vacuum drift.
segments_a=positions[:-1,:,:2].reshape(-1,2);segments_b=positions[1:,:,:2].reshape(-1,2)
rectangle_intersections={}
for foil in ['lower_foil','upper_foil']:
    lower=np.array([options[foil+'.left_m'],options[foil+'.bottom_m']]);upper=np.array([options[foil+'.right_m'],options[foil+'.top_m']])
    delta=segments_b-segments_a;enter=np.zeros(len(delta));leave=np.ones(len(delta));possible=np.ones(len(delta),dtype=bool)
    for axis in range(2):
        moving=delta[:,axis]!=0;possible&=moving|((segments_a[:,axis]>lower[axis])&(segments_a[:,axis]<upper[axis]))
        low=np.full(len(delta),-np.inf);high=np.full(len(delta),np.inf)
        a=np.zeros(len(delta));b=np.zeros(len(delta));a[moving]=(lower[axis]-segments_a[moving,axis])/delta[moving,axis];b[moving]=(upper[axis]-segments_a[moving,axis])/delta[moving,axis]
        low[moving]=np.minimum(a[moving],b[moving]);high[moving]=np.maximum(a[moving],b[moving]);enter=np.maximum(enter,low);leave=np.minimum(leave,high)
    midpoint=segments_a+(.5*(enter+leave))[:,None]*delta
    inside=np.all((midpoint>lower)&(midpoint<upper),axis=1)
    hits=possible&(enter<leave)&inside
    rectangle_intersections[foil]=int(np.count_nonzero(hits));assert not np.any(hits)
def cross(a,b):
    return np.stack([a[...,1]*b[...,2]-a[...,2]*b[...,1],a[...,2]*b[...,0]-a[...,0]*b[...,2],a[...,0]*b[...,1]-a[...,1]*b[...,0]],axis=-1)
qm=np.array([seed['selected'][identifier]['q']/seed['selected'][identifier]['m'] for identifier in ids])
scale=.5*step_dt[:,None]*qm[None,:]
minus=u_before+scale[...,None]*electric;t=scale[...,None]*magnetic/gamma(minus)[...,None]
prime=minus+cross(minus,t);s=2*t/(1+norm2(t))[...,None]
expected_after=(minus+cross(prime,s))+scale[...,None]*electric
boris_error=float(np.max(np.abs(expected_after-force_u)))
boris_limit=256*np.finfo(float).eps*max(1.,float(np.max(np.abs(force_u))),float(np.max(np.abs(u_before))))
assert boris_error<=boris_limit
segment_lengths=np.linalg.norm(np.diff(positions,axis=0),axis=2)
path_lengths=np.sum(segment_lengths,axis=0);endpoint_lengths=np.linalg.norm(positions[-1]-positions[0],axis=1)
speed_before=np.sqrt(norm2(u_before))/gamma(u_before);speed_after=np.sqrt(norm2(force_u))/gamma(force_u)
velocity_lengths=np.sum(.5*step_dt[:,None]*(speed_before+speed_after),axis=0)
path_velocity_error=float(np.max(np.abs(path_lengths-velocity_lengths)))
assert np.all(path_lengths+roundoff_limit>=endpoint_lengths)
assert path_velocity_error<1e-11*max(1e-6,float(np.max(path_lengths)))
progress=json.loads((RUN/'progress.json').read_text());origin=progress['start_time_s']
start_absolute_ns=(origin+times[0])*1e9;end_absolute_ns=(origin+times[-1])*1e9
elapsed_ps=(times-times[0])*1e12
assert abs(elapsed_ps[-1]-10)<1e-8
step_blocks=np.arange(0,nblocks,2)
frame_blocks=[]
for goal in np.linspace(times[0],times[-1],56):
    j=bisect.bisect_left(times[step_blocks],goal);candidates=[min(j,len(step_blocks)-1),max(j-1,0)]
    pick=min(candidates,key=lambda k:abs(times[step_blocks[k]]-goal));frame_blocks.append(int(step_blocks[pick]))
assert len(frame_blocks)==len(set(frame_blocks))==56 and frame_blocks[0]==0 and frame_blocks[-1]==nblocks-1
validation={'status':'passed','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'selected_ids':ids,'id_count':32,'native_steps':nsteps,'rows':count,'initial_rows':32,'force_rows':32*nsteps,'step_end_rows':32*nsteps,
    'all_ids_present_at_all_records':True,'finite_positions_velocities_and_force_samples':True,'state_force_columns_blank':True,
    'initial_and_final_records_exactly_match_checkpoints':True,'times_monotonic':bool(np.all(np.diff(times)>0)),
    'maximum_midpoint_clock_error_s':time_error_mid,'maximum_step_clock_error_s':time_error_end,
    'pre_force_equals_previous_end_velocity_exactly':True,'post_force_equals_step_end_velocity_exactly':True,
    'finite_foils_absorb_instead_of_reflect':True,'finite_axial_extent_absorbing':True,'foil_interior_segment_intersections':rectangle_intersections,
    'maximum_first_half_drift_position_error_m':first_error,'maximum_second_half_drift_position_error_m':second_error,
    'straight_drift_roundoff_check_limit_m':roundoff_limit,'maximum_logged_Boris_reconstruction_error_m_s':boris_error,'Boris_reconstruction_roundoff_limit_m_s':boris_limit,
    'all_selected_half_drifts_verified_straight_without_hidden_reflection_or_clamp':True,
    'numerical_pathlength_definition':'Sum of all logged initial/midpoint/end straight-drift segments. This is discrete numerical model distance, not a continuum true orbit.',
    'maximum_pathlength_minus_integrated_logged_speed_m':path_velocity_error,
    'retrospective_sample_median_numerical_path_um':float(np.median(path_lengths))*1e6,
    'retrospective_sample_median_endpoint_displacement_um':float(np.median(endpoint_lengths))*1e6,
    'retrospective_sample_median_path_to_endpoint_ratio':float(np.median(path_lengths/endpoint_lengths)),
    'start_absolute_ns':start_absolute_ns,'end_absolute_ns':end_absolute_ns,'elapsed_ps':float(elapsed_ps[-1]),
    'frame_count':56,'frame_record_block_indices':frame_blocks,'frame_elapsed_ps':[float(elapsed_ps[k]) for k in frame_blocks],
    'full_comparison_status':comparison['status'],'full_comparison_position_linf_m':comparison['particles']['position_m']['linf_absolute'],
    'full_comparison_is_policy_pass_not_bitwise_claim':True,'logical_cpus':[6],'nice':15,'threads':1,'lease_exclusion_confirmed':True,'new_simulation_run':False,'gpu_work':False}
save_json(OUT/'log-validation.json',validation)
selection={int(row['id']):row['selection'] for row in csv.DictReader((AUDIT/'selected-persistent-ids.csv').open())}
colors=list(plt.colormaps['tab20'].colors)+list(plt.colormaps['tab20b'].colors[:12])
colors=[tuple(float(value*.84) if sum(color)/3>.74 else float(value) for value in color) for color in colors]
with (OUT/'selected-pathlengths-and-colors.csv').open('w',newline='') as output:
    fields=['id','selection','weight','color_hex','numerical_path_um','endpoint_displacement_um','path_to_endpoint_ratio','maximum_abs_delta_z_um']
    writer=csv.DictWriter(output,fieldnames=fields);writer.writeheader()
    for column,identifier in enumerate(ids):
        writer.writerow({'id':identifier,'selection':selection[identifier],'weight':seed['selected'][identifier]['weight'],
            'color_hex':matplotlib.colors.to_hex(colors[column]),'numerical_path_um':path_lengths[column]*1e6,
            'endpoint_displacement_um':endpoint_lengths[column]*1e6,'path_to_endpoint_ratio':path_lengths[column]/endpoint_lengths[column],
            'maximum_abs_delta_z_um':float(np.max(np.abs(positions[:,column,2]-positions[0,column,2])))*1e6})
print('Log validation passed:',count,'rows;',nsteps,'steps; all 32 IDs complete',flush=True)
print('Straight drifts verified; sample median numerical path / endpoint displacement (um):',np.median(path_lengths)*1e6,np.median(endpoint_lengths)*1e6,flush=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.titlesize':12.5,'axes.labelsize':11.5,
    'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','text.color':'#172330','axes.labelcolor':'#172330',
    'xtick.color':'#3b4854','ytick.color':'#3b4854','axes.edgecolor':'#a3adb7'})
fig=plt.figure(figsize=(12.8,7.5),dpi=120,facecolor='white')
left=fig.add_axes([.065,.23,.415,.56]);right=fig.add_axes([.565,.23,.395,.56])
fig.text(.5,.965,'Recorded motion of 32 selected surviving particles',ha='center',va='top',fontsize=19,weight='bold')
fig.text(.5,.916,f'Separate copied-state replay: {start_absolute_ns:.5f} → {end_absolute_ns:.5f} ns  |  Full production remains paused at ≈26.48 ns',ha='center',fontsize=11.5)
clock_text=fig.text(.5,.86,'',ha='center',fontsize=12.5,family='DejaVu Sans Mono',color='#263e50')
left.set_title('Absolute x–y positions and actual foils',pad=11);left.set_xlabel('x (mm)');left.set_ylabel('y (mm)')
left.set_xlim(0,options['width_m']*1e3);left.set_ylim(0,options['gap_m']*1e3);left.set_aspect('equal',adjustable='box')
left.set_xticks([0,2,4,6,8,9]);left.set_yticks([0,2,4,6,8,9]);left.grid(color='#edf0f3',linewidth=.7)
for foil,label in [('lower_foil','lower foil'),('upper_foil','upper foil')]:
    x=options[foil+'.left_m']*1e3;y=options[foil+'.bottom_m']*1e3;w=(options[foil+'.right_m']-options[foil+'.left_m'])*1e3;h=(options[foil+'.top_m']-options[foil+'.bottom_m'])*1e3
    left.add_patch(Rectangle((x,y),w,h,facecolor='#d4dbe2',edgecolor='#687686',linewidth=1.0,zorder=1))
    left.text(x+w/2,y+(.36 if foil=='upper_foil' else -.34),label,ha='center',va='center',fontsize=9.5,color='#596879')
right.set_title('Per-particle displacement in y–z',pad=11);right.set_xlabel('Δy from that ID’s start (µm)');right.set_ylabel('Δz from that ID’s start (µm)')
displacement=(positions-positions[0:1])*1e6
for axis,component in [('x',1),('y',2)]:
    low=float(np.min(displacement[...,component]));high=float(np.max(displacement[...,component]));pad=max((high-low)*.085,1.)
    (right.set_xlim if axis=='x' else right.set_ylim)(low-pad,high+pad)
right.grid(color='#edf0f3',linewidth=.7);right.axhline(0,color='#b7c0c9',linewidth=.8);right.axvline(0,color='#b7c0c9',linewidth=.8)
right.text(.02,.97,'Axes use independent scales',transform=right.transAxes,ha='left',va='top',fontsize=9,color='#5b6772')
left_lines=[];right_lines=[];left_heads=[];right_heads=[]
for column,color in enumerate(colors):
    left_lines.append(left.plot([],[],color=color,linewidth=1.05,alpha=.57,zorder=3)[0])
    right_lines.append(right.plot([],[],color=color,linewidth=1.0,alpha=.48,zorder=3)[0])
    left.scatter([positions[0,column,0]*1e3],[positions[0,column,1]*1e3],s=20,facecolors='white',edgecolors=[color],linewidths=.8,zorder=4)
    left_heads.append(left.plot([],[],marker='o',markersize=4.4,markerfacecolor=color,markeredgecolor='#26303a',markeredgewidth=.3,linestyle='None',zorder=5)[0])
    right_heads.append(right.plot([],[],marker='o',markersize=4.3,markerfacecolor=color,markeredgecolor='#26303a',markeredgewidth=.3,linestyle='None',zorder=5)[0])
legend=[Line2D([],[],marker='o',markerfacecolor='white',markeredgecolor='#667789',linestyle='None',markersize=5,label='initial'),
        Line2D([],[],marker='o',markerfacecolor='#405a72',markeredgecolor='#26303a',linestyle='None',markersize=5,label='current')]
left.legend(handles=legend,loc='lower right',frameon=False,fontsize=9,handletextpad=.4,borderpad=.2)
fig.text(.5,.125,'32 survivors selected retrospectively; ID colors stay fixed. This sample is for inspection, not a population estimate.',ha='center',fontsize=10.5)
fig.text(.5,.082,'Heads use exact recorded step states. Trails join recorded initial / native midpoint / end states; no fitted or interpolated samples.',ha='center',fontsize=9.5,color='#465564')
fig.text(.5,.044,'Numerical straight-drift polylines are not continuum true orbits or a convergence claim. No field grid or field-line overlay is used.',ha='center',fontsize=9.5,color='#465564')
png_paths=[]
for frame,block in enumerate(frame_blocks):
    for column in range(32):
        left_lines[column].set_data(positions[:block+1,column,0]*1e3,positions[:block+1,column,1]*1e3)
        right_lines[column].set_data(displacement[:block+1,column,1],displacement[:block+1,column,2])
        left_heads[column].set_data([positions[block,column,0]*1e3],[positions[block,column,1]*1e3])
        right_heads[column].set_data([displacement[block,column,1]],[displacement[block,column,2]])
    clock_text.set_text(f'Elapsed {elapsed_ps[block]:6.3f} ps   •   recorded step {block//2:4d}/{nsteps}')
    path=OUT/'frames'/f'frame_{frame:03d}.png';fig.savefig(path,dpi=120,facecolor='white');png_paths.append(path)
    if frame%10==0 or frame==55:print('Rendered exact frame',frame+1,'/ 56',flush=True)
plt.close(fig)
shutil.copy2(png_paths[-1],OUT/'trajectory-final.png')
# A fixed palette keeps ID colors and a pure-white background stable in GIF.
palette_colors=[(255,255,255),(23,35,48),(212,219,226),(104,118,134)]
base_colors=[tuple(round(value*255) for value in color) for color in colors]
palette_colors.extend(base_colors)
for alpha in [.57,.48,.25]:palette_colors.extend(tuple(round(alpha*value+(1-alpha)*255) for value in color) for color in base_colors)
palette_colors.extend((value,value,value) for value in range(0,256,3))
palette_colors=palette_colors[:256]
while len(palette_colors)<256:palette_colors.append((255,255,255))
palette=Image.new('P',(1,1));palette.putpalette([channel for color in palette_colors for channel in color])
frames=[]
for path in png_paths:
    with Image.open(path) as image:
        rgb=image.convert('RGB');assert rgb.getpixel((0,0))==(255,255,255)
        quantized=rgb.quantize(palette=palette,dither=Image.Dither.NONE)
        indices=np.asarray(quantized).copy();indices[np.all(np.asarray(rgb)==255,axis=2)]=0
        quantized=Image.fromarray(indices);quantized.putpalette(palette.getpalette())
        frames.append(quantized)
durations=[500]+[85]*54+[1300]
gif_path=OUT/'recorded-particle-trajectories-10ps.gif'
frames[0].save(gif_path,save_all=True,append_images=frames[1:],duration=durations,loop=0,disposal=2,optimize=False,background=0)
decoded=[]
with Image.open(gif_path) as gif:
    assert gif.n_frames==56
    for frame in range(56):
        gif.seek(frame);rgb=gif.convert('RGB');expected=frames[frame].convert('RGB')
        assert rgb.size==expected.size and np.array_equal(np.asarray(rgb),np.asarray(expected))
        for corner in [(0,0),(rgb.width-1,0),(0,rgb.height-1),(rgb.width-1,rgb.height-1)]:assert rgb.getpixel(corner)==(255,255,255)
        decoded.append({'frame':frame,'block':frame_blocks[frame],'elapsed_ps':float(elapsed_ps[frame_blocks[frame]]),'decoded_rgb_sha256':hashlib.sha256(rgb.tobytes()).hexdigest(),'duration_ms':gif.info.get('duration')})
for path in input_paths:assert sha(path)==before_hashes[str(path)],'Input changed during rendering'
save_json(OUT/'decoded-gif-validation.json',{'status':'passed','frame_count':56,'size_pixels':list(frames[0].size),'all_frames_equal_expected_fixed_palette_renders':True,'all_frame_corners_pure_white':True,'frames':decoded})
save_json(OUT/'input-manifest.json',{str(path):{'sha256':before_hashes[str(path)],'bytes':path.stat().st_size} for path in input_paths})
save_json(OUT/'render-metadata.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'Pillow':Image.__version__,'logical_cpus':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':1,'lease_exclusion_confirmed':True,'rendered_frames':56,'frames_are_exact_recorded_step_end_samples':True,'trail_vertices_are_all_recorded_native_samples':True,'B_overlay_used':False,'other_time_field_grids_opened':False,'new_simulation_run':False,'gpu_work':False})
print('GIF decode verification passed for all 56 frames; pure-white corners and exact palette renders',flush=True)
