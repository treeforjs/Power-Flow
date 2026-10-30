from pathlib import Path
import bisect
import csv
import datetime
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import sys

sys.dont_write_bytecode=True
ROOT=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile')
OUT=ROOT/'particle-motion-audit'
SOURCE=ROOT/'cpu-independent-review/checkpoint-validation/compare_saved_pic.py'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
with open(ROOT/'testing-core.lock','a') as lease:
    try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:pass
    else:raise RuntimeError('Testing-core flock lease is required')
spec=importlib.util.spec_from_file_location('saved_pic_endpoint_reader',SOURCE)
reader=importlib.util.module_from_spec(spec);sys.modules[spec.name]=reader;spec.loader.exec_module(reader)
schema=reader.schema()
start=reader.read_state(ROOT/'seed-current',schema)
end=reader.read_state(ROOT/'baseline-long-10ps',schema)
assert start.seed_info==end.seed_info
assert reader.sha(start.checkpoint)==start.seed_info['checkpoint_sha256']
dt=end.ledger['time_s']-start.ledger['time_s'];assert dt>0 and abs(dt-10e-12)<1e-20
for key in ['nx','ny','width_m','gap_m','depth_m','finite_foils','finite_axial_extent','short_z_m','feed_z_m']:
    assert start.options[key]==end.options[key]
initial={p['id']:p for p in start.particles};final={p['id']:p for p in end.particles}
common=sorted(initial.keys()&final.keys());lost=sorted(initial.keys()-final.keys());born=sorted(final.keys()-initial.keys())
assert all(identifier>=start.next_id for identifier in born)
assert all(initial[i][key]==final[i][key] for i in common for key in ['weight','charge_C','mass_kg'])
QUANTILES=[0,.1,.25,.5,.75,.9,.95,.99,1]
def stats(values,weights):
    if not values:return {'count':0,'represented_weight':0}
    total=math.fsum(weights);pairs=sorted(zip(values,weights));cumulative=[];running=0.
    for value,weight in pairs:running+=weight;cumulative.append(running)
    weighted={str(q):pairs[min(bisect.bisect_left(cumulative,q*total),len(pairs)-1)][0] for q in QUANTILES}
    ordered=sorted(values)
    plain={str(q):ordered[max(0,min(len(ordered)-1,math.ceil(q*len(ordered))-1))] for q in QUANTILES}
    return {'count':len(values),'represented_weight':total,'weighted_mean':math.fsum(v*w for v,w in zip(values,weights))/total,
            'weighted_rms':math.sqrt(math.fsum(v*v*w for v,w in zip(values,weights))/total),
            'weighted_quantiles':weighted,'unweighted_quantiles':plain,'minimum':ordered[0],'maximum':ordered[-1]}
assert stats([1.,2.,3.],[1.,8.,1.])['weighted_quantiles']['0.5']==2.
def vector(p,name):return [p[name+'['+str(k)+']'] for k in range(3)]
def speed(p):
    u=math.hypot(*vector(p,'proper_velocity_m_s'));return u,u/math.sqrt(1+u*u/(reader.C_LIGHT**2))
def energy(p):return reader.kinetic(p)/reader.E_CHARGE
def rectangle_distance(p,foil):
    x,y,_=vector(p,'position_m');o=start.options
    return math.hypot(max(o[foil+'.left_m']-x,0,x-o[foil+'.right_m']),max(o[foil+'.bottom_m']-y,0,y-o[foil+'.top_m']))
def cell(p):
    x,y,_=vector(p,'position_m');o=start.options
    return (min(int(x/(o['width_m']/o['nx'])),o['nx']-1),min(int(y/(o['gap_m']/o['ny'])),o['ny']-1))
def midplane(p):return abs(p['position_m[2]'])<.0005
rows=[]
for identifier in common:
    a,b=initial[identifier],final[identifier];x0=vector(a,'position_m');x1=vector(b,'position_m')
    delta=[right-left for left,right in zip(x0,x1)];u0=vector(a,'proper_velocity_m_s');u1=vector(b,'proper_velocity_m_s')
    proper0,v0=speed(a);proper1,v1=speed(b);k0=energy(a);k1=energy(b)
    cosine=sum(x*y for x,y in zip(u0,u1))/(proper0*proper1) if proper0*proper1 else 1.
    rows.append({'id':identifier,'weight':a['weight'],'x0_m':x0[0],'y0_m':x0[1],'z0_m':x0[2],
                 'x1_m':x1[0],'y1_m':x1[1],'z1_m':x1[2],'dx_m':delta[0],'dy_m':delta[1],'dz_m':delta[2],
                 'displacement_3d_m':math.hypot(*delta),'displacement_xy_m':math.hypot(*delta[:2]),'abs_dz_m':abs(delta[2]),
                 'proper_speed0_m_s':proper0,'proper_speed1_m_s':proper1,'speed0_m_s':v0,'speed1_m_s':v1,
                 'kinetic0_eV':k0,'kinetic1_eV':k1,'delta_kinetic_eV':k1-k0,
                 'delta_proper_velocity_m_s':math.hypot(*(y-x for x,y in zip(u0,u1))),
                 'endpoint_velocity_angle_deg':math.degrees(math.acos(max(-1.,min(1.,cosine)))),
                 'midplane0':midplane(a),'midplane1':midplane(b),'same_xy_cell':cell(a)==cell(b),
                 'seed_upper_foil_distance_m':rectangle_distance(a,'upper_foil'),
                 'seed_lower_foil_distance_m':rectangle_distance(a,'lower_foil')})
def cohort(selected):
    weights=[r['weight'] for r in selected]
    result={name:stats([r[name] for r in selected],weights) for name in [
        'displacement_3d_m','displacement_xy_m','abs_dz_m','dx_m','dy_m','dz_m',
        'proper_speed0_m_s','proper_speed1_m_s','speed0_m_s','speed1_m_s','kinetic0_eV','kinetic1_eV',
        'delta_kinetic_eV','delta_proper_velocity_m_s','endpoint_velocity_angle_deg']}
    total=math.fsum(weights)
    result['same_xy_cell']={'macro_count':sum(r['same_xy_cell'] for r in selected),
                            'weighted_fraction':math.fsum(r['weight'] for r in selected if r['same_xy_cell'])/total if total else None}
    result['displacement_threshold_fractions']={str(threshold):{'unweighted':sum(r['displacement_3d_m']>threshold for r in selected)/len(selected) if selected else None,
        'weighted':math.fsum(r['weight'] for r in selected if r['displacement_3d_m']>threshold)/total if total else None} for threshold in [1e-9,1e-7,1e-6,25e-6]}
    return result
cohorts={'all_survivors':rows,'seed_midplane_survivors':[r for r in rows if r['midplane0']],
         'both_endpoints_midplane':[r for r in rows if r['midplane0'] and r['midplane1']],
         'seed_within_100um_upper_foil':[r for r in rows if r['seed_upper_foil_distance_m']<=100e-6],
         'seed_within_100um_lower_foil':[r for r in rows if r['seed_lower_foil_distance_m']<=100e-6]}
def population(particles):
    weights=[p['weight'] for p in particles]
    return {'macro_count':len(particles),'represented_electrons_if_single_electrons':math.fsum(weights),
            'proper_speed_m_s':stats([speed(p)[0] for p in particles],weights),
            'physical_speed_m_s':stats([speed(p)[1] for p in particles],weights),
            'kinetic_eV':stats([energy(p) for p in particles],weights),
            'total_kinetic_J':math.fsum(reader.kinetic(p)*p['weight'] for p in particles)}
def density_change(use_midplane):
    maps=[]
    for state in [start,end]:
        bins={}
        for p in state.particles:
            if use_midplane and not midplane(p):continue
            key=cell(p);bins[key]=bins.get(key,0.)+p['weight']
        maps.append(bins)
    a,b=maps;ta,tb=math.fsum(a.values()),math.fsum(b.values());keys=a.keys()|b.keys()
    return {'start_represented_weight':ta,'end_represented_weight':tb,'start_occupied_cells':len(a),'end_occupied_cells':len(b),
            'occupied_cell_jaccard':len(a.keys()&b.keys())/len(keys) if keys else None,
            'normalized_shape_total_variation':.5*math.fsum(abs(a.get(k,0.)/ta-b.get(k,0.)/tb) for k in keys) if ta and tb else None,
            'raw_cell_weight_l1_change':math.fsum(abs(a.get(k,0.)-b.get(k,0.)) for k in keys)}
created=end.next_id-start.next_id;newborn_lost=created-len(born);assert newborn_lost>=0
single_electrons=all(p['charge_C']==-reader.E_CHARGE for p in start.particles+end.particles)
ledger_keys=['lower_injected_charge_C','upper_injected_charge_C','authored_charge_C','lower_absorbed_charge_C',
             'upper_absorbed_charge_C','short_absorbed_charge_C','escaped_charge_C','axial_escaped_charge_C',
             'lower_absorbed_energy_J','upper_absorbed_energy_J','short_absorbed_energy_J','escaped_kinetic_energy_J','injected_energy_J']
ledger_delta={key:end.ledger[key]-start.ledger[key] for key in ledger_keys}
count_summary={'seed_macro_count':len(initial),'end_macro_count':len(final),'surviving_seed_ids':len(common),
               'lost_seed_ids':len(lost),'new_ids_surviving_at_end':len(born),'created_macro_ids':created,
               'created_then_lost_macro_ids':newborn_lost,'total_macro_losses':len(lost)+newborn_lost,
               'macro_balance_check':len(initial)+created-len(lost)-newborn_lost==len(final),
               'survival_macro_fraction':len(common)/len(initial),'survival_weight_fraction':math.fsum(initial[i]['weight'] for i in common)/math.fsum(p['weight'] for p in start.particles),
               'seed_lost_represented_weight':math.fsum(initial[i]['weight'] for i in lost),
               'new_end_represented_weight':math.fsum(final[i]['weight'] for i in born),
               'all_checkpoint_particles_have_single_electron_charge':single_electrons,
               'created_id_interval':[start.next_id,end.next_id],
               'source_and_loss_charge_deltas_C':{key:value for key,value in ledger_delta.items() if key.endswith('charge_C')},
               'ledger_energy_deltas_J':{key:value for key,value in ledger_delta.items() if key.endswith('energy_J')}}
if single_electrons:
    count_summary['source_and_loss_represented_electrons_from_charge']={key:-value/reader.E_CHARGE for key,value in ledger_delta.items() if key.endswith('charge_C')}
    count_summary['source_loss_note']='Axial escape is a subset of escaped charge; do not add both. Endpoint IDs give exact macro counts; charge deltas give represented physical counts, not per-sink macro counts.'
selected=[];selection_labels={}
for label,wanted in [('seed_midplane_survivors',12),('seed_within_100um_upper_foil',12),('all_survivors',8)]:
    pool=[r for r in cohorts[label] if r['id'] not in selected]
    for j in range(min(wanted,len(pool))):
        index=round(j*(len(pool)-1)/max(1,min(wanted,len(pool))-1));identifier=pool[index]['id']
        if identifier not in selected:selected.append(identifier);selection_labels[identifier]=label
for r in rows:
    if len(selected)>=32:break
    if r['id'] not in selected:selected.append(r['id']);selection_labels[r['id']]='fill_from_all_survivors'
assert len(selected)==32
summary={'status':'read_only_endpoint_audit_complete','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'resources':{'logical_cpu':6,'nice':15,'lease_exclusion_confirmed':True,'gpu_work':False},
         'start_directory':str(start.directory),'end_directory':str(end.directory),
         'seed_identity':start.seed_info,'start_step':start.step,'end_step':end.step,'step_count':end.step-start.step,
         'start_relative_time_s':start.ledger['time_s'],'end_relative_time_s':end.ledger['time_s'],'elapsed_s':dt,'elapsed_ps':dt*1e12,
         'cell_size_m':[start.options['width_m']/start.options['nx'],start.options['gap_m']/start.options['ny']],
         'midplane_slab_rule':'abs(z)<0.0005 m, matching native snapshot()',
         'counts':count_summary,'population_start':population(start.particles),'population_end':population(end.particles),
         'cohorts':{name:cohort(selected_rows) for name,selected_rows in cohorts.items()},
         'density_endpoint_comparison':{'height_average_xy':density_change(False),'one_mm_midplane_slab':density_change(True)},
         'selected_tracking_ids':selected,'tracking_selection_labels':selection_labels,
         'quantile_definition':'Weighted inverse empirical CDF by unchanged seed particle weights; unweighted nearest-rank quantiles. Survivor motion is conditional on IDs present at both endpoints.',
         'limitations':['Endpoint displacement is a chord, not path length, gyro radius, orbit, or trajectory.','Particles born and lost between endpoints have no saved positions here.','No force-loop E/B samples are stored in these checkpoints; endpoint speeds or spatial density do not identify a cause of motion.','The plotted density is an aggregate projection/slab population; persistent density shape does not imply stationary individual particles.','Tracking selection is retrospective and survivor-conditioned, intended for trajectory inspection rather than an unbiased population estimate.']}
(OUT/'endpoint-motion.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
with (OUT/'matched-endpoints.csv').open('w',newline='') as output:
    writer=csv.DictWriter(output,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
with (OUT/'selected-persistent-ids.csv').open('w',newline='') as output:
    chosen=[dict(r,selection=selection_labels[r['id']]) for r in rows if r['id'] in selected]
    writer=csv.DictWriter(output,fieldnames=list(chosen[0]));writer.writeheader();writer.writerows(chosen)
(OUT/'selected-ids.txt').write_text(''.join(str(i)+'\n' for i in sorted(selected)))
inputs={str(path):{'sha256':reader.sha(path),'bytes':path.stat().st_size} for path in start.files+end.files+[SOURCE,SOURCE.parent/'abi-schema.json']}
(OUT/'input-manifest.json').write_text(json.dumps(inputs,indent=2)+'\n')
print('Matched',len(common),'of',len(initial),'seed IDs across',dt*1e12,'ps')
print('Source/loss macro counts',count_summary)
for name in ['all_survivors','seed_midplane_survivors','seed_within_100um_upper_foil']:
    data=summary['cohorts'][name]
    print(name,'count',data['displacement_3d_m']['count'],'weighted median/p90 displacement um',
          data['displacement_3d_m'].get('weighted_quantiles',{}).get('0.5',0)*1e6,
          data['displacement_3d_m'].get('weighted_quantiles',{}).get('0.9',0)*1e6,
          'same xy cell',data['same_xy_cell'])
