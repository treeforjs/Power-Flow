import json
import math
from pathlib import Path
import numpy as np
import h5py

p = Path('runs/neutral_desorption_current_17133_perimeter_v2')
m = json.loads((p/'metadata.json').read_text())
l = json.loads((p/'particle_ledgers.json').read_text())
maps = np.load(p/'neutral_maps.npz')
old = np.load('runs/neutral_desorption_current_17133/source_packets.npz')
b = np.load(p/'source_packets.npz')
nold = len(old['birth_time_s'])
controls = {}
for key in old.files:
    controls[key] = bool(np.array_equal(old[key], b[key][:nold]))
history_controls = {}
for label in ('top', 'bottom'):
    with h5py.File(p/(label+'_release.h5')) as current, h5py.File('runs/neutral_desorption_current_17133/'+label+'_release.h5') as previous:
        history_controls[label] = {key:bool(np.array_equal(current[key][...],previous[key][...]))
            for key in current if isinstance(current[key],h5py.Dataset)}
assert all(value for face in history_controls.values() for value in face.values())
w, v, pos = b['particle_weight'], b['velocity_m_s'], b['position_m']
birth, face, sp = b['birth_time_s'], b['face_index'], b['species_index']
mass, temp = b['mass_kg'], b['temperature_k']
assert np.isfinite(w).all() and (w > 0).all()
assert np.isfinite(v).all() and np.isfinite(pos).all()
assert np.isfinite(temp).all() and (temp > 0).all()
assert maps['time_s'][-1] == 60e-9
assert np.isfinite(maps['number_density_m3']).all() and (maps['number_density_m3'] >= 0).all()
kb = 1.380649e-23
rows = []
for f, name in enumerate(m['packet_face_index_labels']):
    axis = 1 if f < 2 else 0
    coordinate = [.0027, .0025, -.0025, .0025][f]
    sign = [1, -1, -1, 1][f]
    assert np.all(pos[face == f, axis] == coordinate)
    assert np.all(sign*v[face == f, axis] > 0)
    for s, species in enumerate(['H','H2','H2O','OH']):
        select = (face == f) & (sp == s)
        weights = w[select]
        vv, tt = v[select], temp[select]
        normalized = np.average(mass[select,None]*vv**2/(kb*tt[:,None]), axis=0, weights=weights)
        expected = np.ones(3); expected[axis] = 2
        assert np.max(np.abs(normalized/expected-1)) < .08
        age = maps['time_s'][-1]-birth[select]
        rms = np.sqrt(np.average((age[:,None]*vv)**2,axis=0,weights=weights))
        rows.append({'face':name,'species':species,'launched_particles':math.fsum(map(float,weights)),
            'mean_age_s':float(np.average(age,weights=weights)),
            'mean_birth_temperature_k':float(np.average(tt,weights=weights)),
            'normalized_velocity_second_moments_xyz':normalized.tolist(),
            'rms_displacement_xyz_m':rms.tolist()})
balance=max(abs(s['particle_balance_residual'])/s['initial_surface_particles'] for t in l for s in t['species'])
hist=max(abs(s['histogram_balance_residual'])/s['initial_surface_particles'] for t in l for s in t['species'])
assert balance < 1e-12 and hist < 1e-12
omission_bound=math.fsum(r['omitted_particles_upper_bound'] for f in m['source_records'] for r in f['packet_numerical_omissions'])
omission_count=sum(r['omitted_component_count'] for f in m['source_records'] for r in f['packet_numerical_omissions'])
side=sum(r['launched_particles'] for r in rows if r['face'] in {'left','right'})
total=sum(r['launched_particles'] for r in rows)
summary={'status':'provisional','packet_count':len(w),'frame_count':len(l),
    'time_support_s':maps['time_s'][[0,-1]].tolist(),'normal_axis_magnification':m['normal_axis_magnification'],
    'broad_face_source_history_numeric_datasets_unchanged':history_controls,
    'same_broad_face_packet_realization':all(controls.values()),
    'packet_realization_note':'Explicit negligible-tail omission changes minimum-component allocation and therefore the Monte Carlo realization; original broad-face source histories are numerically unchanged.',
    'maximum_relative_particle_balance_residual':balance,
    'maximum_relative_histogram_balance_residual':hist,
    'side_fraction_of_total_launched_particles':side/total,
    'side_area_fraction_relative_to_broad_pair':sum(f['emitting_face_area_m2'] for f in m['source_records'][2:])/sum(f['emitting_face_area_m2'] for f in m['source_records'][:2]),
    'omitted_numerical_tail_component_count':omission_count,
    'omitted_particle_count_upper_bound':omission_bound,'face_species_motion':rows,'final_species':l[-1]['species']}
(p/'verification.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in {'face_species_motion','final_species'}},indent=2))
print(json.dumps([r for r in rows if r['species']=='H'],indent=2))
