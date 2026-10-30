import array,bisect,csv,hashlib,json,math,os,pathlib,struct,sys

W=pathlib.Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs')
O=W/'pic-emission-applicability-v2-20261029';O.mkdir(parents=True,exist_ok=True)
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
metadata=json.loads((W/'heating-history-full.json').read_text())
field=json.loads((W/'cathode-emission-audit-20261029/fields/summary.json').read_text())
data=(W/'heating-history-full.bin').read_bytes()
assert hashlib.sha256(data).hexdigest()==metadata['binary_sha256']
magic,count,*sizes=struct.unpack_from('<6Q',data)
assert magic==0x544845524d303031 and sizes==metadata['patches'] and count==metadata['times']
def numbers(offset,count):
 a=array.array('d');a.frombytes(data[offset:offset+8*count])
 if sys.byteorder!='little':a.byteswap()
 return a
times=numbers(48,count);offset=48+8*count
thermal={}
for face,size in zip(metadata['face_order'],sizes):
 thermal[face]=numbers(offset,count*size);offset+=8*count*size
assert offset==len(data)
t=field['time_s'];k=max(0,min(count-2,bisect.bisect_right(times,t)-1));q=max(0.,min(1.,(t-times[k])/(times[k+1]-times[k])))
profiles={}
for face,size in zip(metadata['face_order'],sizes):
 values=thermal[face];profiles[face]=[(1-q)*values[k*size+i]+q*values[(k+1)*size+i] for i in range(size)]
probes=W/'cathode-emission-audit-20261029/fields/probes.csv'
rows=list(csv.DictReader(probes.open()))
for index,row in enumerate(rows):
 upper=row['name'].startswith('upper_');face=row['name'].split('_',1)[1]
 thermal_face={'back':'top','gap':'bottom','left':'left','right':'right'}[face]
 if face in ['back','gap']:fraction=(float(row['x_m'])-.002)/.005
 else:
  surface_y=float(row['y_m']) if upper else .009-float(row['y_m'])
  fraction=(surface_y-.007)/.0002
 assert 0<=fraction<=1
 size=len(profiles[thermal_face]);patch=min(size-1,int(fraction*size))
 row.update({'legacy_row':index,'historical_T_K':profiles[thermal_face][patch],
             'historical_thermal_face':thermal_face,'historical_thermal_patch':patch,
             'frozen_field_time_s':t,'thermal_lower_time_s':times[k],'thermal_upper_time_s':times[k+1],
             'thermal_time_fraction':q,'lower_foil_temperature_mirrored':int(not upper)})
with (O/'frozen-probes-with-historical-temperature.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
summary={'schema':'frozen-field-instantaneous-historical-temperature-binding-v1','field_time_s':t,
 'probes':len(rows),'probe_sha256':hashlib.sha256(probes.read_bytes()).hexdigest(),
 'thermal_sha256':metadata['binary_sha256'],'thermal_clock_bracket_s':[times[k],times[k+1]],'time_fraction':q,
 'temperature_scope':'time-interpolated historical material-boundary patch-average; lower foil is the established mirrored reconstruction; not accepted heat from new coupled driver',
 'spatial_mapping':'Native uniform exported patch bins, no J-weighted reservoir or multi-time averaging; probe normal offsets do not change tangential patch coordinate',
 'probes_remain_finite_offset_not_actual_one_sided_boundary':True,
 'profiles':{face:{'patches':len(values),'minimum_K':min(values),'maximum_K':max(values)} for face,values in profiles.items()},
 'all_probe_temperature_range_K':[min(r['historical_T_K'] for r in rows),max(r['historical_T_K'] for r in rows)],
 'v1_numerical_T_range_violations':sum(not (1<=r['historical_T_K']<=3000) for r in rows),
 'new_physical_source_regime_admitted':False,'production_resumed':False,'GPU_executed':False}
(O/'temperature-binding.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
