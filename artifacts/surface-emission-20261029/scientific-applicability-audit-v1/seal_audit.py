from pathlib import Path
import json,hashlib,os,datetime,re
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
out=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-emission-scientific-audit-20261029')
W=out.parent
roots={
 'source-v1':Path('/home/treefor/projects/runtimes/pic-thermal-field-emission-20261029/snapshots/usable-api-v1'),
 'phase-guard':W/'pic-emission-applicability-v2-20261029/source-backup-applicability-v2',
 'tangent':W/'pic-magnetized-emission-20261029/source-backup-tangential-onepass-v1/tangential-onepass',
 'surface-v3':W/'pic-electrode-heat-source-v3-20261029/source-snapshot-v3',
 'oblique':W/'pic-oblique-phase-map-20261029/source-backup-full-SN-coordinates-and-classical-bound-v1'}
spec=[
 ('P1','Original v1 leaves unsupported outgoing phases',[('source-v1','src/micro_layer.cpp',68),('source-v1','src/micro_layer.cpp',77),('phase-guard','phase_guard.hpp',12)],['S1']),
 ('P1','Return accommodation does not qualify reciprocal quantum scattering',[('source-v1','src/micro_layer.cpp',127),('surface-v3','trial_source_binding.hpp',98),('tangent','tangential_joint.cpp',176),('tangent','tangential_joint.cpp',200),('tangent','tangential_joint.cpp',203)],['S5']),
 ('P1','Strict tangent and oblique diagnostics do not qualify a full oblique emitter',[('tangent','tangential_joint.cpp',106),('tangent','tangential_joint.cpp',113),('tangent','tangential_joint.cpp',115),('tangent','tangential_joint.cpp',151),('oblique','phase_map.py',55),('oblique','phase_map.py',68),('oblique','phase_map.py',74),('oblique','phase_map.py',86),('oblique','dividing_surface_bound.py',68),('oblique','dividing_surface_bound.py',78),('oblique','dividing_surface_bound.py',81)],['S3','S4']),
 ('P1','Captured finite-offset and reconstructed-temperature data are not coupled boundary inputs',[('source-v1','src/surface_emission.cpp',163),('surface-v3','axial_trial_field.hpp',12),('surface-v3','axial_trial_field.hpp',17),('surface-v3','axial_trial_field.hpp',18),('surface-v3','trial_source_binding.hpp',58),('surface-v3','trial_source_binding.hpp',91),('phase-guard','extract_frozen_temperature.py',21),('phase-guard','extract_frozen_temperature.py',42)],['S1','S2']),
 ('P1','Numerical domain is not material, surface or reservoir admission',[('source-v1','src/surface_emission.cpp',21),('source-v1','src/surface_emission.cpp',106),('source-v1','src/surface_emission.cpp',168),('source-v1','src/surface_emission.cpp',169),('surface-v3','trial_source_binding.hpp',56)],['S6']),
 ('P1','Held accepted temperature and component closure are not complete source midpoint qualification',[('surface-v3','trial_source_binding.hpp',53),('surface-v3','trial_source_binding.hpp',59),('surface-v3','trial_source_binding.hpp',63),('surface-v3','trial_source_binding.hpp',88),('surface-v3','trial_source_binding.hpp',98),('surface-v3','axial_trial_field.hpp',19),('surface-v3','axial_trial_field.hpp',27),('surface-v3','axial_trial_field.hpp',30)],[]),
 ('P2','Preserve unified source, thermal limit and exactly-once energy/current owners',[('source-v1','src/surface_emission.cpp',105),('source-v1','src/surface_emission.cpp',132),('source-v1','src/surface_emission.cpp',144),('source-v1','src/surface_emission.cpp',153),('source-v1','src/surface_emission.cpp',195),('tangent','tangential_joint.cpp',151),('tangent','tangential_joint.cpp',157),('surface-v3','trial_source_binding.hpp',63),('surface-v3','trial_source_binding.hpp',86),('surface-v3','trial_source_binding.hpp',103)],['S1','S2'])]
findings=[]
for i,(severity,title,refs,sources) in enumerate(spec,1):
    code=[]
    for alias,rel,line in refs:
        p=roots[alias]/rel;lines=p.read_text().splitlines();assert line<=len(lines)
        code.append({'path':str(p),'line':line,'file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'line_text':lines[line-1]})
    findings.append({'id':f'E{i}','severity':severity,'title':title,'code_references':code,'primary_source_ids':sources,'remediation':f'See finding {i} in scientific-applicability-audit.txt; no model or full coupled test authorized by this audit'})
verified={}
for name,root in [('source-v1',roots['source-v1']),('phase-guard',roots['phase-guard']),('tangent',roots['tangent'].parent),('oblique',roots['oblique'])]:
    p=root/'SOURCE-SHA256.json';manifest=json.loads(p.read_text())
    for rel,sha in manifest.items():assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==sha
    verified[name]={'manifest_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'verified_files':len(manifest)}
v3seal=W/'pic-electrode-heat-source-v3-20261029/source-seal-v3.json';v3=json.loads(v3seal.read_text())
for f in v3['files']:assert hashlib.sha256((roots['surface-v3']/f['path']).read_bytes()).hexdigest()==f['sha256']
verified['surface-v3']={'manifest_sha256':hashlib.sha256(v3seal.read_bytes()).hexdigest(),'verified_files':len(v3['files'])}
# Frozen userspace owns a private PID namespace. Read the separately recorded
# host observation; never substitute a different process with namespace PID 3.
production=json.loads((out/'production-observation.json').read_text())
events=Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text()
assert production['state']=='T' and production['startticks']==5288339 and 'frozen 1' in events and 'populated 1' in events
primary={
 'S1':{'title':'Murphy and Good, Thermionic Emission, Field Emission, and the Transition Region (1956)','url':'https://doi.org/10.1103/PhysRev.102.1464','access':'publisher abstract and primary-paper excerpt'},
 'S2':{'title':'Kyritsakis and Djurabekova, A general computational method for electron emission and thermal effects in field emitting nanotips (2017)','url':'https://arxiv.org/html/1609.02364','access':'author full HTML, equations 7 and 8 and barrier discussion'},
 'S3':{'title':'Fertig and Halperin, Transmission coefficient through a saddle-point potential in a magnetic field (1987)','url':'https://doi.org/10.1103/PhysRevB.36.7969','access':'publisher abstract, quadratic scope only'},
 'S4':{'title':'Buttiker, Quantized transmission of a saddle-point constriction (1990)','url':'https://doi.org/10.1103/PhysRevB.41.7906','access':'publisher abstract'},
 'S5':{'title':'Buttiker, Four-Terminal Phase-Coherent Conductance (1986)','url':'https://doi.org/10.1103/PhysRevLett.57.1761','access':'publisher abstract, Onsager-Casimir principle'},
 'S6':{'title':'Forbes, Refining the application of Fowler-Nordheim theory (1999)','url':'https://doi.org/10.1016/S0304-3991(99)00097-2','access':'primary publisher search excerpt; full page blocked'}}
record={'schema':'emission-scientific-applicability-audit-v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'outcome':'NOT_QUALIFIED_FOR_FULL_PHYSICAL_SOURCE_INTEGRATION','findings':findings,'severity_counts':{'P0':0,'P1':6,'P2':1},'primary_sources':primary,'verified_sealed_inputs':verified,'pushed_captured_map_commit':'22f32422dbefc2d35710a1d2f8afb3e44cfa75ad','saved_oblique_manifest':verified['oblique']['manifest_sha256'],'native_model_changed':False,'new_full_coupled_test_run':False,'GPU_work':False,'production':production,'report_sha256':hashlib.sha256((out/'scientific-applicability-audit.txt').read_bytes()).hexdigest()}
(out/'audit.json').write_text(json.dumps(record,indent=2)+'\n')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()}
(out/'SOURCE-SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'outcome':record['outcome'],'severity_counts':record['severity_counts'],'code_references':sum(len(f['code_references']) for f in findings),'production_state':production['state'],'production_frozen':True,'all_input_seals_verified':True,'audit_manifest_sha256':hashlib.sha256((out/'SOURCE-SHA256.json').read_bytes()).hexdigest()}))
