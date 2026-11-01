import hashlib,json,os,pathlib,shutil,time
C=pathlib.Path('/home/treefor/projects/runtimes/pic-magnetized-emission-20261029')
O=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs/pic-magnetized-emission-20261029')
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
P=O/'source-backup-guard-map-v1'
assert not P.exists()
shutil.copytree(C,P,ignore=shutil.ignore_patterns('build','snapshots','__pycache__'))
shutil.copytree(O/'rate-map',P/'evidence/rate-map')
shutil.copytree(O/'interfaces',P/'interfaces')
for name in ['guard-assessment-qualification.json','analytic-guard-vs-native-predicate.csv','guard-build-and-check-v1.log','rate-map-batches.log','rate-weighted-admission-map.json','rate-weighted-admission-map.csv']:
 shutil.copy(O/name,P/'evidence'/name)
api=pathlib.Path('/home/treefor/projects/runtimes/pic-thermal-field-emission-20261029/snapshots/usable-api-v1/SOURCE-SHA256.json')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(api)=='9bd0fa5771b8cba405cc0b63b99c0e32a24da8d735eef240606c36b1fb49f0ca'
(P/'PROVENANCE.json').write_text(json.dumps({'API_v1_manifest_sha256':sha(api),'created_unix_s':time.time(),'states':3168,'nominal_guard':.2,'physical_source_admitted':False,'GPU_executed':False,'production_resumed':False},indent=2)+'\n')
(P/'PUBLICATION.txt').write_text('This source-only publication includes the complete rate map in bounded CSV shards. The redundant combined CSV is preserved in the complete outside-project backup and can be regenerated with summarize_map.py. The local hash manifest includes that combined CSV. No physical magnetized current or foil model is admitted.\n')
manifest={str(p.relative_to(P)):sha(p) for p in sorted(P.rglob('*')) if p.is_file()}
(P/'SOURCE-SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
seal={'status':'SEALED_GUARD_RATE_MAP_V1','manifest_sha256':sha(P/'SOURCE-SHA256.json'),'source_snapshot':str(P),'nominal_guard':.2,'actual_magnetized_current_known':False}
(P/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n')
shutil.make_archive(str(P),'gztar',P.parent,P.name)
items=[]
for p in sorted(P.rglob('*')):
 if p.is_file() and p.name!='rate-weighted-admission-map.csv':
  items.append({'path':'artifacts/surface-emission-20261029/magnetic-guard-map-v1/'+str(p.relative_to(P)),'mode':'100644','type':'blob','content':p.read_text()})
batches=[];current=[];size=0
for item in items:
 n=len(json.dumps(item).encode())
 if current and size+n>180000:batches.append(current);current=[];size=0
 current.append(item);size+=n
if current:batches.append(current)
for i,batch in enumerate(batches):(O/('publish-guard-map-batch-'+str(i)+'.json')).write_text(json.dumps(batch))
(O/'GUARD-MAP-V1.json').write_text(json.dumps({**seal,'publication_batches':len(batches),'publication_files':len(items)},indent=2)+'\n')
print(json.dumps({**seal,'publication_batches':len(batches),'publication_files':len(items)}))
