from pathlib import Path
import json,os,hashlib,subprocess,resource,time
P=Path(__file__).parent;O=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-streaming-chemistry-evidence-20261030');G=Path('/home/treefor/projects/runtimes/pic-generic-water-bridge-20261029')
def need(x,s):
 if not x:raise RuntimeError(s)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if (P/'inputs/frozen-source-manifest.json').is_file():
 for name,want in json.loads((P/'inputs/frozen-source-manifest.json').read_text()).items():need(sha(P/name)==want,'stream source/input changed: '+name)
need(os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15,'CPU7/nice15 required')
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(180,180))
need(sha(G/'source/generic_water_bridge.hpp')=='b16d517d0d2cb9aa5da70dde60ddf7a660804ca62018058a3295d8f99189b191','generic-v6 changed')
for name,want in json.loads((G/'inputs/frozen-source-manifest.json').read_text()).items():need(sha(G/name)==want,'generic frozen input changed: '+name)
provenance=json.loads((P/'inputs/provenance.json').read_text());need(sha(P/'inputs/exact_particle_count.h')==provenance['exact_count_header_sha256'],'exact count header changed')
base=json.loads((P/'inputs/run-binding.json').read_text());missing=G/'products/data/collisions/h2o-missing-rate-ledger-v6.json'
bindings={'dataset_sha256':base['dataset_sha256'],'provider_source_sha256':sha(G/'source/generic_water_bridge.hpp'),'kernel_manifest_sha256':base['kernel_manifest_sha256'],'kernel_binding_file_sha256':base['kernel_binding_file_sha256'],'mechanism_contract_sha256':sha(P/'inputs/mechanism-contract.json'),'missing_rate_contract_sha256':sha(missing),'query_input_sha256':'','evidence_component_sha256':sha(P/'source/streaming_coverage.hpp'),'exact_arithmetic_source_sha256':sha(P/'inputs/exact_particle_count.h'),'weight_measure':'represented-electron query weight; repeated queries are not a unique census'}
command=['/usr/bin/g++','-std=c++20','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(P/'source'),'-I'+str(P/'inputs'),'-I'+str(G/'source'),'-I'+str(G/'inputs'),'-isystem',str(G/'products/include'),str(P/'tests/test_stream.cpp'),'-lcrypto','-o',str(P/'build/test_stream')]
commands=[command]
with (O/'build.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=90)
for mode,input_file in [('small',P/'inputs/coverage.json'),('capture',Path(provenance['csv']))]:
 manifest={'schema':'stream-query-input-manifest-v1','mode':mode,'primary_input':str(input_file),'primary_input_sha256':sha(input_file),'query_generator_sha256':sha(P/'tests/test_stream.cpp'),'provider_source_sha256':bindings['provider_source_sha256'],'target_assumption':'stationary ground H2O for captured pre-kick relative energy only','weight_source':'small fixed fixture weights or captured CSV weight; no normalization','physical_evolution':False}
 manifest_path=O/(mode+'-input-manifest.json');manifest_path.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
 selected=dict(bindings);selected['query_input_sha256']=sha(manifest_path);run=O/(mode+'-binding.json');run.write_text(json.dumps({'bindings':selected},indent=2)+'\n')
 command=['/usr/bin/time','-v','-o',str(O/(mode+'-resources.txt')),str(P/'build/test_stream'),mode,str(G/'products/data/collisions/h2o-gas-20261029-elastic.json'),str(run),str(input_file),str(O)];commands.append(command)
 (O/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 with (O/(mode+'-run.log')).open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=150)
 print((O/(mode+'-run.log')).read_text(),flush=True)
(O/'build-receipt.json').write_text(json.dumps({'status':'passed','executable_sha256':sha(P/'build/test_stream'),'stream_header_sha256':sha(P/'source/streaming_coverage.hpp'),'test_driver_sha256':sha(P/'tests/test_stream.cpp'),'exact_count_header_sha256':sha(P/'inputs/exact_particle_count.h'),'source_kernel_changed':False,'cpu':7,'nice':15,'threads':1,'GPU':False},indent=2)+'\n')

# Replay is a second scalar/evidence process, never a native or physical step.
replay=O/'capture-replay';replay.mkdir(exist_ok=True)
command=['/usr/bin/time','-v','-o',str(O/'capture-replay-resources.txt'),str(P/'build/test_stream'),'capture',str(G/'products/data/collisions/h2o-gas-20261029-elastic.json'),str(O/'capture-binding.json'),provenance['csv'],str(replay)]
commands.append(command)
with (O/'capture-replay-run.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=150)
need((O/'capture-receipt.json').read_bytes()==(replay/'capture-receipt.json').read_bytes(),'captured receipt is not exactly reproducible')
(O/'replay-results.json').write_text(json.dumps({'status':'passed','capture_queries':46726,'whole_receipt_byte_exact':True,'receipt_sha256':sha(O/'capture-receipt.json'),'physical_evolution':False},indent=2)+'\n')
command=['/usr/bin/python',str(P/'independent_oracle.py')];commands.append(command)
with (O/'independent-oracle.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=60)
(O/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
print((O/'independent-oracle.log').read_text(),flush=True)
