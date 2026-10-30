from pathlib import Path
import os,json,hashlib,subprocess,resource
P=Path(__file__).resolve().parent;O=Path(os.environ.get('PIC_STREAMING_WRAPPER_OUTPUT','/home/treefor/sim-workspace/powerflow/outputs/pic-generic-water-streaming-wrapper-20261030'));G=P/'provider'
O.mkdir(parents=True,exist_ok=True)
(P/'build').mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,s):
 if not x:raise RuntimeError(s)
need(os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15,'CPU7/nice15 required')
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(180,180))
provenance=json.loads((P/'inputs/provenance.json').read_text())
original=json.loads((G/'inputs/frozen-source-manifest.json').read_text())
for n,want in json.loads((P/'inputs/provider-manifest.json').read_text()).items():
 need(sha(G/n)==want,'sealed generic-v6 input changed: '+n)
 if n in original:need(original[n]==want,'provider copy detached from original source manifest: '+n)
for name,key in [('streaming_coverage.hpp','stream_header_sha256'),('exact_particle_count.h','exact_count_header_sha256')]:need(sha(P/'inputs'/name)==provenance[key],'sealed stream dependency changed')
# Strip only the declared observer plumbing and compare the entire arithmetic,
# validation, traversal, work/allocation budget and sorting body byte-for-byte.
s=(P/'source/observed_candidates.hpp').read_text();a=s.index('template<class Observer>');b=s.index('\n    return out;\n}',a)+len('\n    return out;\n}');body=s[a:b]
body=body.replace('template<class Observer>\ninline std::vector<Candidate> observed_candidates(Observer& observer,','inline std::vector<Candidate> candidates(',1).replace('\n            observer.begin_pair(p,n,*local.cell,kin);','',1).replace('\n                observer.channel(ch,kin,sigma);','',1).replace('\n            observer.end_pair();','',1)
need(body==(P/'inputs/builder-body-reference.txt').read_text(),'frozen physical builder body changed')
original=(G/'inputs/collision_scheduler.hpp').read_text();a=original.index('inline std::vector<Candidate> candidates(');b=original.index('\n    return out;\n}',a)+len('\n    return out;\n}');need(original[a:b]==body,'builder reference detached from frozen original')
s=(P/'source/streaming_water_bridge.hpp').read_text();a=s.index('    Json project(double E,Memo&memo)const {');b=s.index('    struct Recorder',a);projection=s[a:b].rstrip().replace('Json project(double E,Memo&memo)const','Json coverage(double E)const',1).replace('physical_.specs()','specs_').replace('memo.value(s.channel)','data_->observable_value_SI(s.channel,E)').replace('memo.value("H2O:conditional_outer3:unmapped_2a1")','data_->observable_value_SI("H2O:conditional_outer3:unmapped_2a1",E)').replace('memo.value(id)','data_->observable_value_SI(id,E)')
original=(G/'source/generic_water_bridge.hpp').read_text();a=original.index('    Json coverage(double E)const {');b=original.index('    std::vector<pc::Channel> comparison_channels',a);need(projection==original[a:b].rstrip(),'coverage projection changed beyond memo source lookups')
predicate='        return E>=s.low_eV&&E<=s.high_eV&&!(s.model=="parent_ground_equal_pair_CM_v1"&&E>=1000);'
need(predicate in s and predicate in original,'frozen support predicate differs')
(O/'body-equivalence.json').write_text(json.dumps({'status':'passed','physical_indexed_builder_body_exact_after_three_nonnumeric_hooks':True,'coverage_projection_exact_after_diagnostic_memo_substitution':True,'frozen_builder_sha256':hashlib.sha256(body.encode()).hexdigest(),'original_generic_sha256':sha(G/'source/generic_water_bridge.hpp')},indent=2)+'\n')
base=json.loads((P/'inputs/run-binding.json').read_text());sourcehash=sha(P/'source/streaming_water_bridge.hpp');builderhash=sha(P/'source/observed_candidates.hpp')
querymanifest={'schema':'bounded-static-wrapper-test-input-v1','test_generator_sha256':sha(P/'tests/test_wrapper.cpp'),'wrapper_source_sha256':sourcehash,'observed_builder_source_sha256':builderhash,'context':'static scalar/candidate/clock/product comparisons; no native, trajectory or field evolution'}
(O/'query-input-manifest.json').write_text(json.dumps(querymanifest,sort_keys=True,indent=2)+'\n')
bindings={'dataset_sha256':base['dataset_sha256'],'provider_source_sha256':sha(G/'source/generic_water_bridge.hpp'),'kernel_manifest_sha256':base['kernel_manifest_sha256'],'kernel_binding_file_sha256':base['kernel_binding_file_sha256'],'mechanism_contract_sha256':sha(P/'inputs/mechanism-contract.json'),'missing_rate_contract_sha256':sha(G/'products/data/collisions/h2o-missing-rate-ledger-v6.json'),'query_input_sha256':sha(O/'query-input-manifest.json'),'evidence_component_sha256':sha(P/'inputs/streaming_coverage.hpp'),'exact_arithmetic_source_sha256':sha(P/'inputs/exact_particle_count.h'),'weight_measure':'projectile-carrier weight per evaluated electron-H2O pair; repeated pairs are not a unique population census'}
run={'stream_bindings':bindings,'scope_context_sha256':hashlib.sha256(b'static-no-native-state-20261030').hexdigest(),'wrapper_source_sha256':sourcehash,'observed_builder_source_sha256':builderhash};(O/'run-binding.json').write_text(json.dumps(run,indent=2)+'\n')
commands=[]
command=['/usr/bin/g++','-std=c++20','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(P/'source'),'-I'+str(P/'inputs'),'-I'+str(G/'source'),'-I'+str(G/'inputs'),'-isystem',str(G/'products/include'),str(P/'tests/test_wrapper.cpp'),'-lcrypto','-o',str(P/'build/test_wrapper')];commands.append(command)
with (O/'build.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=90)
for mode in ['small','broad']:
 command=['/usr/bin/time','-v','-o',str(O/(mode+'-resources.txt')),str(P/'build/test_wrapper'),str(G/'products/data/collisions/h2o-gas-20261029-elastic.json'),str(O/'run-binding.json'),str(O),mode];commands.append(command)
 (O/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
 with (O/(mode+'-run.log')).open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=150)
 print((O/(mode+'-run.log')).read_text(),flush=True)
(O/'build-receipt.json').write_text(json.dumps({'status':'passed','wrapper_source_sha256':sourcehash,'observed_builder_sha256':builderhash,'stream_header_sha256':sha(P/'inputs/streaming_coverage.hpp'),'executable_sha256':sha(P/'build/test_wrapper'),'cpu':7,'nice':15,'threads':1,'physical_kernel_changed':False,'native_field_or_joint_execution':False},indent=2)+'\n')
