from pathlib import Path
import json,subprocess,os,sys,hashlib,datetime
P=Path(__file__).parent;R=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/reaction-transactions-20261029/native-full-time-candidate003')
def need(ok,message):
 if not ok:raise RuntimeError(message)
need(os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15,'CPU7/nice15 required')
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:need(os.environ.get(name)=='1','one thread required')
manifest=json.loads((P/'SOURCE-SHA256.json').read_text())
for path,digest in manifest['files'].items():need(hashlib.sha256((P/path).read_bytes()).hexdigest()==digest,'source seal changed: '+path)
commands=[]
def run(command,log,timeout=300,environment=None):
 commands.append(command)
 with(R/log).open('w')as out:p=subprocess.run(command,stdout=out,stderr=subprocess.STDOUT,timeout=timeout,env=environment)
 if p.returncode:print((R/log).read_text()[-16000:]);raise RuntimeError('failed '+log+' exit '+str(p.returncode))
run([sys.executable,str(P/'build-components.py')],'final-normal-build.log')
run([sys.executable,str(P/'run-candidate-tests.py')],'final-normal-components.log')
base=json.loads((R/'legacy-commands.json').read_text());C=json.loads((R/'candidate-tests-commands.json').read_text());frozen=base[0][0]
# Only this separately linked adapter has an unresolved genuine Hall receipt function.
for source,name in [('source/src/mfem_pic_hall_motion_adapter.cpp','hall-adapter'),('tests/test_hall_adapter_compile.cpp','hall-header-traits')]:
 command=list(base[0]);command[command.index(str(P/'source/src/mfem_pic_arena.cpp'))]=str(P/source);command[-1]=str(P/'build'/str(name+'.o'));run(command,name+'-compile.log')
run([frozen,'/usr/bin/nm','-u',str(P/'build/arena.o')],'core-undefined-symbols.txt')
run([frozen,'/usr/bin/nm','-u',str(P/'build/hall-adapter.o')],'adapter-undefined-symbols.txt')
need('hall_driver'not in(R/'core-undefined-symbols.txt').read_text(),'core acquired Hall dependency')
need('hall_driver' in(R/'adapter-undefined-symbols.txt').read_text(),'adapter lost genuine Hall binding')
run([frozen,'/usr/bin/ar','rcs',str(P/'build/libpic_native_candidate003.a'),str(P.parent/'native-exact-count/build/native.o'),str(P/'build/arena.o')],'core-archive.log')
run([frozen,'/usr/bin/ar','rcs',str(P/'build/libpic_hall_motion_adapter003.a'),str(P/'build/hall-adapter.o')],'adapter-archive.log')
# Instrument the changed core and both component test programs; native/MFEM stay frozen.
asan=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
for command,name in [(base[0],'arena'),(base[1],'legacy-tests'),(C[0],'candidate-tests')]:
 command=list(command);command=[('-O1'if x in['-O0','-O2']else '-g'if x=='-g0'else x)for x in command];command[command.index('-c'):command.index('-c')]=asan;command[-1]=str(P/'build'/('sanitized-'+name+'.o'));run(command,'sanitized-'+name+'-compile.log',600)
env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
for original,name,test_object,arguments in [(base[2],'legacy','legacy-tests',[str(P/'inputs/h2o-gas-20261029.json')]),(C[1],'candidate','candidate-tests',[])]:
 command=list(original);command[2:2]=asan;command=[str(P/'build/sanitized-arena.o')if x==str(P/'build/arena.o')else str(P/'build'/('sanitized-'+test_object+'.o'))if x in[str(P/'build/tests.o'),str(P/'build/candidate-tests.o')]else x for x in command];command[-1]=str(P/'build'/('test-sanitized-'+name));run(command,'sanitized-'+name+'-link.log')
 outdir=R/('sanitized-'+name);outdir.mkdir(exist_ok=True);run([frozen,command[-1],*arguments,str(outdir)],'sanitized-'+name+'-run.log',180,env)
legacy=json.loads((R/'tests/arena-tests.json').read_text());old=json.loads((R.parent/'native-exact-count/tests/arena-tests.json').read_text());san_legacy=json.loads((R/'sanitized-legacy/arena-tests.json').read_text());need(legacy==old==san_legacy,'legacy proof/assertion/ABI parity changed')
normal=json.loads((R/'component-tests/candidate003-tests.json').read_text());san=json.loads((R/'sanitized-candidate/candidate003-tests.json').read_text());need(normal==san,'normal/sanitized component results differ')
for p,digest in manifest['files'].items():need(hashlib.sha256((P/p).read_bytes()).hexdigest()==digest,'source changed during verification')
result={'status':'passed-bounded-components-only','source_identity':manifest['source_identity'],'legacy_recorded_checks':legacy['check_count'],'legacy_results_and_proofs_identical_to_exact_count':True,'legacy_normal_sanitized_identical':True,'candidate_component_checks':normal['check_count'],'candidate_normal_sanitized_identical':True,'allocation_failure_sweep':normal['allocation_failure_sweep'],'forbidden_native_solver_calls':normal['forbidden_native_solver_calls'],'forbidden_native_advance_calls':normal['forbidden_native_advance_calls'],'sanitizers':'address/undefined/leak clean for changed arena and tests; frozen native/MFEM uninstrumented','native_ABI':legacy['ABI'],'separate_Hall_adapter_header_compile':'passed; core has no Hall symbol dependency','genuine_motion_or_joint_execution':False,'completed_full_time_restart_qualified':False,'GPU':False,'resource':{'cpu':7,'nice':15,'threads':1,'jobs':1,'external_flock':str(R.parent.parent/'testing-core7.lock')},'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'final-commands.json').write_text(json.dumps(commands,indent=2)+'\n');(R/'verification-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
