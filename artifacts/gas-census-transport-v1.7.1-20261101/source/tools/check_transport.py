from pathlib import Path
import os,sys,json,hashlib,subprocess,time
if not __debug__:raise RuntimeError('optimized Python refused')
R=Path(__file__).resolve().parents[1];BASE=R.parent/'gas-target-owner-v1.6-reservoir-release-20261101'
F=R.parent/'arch-base-20261028';FORMS=R.parent/'pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix/forms/include/pf_v7'
mode,label=sys.argv[1:];assert mode in ('normal','sanitized')
D=Path('/home/treefor/projects/.workspace/outputs')/R.name/label;D.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
aff=sorted(os.sched_getaffinity(0));assert len(aff)==1 and os.getenv('SIMQ_CPU_LIST')==str(aff[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow' and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text()
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')=='' and not os.getenv('PF_BUILD_ISA_FLAGS','')
assert int(os.environ['SIMQ_MEM'])>=3*1024**3
rs=Path('/').stat();fs=(F/'root').stat();assert (rs.st_dev,rs.st_ino)==(fs.st_dev,fs.st_ino)
assert sha('/usr/bin/python3')=='375532f83be620283ab780c14ca6bcc5907dc379d2b621de0c5d0a112cf5458b'
def pins():
 return {str(p):sha(p) for base,dirs in ((R,('include','src','device','tests','tools','integration')),(BASE,('include','src')),(FORMS,('',))) for d in dirs for p in (base/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts}
before=pins();receipt={'passed':False,'job':os.environ['SIMQ_JOB_ID'],'mode':mode,'source_inputs':before,'affinity':aff,'frozen_root_identity':[rs.st_dev,rs.st_ino],'production_qualified':False,'device_qualified':False,'runs':[]}
env=dict(os.environ)
if mode=='sanitized':env.update(ASAN_OPTIONS='detect_leaks=1:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
def run(cmd,name):
 cmd=list(map(str,cmd));t=time.monotonic()
 with (D/(name+'.log')).open('w')as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=180)
 receipt['runs'].append({'command':cmd,'rc':p.returncode,'wall_s':time.monotonic()-t,'log_sha256':sha(D/(name+'.log'))})
 print(name,p.returncode,flush=True)
 if p.returncode:raise RuntimeError((D/(name+'.log')).read_text()[-6000:])
try:
 flags=['-std=c++20','-O1','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math']
 if mode=='sanitized':flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-sanitize-recover=all']
 common=['/usr/bin/c++',*flags,'-I'+str(R/'include'),'-I'+str(FORMS)]
 run(common+['-c',R/'src/gas_owner.cpp','-o',D/'gas.o'],'compile-gas')
 run(common+[R/'tests/test_census_transport.cpp',D/'gas.o','-o',D/'controls'],'compile-controls')
 run([D/'controls'],'controls');receipt['result']=json.loads((D/'controls.log').read_text());assert receipt['result']['passed']
 run([D/'controls'],'replay');assert (D/'controls.log').read_bytes()==(D/'replay.log').read_bytes()
 run(common+['-I'+str(R/'integration/include'),'-I'+str(R/'integration/tests'),R/'tests/test_guards.cpp',R/'tests/allocation_probe.cpp',D/'gas.o','-o',D/'guards'],'compile-inherited')
 run([D/'guards'],'inherited-guards')
 if mode=='normal':
  run(['/usr/bin/c++',*flags,'-I'+str(BASE/'include'),'-I'+str(FORMS),R/'tests/test_transport_preimage.cpp',BASE/'src/gas_owner.cpp','-o',D/'preimage'],'compile-preimage')
  run([D/'preimage'],'preimage-v1.6')
 assert before==pins();receipt['inputs_unchanged']=True;receipt['passed']=True
finally:(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
