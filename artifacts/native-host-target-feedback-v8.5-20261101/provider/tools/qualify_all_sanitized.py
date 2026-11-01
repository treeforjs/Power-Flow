from pathlib import Path
import os,sys,json,hashlib,subprocess,time,importlib.util
if not __debug__:raise RuntimeError('optimized Python guards refused')
R=Path(__file__).resolve().parents[1];RT=R.parent
N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.5-host-target-feedback-20261101'
B=RT/'pic-reaction-transactions-20261029/native-full-time-v8.4-water-scope-20261101'
G=RT/'gas-target-owner-v1.4-binding-20261101';F=RT/'arch-base-20261028';A=RT/'native-gas-adapter-v1.2-immutable-export-20261101'
O=Path('/home/treefor/projects/.workspace/outputs/native-gas-target-feedback-v1-20261101')
mode,label=sys.argv[1:];assert mode=='sanitized';D=O/label;D.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
aff=sorted(os.sched_getaffinity(0));assert len(aff)==1 and os.getenv('SIMQ_CPU_LIST')==str(aff[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow' and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text()
assert aff[0] in (6,7) or (aff[0] in (4,5) and os.getenv('SIMQ_BORROW')=='1')
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')=='' and not os.getenv('PF_BUILD_ISA_FLAGS','')
assert int(os.environ['SIMQ_MEM'])>=3*1024**3
rs=Path('/').stat();fs=(F/'root').stat();assert (rs.st_dev,rs.st_ino)==(fs.st_dev,fs.st_ino),'actual frozen filesystem root required'
spec=importlib.util.spec_from_file_location('native_tooling',N/'tooling.py');native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native);native.resource_guard(build=True)
manifest=json.loads((N/'SOURCE-SHA256.json').read_text());assert all(sha(N/n)==v for n,v in manifest['files'].items())
assert json.loads((O.parent/'native-v8.5-host-target-feedback-20261101/build-normal/build-result.json').read_text())['source_identity']==manifest['source_identity'],'successful normal native build for this source required'
assert json.loads((O.parent/'native-v8.5-host-target-feedback-20261101/build-sanitized/build-result.json').read_text())['source_identity']==manifest['source_identity'],'successful sanitized native build required'
base=json.loads((B/'SOURCE-SHA256.json').read_text());assert base['source_identity']=='08ece5073517da1dc0832f8d90f11967c31f9b9e388af409bbd94e05058ee0f2';assert all(sha(B/n)==v for n,v in base['files'].items())
assert sha(B/'build/libpic_native_v8.a')=='de0811fc2c7fcdc4da81ff571fb5cfeb36a99510af241198f7d28eb20a2f5dde'
assert sha(G/'include/pf_gas_v1/gas_owner.hpp')=='208f2ecead8c6b56318bcae6a7162b1592c3db14df45a3988f1ac42937c43e12' and sha(G/'src/gas_owner.cpp')=='9e46282ff3578725d4f2e9454fa884ddaea1f69e63b2b8a7843b94ad3fe8d6d9'
dep=json.loads((R/'DEPENDENCIES.json').read_text());assert sha(dep['binding_header'])==dep['sha256']
def pins():
 result={str(N/n):sha(N/n) for n in manifest['files']}|{str(B/n):sha(B/n) for n in base['files']}
 for root,dirs in ((R,['include','tests','tools']),(G,['include','src'])):
  for d in dirs:
   result.update({str(p):sha(p) for p in (root/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts})
 for p in [R/'DEPENDENCIES.json',A/'include/native_gas_stage.hpp',N/'build/sanitized-libpic_native_v8.a',N/'build/sanitized-kinetic-boundary.o',B/'build/libpic_native_v8.a',B/'build/kinetic-boundary.o',N/'tests/test_v85_registration.cpp',N/'tests/chem_v8_fixture.hpp',B/'tests/chem_v8_fixture.hpp',N/'tooling.py']:result[str(p)]=sha(p)
 return result
before=pins();deps=native.sealed_inputs();receipt={'passed':False,'mode':mode,'job':os.environ['SIMQ_JOB_ID'],'affinity':aff,'source_inputs':before,'dependencies':deps,'native_source_identity':manifest['source_identity'],'frozen_root_identity':[rs.st_dev,rs.st_ino],'runs':[],'production_qualified':False,'joint_energy_closure':False,'sanitizer_scope':'provider/tests/gas plus rebuilt sanitized v8.5 native core; sealed exact-count native.o and frozen external libraries remain prebuilt'}
def run(cmd,name,env=None):
 cmd=list(map(str,cmd));start=time.monotonic()
 with (D/(name+'.log')).open('w')as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=600)
 receipt['runs'].append({'cmd':cmd,'rc':p.returncode,'wall_s':time.monotonic()-start,'log':name+'.log','log_sha256':sha(D/(name+'.log'))});print(name,p.returncode,flush=True)
 if p.returncode:raise RuntimeError((D/(name+'.log')).read_text()[-7000:])
try:
 commands=json.loads((N.parent.parent.parent/'.workspace/outputs/native-v8.5-host-target-feedback-20261101/build-normal/build-commands.json').read_text())
 original=next(x['command'] for x in commands if x['log']=='v85-registration-compile.log');cc=original[original.index('/usr/bin/mpic++'):original.index('-c')]
 san=['-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-sanitize-recover=all'] if mode=='sanitized' else []
 run(cc+san+['-I'+str(R/'include'),'-I'+str(R/'tests'),'-I'+str(A/'include'),'-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),'-MMD','-MF',D/'targets.d','-c',R/'tests/test_targets.cpp','-o',D/'targets.o'],'compile-targets')
 run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),*san,'-c',G/'src/gas_owner.cpp','-o',D/'gas.o'],'compile-gas')
 headers={name:sha(name) for name in (D/'targets.d').read_text().replace('\\\n',' ').split(':',1)[1].split()};receipt['compiler_headers']=headers
 run(['/usr/bin/mpic++',*san,'-Wl,--gc-sections',D/'targets.o',D/'gas.o',N/'build/sanitized-kinetic-boundary.o',N/'build/sanitized-libpic_native_v8.a',*native.LIBS,'-o',D/'test_targets'],'link-targets')
 env=dict(os.environ,PF_V8_CHEM='v3.1');env.update(native.ASAN_ENV if mode=='sanitized' else {})
 run([D/'test_targets',D/'controls'],'controls',env)
 if mode=='normal':
  run([N/'build/test_v85_registration'],'registration-fixed',env)
  oldcc=[str(x).replace(str(N),str(B)) for x in cc]
  run(oldcc+['-c',N/'tests/test_v85_registration.cpp','-o',D/'baseline-registration.o'],'compile-baseline-registration')
  run(['/usr/bin/mpic++','-Wl,--gc-sections',D/'baseline-registration.o',B/'build/kinetic-boundary.o',B/'build/libpic_native_v8.a',*native.LIBS,'-o',D/'baseline-registration'],'link-baseline-registration')
  run([D/'baseline-registration','baseline'],'registration-baseline',env)
 assert pins()==before and native.sealed_inputs()==deps and all(sha(p)==v for p,v in headers.items())
 receipt['result']=json.loads((D/'controls/result.json').read_text());assert receipt['result']['passed'];receipt['binary_sha256']=sha(D/'test_targets');receipt['passed']=True
finally:(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
