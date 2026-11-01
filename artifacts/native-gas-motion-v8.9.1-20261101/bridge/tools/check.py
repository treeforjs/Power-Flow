from pathlib import Path
import os,sys,json,hashlib,subprocess,time,importlib.util
if not __debug__:raise RuntimeError('optimized Python guards refused')
R=Path(__file__).resolve().parents[1];RT=R.parent
N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.9-motion-inventory-20261101'
BASE=RT/'pic-reaction-transactions-20261029/native-full-time-v8.8-provisional-target-owner-20261101'
G=RT/'gas-target-owner-v1.7.1-census-transport-20261101';F=RT/'arch-base-20261028'
P=Path('/home/treefor/projects/Power-Flow/native/coupled_pic/include')
O=Path('/home/treefor/projects/.workspace/outputs/native-gas-common-stage-v1.3.1-motion-unwind-20261101')
mode,label=sys.argv[1:3];reuse=sys.argv[3:]==["--reuse-native"];assert not sys.argv[3:] or reuse;assert mode in ('normal','sanitized');D=O/label;D.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
aff=sorted(os.sched_getaffinity(0));assert len(aff)==1 and os.getenv('SIMQ_CPU_LIST')==str(aff[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow' and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text()
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')=='' and not os.getenv('PF_BUILD_ISA_FLAGS','')
assert int(os.environ['SIMQ_MEM']) >= (8 if mode=='sanitized' else 3)*1024**3
rs=Path('/').stat();fs=(F/'root').stat();assert (rs.st_dev,rs.st_ino)==(fs.st_dev,fs.st_ino),'frozen root required'
spec=importlib.util.spec_from_file_location('native_tooling',N/'tooling.py');native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native);native.resource_guard(build=True)
dep=json.loads((R/'DEPENDENCIES.json').read_text());assert all(sha(p)==v for p,v in dep.items())
def pins():
 result={}
 for root,dirs in ((N,['source','forms','inputs','tests']),(R,['include','tests','tools']),(G,['include','src'])):
  for d in dirs:result.update({str(p):sha(p) for p in (root/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts})
 for p in [R/'DEPENDENCIES.json',R/'NATIVE-BUILD.json',R/'NATIVE-BUILD-SANITIZED.json',*map(Path,dep),N/'tooling.py',N/'build-components.py']:result[str(p)]=sha(p)
 return result
before=pins();deps=native.sealed_inputs();receipt={'passed':False,'mode':mode,'job':os.environ['SIMQ_JOB_ID'],'affinity':aff,'source_inputs':before,'dependencies':deps,'frozen_root_identity':[rs.st_dev,rs.st_ino],'runs':[],'production_qualified':False}
env=dict(os.environ,PF_V8_CHEM='v3.1');
if mode=='sanitized':env.update(native.ASAN_ENV)
def run(cmd,name):
 cmd=list(map(str,cmd));start=time.monotonic()
 with (D/(name+'.log')).open('w')as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=900)
 receipt['runs'].append({'cmd':cmd,'rc':p.returncode,'wall_s':time.monotonic()-start,'log':name+'.log','log_sha256':sha(D/(name+'.log'))});print(name,p.returncode,flush=True)
 if p.returncode:raise RuntimeError((D/(name+'.log')).read_text()[-8000:])
try:
 if reuse:
  build=json.loads((R/('NATIVE-BUILD-SANITIZED.json' if mode=='sanitized' else 'NATIVE-BUILD.json')).read_text());assert build['mode']==mode and all(sha(N/name)==h for name,h in build['source']['files'].items()) and all(sha(p)==h for p,h in build['binaries'].items())
  receipt['reused_native_build']=build
 else:run(['/usr/bin/python3',N/'build-components.py',mode,'core','hall','v84-scope','audit-target-guard','v85-registration','v87-stocks'],'native-build')
 manifest=json.loads((N/'SOURCE-SHA256.json').read_text());assert all(sha(N/n)==v for n,v in manifest['files'].items())
 receipt['native_source']=manifest['source_identity'];tag='sanitized-' if mode=='sanitized' else ''
 commands=json.loads((native.OUT/('build-'+mode)/'build-commands.json').read_text());original=next(x['command'] for x in commands if x['log']==tag+'v87-stocks-compile.log');cc=original[original.index('/usr/bin/mpic++'):original.index('-c')]
 run(cc+['-I'+str(R/'include'),'-I'+str(R/'tests'),'-I'+str(G/'include'),'-I'+str(P),'-I'+str(N/'forms/include/pf_v7'),'-MMD','-MF',D/'common.d','-c',R/'tests/test_common_stage.cpp','-o',D/'common.o'],'compile-common')
 san=['-fsanitize=address,undefined','-fno-omit-frame-pointer'] if mode=='sanitized' else []
 run(['/usr/bin/c++','-std=c++20','-O1' if san else '-O2',*san,'-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),'-c',G/'src/gas_owner.cpp','-o',D/'gas.o'],'compile-gas')
 run(['/usr/bin/c++','-std=c++20','-O1',*san,'-Wall','-Wextra','-Werror','-c',R/'tests/allocation_probe.cpp','-o',D/'allocation.o'],'compile-allocation')
 headers={name:sha(name) for name in (D/'common.d').read_text().replace('\\\n',' ').split(':',1)[1].split()};receipt['compiler_headers']=headers
 run(['/usr/bin/mpic++',*san,'-Wl,--gc-sections',D/'common.o',D/'gas.o',D/'allocation.o',N/'build'/(tag+'kinetic-boundary.o'),N/'build'/(tag+'libpic_native_v8.a'),*native.LIBS,'-o',D/'test_common'],'link-common')
 controls=['normal','motion-empty-flows','allocation-motion-count','allocation-motion-0','allocation-motion-16','allocation-motion-64','allocation-motion-1163','snapshot-controls','empty-stream-caps','correction-exits','correction-exterior','motion-flows','same-digest-motion','empty','abandon','same-digest-stage1','same-digest-stage2','abandon-final','regular-refusal','move-final','cell-crossing','two-step','destroyed-native','allocation-seal-0','allocation-seal-16','allocation-seal-64','allocation-seal-173','allocation-seal-256','allocation-seal-312','allocation-count']
 for control in controls:run([D/'test_common',D/control,control],control)
 for name,args in [('test_v84_scope',[D/'inherited-scope']),('test_audit_target_guard',[D/'inherited-target-energy']),('test_v85_registration',[]),('test_v87_stocks',[D/'inherited-stocks'])]:run([N/'build'/(tag+name),*args],name)
 assert pins()==before and native.sealed_inputs()==deps and all(sha(p)==v for p,v in headers.items())
 receipt['results']={c:json.loads((D/c/'result.json').read_text()) for c in controls};assert all(r['passed'] for r in receipt['results'].values());receipt['binary_sha256']=sha(D/'test_common');receipt['native_archive_sha256']=sha(N/'build'/(tag+'libpic_native_v8.a'));receipt['passed']=True
finally:(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
