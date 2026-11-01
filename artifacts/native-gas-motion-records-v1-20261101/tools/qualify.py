from pathlib import Path
import os,sys,json,hashlib,subprocess,time,importlib.util
if not __debug__:raise RuntimeError('optimized Python refused')
R=Path(__file__).resolve().parents[1];RT=R.parent
N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.7-reservoir-stock-20261101'
G=RT/'gas-target-owner-v1.6-reservoir-release-20261101';F=RT/'arch-base-20261028'
OUT=Path('/home/treefor/projects/.workspace/outputs/native-gas-motion-v1-20261101')
mode,label=sys.argv[1:3];assert mode in ('normal','sanitized') and len(sys.argv)==3
D=OUT/label;D.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
aff=sorted(os.sched_getaffinity(0));assert len(aff)==1 and os.getenv('SIMQ_CPU_LIST')==str(aff[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow' and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text()
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')=='' and not os.getenv('PF_BUILD_ISA_FLAGS','')
assert int(os.environ['SIMQ_MEM'])>=3*1024**3
rs=Path('/').stat();fs=(F/'root').stat();assert (rs.st_dev,rs.st_ino)==(fs.st_dev,fs.st_ino),'frozen root required'
spec=importlib.util.spec_from_file_location('native_tooling',N/'tooling.py');native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
native.resource_guard(build=True)
def pins():
 result={}
 for root,dirs in ((R,['include','tests','tools']),(G,['include','src']),(N,['source/include','source/src','forms/include','inputs'])):
  for d in dirs:result.update({str(p):sha(p) for p in (root/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts})
 for p in [N/'tooling.py',N/'build/libpic_native_v8.a',N/'build/kinetic-boundary.o',N/'SOURCE-SHA256.json',native.OUT/'build-normal/build-commands.json']:
  result[str(p)]=sha(p)
 return result
before=pins();deps=native.sealed_inputs();receipt={'passed':False,'job':os.environ['SIMQ_JOB_ID'],'mode':mode,'affinity':aff,'source_inputs':before,'dependencies':deps,'frozen_root_identity':[rs.st_dev,rs.st_ino],'runs':[],'native_motion_authority_qualified':False,'production_qualified':False}
env=dict(os.environ,PF_V8_CHEM='v3.1')
if mode=='sanitized':env.update(native.ASAN_ENV)
def run(cmd,name):
 cmd=list(map(str,cmd));start=time.monotonic()
 with (D/(name+'.log')).open('w')as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=180)
 receipt['runs'].append({'cmd':cmd,'rc':p.returncode,'wall_s':time.monotonic()-start,'log':name+'.log','log_sha256':sha(D/(name+'.log'))})
 print(name,p.returncode,flush=True)
 if p.returncode:raise RuntimeError((D/(name+'.log')).read_text()[-6000:])
try:
 commands=json.loads((native.OUT/'build-normal/build-commands.json').read_text())
 original=next(x['command'] for x in commands if x['log']=='v87-stocks-compile.log')
 cc=original[original.index('/usr/bin/mpic++'):original.index('-c')]
 san=['-fsanitize=address,undefined','-fno-omit-frame-pointer'] if mode=='sanitized' else []
 run(cc+san+['-I'+str(R/'include'),'-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),'-MMD','-MF',D/'motion.d','-c',R/'tests/test_motion.cpp','-o',D/'motion.o'],'compile-motion')
 run(['/usr/bin/c++','-std=c++20','-O1',*san,'-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),'-c',G/'src/gas_owner.cpp','-o',D/'gas.o'],'compile-gas')
 headers={p:sha(p) for p in (D/'motion.d').read_text().replace('\\\n',' ').split(':',1)[1].split()};receipt['compiler_headers']=headers
 run(['/usr/bin/mpic++',*san,'-Wl,--gc-sections',D/'motion.o',D/'gas.o',N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*native.LIBS,'-o',D/'test_motion'],'link')
 run([D/'test_motion'],'controls')
 receipt['results']=json.loads((D/'controls.log').read_text());assert receipt['results']['passed']
 assert before==pins() and deps==native.sealed_inputs() and all(sha(p)==v for p,v in headers.items())
 receipt['binary_sha256']=sha(D/'test_motion');receipt['passed']=True
 receipt['sanitizer_scope']='adapter and gas owner rebuilt; pinned native archive and external libraries not rebuilt'
finally:(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
