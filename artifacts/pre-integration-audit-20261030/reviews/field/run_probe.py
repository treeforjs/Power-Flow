from pathlib import Path
import fcntl,hashlib,json,os,subprocess,time
here=Path(__file__).resolve().parent
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'])
assert os.getenv('CUDA_VISIBLE_DEVICES')==''
lease='/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock'
with open(lease,'a') as f:
 try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('Acquire CPU6 lease first')
source=here/'inputs/hall-source'
module=source/'development/hall-canonical-driver-20261029'
native=module/'bindings/native-v6.2'
dep=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
build=here/'probe-build';build.mkdir(exist_ok=False)
libs=[dep/'install-mpi/lib/libmfem.a',dep/'tpl/hypre/lib/libHYPRE.a',dep/'tpl/metis/lib/libmetis.a',dep/'tpl/GKlib/lib/libGKlib.a']
flags=['-std=c++20','-O2','-DNDEBUG','-ffp-contract=off','-fno-fast-math','-ffunction-sections','-fdata-sections',
 '-I'+str(source/'development/hall-finite-foil-em-20261029'),'-I'+str(native/'source/include'),
 '-isystem',str(native/'inputs'),'-I'+str(native/'vendor')]
for d in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:flags+=['-isystem',str(dep/d)]
commands=[('native',['/usr/bin/c++',*flags,'-c',str(native/'source/src/mfem_gap_pic.cpp'),'-o',str(build/'native.o')]),
 ('probe',['/usr/bin/c++',*flags,'-c',str(here/'probe_operators.cpp'),'-o',str(build/'probe.o')]),
 ('link',['/usr/bin/c++','-Wl,--gc-sections',str(build/'probe.o'),str(build/'native.o'),*map(str,libs),'/usr/lib/libmpi.so','-lcrypto','-o',str(build/'operators')]),
 ('observe',[str(build/'operators')])]
receipts=[]
try:
 for label,cmd in commands:
  actual=['/home/treefor/projects/runtimes/arch-base-20261028/run-frozen','prlimit','--as=2147483648','--cpu=180','--',*cmd]
  start=time.monotonic()
  with (here/('operator-observations.json' if label=='observe' else label+'.log')).open('w') as log:
   result=subprocess.run(actual,cwd=here,stdout=log,stderr=subprocess.STDOUT,timeout=240)
  receipts.append({'label':label,'command':actual,'returncode':result.returncode,'elapsed_s':time.monotonic()-start,
                   'CPU':6,'nice':15,'threads':1,'lease':lease})
  (here/'probe-commands.json').write_text(json.dumps(receipts,indent=2)+'\n')
  print(label,result.returncode,receipts[-1]['elapsed_s'],flush=True)
  if result.returncode:raise RuntimeError('Probe command failed: '+label)
 manifest=json.loads((here/'input-manifest.json').read_text())
 for entry in manifest['files']:
  assert hashlib.sha256((here/entry['path']).read_bytes()).hexdigest()==entry['sha256']
 (here/'probe-status.json').write_text(json.dumps({'status':'observations_recorded','full_coupled_test':False,
  'sealed_inputs_unchanged':True,'binary_sha256':hashlib.sha256((build/'operators').read_bytes()).hexdigest()},indent=2)+'\n')
except Exception as e:
 (here/'probe-status.json').write_text(json.dumps({'status':'failed','error':str(e)},indent=2)+'\n');raise
