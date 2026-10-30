from pathlib import Path
import os,sys,fcntl,hashlib,json,subprocess,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PREFIX=ROOT.parent
OUT=Path('/home/treefor/sim-workspace/powerflow/outputs/hall-coupled-transport-20261029/finite-foil-em')
launcher='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
lease='/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(v)=='1' for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')) and os.getenv('CUDA_VISIBLE_DEVICES')==''
with open(lease,'a') as f:
 try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('Acquire test lease first')
run=OUT/sys.argv[1];run.mkdir(parents=True,exist_ok=False)
build=PREFIX/'build-finite-foil-em';build.mkdir(exist_ok=True);binary=build/'finite_foil_em'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
dep=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
libs=[dep/'install-mpi/lib/libmfem.a',dep/'tpl/hypre/lib/libHYPRE.a',dep/'tpl/metis/lib/libmetis.a',dep/'tpl/GKlib/lib/libGKlib.a']
files=[ROOT/'src/mfem_gap_pic.cpp',ROOT/'src/mfem_pic_distributed.hpp',ROOT/'include/raikou/mfem_gap_pic.h',ROOT/'include/raikou/q1_scalar_gather.h',HERE/'foil_em.hpp',HERE/'foil_particles.hpp',HERE/'surface_path_port.hpp',HERE/'native_total_field_transport.hpp',HERE/'surface-source-binding.json',HERE/'benchmark.cpp',HERE/'controls.json',HERE/'run.py',HERE/'native-start-binding.json',HERE.parent/'hall-coupled-transport-20261029/coupled_pic.hpp',*libs]
manifest={str(p):sha(p) for p in files};(run/'input-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
flags=['-std=c++20','-O2','-DNDEBUG','-I'+str(ROOT/'include'),'-I'+str(ROOT/'src'),'-isystem',str(dep/'install-mpi/include'),'-isystem',str(dep/'tpl/hypre/include'),'-isystem',str(dep/'tpl/metis/include')]
commands=[('build',['/usr/bin/c++',*flags,str(ROOT/'src/mfem_gap_pic.cpp'),str(HERE/'benchmark.cpp'),'-o',str(binary),*[str(p) for p in libs],'/usr/lib/libmpi.so']),('benchmark',[str(binary),str(run)])]
receipts=[]
for label,cmd in commands:
 actual=[launcher,'prlimit','--as=2147483648','--cpu=180','--',*cmd];start=time.monotonic()
 with (run/(label+'.log')).open('w') as log:result=subprocess.run(actual,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
 receipts.append({'label':label,'command':actual,'elapsed_s':time.monotonic()-start,'returncode':result.returncode,'CPU':[6],'nice':15,'threads':1,'shared_lease':True})
 (run/'commands.json').write_text(json.dumps(receipts,indent=2)+'\n');print(label,result.returncode,receipts[-1]['elapsed_s'],flush=True)
 if result.returncode:raise RuntimeError((run/(label+'.log')).read_text()[-6000:])
assert all(sha(Path(p))==h for p,h in manifest.items())
(run/'binary-sha256.txt').write_text(sha(binary)+'\n');print(run,flush=True)
