from pathlib import Path
import os,sys,fcntl,hashlib,json,subprocess,time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PREFIX=ROOT.parent
OUT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/hall-coupled-transport-20261029')
C=json.loads((HERE/'controls.json').read_text())
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(v)=='1' for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')==''
with open(C['lease'],'a') as probe:
    try:fcntl.flock(probe,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:pass
    else:raise RuntimeError('Acquire shared testing-core.lock before invoking runner')
run=OUT/(sys.argv[1] if len(sys.argv)>1 else 'run-001')
run.mkdir(exist_ok=False)
build=PREFIX/'build';build.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
dep=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
files=[ROOT/'src/mfem_gap_pic.cpp',ROOT/'src/mfem_pic_distributed.hpp',ROOT/'include/raikou/mfem_gap_pic.h',ROOT/'include/raikou/q1_scalar_gather.h',HERE/'benchmark.cpp',HERE/'coupled_pic.hpp',HERE/'controls.json',HERE/'run.py',HERE/'validate.py',HERE/'plot.py',ROOT/'tests/test_mfem_gap_pic.cpp',dep/'install-mpi/lib/libmfem.a',dep/'tpl/hypre/lib/libHYPRE.a',dep/'tpl/metis/lib/libmetis.a',dep/'tpl/GKlib/lib/libGKlib.a']
manifest={str(p):sha(p) for p in files};(run/'input-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
flags=['-std=c++20','-O2','-DNDEBUG','-I'+str(ROOT/'include'),'-I'+str(ROOT/'src'),'-I'+str(HERE),'-isystem',str(dep/'install-mpi/include'),'-isystem',str(dep/'tpl/hypre/include'),'-isystem',str(dep/'tpl/metis/include')]
libs=[str(dep/'install-mpi/lib/libmfem.a'),str(dep/'tpl/hypre/lib/libHYPRE.a'),str(dep/'tpl/metis/lib/libmetis.a'),str(dep/'tpl/GKlib/lib/libGKlib.a'),'/usr/lib/libmpi.so']
commands=[('build',['/usr/bin/c++',*flags,str(ROOT/'src/mfem_gap_pic.cpp'),str(HERE/'benchmark.cpp'),'-o',str(build/'hall_benchmark'),*libs]),('benchmark',[str(build/'hall_benchmark'),str(run)]),('foundation_build',['/usr/bin/c++',*flags,str(ROOT/'src/mfem_gap_pic.cpp'),str(ROOT/'tests/test_mfem_gap_pic.cpp'),'-o',str(build/'pic_foundation_tests'),*libs]),('foundation',[str(build/'pic_foundation_tests')]),('validation',['/usr/bin/python',str(HERE/'validate.py'),str(run)]),('plots',['/home/treefor/projects/Power-Flow/.venv/bin/python',str(HERE/'plot.py'),str(run)])]
receipts=[]
for label,cmd in commands:
    command=[C['frozen_launcher'],'prlimit','--as=2147483648','--cpu=180','--',*cmd]
    start=time.monotonic()
    with (run/(label+'.log')).open('w') as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
    receipts.append({'label':label,'command':command,'elapsed_s':time.monotonic()-start,'returncode':result.returncode,'cpu_affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':1,'shared_lease_exclusion_verified':True})
    (run/'commands.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print(label,result.returncode,receipts[-1]['elapsed_s'],flush=True)
    if result.returncode:raise RuntimeError((run/(label+'.log')).read_text()[-6000:])
assert all(sha(Path(p))==h for p,h in manifest.items())
(run/'binary-sha256.txt').write_text(sha(build/'hall_benchmark')+'\n')
print(run,flush=True)
