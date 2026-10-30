from pathlib import Path
import fcntl, hashlib, json, os, shutil, subprocess, sys, time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/hall-coupled-transport-20261029/canonical-driver')
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LEASE='/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'])
assert os.getenv('CUDA_VISIBLE_DEVICES')==''
with open(LEASE,'a') as f:
    try: fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: pass
    else: raise RuntimeError('Acquire the exclusive CPU6 lease first')
run=OUT/sys.argv[1]
run.mkdir(parents=True,exist_ok=False)
shutil.copytree(HERE,run/'source-snapshot')
build=ROOT.parent/'build-canonical-driver'
build.mkdir(exist_ok=True)
native=HERE/'bindings/native-v6.2'
dep=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
libs=[dep/'install-mpi/lib/libmfem.a',dep/'tpl/hypre/lib/libHYPRE.a',dep/'tpl/metis/lib/libmetis.a',dep/'tpl/GKlib/lib/libGKlib.a']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[p for p in HERE.rglob('*') if p.is_file()]
files += [p for p in (HERE.parent/'hall-finite-foil-em-20261029').iterdir() if p.is_file()]
files += list((HERE.parent/'hall-finite-foil-em-20261029/data').glob('*'))
files += [HERE.parent/'hall-coupled-transport-20261029/coupled_pic.hpp',*libs]
manifest={str(p):sha(p) for p in files}
(run/'input-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
includes=['-I'+str(native/'source/include'),'-isystem',str(native/'inputs'),'-I'+str(native/'vendor')]
for d in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:
    includes += ['-isystem',str(dep/d)]
flags=['-std=c++20','-O2','-DNDEBUG','-ffp-contract=off','-fno-fast-math',*includes]
binary=build/'thermal_v2'
commands=[('build',['/usr/bin/c++',*flags,str(native/'source/src/mfem_gap_pic.cpp'),str(HERE/'test_thermal_v2.cpp'),'-o',str(binary),*[str(p) for p in libs],'/usr/lib/libmpi.so'])]
commands += [(mode,[str(binary),str(run),mode]) for mode in ['benchmark','full','split','resume']]
receipts=[]
try:
    for label,command in commands:
        (run/'run-status.json').write_text(json.dumps({'status':'running','command':label})+'\n')
        actual=[FROZEN,'prlimit','--as=2147483648','--cpu=180','--',*command]
        started=time.monotonic()
        with (run/(label+'.log')).open('w') as log:
            result=subprocess.run(actual,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
        receipts.append({'label':label,'command':actual,'returncode':result.returncode,
            'elapsed_s':time.monotonic()-started,'CPU':[6],'nice':15,'threads':1,'lease':LEASE})
        (run/'commands.json').write_text(json.dumps(receipts,indent=2)+'\n')
        print(label,result.returncode,receipts[-1]['elapsed_s'],flush=True)
        if result.returncode:raise RuntimeError((run/(label+'.log')).read_text()[-6000:])
    assert all(sha(Path(p))==h for p,h in manifest.items())
    checkpoint_files=sorted(p.relative_to(run/'full') for p in (run/'full').rglob('*') if p.is_file())
    assert all((run/'full'/p).read_bytes()==(run/'resume'/p).read_bytes() for p in checkpoint_files)
    validation=json.loads((run/'validation.json').read_text())
    validation['checks'].append({'name':'fresh-process v2 heat and Maxwell fixture restart exact',
        'passed':True,'detail':{'files':len(checkpoint_files),'canonical_arena_checkpoint':False}})
    validation['count']=len(validation['checks'])
    (run/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    (run/'binary-sha256.txt').write_text(sha(binary)+'\n')
    (run/'run-status.json').write_text(json.dumps({'status':'passed','checks':validation['count'],'canonical_arena_integrated':False},indent=2)+'\n')
except Exception as e:
    (run/'run-status.json').write_text(json.dumps({'status':'failed','error':str(e)},indent=2)+'\n')
    raise
print(run,flush=True)
