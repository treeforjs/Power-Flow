from pathlib import Path
import fcntl,hashlib,json,os,shutil,subprocess,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
OUT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/hall-coupled-transport-20261029/canonical-driver')
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LEASE='/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']) and os.getenv('CUDA_VISIBLE_DEVICES')==''
with open(LEASE,'a') as f:
    try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:pass
    else:raise RuntimeError('Acquire the exclusive CPU6 lease first')
run=OUT/sys.argv[1];run.mkdir(parents=True,exist_ok=False);shutil.copytree(HERE,run/'source-snapshot')
build=ROOT.parent/'build-canonical-driver';build.mkdir(exist_ok=True)
native=HERE/'bindings/native-v6.2';dep=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
libs=[dep/'install-mpi/lib/libmfem.a',dep/'tpl/hypre/lib/libHYPRE.a',dep/'tpl/metis/lib/libmetis.a',dep/'tpl/GKlib/lib/libGKlib.a']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[p for p in HERE.rglob('*') if p.is_file()]+[p for p in (HERE.parent/'hall-finite-foil-em-20261029').iterdir() if p.is_file()]
files += [HERE.parent/'hall-coupled-transport-20261029/coupled_pic.hpp',*libs]
manifest={str(p):sha(p) for p in files};(run/'input-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
identity=json.loads((HERE/'bindings.json').read_text())['native_identity']
flags=['-std=c++20','-O2','-DNDEBUG','-ffp-contract=off','-fno-fast-math',
    '-ffunction-sections','-fdata-sections','-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"',
    '-I'+str(native/'source/include'),'-isystem',str(native/'inputs'),'-I'+str(native/'vendor')]
for d in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:flags += ['-isystem',str(dep/d)]
commands=[];objects=[]
for source,label in [(native/'source/src/mfem_gap_pic.cpp','native'),(native/'source/src/mfem_pic_arena.cpp','arena'),(HERE/'kinetic_boundary.cpp','kinetic-boundary'),
        (HERE/'native_arena_motion.cpp','motion'),(HERE/'test_joint_v2.cpp','test-joint-v2')]:
    obj=build/(label+'.o');objects.append(str(obj));commands.append(('build-'+label,['/usr/bin/c++',*flags,'-c',str(source),'-o',str(obj)]))
binary=build/'joint_v2_test'
commands.append(('link',['/usr/bin/c++','-Wl,--gc-sections',
    '-Wl,--wrap=_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv',
    '-Wl,--wrap=_ZN6raikou12mfem_backend6GapPic7advanceEd',*objects,
    *map(str,libs),'/usr/lib/libmpi.so','-lcrypto','-o',str(binary)]))
commands.append(('joint',[str(binary),str(run)]));receipts=[]
try:
    for label,cmd in commands:
        (run/'run-status.json').write_text(json.dumps({'status':'running','command':label})+'\n')
        actual=[FROZEN,'prlimit','--as=2147483648','--cpu=180','--',*cmd];start=time.monotonic()
        with (run/(label+'.log')).open('w') as log:result=subprocess.run(actual,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
        receipts.append({'label':label,'command':actual,'elapsed_s':time.monotonic()-start,'returncode':result.returncode,'CPU':[6],'nice':15,'threads':1,'lease':LEASE})
        (run/'commands.json').write_text(json.dumps(receipts,indent=2)+'\n');print(label,result.returncode,receipts[-1]['elapsed_s'],flush=True)
        if result.returncode:raise RuntimeError((run/(label+'.log')).read_text()[-6000:])
    assert all(sha(Path(p))==h for p,h in manifest.items())
    validation=json.loads((run/'joint-validation.json').read_text())
    (run/'binary-sha256.txt').write_text(sha(binary)+'\n')
    (run/'run-status.json').write_text(json.dumps({'status':'passed','checks':validation['count'],
        'accepted_native_arena_commit':False,'actual_SSA':False},indent=2)+'\n')
except Exception as e:
    (run/'run-status.json').write_text(json.dumps({'status':'failed','error':str(e)},indent=2)+'\n');raise
print(run,flush=True)
