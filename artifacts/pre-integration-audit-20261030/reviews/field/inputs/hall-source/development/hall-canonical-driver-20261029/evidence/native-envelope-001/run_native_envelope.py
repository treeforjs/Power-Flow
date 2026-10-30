from pathlib import Path
import fcntl,hashlib,json,os,shutil,subprocess,sys,time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
OUT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/hall-coupled-transport-20261029/canonical-driver')
NATIVE=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time')
PARENT=Path('/home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029')
DEP=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LEASE='/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']) and os.getenv('CUDA_VISIBLE_DEVICES')==''
with open(LEASE,'a') as f:
    try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:pass
    else:raise RuntimeError('Acquire the exclusive CPU6 lease first')
run=OUT/sys.argv[1];run.mkdir(parents=True,exist_ok=False)
for name in ['native_joint_binding.hpp','root_joint_prepare_example.hpp','test_native_motion_envelope.cpp','run_native_envelope.py']:
    shutil.copyfile(HERE/name,run/name)
shutil.copyfile(NATIVE/'MOTION-RECEIPT-SCHEMA.txt',run/'native-MOTION-RECEIPT-SCHEMA.txt')
shutil.copyfile(PARENT/'source/closed_parent_joint_step.hpp',run/'root-closed_parent_joint_step.hpp')
probe=run/'probe.cpp'
probe.write_text('#include "native_parent_subflow.hpp"\n#include "root_joint_prepare_example.hpp"\n'
    'void probe(raikou::mfem_backend::PicArena&a,raikou::native_water::ParentSubflow&p,const axial::Mesh&g,'
    'const axial::Fields&f,axial::State&s,hall_driver::Thermal&t,const hall_driver::CurrentTrace&c,'
    'const std::string&h){hall_driver::root_two_exposure_example(a,p,g,f,s,t,c,1e-14,h,h,h,0,100000,4096);}\n')
files=[p for directory in [HERE,NATIVE/'source/include',NATIVE/'inputs',NATIVE/'vendor',PARENT/'inputs'] for p in directory.rglob('*') if p.is_file()]
files += [PARENT/'source/native_parent_subflow.hpp',NATIVE/'MOTION-RECEIPT-SCHEMA.txt',PARENT/'source/closed_parent_joint_step.hpp']
files += [p for p in (HERE.parent/'hall-finite-foil-em-20261029').iterdir() if p.is_file()]
files += [HERE.parent/'hall-coupled-transport-20261029/coupled_pic.hpp']
libs=[DEP/x for x in ['install-mpi/lib/libmfem.a','tpl/hypre/lib/libHYPRE.a','tpl/metis/lib/libmetis.a','tpl/GKlib/lib/libGKlib.a']]
files += libs
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest={str(p):sha(p) for p in files};(run/'input-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
includes=['-I'+str(HERE),'-I'+str(PARENT/'source'),'-I'+str(NATIVE/'source/include'),
    '-isystem',str(PARENT/'inputs'),'-isystem',str(PARENT/'inputs/include'),
    '-isystem',str(NATIVE/'inputs'),'-isystem',str(NATIVE/'vendor')]
for d in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:includes+=['-isystem',str(DEP/d)]
flags=['-std=c++20','-O2','-ffp-contract=off','-fno-fast-math','-ffunction-sections','-fdata-sections',*includes]
binary=ROOT.parent/'build-canonical-driver/native_envelope_test'
commands=[('build-envelope',['/usr/bin/c++',*flags,
    str(HERE/'bindings/native-v6.2/source/src/mfem_gap_pic.cpp'),str(HERE/'kinetic_boundary.cpp'),
    str(HERE/'native_arena_motion.cpp'),str(HERE/'test_native_motion_envelope.cpp'),
    '-Wl,--gc-sections',*map(str,libs),'/usr/lib/libmpi.so','-lcrypto','-o',str(binary)]),
    ('envelope',[str(binary),str(run)]),('root-example-syntax',['/usr/bin/c++',*flags,'-fsyntax-only',str(probe)])]
receipts=[]
try:
    for label,cmd in commands:
        actual=[FROZEN,'prlimit','--as=2147483648','--cpu=180','--',*cmd];start=time.monotonic()
        (run/'run-status.json').write_text(json.dumps({'status':'running','command':label})+'\n')
        with (run/(label+'.log')).open('w') as log:result=subprocess.run(actual,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
        receipts.append({'label':label,'command':actual,'returncode':result.returncode,'elapsed_s':time.monotonic()-start,'CPU':[6],'nice':15,'threads':1,'lease':LEASE})
        (run/'commands.json').write_text(json.dumps(receipts,indent=2)+'\n');print(label,result.returncode,flush=True)
        if result.returncode:raise RuntimeError((run/(label+'.log')).read_text()[-5000:])
    assert all(sha(Path(p))==s for p,s in manifest.items()),'Input changed during focused qualification'
    (run/'binary-sha256.txt').write_text(sha(binary)+'\n')
    validation=json.loads((run/'validation.json').read_text())
    (run/'run-status.json').write_text(json.dumps({'status':'passed','envelope_checks':validation['count'],
        'root_example_syntax':'passed','component_matrices_repeated':False,'full_native_link':'pending'},indent=2)+'\n')
except Exception as e:
    (run/'run-status.json').write_text(json.dumps({'status':'failed','error':str(e)},indent=2)+'\n');raise
print(run,flush=True)
