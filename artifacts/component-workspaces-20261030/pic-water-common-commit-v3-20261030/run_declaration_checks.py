#!/usr/bin/env python3
import fcntl, hashlib, json, os, resource, subprocess, time
from pathlib import Path

P=Path(__file__).resolve().parent
R=Path('/home/treefor/projects/runtimes')
O=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-common-commit-v3-20261030')
N=R/'pic-reaction-transactions-20261029/native-full-time-candidate003'
H=R/'hall-coupled-transport-20261029/source'
I=R/'pic-water-full-time-subflow-20261029/inputs'
M=R/'mfem-4.10-pic-20261029'
F=R/'arch-base-20261028/run-frozen'
if not __debug__ or os.getenv('PYTHONOPTIMIZE'):
    raise SystemExit('Optimized guard execution refused')
if os.sched_getaffinity(0)!={7} or os.getpriority(os.PRIO_PROCESS,0)!=15:
    raise SystemExit('CPU7 nice15 required')
lease=open(O.parent/'pic-profile/testing-core7.lock','a')
try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit(75)
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
env={'PATH':'/usr/bin:/bin','LC_ALL':'C','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
     'MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1',
     'OMP_THREAD_LIMIT':'1','CUDA_VISIBLE_DEVICES':'','PYTHONNOUSERSITE':'1'}
roots=[P/'source',P/'tests',N/'source',N/'inputs',H,I,M/'install-mpi/include',M/'tpl/hypre/include',M/'tpl/metis/include']
def pins():
    files={p for d in roots for p in d.rglob('*') if p.is_file()}
    files.update([Path(__file__).resolve(),F,N/'SOURCE-SHA256.json',M/'source-lock.json',M/'installed-files.json',
                  R/'arch-base-20261028/root/usr/bin/g++',R/'arch-base-20261028/root/usr/bin/mpic++'])
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
before=pins()
source_manifest=json.loads((N/'SOURCE-SHA256.json').read_text())
for name,digest in source_manifest['files'].items():
    if before[str(N/name)]!=digest:raise RuntimeError('Native source differs from declared identity: '+name)
for name in ('motion_receipt.hpp','native_arena_motion.cpp','kinetic_boundary.hpp','kinetic_boundary.cpp'):
    if before[str(N/'inputs/hall-motion-v2'/name)]!=before[str(H/'development/hall-canonical-driver-20261029'/name)]:
        raise RuntimeError('Different Hall motion provider: '+name)
O.mkdir(exist_ok=True);(P/'build').mkdir(exist_ok=True)
A=O/f'declaration-{len(list(O.glob("declaration-*")))+1:03d}';A.mkdir()
(A/'inputs.json').write_text(json.dumps(before,indent=2)+'\n')
flags=['-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math',
       '-ffunction-sections','-fdata-sections','-I'+str(P/'source')]
for d in (N/'source/include',I,I/'include',N/'inputs',H/'development/hall-field-correction-20261030',
          M/'install-mpi/include',M/'tpl/hypre/include',M/'tpl/metis/include'):
    flags+=['-isystem',str(d)]
prefix=[str(F),'/usr/bin/mpic++',*flags]
exe=P/'build/test_closed_configuration'
commands=[prefix+['-fsyntax-only',str(P/'tests/declaration_opaque_joint.cpp')],
          prefix+[str(P/'tests/test_closed_configuration.cpp'),'-Wl,--gc-sections','-lcrypto','-o',str(exe)],
          [str(F),str(exe)]]
(A/'commands.json').write_text(json.dumps({'commands':commands,'environment':env},indent=2)+'\n')
start=time.monotonic();status='failed';result=None;error=None
try:
    for index,cmd in enumerate(commands):
        with (A/f'command-{index}.log').open('w') as log:
            subprocess.run(cmd,cwd=P,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=45 if index<2 else 10,check=True)
    result=json.loads((A/'command-2.log').read_text())
    if result.get('status')!='passed' or result.get('full_joint_test') is not False or result.get('accepted_owners_constructed')!=0:
        raise RuntimeError('Incorrect component check scope')
    if pins()!=before:raise RuntimeError('Dependency content or membership changed during checks')
    status='passed'
except Exception as exc:error=repr(exc)
receipt={'status':status,'result':result,'error':error,'elapsed_s':time.monotonic()-start,
    'native_source_identity':source_manifest['source_identity'],'cpu':7,'nice':15,'threads':1,'core7_lease':True,
    'GPU_used':False,'full_joint_execution':False,'physical_time_advanced_s':0,
    'scope':'root header declaration and pure configuration helper; no accepted owner constructed',
    'dependency_scope':'listed application trees and numerical headers pinned before/after; not full toolchain closure qualification',
    'inputs_file_sha256':hashlib.sha256((A/'inputs.json').read_bytes()).hexdigest()}
(A/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'attempt':str(A),**receipt}))
if status!='passed':raise SystemExit(1)
