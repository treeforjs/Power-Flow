#!/usr/bin/env python3
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time

P = Path(__file__).resolve().parent
O = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-water-common-commit-v2-20261030')
R = Path('/home/treefor/projects/runtimes')
N = R/'pic-water-full-time-subflow-20261029/native'
M = R/'mfem-4.10-pic-20261029'
F = R/'arch-base-20261028/run-frozen'
if os.sched_getaffinity(0)!={7} or os.getpriority(os.PRIO_PROCESS,0)!=15:
    raise SystemExit('CPU7 / nice15 admission required')
if os.environ.get('PYTHONOPTIMIZE') or __debug__ is False:
    raise SystemExit('Optimized guard execution refused')
lease = open(O.parent/'pic-profile/testing-core7.lock','a')
try:
    fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:
    raise SystemExit(75)
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
env={'PATH':'/usr/bin:/bin','LC_ALL':'C','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
     'MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1',
     'VECLIB_MAXIMUM_THREADS':'1','OMP_THREAD_LIMIT':'1','CUDA_VISIBLE_DEVICES':'',
     'PYTHONNOUSERSITE':'1'}
inputs=[P/'source/common_commit_contract.hpp',P/'tests/test_common_commit_contract.cpp',
        N/'source/include/raikou/mfem_pic_full_time.h',
        N/'source/include/raikou/mfem_pic_arena.h',F,
        M/'source-lock.json',M/'installed-files.json',
        R/'arch-base-20261028/root/usr/bin/g++',
        R/'arch-base-20261028/root/usr/bin/mpic++']
def pin(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p):pin(p) for p in inputs}
exe=P/'build/test_common_commit_contract'
command=[str(F),'/usr/bin/mpic++','-std=c++20','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror',
         '-ffp-contract=off','-fno-fast-math','-ffunction-sections','-fdata-sections',
         '-I'+str(P/'source'),'-isystem',str(N/'source/include'),'-isystem',str(N/'inputs'),
         '-isystem',str(N/'vendor'),'-isystem',str(M/'install-mpi/include'),
         '-isystem',str(M/'tpl/hypre/include'),'-isystem',str(M/'tpl/metis/include'),
         str(P/'tests/test_common_commit_contract.cpp'),'-Wl,--gc-sections','-o',str(exe)]
attempts=[p for p in O.glob('attempt-*') if p.is_dir()]
A=O/f'attempt-{len(attempts)+1:03d}';A.mkdir()
(A/'commands.json').write_text(json.dumps({'compile':command,'run':[str(F),str(exe)],'environment':env},indent=2)+'\n')
start=time.monotonic()
with (A/'build.log').open('w') as log:
    subprocess.run(command,env=env,cwd=P,stdout=log,stderr=subprocess.STDOUT,timeout=45,check=True)
result=subprocess.run([str(F),str(exe)],env=env,cwd=P,capture_output=True,text=True,timeout=10)
(A/'run.log').write_text(result.stdout+result.stderr)
result.check_returncode()
record=json.loads(result.stdout)
if record.get('status')!='passed' or record.get('full_joint_test') is not False:
    raise RuntimeError('Wrong predicate-test result scope')
after={str(p):pin(p) for p in inputs}
if before!=after:raise RuntimeError('Test input changed')
receipt={'status':'passed','result':record,'input_sha256':before,'executable_sha256':pin(exe),
         'elapsed_s':time.monotonic()-start,'cpu':7,'nice':15,'threads':1,
         'exclusive_core7_lease':True,'frozen_launcher':str(F),'GPU_used':False,
         'native_field_heat_source_operators_run':False,'common_driver_linked':False,
         'scope':'cross-owner value predicate; actual corrected owner integration pending'}
(A/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
(O/'predicate-check-handoff.json').write_text(json.dumps({'attempt':str(A),'receipt':str(A/'receipt.json'),**receipt},indent=2)+'\n')
print(json.dumps({'status':'passed','checks':record['checks'],'elapsed_s':receipt['elapsed_s'],'attempt':str(A)}))
