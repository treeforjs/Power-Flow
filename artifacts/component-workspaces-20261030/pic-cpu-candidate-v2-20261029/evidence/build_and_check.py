from pathlib import Path
import datetime
import fcntl
import hashlib
import json
import os
import shlex
import subprocess
import time

ROOT=Path('/home/treefor/projects/runtimes/pic-cpu-candidate-v2-20261029')
EVIDENCE=ROOT/'evidence'
LAUNCHER='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCK='/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock'
V1=Path('/home/treefor/projects/runtimes/pic-cpu-candidate-20261029/source')
ACTIVE=Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source')
MFEM=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')

assert os.sched_getaffinity(0)=={6}
assert os.getpriority(os.PRIO_PROCESS,0)==15
with open(LOCK,'a') as lease:
    try:
        fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        pass
    else:
        raise RuntimeError('Run this entry point beneath the testing-core flock lease')
commands=[]
def run(label,args):
    command=[LAUNCHER,*map(str,args)]
    start=time.monotonic()
    result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
    record={'label':label,'command':command,'shell_display':shlex.join(command),
            'elapsed_seconds':time.monotonic()-start,'returncode':result.returncode,
            'cpu_affinity':sorted(os.sched_getaffinity(0)),
            'nice':os.getpriority(os.PRIO_PROCESS,0),'lease_exclusion_confirmed':True}
    commands.append(record)
    (EVIDENCE/(label+'.log')).write_text(result.stdout+result.stderr)
    (EVIDENCE/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    if result.returncode:
        raise RuntimeError(label+' failed: '+result.stdout+result.stderr)
    print(label+' passed',flush=True)
    return result.stdout

run('configure',['cmake','-S','source','-B','build','-DCMAKE_BUILD_TYPE=Release',
                 '-DMFEM_DIR='+str(MFEM/'install-mpi/lib/cmake/mfem')])
run('build',['cmake','--build','build','--parallel','1'])
run('foundation',[ROOT/'build/pic_foundation_tests'])
run('q1-invalidation',[ROOT/'build/q1_gather_tests'])
trace=run('serial-trace',[ROOT/'build/pic_cpu_trace'])
baseline=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-audit/baseline-trace.txt').read_text()
assert trace==baseline,'V2 serial trace differs from the sealed-source baseline'

includes=['-isystem',str(MFEM/'install-mpi/include'),'-isystem',str(MFEM/'tpl/hypre/include'),
          '-isystem',str(MFEM/'tpl/metis/include')]
link_command=shlex.split((ROOT/'build/CMakeFiles/compare_checkpoints.dir/link.txt').read_text())
mfem_libraries=link_command[link_command.index('-o')+2:]
abi={}
for label,source in [('active',ACTIVE),('v1',V1),('v2',ROOT/'source')]:
    executable=EVIDENCE/('abi-'+label)
    run('abi-build-'+label,['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG',
                            '-I'+str(source/'include'),*includes,EVIDENCE/'abi_probe.cpp','-o',executable,*mfem_libraries])
    abi[label]=run('abi-'+label,[executable])
assert abi['active']==abi['v1']==abi['v2'],'Checkpoint structure ABI changed'

run('gpu-target-dynamic-links',['readelf','-d',ROOT/'build/heating_history_pic_gpu'])
links=(EVIDENCE/'gpu-target-dynamic-links.log').read_text()
assert 'libcudss.so.0' in links and 'libcudart.so.12' in links
before=json.loads((EVIDENCE/'v1-source-before.json').read_text())
after={str(p.relative_to(V1)):hashlib.sha256(p.read_bytes()).hexdigest() for p in V1.rglob('*') if p.is_file()}
assert before==after,'V1 source changed during V2 preparation'
summary={'status':'passed','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'cpu_affinity':[6],'nice':15,'compiler_jobs':1,'thread_environment':
         {k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CUDA_VISIBLE_DEVICES']},
         'lease_exclusion_confirmed':True,'gpu_target_built':True,'gpu_target_run':False,
         'foundation_passed':True,'q1_invalidation_passed':True,'serial_trace_byte_identical':True,
         'serial_trace_bytes':len(trace.encode()),'serial_trace_sha256':hashlib.sha256(trace.encode()).hexdigest(),
         'abi_identical_active_v1_v2':True,'abi_text_sha256':hashlib.sha256(abi['v2'].encode()).hexdigest(),
         'v1_source_unchanged':True}
(EVIDENCE/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print('V2 build, tests, serial trace and checkpoint ABI checks passed',flush=True)
