from pathlib import Path
import datetime
import fcntl
import hashlib
import json
import os
import shlex
import subprocess
import time

ROOT=Path('/home/treefor/projects/runtimes/pic-repo-integration-20261029')
REPO=Path('/home/treefor/projects/raikou-mfem')
EVIDENCE=ROOT/'evidence'
BASELINE=ROOT/'baseline-source'
LAUNCHER='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCK='<workspace>/outputs/pic-profile/testing-core.lock'
MFEM=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
with open(LOCK,'a') as lease:
    try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:pass
    else:raise RuntimeError('Testing-core flock lease is required')
commands=[]
def run(label,args,cwd=ROOT,allowed=(0,)):
    command=[LAUNCHER,*map(str,args)];start=time.monotonic()
    result=subprocess.run(command,cwd=cwd,capture_output=True,text=True)
    commands.append({'label':label,'command':command,'shell_display':shlex.join(command),
                     'elapsed_seconds':time.monotonic()-start,'returncode':result.returncode,
                     'logical_cpus':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0),
                     'lease_exclusion_confirmed':True})
    (EVIDENCE/(label+'.log')).write_text(result.stdout+result.stderr)
    (EVIDENCE/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    if result.returncode not in allowed:raise RuntimeError(label+' failed:\n'+result.stdout+result.stderr)
    print(label+' passed',flush=True);return result.stdout

targets=['raikou_test_mfem_electrostatics','raikou_test_mfem_gap_pic','raikou_test_mfem_q1_scalar_gather',
         'raikou_mfem_gap_pic','raikou_mfem_finite_foil_pic']
for mode in ['serial','mpi']:
    build=ROOT/(mode+'-build')
    run(mode+'-configure',['cmake','-S',REPO,'-B',build,'-DCMAKE_BUILD_TYPE=Release',
                          '-DRAIKOU_ENABLE_MFEM=ON','-DRAIKOU_ENABLE_MPI='+('ON' if mode=='mpi' else 'OFF'),
                          '-DRAIKOU_ENABLE_CUDA=OFF','-DRAIKOU_BUILD_RUNTIME=OFF','-DRAIKOU_BUILD_TESTS=ON',
                          '-DRAIKOU_BUILD_PYTHON=OFF','-DMFEM_DIR='+str(MFEM/'install-mpi/lib/cmake/mfem')])
    selected=targets+(['raikou_mfem_heating_history_pic','raikou_mfem_heating_pic'] if mode=='mpi' else [])
    run(mode+'-build',['cmake','--build',build,'--parallel','1','--target',*selected])
    run(mode+'-ctest',['ctest','--test-dir',build,'--output-on-failure','--parallel','1','-R',r'^raikou\.mfem_'])

def flags(build):
    text=(build/'CMakeFiles/raikou_mfem.dir/flags.make').read_text()
    includes=shlex.split(next(line.split('=',1)[1] for line in text.splitlines() if line.startswith('CXX_INCLUDES =')))
    raw=shlex.split((build/'CMakeFiles/raikou_test_mfem_gap_pic.dir/link.txt').read_text())
    libraries=[str(build/value) if value.endswith('.a') and not value.startswith('/') else value for value in raw[raw.index('-o')+2:]]
    return includes,libraries
includes,libraries=flags(ROOT/'serial-build')
baseline_includes=[arg.replace(str(REPO/'include'),str(BASELINE/'include')) for arg in includes]
baseline_object=EVIDENCE/'baseline-mfem-gap-pic.o'
run('baseline-core-build',['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG','-fPIC','-Wall','-Wextra','-Wpedantic','-Werror',
                           *baseline_includes,'-c',BASELINE/'src/mfem_gap_pic.cpp','-o',baseline_object])
baseline_libraries=[baseline_object,*[arg for arg in libraries if not arg.endswith('/libraikou_mfem.a')]]
serial_source=EVIDENCE/'serial_state_trace.cpp'
serial_trace={}
for label,inc,libs in [('baseline',baseline_includes,baseline_libraries),('repository',includes,libraries)]:
    executable=EVIDENCE/(label+'-serial-state')
    run(label+'-serial-state-build',['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG',*inc,serial_source,'-o',executable,*libs])
    serial_trace[label]=run(label+'-serial-state',[executable,EVIDENCE/(label+'-checkpoint')])
assert serial_trace['baseline']==serial_trace['repository'],'Repository serial state differs from the original repository'

mpi_includes,mpi_libraries=flags(ROOT/'mpi-build')
mpi_traces={}
for label,inc,libs in [('baseline',baseline_includes,baseline_libraries),('repository',mpi_includes,mpi_libraries)]:
    executable=EVIDENCE/(label+'-mpi-pic')
    run(label+'-mpi-pic-build',['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG',*inc,EVIDENCE/'mpi_pic_trace.cpp','-o',executable,*libs])
    for ranks in [1,2]:
        key=label+'-mpi-'+str(ranks)
        mpi_traces[key]=run(key,['mpiexec','--bind-to','none','--oversubscribe','-n',str(ranks),executable,EVIDENCE/(key+'-record')])
for ranks in [1,2]:assert mpi_traces['baseline-mpi-'+str(ranks)]==mpi_traces['repository-mpi-'+str(ranks)],'MPI state differs'

native_traces={}
for label,source,inc,libs in [('baseline',BASELINE,baseline_includes,baseline_libraries),('repository',REPO,includes,libraries)]:
    executable=EVIDENCE/(label+'-native-fields')
    run(label+'-native-build',['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG',
                              '-DHISTORY_SOURCE="'+str(source/'examples/mfem_heating_history_pic.cpp')+'"',*inc,
                              EVIDENCE/'field_cache_trace.cpp','-o',executable,*libs])
    native_traces[label]=run(label+'-native-fields',[executable,EVIDENCE/(label+'-native-fixture'),EVIDENCE/(label+'-native-checksum.json'),label])
assert native_traces['baseline']==native_traces['repository'],'Native fields or time cache changed results'
checksums=[json.loads((EVIDENCE/(label+'-native-checksum.json')).read_text())['checksum'] for label in ['baseline','repository']]
assert checksums[0]==checksums[1]

abi={}
for label,inc in [('baseline',baseline_includes),('repository',includes)]:
    executable=EVIDENCE/('abi-'+label)
    run('abi-'+label+'-build',['/usr/sbin/c++','-std=c++20','-O3','-DNDEBUG',*inc,EVIDENCE/'abi_probe.cpp','-o',executable,*libraries])
    abi[label]=run('abi-'+label,[executable])
assert abi['baseline']==abi['repository'],'Checkpoint structure ABI changed'
run('git-diff-check',['git','diff','--check'],cwd=REPO)
summary={'status':'passed','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'logical_cpus':[6],'nice':15,'compiler_jobs':1,'lease_exclusion_confirmed':True,
         'serial_ctest_passed':True,'mpi_configuration_ctest_passed':True,
         'serial_state_byte_identical':True,'serial_state_bytes':len(serial_trace['repository'].encode()),
         'serial_state_sha256':hashlib.sha256(serial_trace['repository'].encode()).hexdigest(),
         'mpi_state_byte_identical_ranks':[1,2],
         'mpi_resources':'One MPI job at a time; both ranks inherit logical CPU 6, nice 15 and one OMP/BLAS thread.',
         'native_fields_byte_identical':True,'native_field_bytes':len(native_traces['repository'].encode()),
         'native_checksum_identical':True,'checkpoint_abi_identical':True,
         'gpu_work':False,'production_signals':False,'committed':False,'pushed':False}
(EVIDENCE/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Repository integration checks passed',flush=True)
