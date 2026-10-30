#!/usr/bin/env python3
import fcntl,hashlib,json,os,resource,subprocess,time
from pathlib import Path
P=Path(__file__).resolve().parent
R=Path('/home/treefor/projects/runtimes')
B=Path('<workspace>/outputs')
O=B/'pic-water-common-commit-v3-20261030'
N=R/'pic-reaction-transactions-20261029/native-full-time-candidate003'
I=R/'pic-water-full-time-subflow-20261029/inputs'
M=R/'mfem-4.10-pic-20261029'
F=R/'arch-base-20261028/run-frozen'
if not __debug__ or os.getenv('PYTHONOPTIMIZE'):raise SystemExit('Optimized guards refused')
if os.sched_getaffinity(0)!={7} or os.getpriority(os.PRIO_PROCESS,0)!=15:raise SystemExit('CPU7 nice15 required')
lease=open(B/'pic-profile/testing-core7.lock','a')
try:fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit(75)
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
env={'PATH':'/usr/bin:/bin','LC_ALL':'C','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
     'BLIS_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','OMP_THREAD_LIMIT':'1',
     'CUDA_VISIBLE_DEVICES':'','PYTHONNOUSERSITE':'1'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
native_evidence=B/'pic-profile/reaction-transactions-20261029/native-full-time-candidate003'
identity='c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'
lib=N/'build/libpic_native_candidate003.a'
if sha(lib)!='d2568af9fac13f8aa056ae4a9a50bb2a07eee05d7218d1c4b79aef21317d8a20':raise RuntimeError('Unqualified native library')
manifest=json.loads((N/'SOURCE-SHA256.json').read_text())
if manifest['source_identity']!=identity:raise RuntimeError('Native source identity differs')
for name,h in manifest['files'].items():
    if sha(N/name)!=h:raise RuntimeError('Native source file changed')
baseline=B/'pic-water-exact-count-subflow-20261029'
data=R/'pic-water-exact-count-subflow-20261029/inputs/h2o-gas.json'
binding=baseline/'run-binding.json'
libs=[M/x for x in ('install-mpi/lib/libmfem.a','tpl/hypre/lib/libHYPRE.a','tpl/metis/lib/libmetis.a','tpl/GKlib/lib/libGKlib.a')]
trees=[P/'source',P/'tests',N/'source/include',N/'inputs',I,M/'install-mpi/include',M/'tpl/hypre/include',M/'tpl/metis/include']
def pins():
    files={p for tree in trees for p in tree.rglob('*') if p.is_file()}
    files.update([Path(__file__).resolve(),F,lib,*libs,data,binding,baseline/'result.json',N/'SOURCE-SHA256.json',
                  native_evidence/'build-artifact-hashes.json',R/'arch-base-20261028/root/usr/bin/mpic++',R/'arch-base-20261028/root/usr/bin/g++'])
    return {str(p):sha(p) for p in sorted(files)}
before=pins();A=O/f'parent-regression-{len(list(O.glob("parent-regression-*")))+1:03d}';A.mkdir()
(A/'inputs.json').write_text(json.dumps(before,indent=2)+'\n')
flags=['-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math',
       '-ffunction-sections','-fdata-sections','-I'+str(P/'source'),'-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"']
for d in (N/'source/include',I,I/'include',N/'inputs',M/'install-mpi/include',M/'tpl/hypre/include',M/'tpl/metis/include'):
    flags+=['-isystem',str(d)]
exe=P/'build/test_native_parent_subflow'
commands=[[str(F),'/usr/bin/mpic++',*flags,str(P/'tests/test_native_parent_subflow.cpp'),str(lib),
           *map(str,libs),'-Wl,--gc-sections','-lcrypto','-o',str(exe)],
          [str(F),str(exe),str(data),str(binding),str(A),identity]]
(A/'commands.json').write_text(json.dumps({'commands':commands,'environment':env},indent=2)+'\n')
start=time.monotonic();status='failed';error=None;comparison=None
try:
    for index,command in enumerate(commands):
        with (A/f'command-{index}.log').open('w') as log:
            subprocess.run(command,env=env,cwd=P,stdout=log,stderr=subprocess.STDOUT,timeout=60,check=True)
    actual=json.loads((A/'result.json').read_text());old=json.loads((baseline/'result.json').read_text())
    def events(j):return [[{k:v for k,v in e.items() if k!='native'} for e in receipt['events']] for receipt in j['receipts']]
    comparison={'checks':actual['check_count'],'first_events':actual['first_events'],
      'all_check_results_equal':actual['checks']==old['checks'],'recorded_event_physics_equal':events(actual)==events(old),
      'recorded_event_counts':[len(x['events']) for x in actual['receipts']],
      'physical_trajectory_advanced':actual['physical_trajectory_advanced'],'full_coupled_run':actual['full_coupled_run']}
    if actual['status']!='passed' or comparison['checks']!=442 or comparison['first_events']!=old['first_events'] or not comparison['all_check_results_equal'] or not comparison['recorded_event_physics_equal'] or comparison['physical_trajectory_advanced'] or comparison['full_coupled_run']:
        raise RuntimeError('Legacy parent regression differs')
    if pins()!=before:raise RuntimeError('Regression dependencies changed')
    status='passed'
except Exception as exc:error=repr(exc)
receipt={'status':status,'error':error,'comparison':comparison,'elapsed_s':time.monotonic()-start,
    'native_source_identity':identity,'cpu':7,'nice':15,'threads':1,'exclusive_lease':True,
    'GPU_used':False,'full_time_joint_executed':False,'physical_time_advanced_s':0,
    'scope':'legacy frozen-stage parent operator regression; includes existing fixture Poisson solves, no Hall coupling',
    'input_manifest_sha256':sha(A/'inputs.json')}
(A/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'attempt':str(A),**receipt}))
if status!='passed':raise SystemExit(1)
