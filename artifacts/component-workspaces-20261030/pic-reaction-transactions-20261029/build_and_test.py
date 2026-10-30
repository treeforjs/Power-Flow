from pathlib import Path
import datetime,fcntl,hashlib,json,os,subprocess
PREFIX=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029')
REPORT=Path('<workspace>/outputs/pic-profile/reaction-transactions-20261029')
LEASE=REPORT.parent/'testing-core7.lock'
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
with LEASE.open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
REPORT.mkdir(parents=True,exist_ok=True);(PREFIX/'build').mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[Path('<workspace>/outputs/supplied-chicago-papers-20261029/H2O-readiness.txt'),
 REPORT.parent/'fragmentation-raikou-plan/physics-data-handoff.txt',REPORT.parent/'energy-work-ownership-audit/DESIGN-BRIEF.txt',REPORT.parent/'energy-work-ownership-audit/ENERGY-WORK-DESIGN.txt',
 Path('/home/treefor/projects/raikou-mfem/include/raikou/mfem_gap_pic.h'),Path('/home/treefor/projects/raikou-mfem/src/mfem_gap_pic.cpp'),
 Path('/home/treefor/projects/raikou/include/raikou/electron_coupling.h'),Path('/home/treefor/projects/raikou/include/raikou/json.hpp'),
 Path('/home/treefor/projects/raikou/mechanisms/h2o_research.json'),
 Path('/home/treefor/projects/Power-Flow/mhdlab/coverage_inventory.py'),Path('/home/treefor/projects/Power-Flow/mhdlab/neutral_cr_water_source.py'),
 Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include/raikou/mfem_gap_pic.h'),
 Path('/home/treefor/projects/runtimes/arch-base-20261028/root/usr/lib/libcrypto.so.3')]
before={str(p):sha(p)for p in inputs}
assert sha(PREFIX/'source/vendor/raikou/json.hpp')==sha(Path('/home/treefor/projects/raikou/include/raikou/json.hpp'))
(REPORT/'input-manifest.json').write_text(json.dumps({'before':before},indent=2)+'\n')
frozen='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
base=[frozen,'/usr/bin/mpic++','-std=c++17','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror','-ffunction-sections','-fdata-sections','-Wl,--gc-sections',
 '-I'+str(PREFIX/'source/include'),'-I'+str(PREFIX/'source/vendor'),
 '-I/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include',
 '-isystem','/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/include',
 '-isystem','/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/hypre/include',
 '-isystem','/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/tpl/metis/include',
 str(PREFIX/'source/src/reaction_transactions.cpp'),str(PREFIX/'source/tests/test_reaction_transactions.cpp'),'-lcrypto']
binary=PREFIX/'build/test_reaction_transactions';command=base+['-o',str(binary)]
controls={'cpu':7,'nice':15,'threads':1,'lease':str(LEASE),'jobs':1,'ledger_roundoff_factor':256,'manufactured_only':True,'physical_rates_or_thresholds_bound':False,'build':command}
(REPORT/'controls.json').write_text(json.dumps(controls,indent=2)+'\n')
with(REPORT/'build.log').open('w')as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
run=[frozen,str(binary),str(REPORT/'tests')]
with(REPORT/'tests.log').open('w')as log:subprocess.run(run,stdout=log,stderr=subprocess.STDOUT,check=True)
result=json.loads((REPORT/'tests/test-results.json').read_text());assert result['status']=='passed'and all(c['passed']for c in result['checks'])
after={str(p):sha(p)for p in inputs};assert after==before
(REPORT/'input-manifest.json').write_text(json.dumps({'before':before,'after':after,'unchanged':True},indent=2)+'\n')
(REPORT/'build-test-evidence.json').write_text(json.dumps({'status':'passed','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'build':command,'run':run,'binary_sha256':sha(binary),'checks':result['check_count'],'events':result['event_count'],'inputs_unchanged':True,'native_pic_evolution':False,'gpu_work':False},indent=2)+'\n')
print('Passed',result['check_count'],'manufactured checks,',result['event_count'],'transactions',flush=True)
