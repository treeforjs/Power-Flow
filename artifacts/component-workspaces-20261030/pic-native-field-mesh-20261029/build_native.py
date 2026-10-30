from pathlib import Path
import datetime,fcntl,hashlib,json,os,subprocess
P=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-integration')
R=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/reaction-transactions-20261029/native-integration')
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(R.parent.parent/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
B=P/'build';B.mkdir(exist_ok=True);T=R/'tests';T.mkdir(exist_ok=True)
frozen='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen';base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029'
includes=['-I'+str(P/'source/include'),'-I'+str(P/'reaction-core/include'),'-I'+str(P/'tests'),'-I/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/source/vendor']
for path in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:includes+=['-isystem',base+'/'+path]
compiler=[frozen,'/usr/bin/mpic++','-std=c++17','-O2','-g','-ffunction-sections','-fdata-sections']+includes
source_hashes={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest()for directory in ['source','reaction-core']for f in sorted((P/directory).rglob('*'))if f.is_file()}
source_identity=hashlib.sha256(json.dumps(source_hashes,sort_keys=True).encode()).hexdigest()
(R/'native-source-identity.json').write_text(json.dumps({'sha256':source_identity,'files':source_hashes},indent=2)+'\n')
compiler+=['-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+source_identity+'"']
commands=[];objects=[]
for source,name,strict in [(P/'source/src/mfem_gap_pic.cpp','native',False),(P/'reaction-core/src/reaction_transactions.cpp','transactions',True),(P/'source/src/mfem_pic_reaction_batch.cpp','adapter',True),(P/'tests/test_native_reactions.cpp','tests',True)]:
 command=compiler+(['-Wall','-Wextra','-Wpedantic','-Werror']if strict else[])+['-c',str(source),'-o',str(B/(name+'.o'))];commands.append(command)
 with(R/('build-'+name+'.log')).open('w')as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
 objects.append(str(B/(name+'.o')))
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto']
binary=B/'test_native_reactions';command=[frozen,'/usr/bin/mpic++','-Wl,--gc-sections']+objects+libs+['-o',str(binary)];commands.append(command)
with(R/'link.log').open('w')as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
run=[frozen,str(binary),str(T)]
with(R/'native-tests.log').open('w')as log:subprocess.run(run,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
results=json.loads((T/'native-tests.json').read_text());assert results['status']=='passed'and all(c['passed']for c in results['checks'])
(R/'commands.json').write_text(json.dumps({'build':commands,'run':run,'resource':{'cpu':7,'nice':15,'threads':1,'jobs':1,'lease':str(R.parent.parent/'testing-core7.lock')}},indent=2)+'\n')
(R/'build-result.json').write_text(json.dumps({'status':'passed','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':results['check_count'],'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'GPU':False,'native_CPU_field_and_tiny_step_tests':True,'physical_H2O_rates':False},indent=2)+'\n')
print('Native reaction/Poisson checks passed:',results['check_count'])
