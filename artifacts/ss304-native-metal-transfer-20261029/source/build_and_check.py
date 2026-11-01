from pathlib import Path
import fcntl,hashlib,json,os,subprocess,time
S=Path(__file__).resolve().parent;P=S.parent;O=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-native-metal-transfer-20261029')
lock=O.parent/'pic-profile/testing-core7.lock'
os.sched_setaffinity(0,{7});priority=os.getpriority(os.PRIO_PROCESS,0)
if priority<15:os.nice(15-priority)
assert os.getpriority(os.PRIO_PROCESS,0)==15
for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[n]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hold():
 v=Path('/proc/278902/stat').read_text().split(') ',1)[1].split();e=Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text();assert v[0]=='T' and v[19]=='5288339' and 'frozen 1' in e
 return {'pid':278902,'state':v[0],'start_ticks':v[19],'frozen':True}
with lock.open('a')as lease:
 fcntl.flock(lease,fcntl.LOCK_EX);before=hold();start=time.time();B=P/'build';B.mkdir(exist_ok=True)
 for row in json.loads((O/'input-manifest.json').read_text()):assert sha(row['path'])==row['sha256']
 material_seal=json.loads(Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029/source-manifest.json').read_text())
 for p in (S/'material').iterdir():assert sha(p)==material_seal['files']['include/'+p.name]
 files={str(p.relative_to(S)):sha(p) for p in S.rglob('*') if p.is_file()};identity=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()
 (O/'source-identity.json').write_text(json.dumps({'sha256':identity,'files':files},indent=2)+'\n')
 frozen='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen';base=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
 inc=['-I'+str(S),'-I'+str(S/'arena/include'),'-isystem',str(S/'inputs'),'-isystem',str(S/'vendor')]
 for x in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:inc+=['-isystem',str(base/x)]
 compiler=[frozen,'/usr/bin/mpic++','-std=c++20','-O2','-fno-fast-math','-ffunction-sections','-fdata-sections']+inc+['-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"']
 commands=[]
 for src,name in [(S/'arena/src/mfem_gap_pic.cpp','native'),(S/'arena/src/mfem_pic_arena.cpp','arena'),(S/'test_transfer.cpp','transfer'),(S/'test_existing_arena.cpp','existing')]:
  cmd=compiler+(['-Wall','-Wextra','-Wpedantic','-Werror']if name!='native'else[])+['-c',str(src),'-o',str(B/(name+'.o'))];commands.append(cmd);print('compile',name,flush=True)
  with (O/('compile-'+name+'.log')).open('w')as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
 libs=[str(base/x) for x in ['install-mpi/lib/libmfem.a','tpl/hypre/lib/libHYPRE.a','tpl/metis/lib/libmetis.a','tpl/GKlib/lib/libGKlib.a']]+['-lcrypto']
 for name in ['transfer','existing']:
  binary=B/('test_'+name);cmd=[frozen,'/usr/bin/mpic++','-Wl,--gc-sections',str(B/'native.o'),str(B/'arena.o'),str(B/(name+'.o'))]+libs+['-o',str(binary)];commands.append(cmd)
  with (O/('link-'+name+'.log')).open('w')as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
  out=O/name;out.mkdir(exist_ok=True);cmd=[frozen,str(binary)]+([str(S/'inputs/h2o-gas-20261029.json')]if name=='existing'else[])+[str(out)];commands.append(cmd);print('run',name,flush=True)
  with(O/('run-'+name+'.log')).open('w')as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
 after=hold()
 for row in json.loads((O/'input-manifest.json').read_text()):assert sha(row['path'])==row['sha256']
 assert files=={str(p.relative_to(S)):sha(p) for p in S.rglob('*') if p.is_file()}
 result={'status':'native_tests_passed_pending_independent_highprecision_audit','commands':commands,'source_sha256':identity,'executables':{x:sha(B/('test_'+x))for x in ['transfer','existing']},'production_before':before,'production_after':after,'CPU':[7],'nice':15,'threads':1,'GPU':False,'wall_s':time.time()-start}
 (O/'build-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['status','source_sha256','wall_s']}))
