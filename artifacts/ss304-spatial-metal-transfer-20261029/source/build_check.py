from pathlib import Path
import fcntl,hashlib,json,os,subprocess,sys,time
S=Path(__file__).resolve().parent
P=S.parent
N=P/'prerequisites/native-v6_2'
O=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-spatial-metal-transfer-20261029')
B=P/'build'
B.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hold():
 p=Path('/proc/278902');s=(p/'stat').read_text().split(') ',1)[1].split()
 assert s[0]=='T' and s[19]=='5288339'
 assert os.readlink(p/'exe')=='/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/heating-history-pic-gpu'
 assert 'frozen 1' in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text()
 return {'pid':278902,'start_ticks':s[19],'state':s[0],'frozen':True}
os.sched_setaffinity(0,{7})
priority=os.getpriority(os.PRIO_PROCESS,0)
if priority<15:os.nice(15-priority)
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
stage=sys.argv[1]
frozen='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
base=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
inc=['-I'+str(S),'-I'+str(N/'source/include'),'-isystem',str(N/'inputs'),'-isystem',str(S/'vendor')]
for x in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:inc+=['-isystem',str(base/x)]
files={str(p.relative_to(S)):sha(p)for p in S.rglob('*')if p.is_file() and '__pycache__'not in p.parts}
identity=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()
(O/'source-identity.json').write_text(json.dumps({'sha256':identity,'files':files},indent=2)+'\n')
native_identity='5fa3a59f863e98f21e20af698b9c16417f3749344b555f18e3d2cba1e4f5a705'
compiler=[frozen,'/usr/bin/mpic++','-std=c++20','-O2','-fno-fast-math','-ffunction-sections','-fdata-sections']+inc+['-DSPATIAL_MATERIAL_SOURCE_SHA256="'+identity+'"','-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+native_identity+'"']
libs=[str(base/x)for x in ['install-mpi/lib/libmfem.a','tpl/hypre/lib/libHYPRE.a','tpl/metis/lib/libmetis.a','tpl/GKlib/lib/libGKlib.a']]+['-lcrypto']
sources={'native':N/'source/src/mfem_gap_pic.cpp','arena':N/'source/src/mfem_pic_arena.cpp','spatial':S/'test_spatial.cpp','regression':N/'tests/test_arena.cpp'}
if stage in sources:
 cmd=compiler+(['-Wall','-Wextra','-Wpedantic','-Werror','-Wno-misleading-indentation']if stage=='spatial'else[])+['-c',str(sources[stage]),'-o',str(B/(stage+'.o'))]
elif stage in ['link-spatial','link-regression']:
 name=stage.removeprefix('link-');cmd=[frozen,'/usr/bin/mpic++','-Wl,--gc-sections',str(B/'native.o'),str(B/'arena.o'),str(B/(name+'.o'))]+libs+['-o',str(B/('test_'+name))]
elif stage in ['run','restore']:
 (O/'native').mkdir(exist_ok=True);cmd=[frozen,str(B/'test_spatial'),str(O/'native')]+(['restore']if stage=='restore'else[])
elif stage=='regression-run':
 (O/'regression').mkdir(exist_ok=True);cmd=[frozen,str(B/'test_regression'),str(N/'inputs/h2o-gas-20261029.json'),str(O/'regression')]
elif stage=='highprecision':cmd=['python3',str(S/'audit.py'),str(O/'native/tests.json'),str(O/'highprecision.json')]
else:raise ValueError(stage)
with (O.parent/'pic-profile/testing-core7.lock').open('a')as lease:
 fcntl.flock(lease,fcntl.LOCK_EX)
 before=hold();start=time.time()
 pins=json.loads((S/'prerequisites.json').read_text())
 assert pins['native_source_identity']==native_identity
 assert pins['native_files']=={str(p.relative_to(N)):sha(p)for p in N.rglob('*')if p.is_file()}
 assert pins['material_files']=={p.name:sha(p)for p in (S/'material').iterdir()if p.is_file()}
 assert pins['copied_exact_currency_sha256']==sha(S/'dyadic.hpp')
 with (O/(stage+'.log')).open('w')as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=55)
 after=hold()
 assert files=={str(p.relative_to(S)):sha(p)for p in S.rglob('*')if p.is_file() and '__pycache__'not in p.parts}
 receipt={'stage':stage,'returncode':r.returncode,'command':cmd,'wall_s':time.time()-start,'cpu':[7],'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':1,'GPU':False,'production_before':before,'production_after':after,'source_identity':identity}
 (O/(stage+'-receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({k:receipt[k]for k in ['stage','returncode','wall_s','source_identity']}))
 if r.returncode:print((O/(stage+'.log')).read_text()[-12000:])
 sys.exit(r.returncode)
