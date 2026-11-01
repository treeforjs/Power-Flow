#!/usr/bin/python -I
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
R=Path(__file__).resolve().parents[1];OUT=Path('/home/treefor/projects/.workspace/outputs')/R.name
mode,label=sys.argv[1:];assert mode in ('normal','sanitized')
F=R.parent/'arch-base-20261028';a=Path('/').stat();b=(F/'root').stat();assert (a.st_dev,a.st_ino)==(b.st_dev,b.st_ino)
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_UNIT') in Path('/proc/self/cgroup').read_text()
cpus=sorted(os.sched_getaffinity(0));assert len(cpus)==1 and os.getenv('SIMQ_CPU_LIST')==str(cpus[0]) and (cpus[0] in (6,7) or cpus[0] in (4,5) and os.getenv('SIMQ_BORROW')=='1')
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert Path('/proc/sys/kernel/osrelease').read_text().strip()==(F/'kernel.release').read_text().strip()
assert os.getenv('SIMQ_PARTITION')=='powerflow' and int(os.environ['SIMQ_MEM'])>=2*1024**3
D=OUT/label;D.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def hashes():return {str(p):sha(p) for name in ('fast','baseline-fast','reference-kernel','tests','tools') for p in sorted((R/name).rglob('*')) if p.is_file()}
before=hashes();rec={'passed':False,'mode':mode,'job':os.environ['SIMQ_JOB_ID'],'affinity':cpus,'input_hashes':before,'frozen_root_identity':[a.st_dev,a.st_ino],'commands':[],'binaries':{},'speedup_qualified':False,'runtime_manifest_sha256':sha(F/'root-files.sha256'),'launcher_sha256':sha(F/'run-frozen'),'compiler_sha256':sha(Path('/usr/bin/c++'))}
def run(name,cmd,binary=False):
 t=time.monotonic();p=subprocess.run(list(map(str,cmd)),stdout=subprocess.PIPE,stderr=subprocess.PIPE);(D/(name+'.stdout')).write_bytes(p.stdout);(D/(name+'.stderr')).write_bytes(p.stderr);rec['commands'].append({'name':name,'argv':list(map(str,cmd)),'returncode':p.returncode,'wall_s':time.monotonic()-t});print(name,p.returncode,flush=True)
 if p.returncode:raise RuntimeError((p.stderr or p.stdout).decode(errors='replace')[-4000:])
 return p.stdout
flags=['/usr/bin/c++','-std=c++20','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-ffunction-sections','-fdata-sections','-Wl,--gc-sections','-I',R/'reference-kernel/include','-I',R/'tests']
flags+=['-O2'] if mode=='normal' else ['-O1','-g','-fsanitize=address,undefined','-fno-sanitize-recover=all','-fno-omit-frame-pointer']
try:
 run('build-tests',[*flags,'-I',R/'fast',R/'tests/test_mass_cache.cpp',R/'fast/fast_emission.cpp',R/'reference-kernel/src/surface_emission.cpp','-o',D/'test_mass_cache'])
 rec['binaries']['test_mass_cache']=sha(D/'test_mass_cache');run('tests',[D/'test_mass_cache'])
 if mode=='normal':
  traces=[]
  for name,path in [('old',R/'baseline-fast'),('new',R/'fast')]:
   run('build-legacy-'+name,[*flags,'-I',path,R/'tests/legacy_trace.cpp',path/'fast_emission.cpp',R/'reference-kernel/src/surface_emission.cpp','-o',D/('legacy-'+name)])
   rec['binaries']['legacy-'+name]=sha(D/('legacy-'+name));traces.append(run('legacy-'+name,[D/('legacy-'+name)]))
  assert traces[0]==traces[1];rec['legacy_trace_sha256']=hashlib.sha256(traces[0]).hexdigest()
 rec['inputs_unchanged']=hashes()==before;assert rec['inputs_unchanged'];assert all(sha(D/p)==h for p,h in rec['binaries'].items());rec['passed']=True
finally:(D/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
