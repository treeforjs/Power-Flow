from pathlib import Path
import hashlib,json,os,resource,subprocess,sys,time
R=Path(__file__).resolve().parents[1]
G=R.parent/'gas-target-owner-v1.7.1-census-transport-20261101'
F=R.parent/'arch-base-20261028'
FORMS=R.parent/'pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix/forms/include/pf_v7'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(x,s):
 if not x:raise RuntimeError(s)
def verify(pins):
 own={str(p.relative_to(R)) for d in ('include','tests','tools') for p in (R/d).rglob('*') if p.is_file() and p.suffix in ('.hpp','.cpp','.py')}
 need(own==set(pins['source_files']),'source inventory')
 for section in ('source_files','dependencies'):
  for p,h in pins[section].items():need(sha(R/p)==h,'changed input '+p)
 need(sha('/usr/bin/python3')==pins['frozen_python_sha256'],'frozen Python')
 need(sha('/usr/bin/c++')==pins['frozen_compiler_sha256'],'frozen compiler')
def main():
 job=os.getenv('SIMQ_JOB_ID');need(job and os.getenv('SIMQ_LAUNCH_TOKEN') and os.getenv('SIMQ_PARTITION')=='powerflow','SimQueue powerflow required')
 need(len(os.sched_getaffinity(0))==1 and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text(),'one allocated CPU/cgroup')
 need(os.getpriority(os.PRIO_PROCESS,0)>=15 and os.getenv('CUDA_VISIBLE_DEVICES')=='','CPU nice policy')
 need(not os.getenv('PF_BUILD_ISA_FLAGS') and all(os.getenv(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']),'build/thread policy')
 root=Path('/').stat();frozen=(F/'root').stat();need((root.st_dev,root.st_ino)==(frozen.st_dev,frozen.st_ino),'frozen userspace')
 pins=json.loads((R/'INPUTS.json').read_text());verify(pins);out=R/'evidence'/sys.argv[1];out.mkdir(exist_ok=False);commands=[];start=time.monotonic()
 def run(args,name):
  p=subprocess.run([str(x) for x in args],cwd=R,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=120)
  (out/(name+'.stdout')).write_bytes(p.stdout);(out/(name+'.stderr')).write_bytes(p.stderr);commands.append({'argv':[str(x) for x in args],'returncode':p.returncode})
  need(p.returncode==0,name+': '+p.stderr.decode()[-4000:]);return p.stdout
 flags=['-std=c++20','-O1','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(R/'include'),'-I'+str(G/'include'),'-I'+str(FORMS)]
 cmd=['/usr/bin/c++',*flags,R/'tests/test_guard.cpp',G/'src/gas_owner.cpp','-o',out/'controls']
 run(cmd,'compile');first=run([out/'controls'],'controls');need(run([out/'controls'],'replay')==first,'deterministic controls');parsed=json.loads(first);need(parsed['passed'],'controls passed')
 san=cmd.copy();san[2:2]=['-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-sanitize-recover=all'];san[-1]=out/'sanitized'
 run(san,'compile-sanitized');os.environ['ASAN_OPTIONS']='detect_leaks=1:halt_on_error=1';os.environ['UBSAN_OPTIONS']='halt_on_error=1:print_stacktrace=1'
 need(run([out/'sanitized'],'sanitized')==first,'sanitizer parity');verify(pins)
 receipt={'passed':True,'job':job,'source_identity':pins['source_identity'],'inputs_sha256':sha(R/'INPUTS.json'),'allocated_cpu':sorted(os.sched_getaffinity(0)),'elapsed_seconds':time.monotonic()-start,'maximum_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'controls':parsed,'commands':commands,'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},'native_integrated':False,'production_qualified':False}
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({'passed':True,'job':job,'controls':parsed}))
if __name__=='__main__':main()
