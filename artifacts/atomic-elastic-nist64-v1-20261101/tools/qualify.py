from pathlib import Path
import hashlib,json,os,resource,shutil,subprocess,sys,time
R=Path(__file__).resolve().parents[1]
CHEM=Path('/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
LEGACY=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031/data/tables/chemistry-v3.1.blob')
SHA='8435e873b7c68146fb4c429d96c3e8051334c22525fd508729299e8a2d65ea0f'
LSHA='70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(x,s):
 if not x:raise RuntimeError(s)
def verify(pins):
 names={str(p.relative_to(R)) for d in ('include','tests','tools') for p in (R/d).rglob('*') if p.is_file() and p.suffix in ('.py','.hpp','.cpp')}
 need(names==set(pins['source_files']),'source set')
 for section in ('source_files','dependencies'):
  for p,h in pins[section].items():need(sha(R/p)==h,'changed input '+p)
 need(sha('/usr/bin/python3')==pins['frozen_python_sha256'],'frozen Python')
 need(sha('/usr/bin/c++')==pins['frozen_compiler_sha256'],'frozen compiler')
def main():
 job=os.environ.get('SIMQ_JOB_ID');need(job and os.environ.get('SIMQ_LAUNCH_TOKEN') and os.environ.get('SIMQ_PARTITION')=='powerflow','SimQueue powerflow required')
 need(len(os.sched_getaffinity(0))==1 and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text(),'allocated CPU and cgroup')
 need(os.getpriority(os.PRIO_PROCESS,0)>=15,'nice priority');need(os.environ.get('CUDA_VISIBLE_DEVICES')=='','CPU-only')
 need(not os.environ.get('PF_BUILD_ISA_FLAGS'),'unpinned build flags refused')
 need(all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']),'one thread')
 pins=json.loads((R/'INPUTS.json').read_text());verify(pins)
 out=R/'evidence'/sys.argv[1];out.mkdir(exist_ok=False);build=out/'build';build.mkdir();start=time.monotonic();commands=[]
 def run(cmd,name,expected=0,timeout=120):
  p=subprocess.run([str(x) for x in cmd],cwd=R,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=timeout)
  (out/(name+'.stdout')).write_bytes(p.stdout);(out/(name+'.stderr')).write_bytes(p.stderr);commands.append({'argv':[str(x) for x in cmd],'returncode':p.returncode})
  need(p.returncode==expected if expected==0 else p.returncode!=0,name+' unexpected result: '+p.stderr.decode()[-2000:]);return p.stdout
 args=['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(R/'include'),'-I'+str(CHEM/'include/pf_chem_v3'),'-I'+str(ARR),R/'tests/fixture.cpp','-lcrypto','-o',build/'fixture']
 run(args,'compile')
 inputs=json.loads((R/'source-pins.json').read_text());source=Path(inputs['source_dir'])
 invocation=[build/'fixture',source,CHEM/'data/chemistry-cdf.blob',SHA,'events']
 first=run(invocation,'events');need(run(invocation,'repeat')==first,'fresh replay')
 independent=json.loads(run([sys.executable,'-I',R/'tools/verify.py',out/'events.stdout'],'independent'));need(independent['passed'],'independent controls')
 legacy=run([build/'fixture',source,LEGACY,LSHA,'events'],'legacy');old=json.loads(legacy.splitlines()[0]);new=json.loads(first.splitlines()[0]);need(old['model_sha256']!=new['model_sha256'],'cross-blob provenance')
 need(legacy.splitlines()[1:]==first.splitlines()[1:],'identical registry physical operands')
 run([build/'fixture',source,CHEM/'data/chemistry-cdf.blob',LSHA,'load'],'wrong-blob',1)
 run([build/'fixture',source,CHEM/'data/chemistry-cdf.blob','0'*64,'load'],'unknown-blob',1)
 bad=build/'mutated-source';shutil.copytree(source,bad,ignore=shutil.ignore_patterns('acquisition.json'))
 fname=next(iter(inputs['files']));original=(bad/fname).read_bytes()
 for label,mutation in [('tamper',original.replace(b'1000 eV',b'1001 eV')),('oversized',b'x'*32769),('missing',None)]:
  if mutation is None:(bad/fname).unlink()
  else:(bad/fname).write_bytes(mutation)
  run([build/'fixture',bad,CHEM/'data/chemistry-cdf.blob',SHA,'load'],label,1)
 san=args.copy();san[san.index('-O2')]='-O1';san[2:2]=['-fsanitize=address,undefined','-fno-sanitize-recover=all','-fno-omit-frame-pointer'];san[-1]=build/'fixture-sanitized'
 run(san,'compile-sanitized',timeout=180)
 os.environ['ASAN_OPTIONS']='detect_leaks=1:halt_on_error=1';os.environ['UBSAN_OPTIONS']='halt_on_error=1:print_stacktrace=1'
 need(run([build/'fixture-sanitized',*invocation[1:]],'sanitized')==first,'sanitizer parity')
 verify(pins)
 receipt={'passed':True,'job':job,'source_identity':pins['source_identity'],'input_sha256':sha(R/'INPUTS.json'),'cpus':sorted(os.sched_getaffinity(0)),'elapsed_seconds':time.monotonic()-start,'maximum_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'independent':independent,'provenance_controls':6,'commands':commands,'binary_sha256':sha(build/'fixture'),'sanitized_binary_sha256':sha(build/'fixture-sanitized'),'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},'production_qualified':False,'network_complete':False,'physical_model_validated':False,'native_adopted':False}
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({'passed':True,'job':job,'independent':independent}))
if __name__=='__main__':main()
