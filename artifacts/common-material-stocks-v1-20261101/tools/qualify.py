from pathlib import Path
import hashlib,json,os,resource,subprocess,sys,time,runpy
P=Path(__file__).resolve().parents[1];R=P.parent;F=R/'arch-base-20261028';H=R/'hall-corrective-v3.3-b3-20261101';G=R/'gas-target-owner-v1.7.1-census-transport-20261101';N=R/'pic-reaction-transactions-20261029/native-full-time-v8.7-reservoir-stock-20261101';M=R/'mfem-4.10-pic-20261029';TRACE=R/'hall-coupled-transport-20261029/source/development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(x,s):
 if not x:raise RuntimeError(s)
def verify(pins):
 own={str(p.relative_to(P)) for d in ('include','tests','tools') for p in (P/d).rglob('*') if p.is_file() and p.suffix in ('.hpp','.cpp','.py')}
 need(own==set(pins['source_files']),'source inventory')
 for section in ('source_files','dependencies'):
  for p,h in pins[section].items():need(sha(P/p)==h,'changed input '+p)
 need(sha('/usr/bin/python3')==pins['frozen_python_sha256'],'frozen Python');need(sha('/usr/bin/c++')==pins['frozen_compiler_sha256'],'frozen compiler')
def main():
 job=os.getenv('SIMQ_JOB_ID');need(job and os.getenv('SIMQ_LAUNCH_TOKEN') and os.getenv('SIMQ_PARTITION')=='powerflow','SimQueue powerflow required');need(len(os.sched_getaffinity(0))==1 and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text(),'one allocated CPU/cgroup');need(os.getpriority(os.PRIO_PROCESS,0)>=15 and os.getenv('CUDA_VISIBLE_DEVICES')=='','CPU nice policy');need(not os.getenv('PF_BUILD_ISA_FLAGS') and all(os.getenv(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']),'build/thread policy')
 root=Path('/').stat();frozen=(F/'root').stat();need((root.st_dev,root.st_ino)==(frozen.st_dev,frozen.st_ino),'frozen userspace')
 pins=json.loads((P/'INPUTS.json').read_text());verify(pins);out=P/'evidence'/sys.argv[1];out.mkdir(exist_ok=False);commands=[];start=time.monotonic()
 def run(args,name):
  proc=subprocess.run([str(x) for x in args],cwd=P,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=480);(out/(name+'.stdout')).write_bytes(proc.stdout);(out/(name+'.stderr')).write_bytes(proc.stderr);commands.append({'argv':[str(x) for x in args],'returncode':proc.returncode});need(proc.returncode==0,name+': '+proc.stderr.decode()[-5000:]);return proc.stdout
 SOURCE=runpy.run_path(str(R/'pic-electrode-heat-source-v7.4-mass-cache-20261101/run_source_v7.py'),run_name='metadata');os.environ['PF_V8_CHEM']='v3.1'
 base=json.loads((R/'pic-water-common-commit-v8.0.4-b3-20261101/NATIVE-FLAGS.json').read_text())['flags'];base=[base[0],'-I'+str(P/'include'),'-I'+str(R/'common-conservation-v1-20261101/include'),'-I'+str(R/'conductor-correction-guard-v1-20261101/include'),'-I'+str(R/'pic-water-common-commit-v8.0.4-b3-20261101/source'),'-I'+str(R/'pic-water-common-commit-v8.0.4-b3-20261101/vendor'),'-I'+str(R/'pic-arrival-owner-v1.2-publication-20261101/include'),'-I'+str(R/'pic-source-stock-reference-v2-surface-birth-20261101/include'),'-I'+str(H/'source'),'-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),*base[1:],*SOURCE['owner_flags']([])[1:]]
 libs=[M/'install-mpi/lib/libmfem.a',M/'tpl/hypre/lib/libHYPRE.a',M/'tpl/metis/lib/libmetis.a',M/'tpl/GKlib/lib/libGKlib.a'];results={}
 for variant in ('normal','ubsan'):
  flags=base+(['-fsanitize=undefined','-fno-sanitize-recover=all','-fno-omit-frame-pointer'] if variant=='ubsan' else []);run([*flags,'-MD','-MF',out/(variant+'.d'),'-c',P/'tests/test_material.cpp','-o',out/(variant+'.o')],variant+'-compile')
  run([*flags,'-MD','-MF',out/(variant+'-motion.d'),'-c',H/'source/motion_provider.cpp','-o',out/(variant+'-motion.o')],variant+'-motion-compile')
  run([base[0],'-Wl,--gc-sections',*(['-fsanitize=undefined'] if variant=='ubsan' else []),out/(variant+'.o'),out/(variant+'-motion.o'),N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*libs,'-lcrypto','-ldl','-o',out/variant],variant+'-link')
  raw=run([out/variant,TRACE],variant+'-controls');need(run([out/variant,TRACE],variant+'-replay')==raw,'deterministic replay');results[variant]=json.loads(raw);need(results[variant]['passed'],'controls');
 need(results['normal']==results['ubsan'],'UBSan exact result parity');verify(pins)
 receipt={'passed':True,'job':job,'source_identity':pins['source_identity'],'inputs_sha256':sha(P/'INPUTS.json'),'allocated_cpu':sorted(os.sched_getaffinity(0)),'elapsed_seconds':time.monotonic()-start,'maximum_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'variants':results,'commands':commands,'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},'full_system_qualified':False,'production_qualified':False};(out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({'passed':True,'job':job,'variants':results}))
if __name__=='__main__':main()
