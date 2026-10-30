from pathlib import Path
import os,sys,fcntl,hashlib,json,subprocess,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PREFIX=ROOT.parent
OUT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/hall-coupled-transport-20261029/finite-foil-em')
launcher='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
lease='/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core.lock'
assert os.sched_getaffinity(0)=={6} and os.getpriority(os.PRIO_PROCESS,0)==15
assert all(os.getenv(v)=='1' for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')) and os.getenv('CUDA_VISIBLE_DEVICES')==''
with open(lease,'a') as f:
 try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('Acquire test lease first')
run=OUT/sys.argv[1];run.mkdir(parents=True,exist_ok=False)
build=PREFIX/'build-axial-foil-em';build.mkdir(exist_ok=True);binary=build/'axial_foil_em'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
dep=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
libs=[dep/'install-mpi/lib/libmfem.a',dep/'tpl/hypre/lib/libHYPRE.a',dep/'tpl/metis/lib/libmetis.a',dep/'tpl/GKlib/lib/libGKlib.a']
files=[ROOT/'src/mfem_gap_pic.cpp',ROOT/'src/mfem_pic_distributed.hpp',ROOT/'include/raikou/mfem_gap_pic.h',ROOT/'include/raikou/q1_scalar_gather.h',HERE/'foil_em.hpp',HERE/'foil_particles.hpp',HERE/'surface_path_port.hpp',HERE/'native_total_field_transport.hpp',HERE/'surface-source-binding.json',HERE/'axial_benchmark.cpp',HERE/'controls.json',HERE/'run_axial.py',HERE/'axial_mesh.hpp',HERE/'axial_fields.hpp',HERE/'axial-contract.json',HERE/'native-start-binding.json',HERE.parent/'hall-coupled-transport-20261029/coupled_pic.hpp',*libs]
files += [HERE/'axial_particles.hpp',HERE/'axial_restart.hpp',HERE/'axial_refinement.hpp',HERE/'verify_axial_energy.py',HERE/'axial_surface_path_port.hpp',HERE/'axial_thermal.hpp',HERE/'axial_system.hpp',HERE/'axial_split_stage.hpp',*sorted((HERE/'thermal-owner').rglob('*.*')),*sorted((HERE/'data').glob('*'))]
manifest={str(p):sha(p) for p in files};(run/'input-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
flags=['-std=c++20','-O2','-DNDEBUG','-I'+str(ROOT/'include'),'-I'+str(ROOT/'src'),'-isystem',str(dep/'install-mpi/include'),'-isystem',str(dep/'tpl/hypre/include'),'-isystem',str(dep/'tpl/metis/include')]
commands=[('build',['/usr/bin/c++',*flags,str(ROOT/'src/mfem_gap_pic.cpp'),str(HERE/'axial_benchmark.cpp'),'-o',str(binary),*[str(p) for p in libs],'/usr/lib/libmpi.so']),('benchmark',[str(binary),str(run)])]
commands += [(label,[str(binary),str(run),label]) for label in ('refinement','full','split','resume')]
commands.append(('independent-energy',['/usr/bin/python3',str(HERE/'verify_axial_energy.py'),str(run)]))
receipts=[]
(run/'run-status.json').write_text(json.dumps({'status':'running'})+'\n')
for label,cmd in commands:
 actual=[launcher,'prlimit','--as=2147483648','--cpu=180','--',*cmd];start=time.monotonic()
 with (run/(label+'.log')).open('w') as log:result=subprocess.run(actual,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=240)
 receipts.append({'label':label,'command':actual,'elapsed_s':time.monotonic()-start,'returncode':result.returncode,'CPU':[6],'nice':15,'threads':1,'shared_lease':True})
 (run/'commands.json').write_text(json.dumps(receipts,indent=2)+'\n');print(label,result.returncode,receipts[-1]['elapsed_s'],flush=True)
 (run/'run-status.json').write_text(json.dumps({'status':'failed' if result.returncode else 'running','last_command':label,'returncode':result.returncode})+'\n')
 if result.returncode:raise RuntimeError((run/(label+'.log')).read_text()[-6000:])
assert all(sha(Path(p))==h for p,h in manifest.items())
relative = sorted(p.relative_to(run/'full') for p in (run/'full').rglob('*') if p.is_file())
assert all((run/'full'/p).read_bytes()==(run/'resume'/p).read_bytes() for p in relative), 'fresh-process restart mismatch'
validation=json.loads((run/'validation.json').read_text())
validation['checks'].append({'name':'fresh-process field arena current work and thermal restart exact','passed':True,'detail':{'files':len(relative)}})
validation['checks'].append({'name':'coupled time differences decrease and seven DC mesh/exterior sensitivities retain circuit reference','passed':True,'detail':json.loads((run/'refinement.json').read_text())['time']})
validation['checks'].append({'name':'independent 80-digit energy from checkpoint field state','passed':True,'detail':json.loads((run/'independent-energy.json').read_text())})
validation['count']+=3
(run/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
(run/'binary-sha256.txt').write_text(sha(binary)+'\n');print(run,flush=True)

(run/'run-status.json').write_text(json.dumps({'status':'passed','checks':validation['count'],'commands':len(commands)},indent=2)+'\n')
