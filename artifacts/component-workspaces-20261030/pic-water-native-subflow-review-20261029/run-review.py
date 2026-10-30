from pathlib import Path
import hashlib,json,os,subprocess,shutil
p=Path(__file__).parent;s=p/'snapshot';o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-native-subflow-review-20261029')
assert sorted(os.sched_getaffinity(0))==[7] and os.getpriority(os.PRIO_PROCESS,0)==15
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for rec in json.loads((o/'input-bindings.json').read_text()):assert sha(Path(rec['snapshot']))==rec['sha256']
fixture=(s/'tests/test_native_parent_subflow.cpp').read_text();(p/'fixture.hpp').write_text(fixture[:fixture.index('int main(')])
v5=Path('/home/treefor/projects/runtimes/pic-native-arena-review-v5-20261029')
for name in ('source/src/mfem_gap_pic.cpp','source/include/raikou/mfem_gap_pic.h','source/src/mfem_pic_distributed.hpp','source/include/raikou/q1_scalar_gather.h'):
 assert sha(s/'native'/name)==sha(v5/'snapshot'/name),name
shutil.copy2(v5/'build/native.o',p/'build/native.o')
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029';identity=json.loads((s/'native-handoff.json').read_text())['source_identity']
includes=['-I'+str(s/'source'),'-I'+str(s/'native/source/include'),'-isystem',str(s/'inputs'),'-isystem',str(s/'inputs/include'),'-isystem',str(s/'native/inputs'),'-I'+str(s/'native/vendor')]
for n in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+n]
compiler=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"']
commands=[]
for path,name in ((s/'native/source/src/mfem_pic_arena.cpp','arena'),(p/'probe.cpp','probe')):
 command=compiler+['-c',str(path),'-o',str(p/'build'/f'{name}.o')];commands.append(command)
 with (o/f'build-{name}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto']
command=['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/f'{n}.o') for n in ('native','arena','probe')]+libs+['-o',str(p/'build/probes')];commands.append(command);subprocess.run(command,check=True,stdout=(o/'link.log').open('w'),stderr=subprocess.STDOUT)
command=['/usr/bin/time','-v','-o',str(o/'resources.txt'),str(p/'build/probes'),str(s/'inputs/h2o-gas.json'),str(o/'owner-receipts/run-binding.json'),str(o/'probes.json')];commands.append(command);subprocess.run(command,check=True,stdout=(o/'probes.txt').open('w'),stderr=subprocess.STDOUT)
(o/'commands-independent.json').write_text(json.dumps({'commands':commands,'reused_native_object_sha256':sha(p/'build/native.o'),'reused_native_object_source':'independently compiled v5 foundation, four source/header dependency hashes identical to reviewed v4 foundation'},indent=2)+'\n')
print((o/'probes.txt').read_text())
