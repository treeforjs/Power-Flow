from pathlib import Path
import os,json,hashlib,subprocess,shutil,sys
p=Path(__file__).parent;s=p/'snapshot';o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-native-subflow-v6_1-review-20261029')
assert sorted(os.sched_getaffinity(0))==[7] and os.getpriority(os.PRIO_PROCESS,0)==15
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for rec in json.loads((o/'input-bindings.json').read_text()):assert sha(Path(rec['snapshot']))==rec['sha256']
text=(s/'tests/test_native_parent_subflow.cpp').read_text();(p/'fixture.hpp').write_text(text[:text.index('int main(')])
v5=Path('/home/treefor/projects/runtimes/pic-native-arena-review-v5-20261029')
for name in ('source/src/mfem_gap_pic.cpp','source/include/raikou/mfem_gap_pic.h','source/src/mfem_pic_distributed.hpp','source/include/raikou/q1_scalar_gather.h'):assert sha(s/'native'/name)==sha(v5/'snapshot'/name),name
shutil.copy2(v5/'build/native.o',p/'build/native.o')
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029';identity=json.loads((s/'native-handoff.json').read_text())['source_identity'];assert identity=='fd318a6a9ae081bdaf57656c7e47da3e528c9ea2b532224de15bab37b56f3dc6'
includes=['-I'+str(s/'source'),'-I'+str(s/'native/source/include'),'-isystem',str(s/'inputs'),'-isystem',str(s/'inputs/include'),'-isystem',str(s/'native/inputs'),'-I'+str(s/'native/vendor')]
for n in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+n]
compiler=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"']
commands=[]
for source,name in ((s/'native/source/src/mfem_pic_arena.cpp','arena'),(p/'probe.cpp','probes')):
 if name=='arena' and '--reuse-arena' in sys.argv:
  assert (p/'build/arena.o').exists()
  commands.append({'reuse_after_review_wrapper_compile_failure':str(p/'build/arena.o'),'object_sha256':sha(p/'build/arena.o'),'unchanged_bound_source_sha256':sha(source)})
  continue
 command=compiler+['-c',str(source),'-o',str(p/'build'/f'{name}.o')];commands.append(command)
 with (o/f'build-{name}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto']
command=['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/f'{n}.o') for n in ('native','arena','probes')]+libs+['-o',str(p/'build/probes')];commands.append(command);subprocess.run(command,check=True,stdout=(o/'link.log').open('w'),stderr=subprocess.STDOUT)
command=['/usr/bin/time','-v','-o',str(o/'resources.txt'),str(p/'build/probes'),str(s/'inputs/h2o-gas.json'),str(o/'owner-run-binding.json'),str(o)];commands.append(command);subprocess.run(command,check=True,stdout=(o/'probes.txt').open('w'),stderr=subprocess.STDOUT)
(o/'commands-independent.json').write_text(json.dumps({'commands':commands,'reused_native_object_sha256':sha(p/'build/native.o'),'foundation_source_four_hashes_identical_to_independently_compiled_v5':True},indent=2)+'\n');print((o/'probes.txt').read_text())
