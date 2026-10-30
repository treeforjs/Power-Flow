from pathlib import Path
import os,hashlib,json,subprocess
p=Path(__file__).parent;o=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-native-subflow-v6_2-review-20261029');assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();assert sha(p/'driver/source/native_parent_subflow.hpp')=='ad499c5549a78bae371c5b2f72ee66085c3b141126d1931d1ec654d8a9e25675'
binding=json.loads((o/'arena-build-binding.json').read_text());assert sha(p/'build/arena.o')==binding['arena_object_sha256'] and sha(p/'build/native.o')==binding['native_object_sha256']
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029';includes=['-I'+str(p/'driver/source'),'-I'+str(p/'native/source/include'),'-isystem',str(p/'driver/inputs'),'-isystem',str(p/'driver/inputs/include'),'-isystem',str(p/'native/inputs'),'-I'+str(p/'native/vendor')]
for n in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+n]
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto'];commands=[]
for source,name in ((p/'port-elapsed.cpp','port'),(p/'focused.cpp','focused')):
 command=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-c',str(source),'-o',str(p/'build'/f'{name}.o')];commands.append(command)
 with (o/f'build-{name}.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
 command=['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/f'{n}.o') for n in ('native','arena',name)]+libs+['-o',str(p/'build'/name)];commands.append(command);subprocess.run(command,check=True)
 command=['/usr/bin/time','-v','-o',str(o/f'{name}-resources.txt'),str(p/'build'/name),str(p/'driver/inputs/h2o-gas.json'),str(o/'owner-run-binding.json'),str(o)];commands.append(command)
 with (o/f'{name}-run.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
(o/'commands-focused.json').write_text(json.dumps(commands,indent=2)+'\n');print((o/'focused-run.log').read_text())
