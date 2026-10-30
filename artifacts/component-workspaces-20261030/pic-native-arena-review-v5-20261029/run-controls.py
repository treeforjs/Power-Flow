from pathlib import Path
import os,subprocess,json
p=Path(__file__).parent;s=p/'snapshot';o=Path('<workspace>/outputs/pic-native-arena-review-v5-20261029')
assert sorted(os.sched_getaffinity(0))==[7]
text=(s/'tests/test_arena.cpp').read_text();(p/'fixture.hpp').write_text(text[:text.index('int main(')])
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029';includes=['-I'+str(s/'source/include'),'-isystem',str(s/'inputs'),'-I'+str(s/'vendor')]
for n in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+n]
command=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-c',str(p/'restore-controls.cpp'),'-o',str(p/'build/restore-controls.o')]
subprocess.run(command,check=True,stdout=(o/'build-controls.log').open('w'),stderr=subprocess.STDOUT)
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto']
subprocess.run(['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/f'{n}.o') for n in ('native','arena','restore-controls')]+libs+['-o',str(p/'build/restore-controls')],check=True)
subprocess.run([str(p/'build/restore-controls'),str(o/'owner-tests'),str(o/'restore-controls.json')],check=True)
print((o/'restore-controls.json').read_text())
