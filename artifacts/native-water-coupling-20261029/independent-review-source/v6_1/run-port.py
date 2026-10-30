from pathlib import Path
import os,subprocess,json
p=Path(__file__).parent;s=p/'snapshot';o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-native-subflow-v6_1-review-20261029');assert sorted(os.sched_getaffinity(0))==[7]
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029'
includes=['-I'+str(s/'source'),'-I'+str(s/'native/source/include'),'-isystem',str(s/'inputs'),'-isystem',str(s/'inputs/include'),'-isystem',str(s/'native/inputs'),'-I'+str(s/'native/vendor')]
for n in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+n]
commands=[];command=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-c',str(p/'port-elapsed.cpp'),'-o',str(p/'build/port.o')];commands.append(command)
subprocess.run(command,check=True,stdout=(o/'build-port.log').open('w'),stderr=subprocess.STDOUT)
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto'];command=['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/f'{n}.o') for n in ('native','arena','port')]+libs+['-o',str(p/'build/port')];commands.append(command);subprocess.run(command,check=True)
command=[str(p/'build/port'),str(s/'inputs/h2o-gas.json'),str(o/'owner-run-binding.json'),str(o)];commands.append(command);subprocess.run(command,check=True,stdout=(o/'port-bounds-probes.txt').open('w'),stderr=subprocess.STDOUT)
(o/'commands-port.json').write_text(json.dumps(commands,indent=2)+'\n');print((o/'port-bounds-probes.txt').read_text())
