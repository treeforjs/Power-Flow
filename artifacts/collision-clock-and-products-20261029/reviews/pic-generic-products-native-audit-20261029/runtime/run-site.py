from pathlib import Path
import os,subprocess,json
p=Path(__file__).parent;o=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-generic-products-native-audit-20261029');assert sorted(os.sched_getaffinity(0))==[7]
n=Path('/home/treefor/projects/runtimes/pic-native-subflow-v6_2-review-20261029/native');base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029';includes=['-I'+str(n/'source/include'),'-isystem',str(p/'products/include'),'-isystem',str(n/'inputs'),'-I'+str(n/'vendor')]
for s in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+s]
commands=[];command=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-c',str(p/'site-capability.cpp'),'-o',str(p/'build/site.o')];commands.append(command)
subprocess.run(command,check=True,stdout=(o/'build-site.log').open('w'),stderr=subprocess.STDOUT)
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto'];command=['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/s) for s in ('native.o','arena.o','site.o')]+libs+['-o',str(p/'build/site')];commands.append(command);subprocess.run(command,check=True)
command=[str(p/'build/site'),str(p/'products/data/collisions/h2o-gas-20261029-outer3.json'),str(o/'run-binding.json'),str(o)];commands.append(command);subprocess.run(command,check=True,stdout=(o/'site-run.log').open('w'),stderr=subprocess.STDOUT);(o/'commands-site.json').write_text(json.dumps(commands,indent=2)+'\n')
