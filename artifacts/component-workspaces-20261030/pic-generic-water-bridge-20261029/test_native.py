from pathlib import Path
import os,json,hashlib,subprocess,shutil,resource
p=Path(__file__).parent;o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-generic-water-bridge-20261029');assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(180,180))
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
manifest=json.loads((o/'kernel-binding.json').read_text())
for n,v in manifest['sources'].items():assert sha(p/n)==v,n
native=p;base=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
# Verified immutable cached objects; a clean prefix rebuilds the same source.
reference=json.loads((p/'inputs/native-object-reference.json').read_text())
if (p/'build/rebuilt-native-object-reference.json').is_file():
 reference=json.loads((p/'build/rebuilt-native-object-reference.json').read_text())
 assert reference['source_archive_sha256']=='d14f1e229bbd8a84cf75b782f49562217972d4cd55b086fd8f932cc4f0e37db4'
for name,key in [('native.o','native_object_sha256'),('arena.o','arena_object_sha256')]:
 if not (p/'build'/name).is_file() or sha(p/'build'/name)!=reference[key]:
  raise RuntimeError('Pinned foundation cache absent; run rebuild_foundation.py under this same lease, then test with its recorded hashes.')
for filename,want in json.loads((p/'inputs/frozen-library-hashes.json').read_text()).items():assert sha(Path(filename))==want,filename
for name in ['json.hpp','electron_boltzmann.h']:assert sha(p/'products/include/raikou'/name)==sha(native/'native/inputs/raikou'/name)
cc=['mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(p/'source'),'-I'+str(p/'inputs'),'-I'+str(native/'native/source/include'),'-isystem',str(p/'products/include'),'-isystem',str(native/'native/inputs'),'-isystem',str(native/'native/vendor')]
for name in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:cc+=['-isystem',str(base/name)]
command=cc+['-c',str(p/'tests/test_native_admission.cpp'),'-o',str(p/'build/native_test.o')];commands=[command]
with (o/'build-native.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=90)
libs=[str(base/n) for n in ['install-mpi/lib/libmfem.a','tpl/hypre/lib/libHYPRE.a','tpl/metis/lib/libmetis.a','tpl/GKlib/lib/libGKlib.a']]+['-lcrypto']
command=['mpic++','-Wl,--gc-sections']+[str(p/'build'/n) for n in ['native.o','arena.o','native_test.o']]+libs+['-o',str(p/'build/native_test')];commands.append(command)
with (o/'link-native.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=90)
command=['/usr/bin/time','-v','-o',str(o/'native-resources.txt'),str(p/'build/native_test'),str(p/'products/data/collisions/h2o-gas-20261029-elastic.json'),str(o/'run-binding.json'),str(o)];commands.append(command)
(o/'commands-native.json').write_text(json.dumps({'commands':commands,'native_object_sha256':sha(p/'build/native.o'),'arena_object_sha256':sha(p/'build/arena.o'),'native_source_identity':'5fa3a59f863e98f21e20af698b9c16417f3749344b555f18e3d2cba1e4f5a705'},indent=2)+'\n')
with (o/'native-run.log').open('w') as f:subprocess.run(command,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=90)
print((o/'native-run.log').read_text())
