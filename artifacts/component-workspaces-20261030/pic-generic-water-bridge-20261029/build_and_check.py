from pathlib import Path
import os,json,hashlib,subprocess,time,resource,tarfile
p=Path(__file__).parent;o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-generic-water-bridge-20261029')
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(180,180))
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(p/'products-source.tar.gz')=='f264547006a6c88f8dbb47538906efd66cc6561e35ac9155491b514d2a9d7a0b'
if (p/'inputs/frozen-source-manifest.json').is_file():
 for filename,want in json.loads((p/'inputs/frozen-source-manifest.json').read_text()).items():assert sha(p/filename)==want,filename
with tarfile.open(p/'products-source.tar.gz') as tar:
 for entry in tar.getmembers():
  if entry.isfile():assert sha(p/'products'/entry.name)==hashlib.sha256(tar.extractfile(entry).read()).hexdigest(),entry.name

data=p/'products/data/collisions/h2o-gas-20261029-elastic.json';assert sha(data)=='6ea91f535096beafc19b3d69029597addb10909a60c1fae5064417f0df38b140'
headers=['collision_data.h','h2o_conditional_parent.h','h2o_dea_kernel.h','h2o_two_body_attachment.h','h2o_neutral_dissociation.h','h2o_common_mass_basis.h','h2o_outer_valence.h','h2o_elastic_angular.h']
sources={str(f.relative_to(p)):sha(f) for f in [p/'products/include/raikou'/n for n in headers]+[p/'inputs/collision_scheduler.hpp',p/'inputs/water_parent_bridge.hpp',p/'source/generic_water_bridge.hpp',p/'source/generic_water_native_adapter.hpp']}
kernel=hashlib.sha256(json.dumps(sources,sort_keys=True,separators=(',',':')).encode()).hexdigest()
manifest={'scope':'supported conditional comparison subset; not complete water chemistry','sources':sources,'dataset_sha256':sha(data),'canonical_kernel_manifest_sha256':kernel,'one_keV_policy':'outer3 at exact1000eV','collision_site_policy':'preserve reactant sites; capture requires exact3Dcoincidence'}
(o/'kernel-binding.json').write_text(json.dumps(manifest,indent=2)+'\n');binding={'dataset_sha256':sha(data),'kernel_manifest_sha256':kernel,'kernel_binding_file_sha256':sha(o/'kernel-binding.json')};(o/'run-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
command=['g++','-std=c++20','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(p/'source'),'-I'+str(p/'inputs'),'-isystem',str(p/'products/include'),str(p/'tests/test_bridge.cpp'),'-o',str(p/'build/test_bridge')]
commands=[command]
with (o/'build-bridge.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=90)
command=['/usr/bin/time','-v','-o',str(o/'bridge-resources.txt'),str(p/'build/test_bridge'),str(data),str(o/'run-binding.json'),str(o)];commands.append(command)
(o/'commands-bridge.json').write_text(json.dumps(commands,indent=2)+'\n')
with (o/'bridge-run.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=90)
print((o/'bridge-run.log').read_text())
