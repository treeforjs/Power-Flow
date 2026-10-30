from pathlib import Path
import os,json,hashlib,subprocess,shutil
p=Path(__file__).parent;o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-generic-products-native-audit-20261029');assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();data=p/'products/data/collisions/h2o-gas-20261029-outer3.json';assert sha(data)=='c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551'
headers=['collision_data.h','h2o_conditional_parent.h','h2o_dea_kernel.h','h2o_two_body_attachment.h','h2o_neutral_dissociation.h','h2o_common_mass_basis.h','h2o_outer_valence.h'];sources={n:sha(p/'products/include/raikou'/n) for n in headers};kernel=hashlib.sha256(json.dumps(sources,sort_keys=True,separators=(',',':')).encode()).hexdigest();manifest={'sources':sources,'canonical_kernel_manifest_sha256':kernel,'dataset_sha256':sha(data),'native_source_identity':'5fa3a59f863e98f21e20af698b9c16417f3749344b555f18e3d2cba1e4f5a705','catalog_mass_inputs_kg':{'H2O':18.01056468403*1.66053906660e-27,'H':1.00782503223*1.66053906660e-27}}
(o/'kernel-binding.json').write_text(json.dumps(manifest,indent=2)+'\n');binding={'dataset_sha256':sha(data),'kernel_manifest_sha256':kernel,'kernel_binding_file_sha256':sha(o/'kernel-binding.json'),'mesh_sha256':hashlib.sha256(b'native-45x45-9mm-height25mm').hexdigest(),'clock_sha256':hashlib.sha256(b'forced-v5-native-t0-proof').hexdigest()};(o/'run-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
native=Path('/home/treefor/projects/runtimes/pic-native-subflow-v6_2-review-20261029');(p/'build').mkdir(exist_ok=True)
for n in ('native.o','arena.o'):shutil.copy2(native/'build'/n,p/'build'/n)
for n in ('json.hpp','electron_boltzmann.h'):assert sha(p/'products/include/raikou'/n)==sha(native/'native/inputs/raikou'/n)
base='/home/treefor/projects/runtimes/mfem-4.10-pic-20261029';includes=['-I'+str(native/'native/source/include'),'-isystem',str(p/'products/include'),'-isystem',str(native/'native/inputs'),'-I'+str(native/'native/vendor')]
for n in ('install-mpi/include','tpl/hypre/include','tpl/metis/include'):includes+=['-isystem',base+'/'+n]
commands=[];command=['/usr/bin/mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-Wall','-Wextra','-Wpedantic','-Werror']+includes+['-c',str(p/'probe.cpp'),'-o',str(p/'build/probes.o')];commands.append(command)
with (o/'build.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
libs=[base+'/install-mpi/lib/libmfem.a',base+'/tpl/hypre/lib/libHYPRE.a',base+'/tpl/metis/lib/libmetis.a',base+'/tpl/GKlib/lib/libGKlib.a','-lcrypto'];command=['/usr/bin/mpic++','-Wl,--gc-sections']+[str(p/'build'/n) for n in ('native.o','arena.o','probes.o')]+libs+['-o',str(p/'build/probes')];commands.append(command);subprocess.run(command,check=True)
command=['/usr/bin/time','-v','-o',str(o/'resources.txt'),str(p/'build/probes'),str(data),str(o/'run-binding.json'),str(o)];commands.append(command)
with (o/'run.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
(o/'commands.json').write_text(json.dumps({'commands':commands,'native_object_sha256':sha(p/'build/native.o'),'arena_object_sha256':sha(p/'build/arena.o'),'reused_from_qualified_v6_2_source_identity':manifest['native_source_identity']},indent=2)+'\n');print((o/'run.log').read_text())
