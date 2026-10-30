from pathlib import Path
import os,json,hashlib,subprocess,tarfile
p=Path(__file__).parent
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(p/'inputs/native-source.tar.gz')=='d14f1e229bbd8a84cf75b782f49562217972d4cd55b086fd8f932cc4f0e37db4'
with tarfile.open(p/'inputs/native-source.tar.gz') as t:
 for x in t.getmembers():
  if x.isfile():assert sha(p/'native'/x.name)==hashlib.sha256(t.extractfile(x).read()).hexdigest(),x.name
m=Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
cc=['mpic++','-std=c++20','-O2','-g','-ffunction-sections','-fdata-sections','-DPIC_REACTION_NATIVE_SOURCE_SHA256="5fa3a59f863e98f21e20af698b9c16417f3749344b555f18e3d2cba1e4f5a705"','-I'+str(p/'native/source/include'),'-isystem',str(p/'native/inputs'),'-isystem',str(p/'native/vendor')]
for n in ['install-mpi/include','tpl/hypre/include','tpl/metis/include']:cc+=['-isystem',str(m/n)]
commands=[]
for source,name in [('mfem_gap_pic.cpp','native.o'),('mfem_pic_arena.cpp','arena.o')]:
 command=cc+['-c',str(p/'native/source/src'/source),'-o',str(p/'build'/name)];commands.append(command);subprocess.run(command,check=True,timeout=120)
# This local cache is now bound to a source-verified rebuild, not an old hash.
(p/'build/rebuilt-native-object-reference.json').write_text(json.dumps({'native_object_sha256':sha(p/'build/native.o'),'arena_object_sha256':sha(p/'build/arena.o'),'rebuild_commands':commands,'source_archive_sha256':sha(p/'inputs/native-source.tar.gz')},indent=2)+'\n')
