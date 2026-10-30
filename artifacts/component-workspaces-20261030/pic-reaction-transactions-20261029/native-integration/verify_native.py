from pathlib import Path
import datetime,fcntl,hashlib,json,os,subprocess
P=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-integration')
R=Path('<workspace>/outputs/pic-profile/reaction-transactions-20261029/native-integration')
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(R.parent.parent/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
C=json.loads((R/'commands.json').read_text());normal=json.loads((R/'tests/native-tests.json').read_text());assert normal['check_count']==114
frozen=C['run'][0];build=P/'build';san_objects=[str(build/'native.o')]
commands=[]
for command in C['build'][1:4]:
 command=[('-O1'if x=='-O2'else x)for x in command];command[-1]=str(build/('sanitized-'+Path(command[-1]).name));command[command.index('-c'):command.index('-c')]=['-fsanitize=address,undefined','-fno-omit-frame-pointer'];commands.append(command)
 with(R/('sanitizer-'+Path(command[-1]).stem+'.log')).open('w')as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
 san_objects.append(command[-1])
link=C['build'][-1];start=link.index(str(build/'native.o'));stop=start+4
link=link[:start]+san_objects+link[stop:];link[-1]=str(build/'test_native_reactions_sanitized');link[2:2]=['-fsanitize=address,undefined'];commands.append(link)
with(R/'sanitizer-link.log').open('w')as log:subprocess.run(link,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
out=R/'sanitized-tests';out.mkdir(exist_ok=True);env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
with(R/'sanitizer-tests.log').open('w')as log:subprocess.run([frozen,link[-1],str(out)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,timeout=60)
san=json.loads((out/'native-tests.json').read_text());assert san==normal
# Compile/run the unchanged committed native foundation tests against the copied
# implementation. Only the new friend declaration differs in its existing API.
compile=[('-std=c++20' if x=='-std=c++17' else x) for x in C['build'][3]];compile[compile.index(str(P/'tests/test_native_reactions.cpp'))]=str(P/'source/tests/test_mfem_gap_pic.cpp');compile[-1]=str(build/'foundation.o');commands.append(compile)
with(R/'foundation-build.log').open('w')as log:subprocess.run(compile,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
foundation=list(C['build'][-1]);i=foundation.index(str(build/'native.o'));j=i+4;foundation=foundation[:i]+[str(build/'native.o'),str(build/'foundation.o')]+foundation[j:];foundation[-1]=str(build/'test_native_foundation');commands.append(foundation)
with(R/'foundation-link.log').open('w')as log:subprocess.run(foundation,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
with(R/'foundation-tests.log').open('w')as log:subprocess.run([frozen,foundation[-1]],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
assert 'foundation tests passed' in(R/'foundation-tests.log').read_text()
# ABI probes use actual baseline/current/sealed headers, not a copied layout.
headers={'hall_baseline':'/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source/include/raikou/mfem_gap_pic.h','native_candidate':str(P/'source/include/raikou/mfem_gap_pic.h'),'sealed_release':'/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include/raikou/mfem_gap_pic.h'}
abi={}
for label,header in headers.items():
 source=build/(label+'-abi.cpp');source.write_text('#include "'+header+'"\n#include <iostream>\nint main(){using namespace raikou::mfem_backend;std::cout<<sizeof(GapPic)<<" "<<sizeof(GapParticle)<<" "<<sizeof(GapPicOptions)<<" "<<sizeof(GapPicDiagnostics)<<"\\n";}\n')
 probe=list(C['build'][3]);probe[probe.index(str(P/'tests/test_native_reactions.cpp'))]=str(source);probe[-1]=str(build/(label+'-abi.o'));commands.append(probe)
 with(R/(label+'-abi-build.log')).open('w')as log:subprocess.run(probe,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
 executable=build/(label+'-abi');cmd=list(C['build'][-1]);i=cmd.index(str(build/'native.o'));cmd=cmd[:i]+[probe[-1]]+cmd[i+4:];cmd[-1]=str(executable);commands.append(cmd)
 with(R/(label+'-abi-link.log')).open('w')as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
 abi[label]=list(map(int,subprocess.check_output([frozen,str(executable)],text=True).split()))
assert abi['hall_baseline']==abi['native_candidate']and abi['sealed_release'][1:]==abi['native_candidate'][1:]==[80,168,200]
(R/'abi-comparison.json').write_text(json.dumps({'sizes_order':['GapPic','Particle','Options','Diagnostics'],'sizes':abi,'candidate_matches_committed_Hall_class_layout':True,'sealed_raw_checkpoint_structs_unchanged':True,'warning':'Hall observer baseline can differ in GapPic C++ class size from sealed runtime; rebuild all consumers together. Never mix C++ binaries.'},indent=2)+'\n')
(R/'verification-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
(R/'verification-result.json').write_text(json.dumps({'status':'passed','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'native_checks':114,'sanitizers':'address/undefined/leak clean on adapter/core/tests; frozen native/MFEM libraries uninstrumented','normal_sanitized_results_identical':True,'unchanged_native_foundation_suite':'passed','ABI':abi,'GPU':False},indent=2)+'\n')
print('Native verification passed:114 checks, sanitizers, unchanged foundation and actual ABI probes')
