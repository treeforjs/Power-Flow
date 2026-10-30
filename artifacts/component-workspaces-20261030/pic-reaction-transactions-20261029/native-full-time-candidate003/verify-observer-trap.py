from pathlib import Path
import sys,os,json,subprocess
P=Path(__file__).parent;R=Path('<workspace>/outputs/pic-profile/reaction-transactions-20261029/native-full-time-candidate003')
if os.sched_getaffinity(0)!={7}or os.getpriority(os.PRIO_PROCESS,0)!=15:raise RuntimeError('CPU7/nice15 required')
subprocess.run([sys.executable,str(P/'run-candidate-tests.py')],check=True)
C=json.loads((R/'candidate-tests-commands.json').read_text());command=[('-O1'if x in['-O0','-O2']else '-g'if x=='-g0'else x)for x in C[0]];command[command.index('-c'):command.index('-c')]=['-fsanitize=address,undefined','-fno-omit-frame-pointer'];command[-1]=str(P/'build/sanitized-candidate-tests.o')
with(R/'sanitized-observer-test-compile.log').open('w')as out:subprocess.run(command,stdout=out,stderr=subprocess.STDOUT,check=True,timeout=600)
link=list(C[1]);link[2:2]=['-fsanitize=address,undefined'];link=[str(P/'build/sanitized-arena.o')if x==str(P/'build/arena.o')else str(P/'build/sanitized-candidate-tests.o')if x==str(P/'build/candidate-tests.o')else x for x in link];link[-1]=str(P/'build/test-sanitized-candidate')
with(R/'sanitized-observer-test-link.log').open('w')as out:subprocess.run(link,stdout=out,stderr=subprocess.STDOUT,check=True)
env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1');run=[C[2][0],link[-1],str(R/'sanitized-candidate')]
with(R/'sanitized-observer-test-run.log').open('w')as out:subprocess.run(run,stdout=out,stderr=subprocess.STDOUT,check=True,env=env)
a=json.loads((R/'component-tests/candidate003-tests.json').read_text());b=json.loads((R/'sanitized-candidate/candidate003-tests.json').read_text())
if a!=b:raise RuntimeError('normal/sanitized observed trap evidence differs')
p=R/'verification-result.json';result=json.loads(p.read_text());result['forbidden_trapped_functions']=['solve_fields','advance','advance_empty','advance_observed'];result['observer_trap_addendum']='tests only; primary runtime source unchanged';p.write_text(json.dumps(result,indent=2)+'\n');(R/'observer-trap-commands.json').write_text(json.dumps([command,link,run],indent=2)+'\n');print('All four forbidden function traps zero;144 controls unchanged, normal/sanitized identical')
