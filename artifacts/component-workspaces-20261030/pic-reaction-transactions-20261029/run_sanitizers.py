from pathlib import Path
import fcntl,json,os,subprocess
PREFIX=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029')
REPORT=Path('<workspace>/outputs/pic-profile/reaction-transactions-20261029')
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(REPORT.parent/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
command=json.loads((REPORT/'controls.json').read_text())['build']
command=[('-O1'if argument=='-O2'else argument)for argument in command]
command[-1]=str(PREFIX/'build/test_reaction_transactions_sanitized')
index=command.index('-lcrypto');command[index:index]=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1:abort_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
run=[command[0],command[-1],str(REPORT/'sanitized-tests')]
(REPORT/'sanitizer-command.json').write_text(json.dumps({'build':command,'run':run,'ASAN_OPTIONS':env['ASAN_OPTIONS'],'UBSAN_OPTIONS':env['UBSAN_OPTIONS']},indent=2)+'\n')
with(REPORT/'sanitizer-build.log').open('w')as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
with(REPORT/'sanitizer-run.log').open('w')as log:subprocess.run(run,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
normal=json.loads((REPORT/'tests/test-results.json').read_text());sanitized=json.loads((REPORT/'sanitized-tests/test-results.json').read_text())
assert normal==sanitized
(REPORT/'sanitizer-result.json').write_text(json.dumps({'status':'passed','address_sanitizer':True,'undefined_behavior_sanitizer':True,'leak_detection':True,'checks':normal['check_count'],'normal_and_sanitized_results_identical':True},indent=2)+'\n')
print('Address/undefined/leak sanitizer checks passed:',normal['check_count'])
