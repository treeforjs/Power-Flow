from pathlib import Path
import json,subprocess,os
P=Path(__file__).parent;R=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/reaction-transactions-20261029/native-full-time');assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
C=json.loads((R/'legacy-commands.json').read_text());cmd=[x.replace('tests/test_arena.cpp','tests/test_full_time_storage.cpp').replace('build/tests.o','build/fulltime-tests.o')for x in C[1]]
with(R/'storage-compile.log').open('w')as out:p=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
if p.returncode:print((R/'storage-compile.log').read_text()[-8000:]);raise SystemExit(p.returncode)
link=[x.replace('build/tests.o','build/fulltime-tests.o').replace('build/test_arena','build/test_fulltime')for x in C[2]]
with(R/'storage-link.log').open('w')as out:subprocess.run(link,stdout=out,stderr=subprocess.STDOUT,check=True)
S=R/'storage-tests';S.mkdir(exist_ok=True);run=[C[3][0],str(P/'build/test_fulltime'),str(S)]
with(R/'storage-tests.log').open('w')as out:p=subprocess.run(run,stdout=out,stderr=subprocess.STDOUT)
(R/'storage-commands.json').write_text(json.dumps([cmd,link,run],indent=2)+'\n');print('storage exit',p.returncode);print((R/'storage-tests.log').read_text()[-5000:]);raise SystemExit(p.returncode)
