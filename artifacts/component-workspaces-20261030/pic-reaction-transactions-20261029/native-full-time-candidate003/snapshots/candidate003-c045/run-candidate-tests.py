from pathlib import Path
import json,subprocess,os
P=Path(__file__).parent;R=Path('<workspace>/outputs/pic-profile/reaction-transactions-20261029/native-full-time-candidate003')
if os.sched_getaffinity(0)!={7}or os.getpriority(os.PRIO_PROCESS,0)!=15:raise RuntimeError('CPU7/nice15 required')
C=json.loads((R/'legacy-commands.json').read_text());commands=[]
cmd=[x.replace('tests/test_arena.cpp','tests/test_candidate003.cpp').replace('build/tests.o','build/candidate-tests.o')for x in C[1]];commands.append(cmd)
with(R/'candidate-tests-compile.log').open('w')as out:p=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
if p.returncode:print((R/'candidate-tests-compile.log').read_text()[-10000:]);raise SystemExit(p.returncode)
link=[x.replace('build/tests.o','build/candidate-tests.o').replace('build/test_arena','build/test_candidate003')for x in C[2]]
wrappers=['_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv','_ZN6raikou12mfem_backend6GapPic7advanceEd','_ZN6raikou12mfem_backend6GapPic13advance_emptyEd','_ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE'];link[2:2]=['-Wl,--wrap='+w for w in wrappers];commands.append(link)
with(R/'candidate-tests-link.log').open('w')as out:subprocess.run(link,stdout=out,stderr=subprocess.STDOUT,check=True)
S=R/'component-tests';S.mkdir(exist_ok=True);run=[C[3][0],str(P/'build/test_candidate003'),str(S)];commands.append(run)
with(R/'candidate-tests-run.log').open('w')as out:p=subprocess.run(run,stdout=out,stderr=subprocess.STDOUT)
(R/'candidate-tests-commands.json').write_text(json.dumps(commands,indent=2)+'\n');print('component exit',p.returncode);print((R/'candidate-tests-run.log').read_text()[-7000:]);raise SystemExit(p.returncode)
