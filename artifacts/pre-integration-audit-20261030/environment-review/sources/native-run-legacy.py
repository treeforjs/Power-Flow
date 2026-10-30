from pathlib import Path
import json,subprocess,os
P=Path(__file__).parent;R=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/reaction-transactions-20261029/native-full-time');assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
C=json.loads((R.parent/'native-exact-count/commands.json').read_text());commands=[]
for index,name in [(1,'arena'),(2,'tests')]:
 cmd=[x.replace('native-exact-count','native-full-time')for x in C['build'][index]];cmd=[('-O0'if x=='-O2'else '-g0'if x=='-g'else x)for x in cmd];cmd=[x.replace('b8a9a27e6a0fcd8c54e49dd341503fed165da460180b608aa0ae024b0436fb0e','draft-native-full-time')for x in cmd];commands.append(cmd)
 with(R/(name+'-compile.log')).open('w')as out:p=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
 if p.returncode:print((R/(name+'-compile.log')).read_text()[-8000:]);raise SystemExit(p.returncode)
link=[x.replace('native-exact-count','native-full-time')for x in C['build'][-1]];link[link.index(str(P/'build/native.o'))]=str(P.parent/'native-exact-count/build/native.o');commands.append(link)
with(R/'legacy-link.log').open('w')as out:subprocess.run(link,stdout=out,stderr=subprocess.STDOUT,check=True)
run=[x.replace('native-exact-count','native-full-time')for x in C['run']];Path(run[-1]).mkdir(exist_ok=True);commands.append(run)
with(R/'legacy-tests.log').open('w')as out:p=subprocess.run(run,stdout=out,stderr=subprocess.STDOUT)
(R/'legacy-commands.json').write_text(json.dumps(commands,indent=2)+'\n');print('legacy return',p.returncode);print((R/'legacy-tests.log').read_text()[-4000:]);raise SystemExit(p.returncode)
