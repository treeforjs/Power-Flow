from pathlib import Path
import json,subprocess,os
P=Path(__file__).parent;R=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/reaction-transactions-20261029/native-full-time');R.mkdir(exist_ok=True);(P/'build').mkdir(exist_ok=True)
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
C=json.loads((R.parent/'native-exact-count/commands.json').read_text());cmd=[x.replace('native-exact-count','native-full-time')for x in C['build'][1]];cmd=[x for x in cmd if x!='-Werror'];cmd[cmd.index('-O2')]='-O0';cmd[cmd.index('-g')]='-g0';cmd[cmd.index('-DPIC_REACTION_NATIVE_SOURCE_SHA256="b8a9a27e6a0fcd8c54e49dd341503fed165da460180b608aa0ae024b0436fb0e"')]='-DPIC_REACTION_NATIVE_SOURCE_SHA256="draft-native-full-time"'
(R/'compile-command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with(R/'compile-intermediate.log').open('w')as out:p=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
print('compile exit',p.returncode);print((R/'compile-intermediate.log').read_text()[-24000:]);raise SystemExit(p.returncode)
