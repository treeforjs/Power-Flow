from pathlib import Path
import os,sys,runpy,subprocess,json,hashlib
P=Path(__file__).resolve().parent;R=P.parent
S=R/'pic-electrode-heat-source-v7.4-mass-cache-20261101';H=R/'hall-corrective-v3.3-b3-20261101'
N=R/'pic-reaction-transactions-20261029/native-full-time-v8.7-reservoir-stock-20261101'
M=runpy.run_path(str(S/'run_source_v7.py'),run_name='metadata');adm=M['admission']('run')
assert len(sys.argv)==2 and '/' not in sys.argv[1]
out=R.parent/'.workspace/outputs/reaction-local-charge-audit-v1-20261101'/sys.argv[1];out.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
inputs={str(p):sha(p) for p in P.iterdir() if p.is_file()};rec={'passed':False,'admission':adm,'inputs':inputs,'runs':[]}
def run(cmd,name):
 with (out/(name+'.log')).open('w')as f:p=subprocess.run([str(R/'arch-base-20261028/run-frozen'),*map(str,cmd)],env=M['ENV'],stdout=f,stderr=subprocess.STDOUT,timeout=180)
 rec['runs'].append({'command':list(map(str,cmd)),'returncode':p.returncode,'log_sha256':sha(out/(name+'.log'))})
 if p.returncode:raise RuntimeError((out/(name+'.log')).read_text()[-5000:])
try:
 flags=json.loads((R/'pic-water-common-genesis-v9-20261101/stage-snapshot/NATIVE-FLAGS.json').read_text())['flags'];flags=[flags[0],'-I'+str(H/'source'),*flags[1:],*M['owner_flags']([])[1:]]
 flags=[x.replace(str(M['N']/'source/include'),str(N/'source/include')).replace(str(M['N']/'inputs'),str(N/'inputs')) for x in flags]
 archive=N/'build/libpic_native_v8.a';assert sha(archive)=='1778aaecc0d5b9ef2c7b296a63347de44d6df0fccac28fc0f74febdc3736cb1e'
 run(flags+['-MMD','-MF',out/'audit.d','-c',P/'audit.cpp','-o',out/'audit.o'],'compile')
 run(flags+['-MMD','-MF',out/'motion.d','-c',H/'source/motion_provider.cpp','-o',out/'motion.o'],'compile-motion')
 deps={p:sha(p) for p in (out/'audit.d').read_text().replace('\\\n',' ').split(':',1)[1].split()};deps.update({p:sha(p) for p in (out/'motion.d').read_text().replace('\\\n',' ').split(':',1)[1].split()});rec['compiler_headers']=deps
 run(['/usr/bin/mpic++','-Wl,--gc-sections',out/'audit.o',out/'motion.o',N/'build/kinetic-boundary.o',archive,*M['LIBS'],'-lcrypto','-ldl','-o',out/'audit'],'link')
 run([out/'audit',out/'result.json'],'audit')
 rec['result']=json.loads((out/'result.json').read_text());assert rec['result']['passed']
 assert all(sha(p)==h for p,h in inputs.items()) and all(sha(p)==h for p,h in deps.items())
 rec['binary_sha256']=sha(out/'audit');rec['passed']=True
except Exception as e:rec['error']=str(e);print(str(e))
(out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
raise SystemExit(0 if rec['passed'] else 1)
