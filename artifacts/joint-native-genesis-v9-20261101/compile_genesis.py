from pathlib import Path
import os,sys,runpy,subprocess,json,hashlib
R=Path(__file__).resolve().parent;RT=R.parent
S=RT/'pic-electrode-heat-source-v7.4-mass-cache-20261101';H=RT/'hall-corrective-v3.3-b3-20261101'
A=RT/'pic-arrival-owner-v1.2-publication-20261101';N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.7-reservoir-stock-20261101'
M=runpy.run_path(str(S/'run_source_v7.py'),run_name='metadata');admission=M['admission']('run')
assert len(sys.argv)==2 and '/' not in sys.argv[1]
out=RT.parent/'.workspace/outputs/joint-genesis-v9-20261101'/sys.argv[1];out.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def own():return {str(p):sha(p) for p in R.rglob('*') if p.is_file()}
before=own();rec={'passed':False,'scope':'C++ syntax/API check only; SourceStage snapshot is unqualified development','admission':admission,'source':before}
try:
 base=json.loads((R/'stage-snapshot/NATIVE-FLAGS.json').read_text())['flags']
 base=[base[0],'-I'+str(R/'include'),'-I'+str(R/'stage-snapshot'),'-I'+str(H/'source'),'-I'+str(A/'include'),'-I'+str(RT/'gas-target-owner-v1.7.1-census-transport-20261101/include'),'-I'+str(N/'forms/include/pf_v7'),*base[1:]]
 base += M['owner_flags']([])[1:]
 base=[x.replace(str(M['N']/'source/include'),str(N/'source/include')).replace(str(M['N']/'inputs'),str(N/'inputs')) for x in base]
 cmd=base+M['VARIANTS']['normal']+['-MMD','-MF',str(out/'genesis.d'),'-fsyntax-only',str(R/'check.cpp')]
 rec['command']=cmd
 with (out/'compile.log').open('w')as f:p=subprocess.run([str(RT/'arch-base-20261028/run-frozen'),*cmd],env=M['ENV'],stdout=f,stderr=subprocess.STDOUT,timeout=120)
 rec['returncode']=p.returncode
 if p.returncode:raise RuntimeError((out/'compile.log').read_text()[-7000:])
 assert own()==before
 rec['passed']=True
except Exception as e:rec['error']=str(e);print(str(e))
(out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
raise SystemExit(0 if rec['passed'] else 1)
