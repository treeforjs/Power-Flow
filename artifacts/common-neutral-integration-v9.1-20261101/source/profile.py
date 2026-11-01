from pathlib import Path
import hashlib,json,runpy,subprocess,sys
P=Path(__file__).resolve().parent;R=P.parent;F=R/'arch-base-20261028/run-frozen'
M=runpy.run_path(str(R/'pic-electrode-heat-source-v7.4-mass-cache-20261101/run_source_v7.py'),run_name='metadata')
TOOLS=R/'cpu-isa-study-20261031/tools'
BASE=R/'pic-water-common-test-v9.1-profile-20261101/build/c1/attempt-001/normal/test_stage'
QUAL=R.parent/'.workspace/outputs/common-test-v9.1-profile-20261101/c1/attempt-001/qualification.json'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert len(sys.argv)==2 and '/' not in sys.argv[1]
out=R.parent/'.workspace/outputs/pic-common-profile-20261101'/sys.argv[1];out.mkdir(parents=True,exist_ok=False)
record={'passed':False,'scope':'diagnostic sampling of manufactured coupled neutral fixture; not physical full-run throughput or comparative speedup','admission':M['admission']('run')}
try:
 q=json.loads(QUAL.read_text());assert q['all_passed'];assert q['variants']['normal']['binary_sha256'][str(BASE)]==sha(BASE)
 files=[Path(__file__),BASE,QUAL,TOOLS/'pfsample',TOOLS/'pfsample.c',TOOLS/'pfreport.py',F]
 record['inputs']={str(p):sha(p) for p in files}
 runner=out/'repeat.py';runner.write_text('import subprocess,sys\nfrom pathlib import Path\nexe,trace,out=sys.argv[1:]\nfor i in range(5):\n subprocess.run([exe,trace,str(Path(out)/str(i))],check=True)\n')
 command=[str(TOOLS/'pfsample'),str(out/'samples.txt'),'200','--','/usr/bin/python','-I',str(runner),str(BASE),str(M['TRACE']),str(out/'runs')]
 record['command']=command
 with (out/'sample.log').open('w')as log:r=subprocess.run([str(F),*command],env=M['ENV'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
 if r.returncode:raise RuntimeError((out/'sample.log').read_text()[-5000:])
 for i in range(5):assert json.loads((out/'runs'/str(i)/'common.json').read_text())['passed']
 with (out/'report.txt').open('w')as log:r=subprocess.run([str(F),'/usr/bin/python','-I',str(TOOLS/'pfreport.py'),str(out/'samples.txt'),'40','--json',str(out/'report.json')],env=M['ENV'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
 if r.returncode:raise RuntimeError((out/'report.txt').read_text()[-5000:])
 assert all(sha(p)==v for p,v in record['inputs'].items());record['passed']=True
except Exception as e:record['error']=str(e);print(str(e),flush=True)
(out/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
raise SystemExit(0 if record['passed'] else 1)
