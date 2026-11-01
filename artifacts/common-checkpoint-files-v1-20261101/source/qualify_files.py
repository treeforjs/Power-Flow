from pathlib import Path
import hashlib,json,runpy,sys,subprocess
P=Path(__file__).resolve().parent;R=P.parent;F=R/'arch-base-20261028/run-frozen'
M=runpy.run_path(str(R/'pic-electrode-heat-source-v7.4-mass-cache-20261101/run_source_v7.py'),run_name='metadata')
assert len(sys.argv)==2 and '/' not in sys.argv[1]
out=R.parent/'.workspace/outputs/common-checkpoint-files-v1-20261101'/sys.argv[1];out.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
own={str(p):sha(p) for p in P.rglob('*') if p.is_file()}
receipt={'passed':False,'source':own,'admission':M['admission']('run'),'variants':{}}
try:
 for name,san in [('normal',[]),('sanitizer',['-fsanitize=address,undefined','-fno-sanitize-recover=all','-fno-omit-frame-pointer'])]:
  dest=out/name;dest.mkdir();exe=dest/'test';log=dest/'compile.log'
  cmd=['/usr/bin/c++','-std=c++20','-O1','-g0','-Wall','-Wextra','-Werror',*san,'-I'+str(P/'include'),'-I'+str(R/'pic-reaction-transactions-20261029/native-full-time-v8.9-motion-inventory-20261101/inputs'),'-I'+str(R/'pic-reaction-transactions-20261029/source/vendor'),str(P/'test_files.cpp'),'-lcrypto','-o',str(exe)]
  with log.open('w')as f:r=subprocess.run([str(F),*cmd],env=M['ENV'],stdout=f,stderr=subprocess.STDOUT,timeout=120)
  if r.returncode:raise RuntimeError(log.read_text()[-6000:])
  with (dest/'test.log').open('w')as f:r=subprocess.run([str(F),str(exe),str(dest/'test-output')],env=M['ENV'],stdout=f,stderr=subprocess.STDOUT,timeout=120)
  if r.returncode:raise RuntimeError((dest/'test.log').read_text()[-6000:])
  result=json.loads((dest/'test-output/result.json').read_text());assert result['passed'];receipt['variants'][name]={'checks':len(result['cases']),'result_sha256':sha(dest/'test-output/result.json'),'binary_sha256':sha(exe),'command':cmd}
 assert all(sha(p)==v for p,v in own.items());receipt['passed']=True
except Exception as e:receipt['error']=str(e);print(str(e),flush=True)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
raise SystemExit(0 if receipt['passed'] else 1)
