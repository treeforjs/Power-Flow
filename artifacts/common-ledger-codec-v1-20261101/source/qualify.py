#!/usr/bin/python -I
from pathlib import Path
import hashlib,json,os,runpy,subprocess,sys
P=Path(__file__).resolve().parent;R=P.parent;C=R/'pic-water-common-v9.2-conservation-20261101';B=R/'pic-water-common-commit-v8.0.4-b3-20261101';N=R/'pic-reaction-transactions-20261029/native-full-time-v8.9-motion-inventory-20261101';S=R/'pic-electrode-heat-source-v7.4-mass-cache-20261101';H=R/'hall-corrective-v3.3-b3-20261101';Q=R.parent/'.workspace/outputs/common-v9.2-conservation-20261101/c1/attempt-001/qualification.json'
SOURCE=runpy.run_path(str(S/'run_source_v7.py'),run_name='metadata');F=R/'arch-base-20261028/run-frozen'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def own():return {str(p.relative_to(P)):sha(p) for p in P.rglob('*') if p.is_file() and 'build' not in p.parts and 'evidence' not in p.parts}
def qualify(label):
 admission=SOURCE['admission']('run');out=P/'evidence'/label;out.mkdir(parents=True,exist_ok=False);build=P/'build'/label;build.mkdir(parents=True,exist_ok=False)
 base_receipt=json.loads(Q.read_text());assert base_receipt['all_passed'];pins={str(Q):sha(Q),**base_receipt['own_sha256']}
 for variant in base_receipt['variants'].values():pins.update(variant['dependency_sha256']);pins.update(variant['binary_sha256'])
 assert all(sha(k)==v for k,v in pins.items()),'qualified common dependency changed'
 before=own();result={'admission':admission,'source':before,'pins':pins,'passed':False,'variants':{},'commands':[]}
 def execute(cmd,name):
  proc=subprocess.run([str(F),*map(str,cmd)],cwd=P,env=SOURCE['ENV'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=480)
  (out/(name+'.stdout')).write_bytes(proc.stdout);(out/(name+'.stderr')).write_bytes(proc.stderr);result['commands'].append({'argv':list(map(str,cmd)),'returncode':proc.returncode});print(name,proc.returncode,flush=True)
  if proc.returncode:raise RuntimeError(proc.stderr.decode()[-5000:]+proc.stdout.decode()[-5000:])
  return proc.stdout
 try:
  include=[P/'include',C/'include',C/'tests',R/'native-gas-common-stage-v1.3.2-retained-proofs-20261101/include',R/'pic-water-common-genesis-v9-20261101/include',R/'reaction-final-charge-guard-v1-20261101',R/'common-material-stocks-v1-20261101/include',R/'common-conservation-v1-20261101/include',R/'conductor-correction-guard-v1-20261101/include',B/'source',R/'native-gas-adapter-v1.2-immutable-export-20261101/include',R/'gas-target-owner-v1.7.1-census-transport-20261101/include',N/'forms/include/pf_v7',H/'source',B/'vendor',R/'pic-arrival-owner-v1.2-publication-20261101/include',R/'pic-source-stock-reference-v2-surface-birth-20261101/include',N/'tests',R.parent/'Power-Flow/native/coupled_pic/include']
  flags=json.loads((B/'NATIVE-FLAGS.json').read_text())['flags'];flags=[flags[0],*['-I'+str(x) for x in include],*flags[1:],*SOURCE['owner_flags']([])[1:]]
  flags=[x.replace('native-full-time-v8.7-reservoir-stock-20261101','native-full-time-v8.9-motion-inventory-20261101').replace('d795695aca20ddef0efadedaebea8baa6c994e819514c0258bbac290fb66103d','8e8d9553ccac604e7e69ecc3354f124e71643a17e8fbc85b9883dd4f798897db').replace(str(SOURCE['N']/'source/include'),str(N/'source/include')).replace(str(SOURCE['N']/'inputs'),str(N/'inputs')) for x in flags]
  for variant in ('normal','ubsan'):
   dest=build/variant;dest.mkdir();base=C/'build/c1/attempt-001'/variant
   execute([*flags,*SOURCE['VARIANTS'][variant],'-MD','-MF',dest/'test.d','-c',P/'tests/test_codec.cpp','-o',dest/'test.o'],variant+'-compile')
   objects=[p for p in sorted(base.glob('*.o')) if p.name!='test_stage.o'];assert len(objects)==6
   exe=dest/'codec';execute(['/usr/bin/mpic++','-Wl,--gc-sections',*SOURCE['LINK_SAN'][variant],dest/'test.o',*objects,N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*SOURCE['LIBS'],'-lcrypto','-ldl','-o',exe],variant+'-link')
   raw=execute([exe,SOURCE['TRACE']],variant+'-controls');data=json.loads(raw);assert data['passed'];assert execute([exe,SOURCE['TRACE']],variant+'-replay')==raw,'deterministic replay differs'
   result['variants'][variant]=data
   paths=(dest/'test.d').read_text().replace('\\\n',' ').split(':',1)[1].split();result.setdefault('compiled_inputs',{}).update({str(SOURCE['dependency_path'](x)):sha(SOURCE['dependency_path'](x)) for x in paths});result.setdefault('binaries',{}).update({str(x):sha(x) for x in (dest/'test.o',exe)})
  assert result['variants']['normal']==result['variants']['ubsan'],'normal/sanitizer mismatch';assert own()==before and all(sha(k)==v for k,v in pins.items()) and all(sha(k)==v for k,v in result['compiled_inputs'].items()),'input changed'
  result['passed']=True
 except Exception as e:result['error']=str(e);print(str(e),flush=True)
 (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(qualify(sys.argv[1]))
