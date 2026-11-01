#!/usr/bin/python -I
"""Guarded component qualification; one explicit-memory SimQueue CPU, frozen userspace."""
from pathlib import Path
import hashlib,json,os,runpy,subprocess,sys,time
P=Path(__file__).resolve().parent
R=P.parent
S=R/'pic-electrode-heat-source-v7.4-mass-cache-20261101'
HALL=R/'hall-corrective-v3.3-b3-20261101'
ARRIVAL=R/'pic-arrival-owner-v1.2-publication-20261101'
N=R/'pic-reaction-transactions-20261029/native-full-time-v8.9-motion-inventory-20261101'
OUT=R.parent/'.workspace/outputs/common-test-v9.1-profile-20261101'
BASELINE=R/'hall-corrective-v3.3-20261031'
ADAPTER=R/'pic-source-stock-reference-v2-surface-birth-20261101'
F=R/'arch-base-20261028/run-frozen'
SOURCE=runpy.run_path(str(S/'run_source_v7.py'),run_name='source_metadata')
def sha(p):return SOURCE['sha'](p)
def own():return {str(p):sha(p) for p in sorted(P.rglob('*')) if p.is_file() and 'build' not in p.parts and p.suffix in ('.py','.hpp','.cpp','.inc','.json','.txt')}
def sealed_source(runtime,artifact,expected):
 manifest=R.parent/'Power-Flow/artifacts'/artifact/'SOURCE-SHA256.json';data=json.loads(manifest.read_text())
 assert data['identity_sha256']==expected and hashlib.sha256(json.dumps(data['files'],sort_keys=True,separators=(',',':')).encode()).hexdigest()==expected,'sealed source manifest differs'
 assert all(sha(runtime/k)==h for k,h in data['files'].items()),'sealed runtime source differs'
 return {'manifest':str(manifest),'manifest_sha256':sha(manifest),'identity_sha256':expected}
def pins():
 seals=[sealed_source(HALL,'hall-corrective-v3.3-b3-20261101','231b92ce701798fd93df348c1d48c3ebfd7455de0db8fe1d3d0b6a2e5236108b'),sealed_source(S,'source-owner-v7.4-mass-cache-20261101','33e76e550a1eea93f7650fb7a5145baea8ee195a909b29ca2d85f26cb578c602'),sealed_source(ADAPTER,'source-stock-reference-v2-surface-birth-20261101','7d015d1df6ab7c0d1dfa5a3382c27d071d21b916930338b48e484b8ddcea1b3a')]
 b=SOURCE['usable']('c1');u=SOURCE['usable']('u1')
 runs=sorted((SOURCE['OUT']/'c1').glob('run-*/run-receipt.json'));assert runs and json.loads(runs[-1].read_text())['all_passed'],'source v7.4 final qualification required'
 a=N/'build/libpic_native_v8.a';assert sha(a)=='bb3def231e5d6d1506016238ff7261d4de284a47485863381be2cc754ebe84b8','native v8.7 archive changed'
 files=[*[x for x in sorted((R/'gas-target-owner-v1.7.1-census-transport-20261101/include').rglob('*')) if x.is_file()],R/'gas-target-owner-v1.7.1-census-transport-20261101/src/gas_owner.cpp',R/'native-gas-adapter-v1.2-immutable-export-20261101/include/native_gas_stage.hpp',*[x for x in sorted((HALL/'source').rglob('*')) if x.is_file()],*[x for x in sorted((ARRIVAL/'include').rglob('*')) if x.is_file()],ADAPTER/'include/source_stock_reference.hpp',ADAPTER/'tests/test_reference.cpp',ADAPTER/'tests/surface_mapping_checks.inc',a,N/'build/kinetic-boundary.o',F,R/'arch-base-20261028/root-files.json',R/'arch-base-20261028/root-files.sha256',*SOURCE['LIBS'],S/'run_source_v7.py',R/'pic-water-common-commit-v8.0.4-b3-20261101/NATIVE-FLAGS.json']
 for dep in ('native-gas-common-stage-v1.3.1-motion-unwind-20261101/include','pic-water-common-genesis-v9-20261101/include','pic-source-hall-iteration-v1-20261101/include','pic-positive-source-fixture-v1-20261101/include'):
  files.extend(q for q in (R/dep).rglob('*') if q.is_file())
 files.append(R/'reaction-final-charge-guard-v1-20261101/final_charge_guard.hpp')
 return {'sealed_components':seals,'source_build_c1_sha256':sha(SOURCE['OUT']/'c1/build/build-receipt.json'),'source_build_u1_sha256':sha(SOURCE['OUT']/'u1/build/build-receipt.json'),'source_qualification_sha256':sha(runs[-1]),'source_files':SOURCE['own'](),'files':{str(p):sha(p) for p in files}}
def execute(cmd,log):
 with log.open('w') as f:r=subprocess.run([str(F),*map(str,cmd)],cwd=P,env=SOURCE['ENV'],stdout=f,stderr=subprocess.STDOUT,timeout=900)
 print(log.name,r.returncode,flush=True)
 if r.returncode:raise RuntimeError(log.read_text()[-6000:])
def qualify(label):
 adm=SOURCE['admission']('run');series=OUT/label;series.mkdir(parents=True,exist_ok=True)
 attempt=1
 while (series/('attempt-%03d'%attempt)).exists():attempt+=1
 actual=label+'/attempt-%03d'%attempt;out=OUT/actual;out.mkdir(parents=True,exist_ok=False);build=P/'build'/actual;build.mkdir(parents=True,exist_ok=False)
 print('qualification_attempt',actual,flush=True)
 before=own();receipt={'admission':adm,'own_sha256':before,'variants':{},'all_passed':False}
 try:
  deps=pins();receipt['pins']=deps
  base=json.loads((R/'pic-water-common-commit-v8.0.4-b3-20261101/NATIVE-FLAGS.json').read_text())['flags'];base=[base[0],'-I'+str(P/'include'),'-I'+str(R/'pic-water-common-commit-v8.0.4-b3-20261101/source'),'-I'+str(R/'native-gas-adapter-v1.2-immutable-export-20261101/include'),'-I'+str(R/'gas-target-owner-v1.7.1-census-transport-20261101/include'),'-I'+str(N/'forms/include/pf_v7'),'-I'+str(HALL/'source'),'-I'+str(R/'pic-water-common-commit-v8.0.4-b3-20261101/vendor'),'-I'+str(ARRIVAL/'include'),*base[1:]];base+=[*SOURCE['owner_flags']([])[1:],'-I'+str(ADAPTER/'include'),'-I'+str(N/'tests')]
  includes=[R/'native-gas-common-stage-v1.3.1-motion-unwind-20261101/include',R/'pic-water-common-genesis-v9-20261101/include',R/'pic-source-hall-iteration-v1-20261101/include',R/'reaction-final-charge-guard-v1-20261101',R/'pic-positive-source-fixture-v1-20261101/include',R.parent/'Power-Flow/native/coupled_pic/include']
  base=[base[0],*['-I'+str(q) for q in includes],*base[1:]]
  base=[x.replace('native-full-time-v8.7-reservoir-stock-20261101','native-full-time-v8.9-motion-inventory-20261101').replace('d795695aca20ddef0efadedaebea8baa6c994e819514c0258bbac290fb66103d','8e8d9553ccac604e7e69ecc3354f124e71643a17e8fbc85b9883dd4f798897db') for x in base]
  # Native v8.7 include comes first; shared native types never use candidate003.
  base=[x.replace(str(SOURCE['N']/'source/include'),str(N/'source/include')).replace(str(SOURCE['N']/'inputs'),str(N/'inputs')) for x in base]
  for variant in ('normal','ubsan'):
   dest=build/variant;dest.mkdir();logs=out/variant;logs.mkdir();flags=base+SOURCE['VARIANTS'][variant];link=SOURCE['LINK_SAN'][variant];objects={
    'test_stage':P/'tests/test_common.cpp','fast_emission':S/'vendor/fast/fast_emission.cpp','emission_hot_v7':S/'native/mfem_electrode_heat_source_v7/emission_hot_v7.cpp','surface_emission':S/'vendor/heat-v2/source-kernel/src/surface_emission.cpp','micro_layer':S/'vendor/heat-v2/source-kernel/src/micro_layer.cpp','motion_provider_v3':HALL/'source/motion_provider.cpp','gas_owner':R/'gas-target-owner-v1.7.1-census-transport-20261101/src/gas_owner.cpp'}
   for name,src in objects.items():execute([*flags,'-MD','-MF',dest/(name+'.d'),'-c',src,'-o',dest/(name+'.o')],logs/(name+'.log'))
   executables={name:dest/name for name in ('test_stage',)}
   for name,exe in executables.items():
    execute(['/usr/bin/mpic++','-Wl,--gc-sections',*link,*[dest/(key+'.o') for key in objects if key not in ('test_stage',) or key==name],N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*SOURCE['LIBS'],'-lcrypto','-ldl','-o',exe],logs/(name+'-link.log'))
   paths=set()
   for d in dest.glob('*.d'):
    for line in d.read_text().replace('\\\n',' ').splitlines():
     if ':' in line:paths.update(line.split(':',1)[1].split())
   compiled={str(SOURCE['dependency_path'](x)):sha(SOURCE['dependency_path'](x)) for x in paths};binaries={str(x):sha(x) for x in [*executables.values(),*dest.glob('*.o')]}
   result=logs/'result';result.mkdir();execute([executables['test_stage'],SOURCE['TRACE'],result],logs/'test.log');data=json.loads((result/'common.json').read_text());assert data['passed']
   baseline_path=R.parent/'.workspace/outputs/common-test-v9-20261101/c3/attempt-001/normal/result/common.json'
   baseline=json.loads(baseline_path.read_text())
   compared=('step','begin_s','end_s','final_carriers','sha256','final_census','stage1_events','stage2_events')
   assert [{k:x[k] for k in compared} for x in data['steps']]==[{k:x[k] for k in compared} for x in baseline['steps']],'actual common roots/census changed'
   receipt['baseline_common_parity']={'path':str(baseline_path),'sha256':sha(baseline_path),'fields':compared,'identical':True}
   assert all(sha(p)==h for p,h in compiled.items()) and all(sha(p)==h for p,h in binaries.items()),'compiled input/binary changed'
   receipt['variants'][variant]={'stage_passed':len(data['cases']),'dependency_sha256':compiled,'binary_sha256':binaries,'test_sha256':sha(result/'common.json'),'sanitizer_scope':'adapter/source/Hall translation units; shared native archive separately qualified'}
  assert own()==before and pins()==deps,'source/dependency changed'
  receipt['all_passed']=True
 except Exception as e:receipt['error']=str(e);print(str(e),flush=True)
 (out/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');return 0 if receipt['all_passed'] else 1
if __name__=='__main__':
 if len(sys.argv)!=2 or '/' in sys.argv[1]:raise SystemExit('run_source_stage.py LABEL')
 raise SystemExit(qualify(sys.argv[1]))
