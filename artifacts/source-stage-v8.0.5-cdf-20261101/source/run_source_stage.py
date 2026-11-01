#!/usr/bin/python -I
"""Guarded component qualification; one explicit-memory SimQueue CPU, frozen userspace."""
from pathlib import Path
import hashlib,json,os,runpy,subprocess,sys,time
P=Path(__file__).resolve().parent
R=P.parent
S=R/'pic-electrode-heat-source-v7.5-cdf-registry-20261101'
HALL=R/'hall-corrective-v3.3-b3-20261101'
ARRIVAL=R/'pic-arrival-owner-v1.2-publication-20261101'
N=R/'pic-reaction-transactions-20261029/native-full-time-v8.7-cdf-v1-20261101'
OUT=R.parent/'.workspace/outputs/source-stage-v8.0.5-cdf-20261101'
BASELINE=R/'hall-corrective-v3.3-20261031'
ADAPTER=R/'pic-source-stock-reference-v3-cdf-registry-20261101'
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
 seals=[sealed_source(HALL,'hall-corrective-v3.3-b3-20261101','231b92ce701798fd93df348c1d48c3ebfd7455de0db8fe1d3d0b6a2e5236108b'),sealed_source(S,'source-owner-v7.5-cdf-registry-20261101','9a59cda7f6f1779d1d9cac16ddcebf48f870aaf97f8032e42f574404719d11bf'),sealed_source(ADAPTER,'source-stock-reference-v3-cdf-registry-20261101','06101849b9ae9d99a16fe58f3cf0676eacdea293fb5060c61b6eac8118fe306c')]
 b=SOURCE['usable']('c1');u=SOURCE['usable']('u1')
 runs=sorted((SOURCE['OUT']/'c1').glob('run-*/run-receipt.json'));assert runs and json.loads(runs[-1].read_text())['all_passed'],'source v7.5 final qualification required'
 native_pin=json.loads((P/'NATIVE-CDF-PINS.json').read_text());assert native_pin['qualified'] and all(sha(Path(k))==v for k,v in native_pin['files'].items()),'qualified native CDF pins changed'
 a=N/'build/libpic_native_v8.a';assert sha(a)=='aceddaf950dd43dd5e88e3b4e8d59906aa4cc9f0e266d4b68432453fcc155dfe','native v8.7 archive changed'
 files=[P/'NATIVE-CDF-PINS.json',*[Path(x) for x in native_pin['files']],R/'chemistry-dynamic-cdf-v1.1-20261101/data/chemistry-cdf.blob',*[x for x in sorted((HALL/'source').rglob('*')) if x.is_file()],*[x for x in sorted((ARRIVAL/'include').rglob('*')) if x.is_file()],ADAPTER/'include/source_stock_reference.hpp',ADAPTER/'tests/test_reference.cpp',ADAPTER/'tests/surface_mapping_checks.inc',a,N/'build/kinetic-boundary.o',F,R/'arch-base-20261028/root-files.json',R/'arch-base-20261028/root-files.sha256',*SOURCE['LIBS'],S/'run_source_v7.py',P/'NATIVE-FLAGS.json']
 return {'sealed_components':seals,'source_build_c1_sha256':sha(SOURCE['OUT']/'c1/build/build-receipt.json'),'source_build_u1_sha256':sha(SOURCE['OUT']/'u1/build/build-receipt.json'),'source_qualification_sha256':sha(runs[-1]),'source_files':SOURCE['own'](),'files':{str(p):sha(p) for p in files}}
def execute(cmd,log):
 with log.open('w') as f:r=subprocess.run([str(F),*map(str,cmd)],cwd=P,env={**SOURCE['ENV'],'PF_SOURCE_REGISTRY':'cdf'},stdout=f,stderr=subprocess.STDOUT,timeout=900)
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
  base=json.loads((P/'NATIVE-FLAGS.json').read_text())['flags'];base=[base[0],'-I'+str(P/'tests/fixtures'),'-I'+str(P/'source'),'-I'+str(HALL/'source'),'-I'+str(P/'vendor'),'-I'+str(ARRIVAL/'include'),*base[1:]];base+=[*SOURCE['owner_flags']([])[1:],'-I'+str(ADAPTER/'include'),'-I'+str(N/'tests')]
  # Native v8.7 include comes first; shared native types never use candidate003.
  base=[x.replace(str(SOURCE['N']/'source/include'),str(N/'source/include')).replace(str(SOURCE['N']/'inputs'),str(N/'inputs')) for x in base]
  for variant in ('normal','ubsan'):
   dest=build/variant;dest.mkdir();logs=out/variant;logs.mkdir();flags=base+SOURCE['VARIANTS'][variant];link=SOURCE['LINK_SAN'][variant];objects={
    'test_stage':P/'tests/test_source_stage.cpp','test_reference':ADAPTER/'tests/test_reference.cpp','fast_emission':S/'vendor/fast/fast_emission.cpp','emission_hot_v7':S/'native/mfem_electrode_heat_source_v7/emission_hot_v7.cpp','surface_emission':S/'vendor/heat-v2/source-kernel/src/surface_emission.cpp','micro_layer':S/'vendor/heat-v2/source-kernel/src/micro_layer.cpp','motion_provider_v3':HALL/'source/motion_provider.cpp'}
   for name,src in objects.items():execute([*flags,'-MD','-MF',dest/(name+'.d'),'-c',src,'-o',dest/(name+'.o')],logs/(name+'.log'))
   executables={name:dest/name for name in ('test_stage','test_reference')}
   for name,exe in executables.items():
    execute(['/usr/bin/mpic++','-Wl,--gc-sections',*link,*[dest/(key+'.o') for key in objects if key not in ('test_stage','test_reference') or key==name],N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*SOURCE['LIBS'],'-lcrypto','-ldl','-o',exe],logs/(name+'-link.log'))
   paths=set()
   for d in dest.glob('*.d'):
    for line in d.read_text().replace('\\\n',' ').splitlines():
     if ':' in line:paths.update(line.split(':',1)[1].split())
   compiled={str(SOURCE['dependency_path'](x)):sha(SOURCE['dependency_path'](x)) for x in paths};binaries={str(x):sha(x) for x in [*executables.values(),*dest.glob('*.o')]}
   result=logs/'result';result.mkdir();execute([executables['test_stage'],SOURCE['TRACE'],result],logs/'test.log');data=json.loads((result/'source-stage.json').read_text());assert data['failed']==0
   reference_result=logs/'reference-result';reference_result.mkdir();execute([executables['test_reference'],reference_result],logs/'reference-test.log');reference=json.loads((reference_result/'reference.json').read_text());assert reference['failed']==0 and len(reference['cases'])==85
   assert all(sha(p)==h for p,h in compiled.items()) and all(sha(p)==h for p,h in binaries.items()),'compiled input/binary changed'
   receipt['variants'][variant]={'stage_passed':len(data['cases']),'reference_passed':len(reference['cases']),'dependency_sha256':compiled,'binary_sha256':binaries,'test_sha256':sha(result/'source-stage.json'),'reference_test_sha256':sha(reference_result/'reference.json'),'result_files':{str(f.relative_to(logs)):{'sha256':sha(f),'bytes':f.stat().st_size} for directory in (result,reference_result) for f in directory.rglob('*') if f.is_file()},'sanitizer_scope':'adapter/source/Hall translation units; shared native archive separately qualified'}
  assert own()==before and pins()==deps,'source/dependency changed'
  receipt['all_passed']=True
 except Exception as e:receipt['error']=str(e);print(str(e),flush=True)
 (out/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');return 0 if receipt['all_passed'] else 1
if __name__=='__main__':
 if len(sys.argv)!=2 or '/' in sys.argv[1]:raise SystemExit('run_source_stage.py LABEL')
 raise SystemExit(qualify(sys.argv[1]))
