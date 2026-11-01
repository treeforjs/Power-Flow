#!/usr/bin/python -I
"""Guarded component qualification; one explicit-memory SimQueue CPU, frozen userspace."""
from pathlib import Path
import hashlib,json,os,runpy,subprocess,sys,time
P=Path(__file__).resolve().parent
R=P.parent
S=R/'pic-electrode-heat-source-v7.2-prepared-state-20261101'
N=R/'pic-reaction-transactions-20261029/native-full-time-v8.7-reservoir-stock-20261101'
OUT=R.parent/'.workspace/outputs/source-stock-reference-v1-20261101'
F=R/'arch-base-20261028/run-frozen'
SOURCE=runpy.run_path(str(S/'run_source_v7.py'),run_name='source_metadata')
def sha(p):return SOURCE['sha'](p)
def own():return {str(p):sha(p) for p in sorted(P.rglob('*')) if p.is_file() and 'build' not in p.parts and p.suffix in ('.py','.hpp','.cpp','.inc','.json','.txt')}
def pins():
 b=SOURCE['usable']('c4');u=SOURCE['usable']('u2')
 runs=sorted((SOURCE['OUT']/'c4').glob('run-*/run-receipt.json'));assert runs and json.loads(runs[-1].read_text())['all_passed'],'source v7.2 final qualification required'
 a=N/'build/libpic_native_v8.a';assert sha(a)=='1778aaecc0d5b9ef2c7b296a63347de44d6df0fccac28fc0f74febdc3736cb1e','native v8.7 archive changed'
 files=[a,N/'build/kinetic-boundary.o',F,R/'arch-base-20261028/root-files.json',R/'arch-base-20261028/root-files.sha256',*SOURCE['LIBS'],S/'run_source_v7.py',P/'NATIVE-FLAGS.json']
 return {'source_build_c4_sha256':sha(SOURCE['OUT']/'c4/build/build-receipt.json'),'source_build_u2_sha256':sha(SOURCE['OUT']/'u2/build/build-receipt.json'),'source_qualification_sha256':sha(runs[-1]),'source_files':SOURCE['own'](),'files':{str(p):sha(p) for p in files}}
def execute(cmd,log):
 with log.open('w') as f:r=subprocess.run([str(F),*map(str,cmd)],cwd=P,env=SOURCE['ENV'],stdout=f,stderr=subprocess.STDOUT,timeout=900)
 print(log.name,r.returncode,flush=True)
 if r.returncode:raise RuntimeError(log.read_text()[-6000:])
def qualify(label):
 adm=SOURCE['admission']('run');out=OUT/label;out.mkdir(parents=True,exist_ok=False);build=P/'build'/label;build.mkdir(parents=True,exist_ok=False)
 before=own();deps=pins();receipt={'admission':adm,'own_sha256':before,'pins':deps,'variants':{},'all_passed':False}
 try:
  base=json.loads((P/'NATIVE-FLAGS.json').read_text())['flags'];base+=[*SOURCE['owner_flags']([])[1:],'-I'+str(P/'include'),'-I'+str(N/'tests')]
  # Native v8.7 include comes first; shared native types never use candidate003.
  base=[x.replace(str(SOURCE['N']/'source/include'),str(N/'source/include')).replace(str(SOURCE['N']/'inputs'),str(N/'inputs')) for x in base]
  for variant in ('normal','ubsan'):
   dest=build/variant;dest.mkdir();logs=out/variant;logs.mkdir();flags=base+SOURCE['VARIANTS'][variant];link=SOURCE['LINK_SAN'][variant];objects={
    'test_reference':P/'tests/test_reference.cpp','fast_emission':S/'vendor/fast/fast_emission.cpp','emission_hot_v7':S/'native/mfem_electrode_heat_source_v7/emission_hot_v7.cpp','surface_emission':S/'vendor/heat-v2/source-kernel/src/surface_emission.cpp','micro_layer':S/'vendor/heat-v2/source-kernel/src/micro_layer.cpp','motion_provider_v3':SOURCE['HV3']/'source/motion_provider.cpp'}
   for name,src in objects.items():execute([*flags,'-MD','-MF',dest/(name+'.d'),'-c',src,'-o',dest/(name+'.o')],logs/(name+'.log'))
   exe=dest/'test_reference';execute(['/usr/bin/mpic++','-Wl,--gc-sections',*link,*[dest/(x+'.o') for x in objects],N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*SOURCE['LIBS'],'-lcrypto','-ldl','-o',exe],logs/'link.log')
   paths=set()
   for d in dest.glob('*.d'):
    for line in d.read_text().replace('\\\n',' ').splitlines():
     if ':' in line:paths.update(line.split(':',1)[1].split())
   compiled={str(SOURCE['dependency_path'](x)):sha(SOURCE['dependency_path'](x)) for x in paths};binaries={str(x):sha(x) for x in [exe,*dest.glob('*.o')]}
   result=logs/'result';result.mkdir();execute([exe,result],logs/'test.log');data=json.loads((result/'reference.json').read_text());assert data['failed']==0
   assert all(sha(p)==h for p,h in compiled.items()) and all(sha(p)==h for p,h in binaries.items()),'compiled input/binary changed'
   receipt['variants'][variant]={'passed':len(data['cases']),'dependency_sha256':compiled,'binary_sha256':binaries,'test_sha256':sha(result/'reference.json'),'sanitizer_scope':'adapter/source translation units; shared native archive separately qualified'}
  assert own()==before and pins()==deps,'source/dependency changed'
  receipt['all_passed']=True
 except Exception as e:receipt['error']=str(e);print(str(e),flush=True)
 (out/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');return 0 if receipt['all_passed'] else 1
if __name__=='__main__':
 if len(sys.argv)!=2 or '/' in sys.argv[1]:raise SystemExit('run_reference.py LABEL')
 raise SystemExit(qualify(sys.argv[1]))
