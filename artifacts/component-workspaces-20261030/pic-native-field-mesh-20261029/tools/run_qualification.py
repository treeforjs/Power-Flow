#!/usr/bin/env python3
"""Run only bounded manufactured CPU controls under an external core7 lease."""
import pathlib,os,json,subprocess,time,hashlib
P=pathlib.Path('/home/treefor/projects/runtimes/pic-native-field-mesh-20261029')
O=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs/pic-native-field-mesh-20261029')
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
def main():
 assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
 attempt='v3';records=[];identity=json.loads((P/'native-source-identity.json').read_text())
 for name,value in identity['files'].items():assert hashlib.sha256((P/name).read_bytes()).hexdigest()==value,'Source identity changed'
 def run(name,argv,limit):
  start=time.monotonic()
  with (O/(name+'-'+attempt+'.log')).open('x') as log:
   result=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=limit,check=False)
  records.append(dict(name=name,argv=argv,returncode=result.returncode,wall_s=time.monotonic()-start));(O/('qualification-commands-'+attempt+'.json')).write_text(json.dumps(records,indent=2)+'\n');assert result.returncode==0,(name,(O/(name+'-'+attempt+'.log')).read_text()[-2000:])
 run('configure-final',[FROZEN,'cmake','-S',str(P),'-B',str(P/'build'),'-G','Ninja','-DCMAKE_BUILD_TYPE=Release','-DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem','-DPIC_REACTION_NATIVE_SOURCE_SHA256='+identity['sha256']],20)
 run('build-final',[FROZEN,'cmake','--build',str(P/'build'),'--parallel','1'],180)
 run('field-tests-final',[FROZEN,str(P/'build/test_native_field_mesh'),'/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture',str(O/('cpu-controls-'+attempt))],60)
 run('foundation-final',[FROZEN,str(P/'build/test_native_foundation')],60)
 run('q1-final',[FROZEN,str(P/'build/test_q1_scalar_gather')],30)
 (O/('reaction-tests-'+attempt)).mkdir()
 run('reactions-final',[FROZEN,str(P/'build/test_native_reactions'),str(O/('reaction-tests-'+attempt))],60)
 reaction=json.loads((O/('reaction-tests-'+attempt)/'native-tests.json').read_text());assert reaction['status']=='passed' and reaction['check_count']==114
 (O/'qualification-result.json').write_text(json.dumps(dict(status='passed',source_sha256=identity['sha256'],native_field_checks=json.loads((O/('cpu-controls-'+attempt)/'native-field-tests.json').read_text()),unchanged_reaction_checks=114,foundation_passed=True,q1_passed=True,GPU_execution=False,production_actions=False,actual_saved_state_evolution=False),indent=2)+'\n')
 print('Native field geometry and unchanged foundation/Q1/114 reaction controls PASS')
if __name__=='__main__':main()
