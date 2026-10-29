import pathlib,json,copy,subprocess,sys,hashlib,ast
import compare_saved_pic as c
OUT=pathlib.Path(__file__).parent
ROOT=OUT.parent.parent

def main():
 cases=[]
 def case(name,fn,accepted):
  try:fn();actual=True;error=None
  except (c.Invalid,ValueError,OSError,KeyError,TypeError)as e:actual=False;error=str(e)
  assert actual==accepted,(name,actual,error)
  cases.append({'name':name,'expected_accept':accepted,'accepted':actual,'reason':error})
 with c.cpu_lease()as execution:
  proof=c.strict_json(c.ALLOCATION)
  receipt=c.strict_json(ROOT/'phase-baseline-late-cpu7/run-receipt.json')
  p6=c.replay_policy();p7=c.replay_policy(7)
  case('CPU6 default with CPU6 receipt',lambda:c.validate_replay_receipt({'profile_cpu':6},p6),True)
  case('CPU7 is not accepted by default6',lambda:c.validate_replay_receipt(receipt,p6),False)
  case('CPU7 explicit with allocation and recorded affinity/stopped identity',lambda:c.validate_replay_receipt(receipt,p7),True)
  case('CPU6 receipt rejected when CPU7 explicitly expected',lambda:c.validate_replay_receipt({'profile_cpu':6},p7),False)
  case('CPU8 requested',lambda:c.replay_policy(8),False)
  case('CPU8 receipt against CPU6 policy',lambda:c.validate_replay_receipt({'profile_cpu':8},p6),False)
  bad=copy.deepcopy(receipt);bad['profile_cpu']=8
  case('CPU8 receipt against CPU7 policy',lambda:c.validate_replay_receipt(bad,p7),False)
  case('CPU7 missing allocation',lambda:c.replay_policy(7,OUT/'missing-allocation.json'),False)
  wrong=copy.deepcopy(proof);wrong['production_cgroup_frozen']=False
  path=OUT/'fixture-allocation-not-frozen.json';path.write_text(json.dumps(wrong))
  case('CPU7 allocation not frozen',lambda:c.replay_policy(7,path),False)
  bad=copy.deepcopy(receipt);bad['production_left_stopped']=False
  case('CPU7 receipt did not preserve stopped production',lambda:c.validate_replay_receipt(bad,p7),False)
  bad=copy.deepcopy(receipt);bad['production_identity']['pid']+=1
  case('CPU7 identity proof mismatched',lambda:c.validate_replay_receipt(bad,p7),False)
  bad=copy.deepcopy(receipt);bad['command']=[]
  case('CPU7 command affinity proof missing',lambda:c.validate_replay_receipt(bad,p7),False)
  bad=copy.deepcopy(receipt);i=bad['command'].index('taskset');bad['command'][i+2]='8'
  case('CPU7 label with CPU8 command rejected',lambda:c.validate_replay_receipt(bad,p7),False)
  r=c.compare(ROOT/'baseline-cpu6-run1',ROOT/'cpu-candidate-cpu6-run1',ROOT/'seed-current')
  assert r['status']=='passed';r['execution']=execution
  (OUT/'comparison-default6.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
  old=c.strict_json(OUT.parent/'checkpoint-validation/comparison.json')
  keys=['policy','options','particle_counts','particles','raw_diagnostics','seed_to_final_public_ledger_increments','history','source_reservoirs','frame','density_analysis','field_error_budget_impulse_m_s']
  unchanged={key:r[key]==old[key]for key in keys};assert all(unchanged.values())
  case('finite checks still reject a v1 nonfinite fixture',lambda:c.compare(ROOT/'baseline-cpu6-run1',c.strict_json(OUT.parent/'checkpoint-validation/fixture-tests.json')['fixtures']+'/particle_position_m_0_nan',ROOT/'seed-current'),False)
  before=c.strict_json(OUT/'v1-preservation-before.json')
  preserved={name:hashlib.sha256((pathlib.Path(before['v1'])/name).read_bytes()).hexdigest()==digest for name,digest in before['files'].items()};assert all(preserved.values())
 # Exercise actual CLI forwarding outside the internally leased API tests.
 command=[sys.executable,str(OUT/'compare_saved_pic.py'),str(ROOT/'phase-baseline-late-cpu7'),str(ROOT/'phase-cpu-v2-late-cpu7'),str(ROOT/'seed-late-phase'),'--expected-replay-cpu','7','--output',str(OUT/'comparison-phase-cpu7.json')]
 proc=subprocess.run(command,capture_output=True,text=True)
 (OUT/'cpu7-cli.log').write_text(proc.stdout+proc.stderr)
 assert proc.returncode==0,(proc.returncode,proc.stdout,proc.stderr)
 cpu7=c.strict_json(OUT/'comparison-phase-cpu7.json');assert cpu7['execution']['cpu']==6 and cpu7['replay_cpu_policy']['expected_replay_cpu']==7
 bad_command=list(command);bad_command[bad_command.index('7')]='8'
 proc8=subprocess.run(bad_command,capture_output=True,text=True)
 assert proc8.returncode==2 and 'invalid choice' in proc8.stderr
 (OUT/'cpu8-cli.log').write_text(proc8.stdout+proc8.stderr)
 report={'status':'passed','policy_cases':len(cases),'cases':cases,'default6_real_comparison_passed':True,'default6_scientific_metrics_identical_to_v1':unchanged,'explicit7_real_phase_comparison_passed':True,'cpu7_cli_command':command,'CPU8_cli_returncode':proc8.returncode,'v1_all_original_files_preserved':preserved,'execution':execution,'allocation_receipt':cpu7['replay_cpu_policy']['allocation']}
 (OUT/'policy-tests.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'status':'passed','policy_cases':len(cases),'default6':'passed; scientific metrics identical','explicit7':'actual phase outputs passed','CPU8':'rejected','v1_preserved':all(preserved.values())},indent=2))
if __name__=='__main__':main()
