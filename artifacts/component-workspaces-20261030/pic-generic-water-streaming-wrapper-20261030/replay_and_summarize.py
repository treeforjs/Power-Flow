from pathlib import Path
import os,json,subprocess,resource
P=Path(__file__).resolve().parent
O=Path(os.environ.get('PIC_STREAMING_WRAPPER_OUTPUT','/home/treefor/sim-workspace/powerflow/outputs/pic-generic-water-streaming-wrapper-20261030'))
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(180,180))
replay=O/'independent-process-replay';replay.mkdir(exist_ok=True)
cmd=[str(P/'build/test_wrapper'),str(P/'provider/products/data/collisions/h2o-gas-20261029-elastic.json'),str(O/'run-binding.json'),str(replay),'broad']
with (replay/'run.log').open('w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=150)
old=json.loads((O/'broad-cases.json').read_text());new=json.loads((replay/'broad-cases.json').read_text())
for records in [old,new]:
 for row in records:
  del row['wrapper_wall_seconds'];del row['peak_RSS_KiB']
assert old==new
small=json.loads((O/'small-results.json').read_text());broad=json.loads((O/'broad-results.json').read_text())
summary={'status':'passed','assertions':small['assertions']+broad['assertions'],'independent_process_case_receipts_identical':len(old),'replay_command':cmd,'physical_kernel_changed':False,'native_field_or_joint_execution':False,'cases':[]}
for row in json.loads((O/'broad-cases.json').read_text()):
 ev=row['evidence'];stream=ev['stream']
 summary['cases'].append({'name':row['name'],'wrapper_wall_seconds':row['wrapper_wall_seconds'],'candidate_rows':row['candidate_rows'],'candidate_bytes_sha256':row['candidate_bytes_sha256'],'pair_visits':row['work']['pair_visits'],'physical_channel_queries':row['work']['channel_queries'],'retained_category_slots':stream['retained_category_slots'],'stream_queries':stream['queries'],'commitment_sha256':stream['content_sha256'],'diagnostic_source_queries':ev['diagnostic_source_queries'],'diagnostic_lookup_reuses':ev['diagnostic_lookup_reuses'],'actual_selected_return_seeds':ev['actual_selected_return_seeds'],'last_energy_record_reuses':ev['last_energy_record_reuses']})
(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
