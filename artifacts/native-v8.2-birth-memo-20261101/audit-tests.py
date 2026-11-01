"""Targeted audit qualification. Submit through SimQueue with one CPU, explicit memory;
GPU mode additionally requires a GPU share. No timing conclusions are drawn."""
import sys
from tooling import *
mode=sys.argv[1]
resource_guard(gpu=mode=='gpu')
label=sys.argv[2] if len(sys.argv)>2 else 'final-'+mode
out=OUT/label
r=Runner(out,'audit')
inputs=sealed_inputs()
B=P/'build'
if mode=='cpu':
    for directory in ['flows-v5','work-gross']:
        (out/directory).mkdir(parents=True,exist_ok=True)
    cmds=[('sinks',[B/'test_v8_sinks',out/'sinks']),
          ('target-energy',[B/'test_audit_target_guard',out/'target-energy']),
          ('parcels',[B/'test_audit_parcels_codec']),
          ('chemistry',[B/'test_v8_chemistry',out/'chemistry','full']),
          ('flows-v5',[B/'test_flows',out/'flows-v5','reduced']),
          ('work-gross',[B/'test_work_gross',out/'work-gross']),
          ('stage-b',[B/'test_v8_stage_b',out/'stage-b'])]
    for k in [2,3]:
        for action in ['unbroken','save','restore']+(['tamper'] if k==3 else []):
            cmds.append((f'restore-{k}-{action}',[B/'test_v8_restore',out/f'restore-{k}',action,k]))
elif mode=='gpu':
    need(bool(os.environ.get('SIMQ_GPUS')),'SimQueue GPU allocation required')
    cmds=[('device-only-verified',[B/'test_step_timing',out/'device-only.json','device_verified',1000,2,45,NVRTC,GPU_UUID]),
          ('flows-gpu',[B/'test_v8_flows_gpu',out/'flows-gpu',NVRTC,GPU_UUID,32]),
          ('chemistry-gpu',[B/'test_v8_chemistry_gpu',out/'chemistry-gpu',NVRTC,GPU_UUID]),
          ('stage-b-gpu',[B/'test_v8_stage_b',out/'stage-b-gpu',NVRTC,GPU_UUID])]
else: raise ValueError(mode)
if len(sys.argv)>3:
    selected=set(sys.argv[3:])
    need(selected <= {name for name, _ in cmds}, 'unknown audit control selection')
    cmds=[(name,argv) for name,argv in cmds if name in selected]
for name,argv in cmds:
    r.run(argv,name+'.log',timeout=240,env={'PF_V8_CHEM':'v3.1'},rlimit=3*1024**3,limit=mode!='gpu')
need(sealed_inputs()==inputs,'dependency inputs changed during audit')
(out/'audit-result.json').write_text(json.dumps({'passed':True,'mode':mode,'source_identity':json.loads((P/'SOURCE-SHA256.json').read_text())['source_identity'],
    'library_sha256':sha256(B/'libpic_native_v8.a'),'commands':r.records,
    'note':'component qualification only; target energy/momentum owner absent; no full-system energy closure'},indent=2)+'\n')
print(mode,'audit controls passed',flush=True)
