"""Inherited native host/device controls against the qualified RNG2 legacy table."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tooling import *
resource_guard(gpu=True)
need(bool(os.environ.get('SIMQ_GPUS')),'scheduled GPU required')
out=OUT/sys.argv[1]
need(not out.exists(),'preserve evidence')
r=Runner(out,'legacy-gpu')
inputs=sealed_inputs()
manifest=json.loads((P/'SOURCE-SHA256.json').read_text())
source={n:sha256(P/n) for n in manifest['files']}
need(source==manifest['files'],'source changed from build')
fixture=tree_hashes(P/'tests')
exe=P/'build/test_v8_chemistry_gpu'
pins={str(f):sha256(f) for f in [exe,P/'build/libpic_native_v8.a',P/'qualify-legacy-gpu.py',P/'tooling.py']}
r.run([exe,out/'native',NVRTC,GPU_UUID],'legacy.log',timeout=240,env={'PF_V8_CHEM':'v3.1'},limit=False)
result=json.loads((out/'native/v8-chemistry-gpu-tests.json').read_text())
need(result['failures']==0 and result['checks']>0,'native legacy checks failed')
need(inputs==sealed_inputs() and source=={n:sha256(P/n) for n in source} and fixture==tree_hashes(P/'tests'),'inputs changed')
need(pins=={f:sha256(f) for f in pins},'binary/harness changed')
(out/'qualification.json').write_text(json.dumps({'passed':True,'source_identity':manifest['source_identity'],'result':result,'pins':pins,'fixture_pins':fixture,'inputs':inputs,'commands':r.records,'production_qualified':False,'network_complete':False,'note':'legacy physical model; RNG2 reference, not old-seed parity'},indent=2)+'\n')
print('native legacy GPU parity passed')
