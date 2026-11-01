import subprocess,sys
from tooling import *
subprocess.run([sys.executable,"build-components.py","sanitized","core","hall","v83-export","v8-restore"],check=True)
resource_guard()
pins=json.loads((P/'SOURCE-SHA256.json').read_text())
archive=sha256(P/'build/sanitized-libpic_native_v8.a')
inputs=sealed_inputs()
fixtures=tree_hashes(P/'tests')
binaries={name:sha256(P/'build'/name) for name in ['sanitized-test_v83_export','sanitized-test_v8_restore']}
r=Runner(OUT/'sanitized-export','sanitized')
r.run([P/'build/sanitized-test_v83_export',OUT/'sanitized-export/controls'],'controls.log',timeout=240,env={'PF_V8_CHEM':'v3.1'},sanitized=True)
for action in ['unbroken','save','restore','tamper']:
    r.run([P/'build/sanitized-test_v8_restore',OUT/'sanitized-export/restore',action,'3'],'restore-'+action+'.log',timeout=240,env={'PF_V8_CHEM':'v3.1'},sanitized=True)
need(inputs==sealed_inputs(),'dependency changed')
need(fixtures==tree_hashes(P/'tests'),'fixtures changed')
need(all(sha256(P/name)==digest for name,digest in pins['files'].items()),'source changed')
need(archive==sha256(P/'build/sanitized-libpic_native_v8.a'),'archive changed')
need(binaries=={name:sha256(P/'build'/name) for name in binaries},'binary changed')
(OUT/'sanitized-export/result.json').write_text(json.dumps({'passed':True,'source_identity':pins['source_identity'],'archive_sha256':archive,'binary_pins':binaries,'fixture_pins':fixtures,'dependencies':inputs,'commands':r.records},indent=2)+'\n')
