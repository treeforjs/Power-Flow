from tooling import *
resource_guard()
pins=json.loads((P/"SOURCE-SHA256.json").read_text())
lib=sha256(P/"build/libpic_native_v8.a")
inputs=sealed_inputs()
fixture_pins=tree_hashes(P/'tests')
binary_pin=sha256(P/'build/test_v83_export')
harness_pins={name:sha256(P/name) for name in ['check-export.py','tooling.py']}
label=__import__("sys").argv[1] if len(__import__("sys").argv)>1 else "first-export"
r=Runner(OUT/label, "export")
r.run([P/"build/test_v83_export", OUT/label/"controls"], "controls.log", timeout=240, env={"PF_V8_CHEM":"v3.1"})
need(sealed_inputs()==inputs,"dependency changed")
need(tree_hashes(P/'tests')==fixture_pins,'fixture changed')
need(sha256(P/'build/test_v83_export')==binary_pin,'test binary changed')
need({name:sha256(P/name) for name in harness_pins}==harness_pins,'harness changed')
need(sha256(P/"build/libpic_native_v8.a")==lib,"archive changed")
for name,digest in pins["files"].items(): need(sha256(P/name)==digest,"source changed")
(OUT/label/"result.json").write_text(json.dumps({"passed":True,"source_identity":pins["source_identity"],"archive_sha256":lib,"commands":r.records,"fixture_pins":fixture_pins,"binary_sha256":binary_pin,"harness_pins":harness_pins,"dependencies":inputs},indent=2)+"\n")
