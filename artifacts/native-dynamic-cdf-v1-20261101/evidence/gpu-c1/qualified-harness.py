"""Bounded native CDF qualification; SimQueue/frozen userspace only."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tooling import *
mode = sys.argv[1]
need(mode in {"cpu", "gpu", "sanitized"}, "mode")
resource_guard(gpu=mode=="gpu")
if mode=="gpu": need(bool(os.environ.get("SIMQ_GPUS")), "scheduled GPU required")
label = sys.argv[2]
out=OUT/label
need(not out.exists(), "do not overwrite evidence")
r=Runner(out,"qualification")
inputs=sealed_inputs()
manifest=json.loads((P/"SOURCE-SHA256.json").read_text())
source={name:sha256(P/name) for name in manifest["files"]}
need(source==manifest["files"], "source manifest mismatch")
fixtures=tree_hashes(P/"tests")
adoption=tree_hashes(P/"adoption")
gas=Path("/home/treefor/projects/runtimes/gas-target-owner-v1.7.1-census-transport-20261101")
gas_pins=tree_hashes(gas/"include")
exe=P/"build"/(("sanitized-" if mode=="sanitized" else "")+"test_native_cdf")
archive=P/"build"/(("sanitized-" if mode=="sanitized" else "")+"libpic_native_v8.a")
pins={str(x):sha256(x) for x in [exe,archive,P/"qualify-native-cdf.py",P/"tooling.py"]}
args=[exe,out/"native"]
if mode=="gpu":args += [NVRTC,GPU_UUID]
r.run(args,"native.log",timeout=240,env={"PF_V8_CHEM":"v3.1"},sanitized=mode=="sanitized",limit=mode!="gpu")
if mode=="cpu":
    inc=["-I"+str(P/"forms/include/pf_v7"),"-I"+str(P/"adoption"),"-I"+str(P/"source/include"),"-I"+str(P/"source/src"),"-isystem",str(CHEM/"include/pf_chem_v3"),"-isystem",str(gas/"include"),"-isystem",str(P/"inputs"),"-I"+str(VENDOR)]+SYSTEM
    r.run(["/usr/bin/c++","-std=c++20","-O1","-ffp-contract=off","-fno-fast-math","-Wall","-Wextra","-Wpedantic","-Werror",*inc,P/"adoption/test_kind32.cpp","-o",out/"test_kind32"],"adapter-compile.log",timeout=60)
    r.run([out/"test_kind32"],"adapter.log",timeout=30)
need(inputs==sealed_inputs(),"dependencies changed")
need(adoption==tree_hashes(P/"adoption") and gas_pins==tree_hashes(gas/"include"),"adapter inputs changed")
need(source=={name:sha256(P/name) for name in source},"source changed")
need(fixtures==tree_hashes(P/"tests"),"fixtures changed")
need(pins=={name:sha256(name) for name in pins},"binary/harness changed")
result=json.loads((out/"native/result.json").read_text())
need(result["passed"] and result["source"]==manifest["source_identity"],"actual result identity")
(out/"qualification.json").write_text(json.dumps({"passed":True,"mode":mode,"source_identity":manifest["source_identity"],"result":result,"inputs":inputs,"fixture_pins":fixtures,"adapter_pins":adoption,"gas_header_pins":gas_pins,"pins":pins,"commands":r.records,"production_qualified":False,"network_complete":False},indent=2)+"\n")
print("native CDF qualification passed",mode)
