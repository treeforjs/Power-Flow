#!/bin/bash
set -euo pipefail
p=/home/treefor/projects/runtimes/pic-water-event-bridge-20261029
o=/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-water-event-bridge-20261029
mkdir -p "$o"
exec flock /home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock \
 taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /bin/bash -c '
  set -euo pipefail
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ulimit -v 1048576
  ulimit -t 60
  python3 - "$1" "$2" <<"PY"
from pathlib import Path
import json,hashlib,sys
p=Path(sys.argv[1])
for name,want in json.loads((p/"input-sha256.json").read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==want,name
names=["inputs/collision_scheduler.hpp","inputs/include/raikou/collision_data.h","inputs/include/raikou/h2o_conditional_parent.h","source/water_parent_bridge.hpp"]
binding={name:hashlib.sha256((p/name).read_bytes()).hexdigest() for name in names}
canonical=json.dumps(binding,sort_keys=True,separators=(",",":"))
digest=hashlib.sha256(canonical.encode()).hexdigest()
out=Path(sys.argv[2]);(out/"kernel-binding.json").write_text(json.dumps({"sources":binding,"canonical_manifest_sha256":digest},indent=2)+"\n")
(out/"kernel.sha256").write_text(digest+"\n")
PY
  g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I"$1/inputs" -I"$1/inputs/include" "$1/source/test_water_parent_bridge.cpp" -o "$1/build/test_water_parent_bridge"
  read -r kernel_sha < "$2/kernel.sha256"
  read -r data_sha _ < <(sha256sum "$1/inputs/h2o-gas.json")
  "$1/build/test_water_parent_bridge" "$1/inputs/h2o-gas.json" "$kernel_sha" "$data_sha" "$2/checks.json"
 ' -- "$p" "$o"
