#!/bin/bash
set -euo pipefail
p=/home/treefor/projects/runtimes/pic-collision-scheduler-indexed-20261029
o=/home/treefor/sim-workspace/powerflow/outputs/pic-collision-scheduler-indexed-20261029
exec flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock \
 taskset -c 7 nice -n 15 /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /bin/bash -c '
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=
ulimit -v 1048576
ulimit -t 60
p=$1
o=$2
python3 - "$p" "$o" <<"PY"
from pathlib import Path
import hashlib,json,sys,difflib,shutil
p,o=map(Path,sys.argv[1:]);sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for rec in json.loads((o/"input-bindings.json").read_text()):assert sha(Path(rec["snapshot"]))==rec["sha256"]
base=(p/"inputs/collision_scheduler.hpp").read_text();start=base.index("inline std::vector<Candidate> candidates(");end=base.index("\nclass Scheduler;",start)
new=base[:start]+(p/"source/indexed-builder.inc").read_text()+base[end:]
new=new.replace("#include <memory>","#include <memory>\n#include <map>").replace("#include <string>","#include <string>\n#include <string_view>")
assert new==(p/"source/collision_scheduler.hpp").read_text(),"unexplained source change"
shutil.copy2(p/"inputs/test_collision_scheduler.cpp",p/"source/test_existing_contract.cpp")
(o/"indexed-builder.patch").write_text("".join(difflib.unified_diff(base.splitlines(True),new.splitlines(True),fromfile="a/source/collision_scheduler.hpp",tofile="b/source/collision_scheduler.hpp")))
PY
g++ --version > "$o/compiler-version.txt"
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I"$p/source" "$p/source/test_existing_contract.cpp" -o "$p/build/test_existing_contract"
"$p/build/test_existing_contract" "$o/existing-contract-checks.json" > "$o/existing-contract-checks.txt"
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I"$p/source" -I"$p/inputs" "$p/source/test_indexed.cpp" -o "$p/build/test_indexed"
/usr/bin/time -v -o "$o/resources.txt" "$p/build/test_indexed" "$o/results.json" > "$o/checks.txt"
cat "$o/checks.txt"
' -- "$p" "$o"
