#!/bin/bash
# Build and run the JSON destructor reproducer inside the frozen launcher.
# Launch: /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait env CUDA_VISIBLE_DEVICES= repro/run_json_repro.sh LABEL
set -u
R=/home/treefor/projects/runtimes; P=$R/pic-water-common-commit-v4-20261030
O=/home/treefor/sim-workspace/powerflow/outputs/pic-water-common-commit-v4-20261030/repro-json-destructor-$1
J=$R/pic-reaction-transactions-20261029/native-full-time-candidate003/inputs
mkdir "$O" || exit 1
ulimit -v 3145728
cd "$P"
B=( /usr/bin/c++ -std=c++20 -O2 -isystem "$J" repro/repro_json_destructor.cpp -o "$O/repro_json_destructor" )
$R/arch-base-20261028/run-frozen "${B[@]}" > "$O/build.log" 2>&1; brc=$?
$R/arch-base-20261028/run-frozen "$O/repro_json_destructor" > "$O/run.log" 2>&1; rrc=$?
python3 - "$O" "$brc" "$rrc" <<'PY'
import hashlib,json,sys,os
from pathlib import Path
o=Path(sys.argv[1]);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
J='/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate003/inputs/raikou/json.hpp'
src='/home/treefor/projects/runtimes/pic-water-common-commit-v4-20261030/repro/repro_json_destructor.cpp'
r={'schema':'json-destructor-reproducer-receipt-v1','build_returncode':int(sys.argv[2]),'run_returncode':int(sys.argv[3]),
   'expected_run_returncode_when_defect_present':3,'cpu_affinity':sorted(os.sched_getaffinity(0)),'PF_CORE':os.getenv('PF_CORE'),
   'inputs_sha256':{src:sha(src),J:sha(J),'/home/treefor/projects/runtimes/arch-base-20261028/run-frozen':sha('/home/treefor/projects/runtimes/arch-base-20261028/run-frozen')},
   'binary_sha256':sha(o/'repro_json_destructor') if (o/'repro_json_destructor').exists() else None,
   'observation':json.loads((o/'run.log').read_text().strip().splitlines()[-1]) if (o/'run.log').read_text().strip() else None}
(o/'receipt.json').write_text(json.dumps(r,indent=1)+'\n');print(json.dumps(r['observation']))
PY
