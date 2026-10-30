#!/bin/bash
set -euo pipefail
observer_prefix=/home/treefor/projects/runtimes/pic-rhs-observer-20261029
observer_frozen=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
"$observer_frozen" /usr/bin/cmake -S "$observer_prefix/source" -B "$observer_prefix/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem \
    > "$observer_prefix/logs/configure.log" 2>&1
"$observer_frozen" /usr/bin/cmake --build "$observer_prefix/build" \
    --target heating_history_pic_gpu --parallel 1 > "$observer_prefix/logs/build.log" 2>&1
"$observer_frozen" /usr/bin/c++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
    "$observer_prefix/tests/observer_contract.cpp" -o "$observer_prefix/tests/observer_contract" \
    > "$observer_prefix/logs/contract-build.log" 2>&1
observer_test_root=$(mktemp -d "$observer_prefix/evidence/contracts-XXXXXX")
"$observer_frozen" /usr/bin/prlimit --cpu=10 --as=268435456 -- \
    "$observer_prefix/tests/observer_contract" "$observer_test_root" \
    > "$observer_prefix/logs/contract-test.log" 2>&1
python "$observer_prefix/validate_capture.py" "$observer_test_root/normal" \
    --allow-contract-fixture --write-receipt > "$observer_prefix/logs/contract-validation.log" 2>&1
"$observer_frozen" /usr/bin/cmake -S "$observer_prefix/tests" -B "$observer_prefix/test-build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DMFEM_DIR=/home/treefor/projects/runtimes/mfem-4.10-pic-20261029/install-mpi/lib/cmake/mfem \
    > "$observer_prefix/logs/native-contract-configure.log" 2>&1
"$observer_frozen" /usr/bin/cmake --build "$observer_prefix/test-build" --parallel 1 \
    > "$observer_prefix/logs/native-contract-build.log" 2>&1
"$observer_frozen" /usr/bin/prlimit --cpu=10 --as=536870912 -- \
    "$observer_prefix/test-build/native_form_contract" "$observer_test_root/native-form" \
    > "$observer_prefix/logs/native-contract-test.log" 2>&1
python "$observer_prefix/validate_capture.py" "$observer_test_root/native-form" \
    --allow-contract-fixture --write-receipt > "$observer_prefix/logs/native-contract-validation.log" 2>&1
printf '%s\n' "$observer_test_root" > "$observer_prefix/evidence/latest-contract-root.txt"
cat "$observer_prefix/logs/contract-test.log"
cat "$observer_prefix/logs/native-contract-test.log"
printf '%s\n' 'Diagnostic candidate and CPU-only capture validation passed.'
