#!/bin/sh
# ASan/UBSan runs (host mode, CPU only), one at a time on one leased core. usage: scripts/sanitized.sh OUTDIR
O=$1; mkdir -p $O
D=/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time-candidate008-auditfix
F=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
B=$D/build
E="env CUDA_VISIBLE_DEVICES= ASAN_OPTIONS=detect_leaks=1:abort_on_error=1:hard_rss_limit_mb=4096 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
run() { name=$1; shift; mkdir -p $O/$name; $D/scripts/lease.sh $E $F "$@" > $O/$name.log 2>&1; echo "$name exit $? $(tail -n 1 $O/$name.log)" >> $O/progress.txt; }
run san-v8-containers $B/sanitized-test_v8_containers
run san-v8-chemistry $B/sanitized-test_v8_chemistry $O/san-v8-chemistry reduced
for m in unbroken save restore tamper; do run san-v8-restore-$m $B/sanitized-test_v8_restore $O/san-v8-restore $m 3; done
run san-legacy $B/sanitized-test_legacy $D/inputs/h2o-gas-20261029.json $O/san-legacy
run san-ported $B/sanitized-test_ported $O/san-ported
run san-history $B/sanitized-test_history $O/san-history reduced
run san-flows $B/sanitized-test_flows $O/san-flows reduced
run san-forms-cpu $B/sanitized-test_forms_cpu $O/san-forms-cpu 100000
run san-forms-native $B/sanitized-test_forms_native $O/san-forms-native reduced
run san-work-gross $B/sanitized-test_work_gross $O/san-work-gross
