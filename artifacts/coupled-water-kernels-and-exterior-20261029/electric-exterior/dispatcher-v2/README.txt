Version 2: additive-evidence followup dispatcher, not executed here

The original exact-output uniform control and version 1 tools remain failed
and unchanged. uniform9-additive-review-v 2/report.txt and disposition.json
separately document the CPU causal proof, source/CSR/RHS equality, unchanged
original guards and observed scale. Version 2 cannot launch uniform 9. It
requires a new explicit approval of that evidence and explicit case ids.

Read-only preparation verification passed 93 pinned records and exact frozen
PID 278902/start 5288339/executable/cgroup identity. Four bounded CPU approval/
corruption tests passed. No new GPU dispatch or production action occurred.

After inspecting the additive evidence and version 2 source, copy
approval.template.json to a fresh approval file. Set both approval booleans
true only when explicitly authorized, and keep only approved case ids.
The shipped template has both booleans false and authorizes nothing.

Execution syntax after that separate approval:

 taskset -c 7 nice -n 15 env PYTHONDONTWRITEBYTECODE=1 python 3 \
 /home/treefor/sim-workspace/powerflow/outputs/pic-electric-exterior-20261029/dispatcher-v2/dispatch.py \
 approved-followup --approval ABSOLUTE_NEW_APPROVAL_FILE \
 --case nested-105-40 \
 --receipt /home/treefor/sim-workspace/powerflow/outputs/pic-electric-exterior-20261029/parent-followup-v2-40.json

Repeat --case only for the explicitly selected cases. No stage or followup
launches automatically. GPU→core 7 locks are acquired internally; do not wrap
this command in an outer core 7 flock. All original candidate/fixture/library/
tool/hold pins and all additive evidence are verified before/after execution.
The version 2 manifest is also required to remain identical during the process.

Caps remain unchanged: CPU 7,nice 15,one thread,120 s timeout,4 GiB sampled owned
RSS,10 GiB sampled global GPU memory,6 GiB sampled GPU increase,64 MiB output.
Exact candidate exe/complete argv/output/start identity and an owned-ancestry
PIDFD are required. Only the newly owned candidate/timeout may be signaled.
Each case requires original source/charge/work and FP 64 residual guards,
default 2 IR,one factorization,two solves,zero advance and byte-identical B.

CPU field comparison uses the original unchanged tools/compare_fields.py
and field-comparison-plan.json. Results must remain labeled frozen-charge
static sensitivity. No finite pair establishes infinite-domain convergence.
