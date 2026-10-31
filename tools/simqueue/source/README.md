# simqueue

## 0.6: node scheduler

With a node configuration (`SIMQUEUE_CONFIG`, else `<state dir>/node.json`,
schema `simqueue-node-config-v1`) SimQueue 0.6 schedules a whole node:

- **Partitions** name physical cores (by their lowest logical CPU), the GPUs
  (by UUID) they may use, default and maximum memory, default and maximum time
  limit, and default core count. `default_partition` receives jobs submitted
  without `-p`.
- **Concurrent jobs with resource requests**:
  `simq submit -p PARTITION [-c N | --cpus LIST] [--gpus N] [--gpu-mem SIZE] [--mem SIZE] [-t TIME] [--priority N] -- cmd`.
  One worker starts every job that fits, in priority/queue order, with EASY
  backfill: once a job is blocked, a later job starts early only if, given all
  time limits, it cannot delay that job. GPUs are counted resources with a
  memory budget: a job is exclusive on its GPU unless `--gpu-mem` is given and
  the claims of all sharing jobs fit the budget. Memory is counted against
  `memory_schedulable`. `simq submit --test-only` shows where and whether a
  request would start now.
- **Jobs are independent of the worker.** Each job runs as its own transient
  systemd user service (`<prefix>-<id>-<token>.service`) with `CPUAffinity`,
  `MemoryMax`, `MemorySwapMax`, `OOMPolicy=kill`, `RuntimeMaxSec`, working
  directory and placement environment on the unit. The unit's main process is
  `simqueue/jobwrap.py`, which imports nothing from the scheduler, holds the
  job's lease files and records the exit status; `ExecStopPost` records
  systemd's result. The worker can be stopped, crash, or be replaced by
  another build: on start it re-adopts running units by name and launch
  token, records jobs that ended meanwhile, and never starts a duplicate.
- **Placement**: one logical CPU per allocated physical core (SMT siblings
  unused), `CUDA_VISIBLE_DEVICES` = allocated GPU UUIDs (empty if none), and
  `SIMQ_JOB_ID`, `SIMQ_LAUNCH_TOKEN`, `SIMQ_PARTITION`, `SIMQ_NODE`,
  `SIMQ_CPU_LIST`, `SIMQ_CPU_COUNT`, `SIMQ_NCPUS`, `SIMQ_GPUS`,
  `SIMQ_GPU_MODE`, `SIMQ_GPU_MEM`, `SIMQ_MEM`, `SIMQ_TIME_LIMIT`, `SIMQ_LOG`,
  `SIMQ_UNIT`. `SIMQ_CPU_LIST`/`SIMQ_CPU_COUNT` keep their 0.5 meaning. No
  `SLURM_*` variable is exported, and inherited ones are removed unless given
  with `-e`, because launchers switch to Slurm code paths on them.
- **Lease files**: with `lock_files.enabled`, a job unit takes the flock(2)
  lease files of its cores and GPUs (the files used by `any-core.sh` and
  pf-launch) for the job's lifetime, and cores or GPUs whose lease is held
  outside the queue are treated as busy (`simq partitions` names the holder).
- **Compatibility**: every 0.5 command and option works; `-c` keeps its
  meaning (cores of the partition, default all). Schema 6 databases migrate
  in place to schema 7 (new columns only). Idempotency keys replayed with 0.5
  syntax keep their 0.5 fingerprint. Ensembles run inside their own
  allocation (`simq ensemble submit -p PARTITION`; a GPU budget becomes a
  shared GPU claim of that size).
- New commands: `simq partitions` (alias `sinfo`), `simq priority ID N`,
  `simq worker --once`.

Without a node configuration SimQueue runs exactly as 0.5 (serial worker,
`SIMQUEUE_PCORES`); the rest of this README describes that mode and the
features shared by both. Operations, the zenobia configuration, the switchover
from 0.5 and the mapping to Slurm are in `../OPERATIONS.md`.

Tests: `./run-tests.py --cores LIST` takes the lease files of the cores it uses,
pins itself there at nice 15 and runs real job units only on those cores.

## 0.5 reference

`simq` is a persistent job queue for simulations under `~/projects`. Outer
queue jobs are FIFO by default and remain serial; `promote` and `demote` allow
an explicit queued-job reorder. A normal submission is one process tree;
an ensemble submission is one outer allocation whose supervisor can run
several independent, one-core simulation members concurrently. It is tuned
for this workstation's eight P-cores (Linux CPU IDs `0-7`). In the
recommended installation, the queue dynamically reserves
`0-7` as an exclusive, load-balanced cgroup-v2 partition. Ordinary
cgroup-managed userspace outside the queue can use only the E-cores
(`8-23`) until the queued batch finishes; then the P-cores are returned to
the rest of the machine.

The command inside the allocation still decides whether to use MPI ranks,
OpenMP threads, CUDA, or a mixture.

Normal queue jobs deliberately do **not** set `OMP_NUM_THREADS`, MPI rank
counts, or BLAS thread counts. Ensemble lanes are the exception: each member
starts with singleton CPU affinity and one-thread settings for OpenMP,
Kokkos, and common BLAS runtimes. It also disables automatic rank binding in
Open MPI/PRRTE, Intel MPI, and MVAPICH so a rank-1 launcher inherits its lane.
The member command still controls MPI and should launch one rank for this
workflow. Explicit launcher options can reset those settings, so they must
honor the lane contract; process diagnostics report any observed affinity
widening.

## Install

```bash
cd ~/projects/simqueue
sudo ./install-exclusive.sh
source ~/.profile  # only needed if simq is not already on PATH
simq doctor
```

Security: run the sudo installer only from this trusted checkout. It promotes
the helper, systemd units, and Polkit rule staged from this tree into root-owned
system paths, so do not run it from an untrusted or unexpectedly modified copy.

Wait for any currently active queue job to finish before installing. The
installer refuses to migrate an active job. It replaces the earlier user
service with a system-managed runner that still executes as your account,
installs a small root-owned shield helper, and grants your local account
permission to start only `simqueue-runner.service`. It does not grant broad
systemd or root access.

Code stays in this directory. Mutable state remains private to your account
and lives outside every science repository:

```text
~/.local/state/simqueue/queue.sqlite3
~/.local/state/simqueue/logs/000001.log
~/.local/state/simqueue/ensembles/000001/
~/.local/state/simqueue/worker.log
```

`./install.sh` remains available as a no-root, cooperative-affinity setup for
development or another machine. It pins queue jobs to P-cores but cannot keep
other work off them; it does not satisfy exclusive mode.

## Everyday use

When asking local to launch an actual Charon, Dakota, Xyce, Zapdos, or
Multibolt simulation anywhere under `~/projects`, append `--simqueue` to the
request:

```text
Run the FZ750 Charon case --simqueue
```

Local treats that token as a workspace directive, removes it from the real
solver command, and submits the run here. It does not apply the directive to
builds, meshing, preprocessing, postprocessing, or tests merely because they
live in one of those repositories. `--simqueue` is not a literal option
understood by the solvers; for direct terminal use, use `simq submit` as shown
below.

Submit from the directory in which the simulation should run:

```bash
cd ~/projects/charon-xyce-dakota
simq submit -n aging -- ./run_aging_pcores.sh
```

Or name the working directory explicitly:

```bash
simq submit \
  -n pore-breakdown \
  -C ~/projects/zapdos-multibolt \
  -e P_CPUS=0-7 \
  -- ./run.sh
```

The `P_CPUS=0-7` override is currently important for
`zapdos-multibolt/run.sh`: that older launcher otherwise applies a stale
WSL-era CPU mask of `0,1,10,11,12,13,22,23` inside the queue's outer
`taskset` command.

Common commands:

```bash
simqueue -active             # only jobs running right now
simqueue -queue              # only jobs waiting to run
simq                         # active jobs and recent history
simq list --all              # longer history
simq show 12                 # exact command, directory, timestamps, result
simq log -f 12               # follow combined stdout/stderr
simq cancel 12               # cancel queued or running work
simq retry 12                # enqueue a new attempt; preserve old history
simq wait 12                 # wait; Ctrl-C does not cancel the job
simq promote 12              # move a queued outer job to the front
simq demote 12               # move a queued outer job to the back
simq pause                   # global: let the current outer job finish, dispatch nothing else
simq resume                  # global: resume outer-job dispatch
```

`simq promote ID` and `simq demote ID` change only the order of queued outer
jobs. They do not preempt or otherwise affect an active job, and they preserve
the selected job's ID, attempt, members, artifacts, and history. Reordering is
also allowed while the global queue is paused, but top-level `simq pause` still
controls dispatch: no queued job starts until top-level `simq resume`. A paused
ensemble participates in the reordered queue like any other outer job. When it
reaches the front it blocks every job behind it until it is resumed or reordered
again.

The submitting shell's environment is captured for the job. Add or replace
individual variables with repeatable `-e KEY=VALUE` options. The queue also
provides these variables to the running process:

```text
SIMQ_JOB_ID
SIMQ_CPU_LIST
SIMQ_CPU_COUNT
SIMQ_LOG
```

Controllers that must recover safely across a crash between queue submission
and their own ledger update should use an idempotency key and structured
output:

```bash
simq submit \
  --idempotency-key msc080-eq:v1:CAMPAIGN:reconciler:0 \
  --json -n reconcile -C ~/projects/my-campaign -- ./reconcile.py
```

The first call creates the job and returns `"created": true`. Replaying the
same key with the same command, directory, core count, and effective
environment returns the original job ID with `"created": false`, regardless
of whether that job is queued, active, or finished. Reusing a key with a
changed request fails closed. Keys are global to this queue database, so bind
them to a durable logical operation such as campaign/chunk or campaign/round;
do not use a display name as an identity. `simq retry ID` remains an explicit
new attempt and does not inherit the original key.

Use explicit shell mode only when shell syntax is needed:

```bash
simq submit \
  -n activated-run \
  -C ~/projects/some-project \
  --shell 'source .venv/bin/activate && exec ./run.sh > extra-output.txt'
```

Normal `-- COMMAND ARG...` submissions do not pass through a shell, so
arguments are preserved literally.

## Independent simulation ensembles

Use an ensemble for a parameter sweep or training-data batch whose members
are separate simulations. SimQueue still claims only one outer FIFO job, then
fills up to `max_parallel` one-core lanes inside that allocation. For MPI
Charon, each member should launch its own rank-1 process (for example,
`mpiexec -np 1 ...`); do not turn eight independent cases into one
`mpiexec -np 8` job. Keep each rank at one host thread and do not let its
launcher widen `SIMQ_CPU_LIST`/process affinity.

For common launchers, SimQueue forces `PRTE_MCA_bindto=none` (Open MPI 5),
`OMPI_MCA_hwloc_base_binding_policy=none` (older Open MPI), `I_MPI_PIN=0`,
and `MV2_ENABLE_AFFINITY=0`, then the outer `taskset` supplies the singleton
lane. Do not add a conflicting command-line binding such as
`--bind-to core --map-by core:PE=8`: command-line launcher policy can override
these variables. Affinity samples make that mistake visible.

An ensemble manifest is JSON with schema `simqueue-ensemble-v1`:

```json
{
  "schema_version": "simqueue-ensemble-v1",
  "name": "charon-training-sweep",
  "max_parallel": 8,
  "gpu": {
    "device": "0",
    "budget": "8GiB",
    "mps": "off"
  },
  "defaults": {
    "cwd": "/home/tjsmith/projects/my-charon-sweep",
    "env": {
      "OMP_NUM_THREADS": "1",
      "OPENBLAS_NUM_THREADS": "1",
      "MSIC_KOKKOS_CUDAUVM_BACKING": "managed-pool",
      "MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES": "536870912"
    },
    "gpu_memory": "768MiB"
  },
  "members": [
    {
      "name": "sample-0000",
      "cwd": "/home/tjsmith/projects/my-charon-sweep/sample-0000",
      "argv": ["/usr/bin/mpiexec", "-np", "1", "/absolute/path/to/charon", "-i", "input.xml"]
    },
    {
      "name": "sample-0001",
      "cwd": "/home/tjsmith/projects/my-charon-sweep/sample-0001",
      "argv": ["/usr/bin/mpiexec", "-np", "1", "/absolute/path/to/charon", "-i", "input.xml"],
      "env": {"TRAINING_SAMPLE": "1"},
      "gpu_memory": "896MiB"
    }
  ]
}
```

Every member needs a unique name and an `argv` array; arguments are executed
literally, without shell interpretation. `cwd`, `env`, and `gpu_memory` may
come from `defaults` and be overridden by a member. Working directories must
remain under `~/projects`, and independent members should not write into the
same case directory.

`gpu_memory` is an admission claim, not an allocation: declare the measured
working set plus headroom for that member. SimQueue launches a member only
when both a CPU lane is free and the sum of active claims stays within
`gpu.budget`. A member's `MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES` arena, when set,
cannot exceed its claim. The claim should also cover CUDA context and
off-pool allocations. Do not treat UVM oversubscription as a swap strategy;
repeated migration can erase the benefit of concurrency.

Submit and inspect an ensemble with:

```bash
simq ensemble validate sweep.json
simq ensemble submit sweep.json
simq ensemble submit sweep.json --idempotency-key CAMPAIGN:chunk:000 --json
simq ensemble status JOB_ID
simq ensemble log JOB_ID sample-0000
simq ensemble log -f JOB_ID 0       # exact name or zero-based index
```

Gracefully drain and later continue one ensemble without discarding completed
members:

```bash
simq ensemble pause JOB_ID
simq ensemble status JOB_ID
simq ensemble resume JOB_ID
```

`simq ensemble pause` is the per-ensemble drain operation. It atomically closes
future member-launch admission; members that have already been admitted are not
signaled and may finish normally, but their freed lanes are not backfilled. Once
the last admitted member and its descendants exit, SimQueue records the drain,
releases the outer CPU allocation, and returns the same parent job ID to the
queue. Completed member states and logs are preserved, while members that never
started remain queued. `simq ensemble resume` reopens launch admission. If the
drain is still in progress, the active supervisor resumes backfill; if the
allocation was already released, the protected runner is woken and continues
the remaining members under the same job ID. Repeated pause or resume requests
are idempotent. Add `--json` to `pause`, `resume`, or `status` for structured
control data. Status reports `running`, `draining`, `paused`, `resuming`, or a
terminal control state, together with whether the allocation is held or
released.

A paused ensemble retains its explicit queue position. If it is the first
queued job, it also holds every outer job behind it until that ensemble is
resumed or reordered; dispatch does not skip around it. The top-level
`simq pause` and `simq resume` commands are a
separate global control: global pause prevents another *outer* job from being
claimed, but does not stop member backfill inside an ensemble that is already
running. Use `simq ensemble pause JOB_ID` when the active ensemble itself must
drain. The controls remain independent: resuming an ensemble while the global
queue is paused leaves it waiting for top-level `simq resume`, and top-level
resume does not clear a per-ensemble pause.

The default is keep-going: a failed member is recorded, another eligible
member can fill its lane, and the outer job ultimately fails if any member
failed. A manifest with `"fail_fast": true` stops launching pending work and
terminates running members after the first failure. `simq cancel JOB_ID`
terminates the ensemble supervisor and all tracked member process groups.

Fence-heavy CUDA setup can be serialized without serializing the whole
ensemble. Add an optional rolling launch gate whose marker means that the
current member has completed and fenced its setup phase:

```json
"launch_gate": {
  "kind": "member-log-substring",
  "contains": "MSIC_CHARON_PHASE_V1 event=begin phase=solver_eval rank=0",
  "timeout_seconds": 300
}
```

With this gate, SimQueue permits at most one pre-marker member on the GPU.
Once that member emits the exact substring, one follower may start setup while
the ready member continues solving; the follower then becomes the gate holder.
This is a rolling setup gate, not a one-time delay and not a global lock on
the already-ready simulations. The marker must be flushed to the member log
only after setup is genuinely complete and all required device work is
fenced. SimQueue scans output written after its own log header, so marker text
inside the stored command or environment cannot open the gate.

The gate fails closed. A timeout, log truncation/read error, or member exit
before its marker stops follower admission, terminates the cohort, and fails
the outer job. Gate waits, releases, and failures are recorded in
`events.jsonl`; `simq ensemble status` shows the immutable gate configuration.
Retries preserve the gate configuration but use fresh logs and runtime state.

`simq ensemble status` reports queue/member state, assigned CPU, PID, elapsed
time, exit result, GPU claim, parsed allocator high-water data, and the exact
artifact paths. The durable artifacts are:

```text
~/.local/state/simqueue/logs/NNNNNN.log
~/.local/state/simqueue/ensembles/NNNNNN/events.jsonl
~/.local/state/simqueue/ensembles/NNNNNN/gpu_samples.jsonl
~/.local/state/simqueue/ensembles/NNNNNN/process_samples.jsonl
~/.local/state/simqueue/ensembles/NNNNNN/members/NNNN-name.log
```

The outer log and event stream record admission waits, launches, exits, and
the final cohort summary. GPU samples capture device memory, utilization,
power, temperature, and per-process residency when `nvidia-smi` exposes it,
while process samples make lane CPU affinity, RSS, thread count, and runtime
visible. Each member log starts with its resolved command, directory, lane,
claim, and relevant CUDA/threading environment. On exit, SimQueue parses the
`MSIC_KOKKOS_CUDAUVM_BACKING_V1 event=finalize` marker and exposes peak live
and reserved bytes, reclaim count, teardown status, and release mode. Missing
GPU tooling or markers are reported as diagnostics rather than invented
measurements.

`gpu.mps` currently accepts only `off`, and every member receives a private,
nonexistent MPS pipe path so an inherited controller cannot capture it.
SimQueue never starts or stops MPS and never changes compute mode, clocks, or
other GPU state. This is deliberate: [NVIDIA requires controller-mediated
client teardown](https://docs.nvidia.com/deploy/mps/latest/when-to-use-mps.html#client-early-termination)
before an MPS process may safely receive cancellation signals, and the queue's
crash-recovery path cannot yet guarantee that protocol. The
non-MPS CUDA scheduler time-slices independent contexts; use the utilization
samples to decide whether a future, fully managed MPS implementation is worth
adding.

Start a real GPU cohort at two members, measure allocator and GPU high-water
behavior, then try four and eight. In particular, do not reuse a 12 GiB
per-process UVM arena eight times on a 16 GiB device. Set `gpu.budget` below
practical free VRAM and leave room for the display, CUDA contexts, libraries,
and non-UVM allocations. Reduce `max_parallel` or member claims if sampling
shows migration, reclaim churn, or sustained memory pressure.

The existing MSC launcher at
`charon-xyce-dakota/charon-inference/MSC080SMA330D/scripts/run_charon_cuda_msc.sh`
is not an ensemble-lane command: it currently selects all eight P-cores,
requests eight host threads, and loads a 12 GiB managed-pool profile. Do not
place that wrapper unchanged into eight members. Use a rank-1, one-host-thread
launcher whose measured arena and total GPU claim have been qualified for the
exact deck first.

This harmless CPU-only manifest exercises lane scheduling and logs without
touching the GPU; a `0B` claim does not consume the GPU budget:

```json
{
  "schema_version": "simqueue-ensemble-v1",
  "name": "ensemble-dry-run",
  "max_parallel": 2,
  "allow_shared_cwd": true,
  "defaults": {
    "cwd": "/home/tjsmith/projects",
    "gpu_memory": "0B"
  },
  "members": [
    {"name": "dry-0", "argv": ["/usr/bin/python3", "-c", "import time; print('dry-0'); time.sleep(1)"]},
    {"name": "dry-1", "argv": ["/usr/bin/python3", "-c", "import time; print('dry-1'); time.sleep(1)"]},
    {"name": "dry-2", "argv": ["/usr/bin/python3", "-c", "import time; print('dry-2'); time.sleep(1)"]}
  ]
}
```

Ensembles inherit the runner's isolation mode; they do not upgrade it. Check
`simq doctor` before interpreting timing. In `cooperative-affinity` mode (the
current local installation when this feature was added), lanes are pinned to
P-cores but unrelated work can still use those cores. Only the exclusive
cgroup-v2 installation shields them from ordinary outside workloads.

## CPU and lifecycle behavior

- Outer jobs are FIFO unless explicitly promoted or demoted, and are strictly
  serial even when `--cores N` selects fewer than all eight P-cores. The default
  is all eight (`0-7`). An ensemble's concurrent members remain inside its one
  outer allocation.
- The exclusive system runner promotes `simqueue.slice` to a `root` cpuset
  partition before claiming work. This makes `0-7` unavailable to ordinary
  userspace outside the queue while retaining scheduler load balancing across
  all eight P-cores. It keeps the shield for a contiguous batch and releases
  it as soon as the queue drains or pauses.
- A normal job and its child processes receive one combined log. An ensemble
  adds the member logs and diagnostic streams listed above. A nonzero result
  marks the outer job failed and the next queued outer job still runs.
- Cancellation sends `SIGTERM` across the job's tracked process groups,
  including adopted descendants that cleared their environment or created new
  sessions. Exclusive mode also verifies the complete runner cgroup. It
  rescans during the ten-second grace period, then uses `SIGKILL` if needed;
  the next serial job is not released until those processes are gone.
- If the worker or machine stops during a job, that job is marked
  `interrupted` on restart. It is never automatically duplicated; use
  `simq retry ID` after checking its outputs.
- Working directories must remain under `~/projects`; executables may live
  elsewhere (for example `/usr/bin/python3` or a Conda environment).
- In exclusive mode, nested launchers cannot widen affinity beyond the
  queue's P-core cpuset. They can still narrow themselves accidentally, so
  correct stale inner CPU masks with `-e`, as in the Zapdos example above.
- The worker refuses to launch jobs if the protected system service, helper,
  cgroup placement, or effective exclusive CPU set cannot be verified. A
  fresh shield check is made before every job in a batch. A failed release
  stops the runner rather than dispatching more work; if job processes survive
  shutdown, the helper refuses to remove their reservation.

The exclusive runner is a system service, so queued work survives terminal
closure and full logout without enabling user linger. It starts at boot to
recover or drain pending work, but it exits and releases the P-core shield
when the queue is empty.

The exclusivity guarantee covers normal processes and services managed by
cgroup v2. It cannot eliminate hard IRQs, firmware activity, or mandatory
per-CPU kernel threads. Those may still execute briefly on a P-core. Moving
that remaining kernel-level noise would require separate IRQ and boot-time
isolation tuning, which this dynamic queue intentionally does not change.

## Development and recovery

Run the test suite without installing anything:

```bash
cd ~/projects/simqueue
python3 -m unittest discover -s tests -v
```

Before exclusive mode is installed, foreground drain mode is useful for
diagnosis or a machine without a systemd manager:

```bash
SIMQUEUE_NO_AUTOSTART=1 simq submit -- ./run.sh
simq worker --drain
```

After exclusive installation, foreground/detached workers are deliberately
rejected because they would bypass the shield. `simq submit`, `simq retry`,
and `simq resume` wake the fixed protected runner instead.

Inspect service diagnostics with `simq doctor`, and detailed worker events
in `~/.local/state/simqueue/worker.log`.
