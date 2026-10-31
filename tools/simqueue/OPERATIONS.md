# SimQueue 0.6.1 on zenobia: operations

SimQueue 0.6 is a single-node, Slurm-shaped scheduler. One worker
(`simqueue.service`) schedules the whole node; every job runs in its own
transient systemd user service, independent of the worker.

0.6.1 adds **opportunistic borrowing**: restartable work submitted with
`simq submit --borrow` (builds, tests) may run on idle cores of any partition
that lends them, at the lowest CPU priority, and is frozen, requeued or
canceled the moment an owner needs those cores. Owners never wait for a
borrower. With no `--borrow` job in the queue, 0.6.1 schedules exactly as
0.6.0 (see "Borrowing idle cores" below).

| Item | Path |
|---|---|
| Runtime (sealed, read-only once installed) | `/home/treefor/projects/runtimes/simqueue-0.6.1-20261031` |
| Previous runtime (kept: running 0.6.0 job units execute its `jobwrap.py`) | `/home/treefor/projects/runtimes/simqueue-0.6.0-20261030` |
| Interpreter (unchanged from 0.5) | `/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen` |
| Code | `source/simqueue/` (`nodeworker.py`, `scheduler.py`, `borrow.py`, `telemetry.py`, `units.py`, `jobwrap.py`, `locks.py`, `nodeconfig.py`, plus the 0.5 modules) |
| Node configuration (template) | `config/zenobia-node.json` |
| Node configuration (live) | `~/.local/state/simqueue/node.json` (`SIMQUEUE_CONFIG`) |
| Database | `~/.local/state/simqueue/queue.sqlite3` (schema 8; 0.6.0 used schema 7) |
| Borrow slice | `~/.config/systemd/user/simqborrow.slice` (`CPUWeight=idle`) |
| Per-job files | `~/.local/state/simqueue/jobs/NNNNNN-<token>/{spec,status,stop-post}.json` |
| Logs | `~/.local/state/simqueue/logs/NNNNNN.log`, `~/.local/state/simqueue/worker.log` |
| Upgrade 0.6.0 to 0.6.1 / rollback | `upgrade.py` (dry run unless `--execute`) |
| Switchover from 0.5 (historical) | `switchover.py` |
| Tests | `source/run-tests.py` (takes the cores' lease files first; stops when other work needs those cores) |

## Partitions on zenobia

| Partition | Cores (logical CPU per core) | GPU | Default / max memory | Default / max time | Default cores | Lending (0.6.1) |
|---|---|---|---|---|---|---|
| `spheral` (default) | 0-3 | none | 18G / 20G | 2-00:00:00 / 14-00:00:00 | all 4 | lends free and idle-in-job cores |
| `raikou-mfem` | 4-5 | RTX 3080 Ti | 4G / 12G | 12:00:00 / 7-00:00:00 | both | lends free and idle-in-job cores |
| `powerflow` | 6-7 | RTX 3080 Ti | 4G / 12G | 12:00:00 / 7-00:00:00 | both | lends free and idle-in-job cores; `--exclusive-timing` jobs are never lent |

SMT siblings 8-15 are never allocated. Swap per job: 1G. Memory is a counted
resource: the sum of `--mem` over running jobs stays within 28G
(`memory_schedulable`). The GPU has an 11G admission budget: a GPU job is
exclusive unless it passes `--gpu-mem`, in which case jobs share the GPU while
their claims fit the budget. Time limits are enforced by systemd
(`RuntimeMaxSec`); an expired job ends `failed` with result `timeout`.

## Using the queue

```sh
simq submit -c 4 -- ./run.sh                          # 0.5 syntax: default partition (spheral), 4 cores
simq submit -p powerflow -c 1 -t 2:00:00 --mem 2G -- ./case.sh
simq submit -p powerflow --cpus 7 -- ./pinned.sh      # explicit core (partition CPU number)
simq submit -p raikou-mfem --gpus 1 -- ./gpu.sh       # exclusive GPU
simq submit -p powerflow --gpus 1 --gpu-mem 3G -- ./x # shared GPU, 3G claim
simq submit --test-only -p raikou-mfem --gpus 1 -- /bin/true   # where/when it would run; queues nothing
simq                     # active and recent jobs (PARTITION, CPU, GPU, reason for waiting)
simq partitions          # free cores, GPU claims, lease files held outside the queue (alias: sinfo)
simq show ID             # request, allocation, unit, result, job files
simq log -f ID
simq cancel ID           # works with or without a running worker
simq priority ID N       # higher runs first; promote/demote still move to the ends
simq ensemble submit -p powerflow sweep.json          # an ensemble runs inside its own allocation
simq submit --borrow -p powerflow --max-borrow-cores 4 -- make -j"$(nproc)"   # 0.6.1: borrow idle cores
simq submit -p powerflow --exclusive-timing -- ./timing-run.sh               # 0.6.1: never lend its cores
```

Options: `-p/--partition`, `-c/--cores/--cpus-per-task N`, `--cpus LIST`,
`--gpus N`, `--gpu-mem SIZE`, `--mem SIZE` (`4G`; bare number = MiB, as
`sbatch --mem`), `-t/--time` (`MM`, `MM:SS`, `HH:MM:SS`, `D-HH[:MM[:SS]]`),
`--priority N`, `-n/-J/--name`, plus every 0.5 option (`-C`, `-e`, `--shell`,
`--idempotency-key`, `--json`). Without `-p` a job goes to the default
partition; without `-c` it gets the partition's default (all of its cores, as
in 0.5). The 0.5 options are accepted unchanged and an idempotency key
replayed with 0.5 syntax keeps its 0.5 fingerprint.

Scheduling: waiting jobs are ordered by priority, then queue order. Each worker
pass starts every job that fits, except that once a job is blocked, a later
job starts only if (given all time limits) it cannot delay that blocked job
(EASY backfill). A paused ensemble holds its position and blocks later jobs of
its own partition only. Borrow jobs (0.6.1) always sort after every owner job,
whatever their priority, and never take part in backfill reservations.

### Environment of every job

| Variable | Meaning |
|---|---|
| `SIMQ_JOB_ID` | job ID (as in 0.5) |
| `SIMQ_LAUNCH_TOKEN` | random per launch (as in 0.5; used by `bounded-job.py`) |
| `SIMQ_CPU_LIST` | allocated logical CPUs, one per physical core (as in 0.5), e.g. `0-3`, `6` |
| `SIMQ_CPU_COUNT` | their number (as in 0.5) |
| `SIMQ_NCPUS` | same as `SIMQ_CPU_COUNT` |
| `SIMQ_PARTITION`, `SIMQ_NODE` | partition, node name |
| `SIMQ_GPUS` | allocated GPU UUIDs, comma separated; empty when none |
| `SIMQ_GPU_MODE` | `none`, `exclusive` or `shared` |
| `SIMQ_GPU_MEM` | GPU claim in bytes (the budget when exclusive; 0 when none) |
| `SIMQ_MEM` | memory limit in bytes |
| `SIMQ_TIME_LIMIT` | time limit in seconds (empty when unlimited) |
| `SIMQ_LOG`, `SIMQ_UNIT` | job log, systemd unit name |
| `CUDA_VISIBLE_DEVICES` | the allocated GPU UUIDs (empty when none); `CUDA_DEVICE_ORDER=PCI_BUS_ID` |

The rest of the environment is the submitting shell's, captured at submit time
(as in 0.5), except that inherited `SLURM_*` variables are removed unless given
with `-e`.

**No `SLURM_*` variables are exported.** `SLURM_JOB_ID`/`SLURM_JOBID` switch
pf-launch, Open MPI/PRRTE and MPICH Hydra to their Slurm code paths (srun,
Slurm host lists), and `SLURM_CPUS_PER_TASK`, `SLURM_NTASKS` or
`SLURM_JOB_NODELIST` without a running slurmd are read by several libraries and
launchers to size thread pools or build host lists. No subset is safe to fake,
so the allocation is described only by `SIMQ_*` and `CUDA_VISIBLE_DEVICES`. The
mapping below shows how a job script can read either.

### Lease files (interoperation with manual and legacy jobs)

The lease files of `any-core.sh`, `both-cores.sh` and pf-launch are used as is:
core K is `/home/treefor/projects/.workspace/outputs/pic-profile/testing-core<K>.lock`
(core 6: `testing-core.lock`), the GPU is `testing-gpu-<UUID>.lock` plus
`testing-gpu.lock`, all flock(2).

- Each job unit's wrapper takes the lease files of its cores (exclusive) and
  GPU (exclusive, or shared for a `--gpu-mem` job) before the command starts,
  in ascending path order, and passes them to the job's processes, so the
  leases stay held while any process of the job lives, including across worker
  restarts.
- The scheduler never takes a lease itself. It reads `/proc/locks`; a core or
  GPU whose lease is held while no queue job owns it is busy, and
  `simq partitions` names the holder.
- If a manual process takes a lease between the scheduler's check and the
  wrapper's start, the wrapper refuses before running anything and the job is
  requeued under the same ID.
- A job that moves its processes into another unit (for example
  `bounded-job.py`, which creates its own scope) keeps the leases through the
  processes that stay in the job unit (the wrapper and `bounded-job.py`).
- Borrow jobs (0.6.1) never take a core's lease file. They hold
  `borrow-core<K>.lock` (same directory) shared instead, so the core's own
  lease stays free for its owner: a tool such as any-core.sh can take it at
  any time, and the scheduler then preempts the borrower within a pass.

## Borrowing idle cores (0.6.1)

### Submitting borrow jobs

```sh
simq submit --borrow -p raikou-mfem -c 2 -- ./run-tests.sh            # alias: --opportunistic
simq submit --borrow-from spheral,powerflow -p raikou-mfem -- make     # only these lenders (implies --borrow)
simq submit --borrow --mem 4G --max-borrow-cores 6 -- make -j"$(nproc)"  # elastic: 1..6 cores, as many as are idle
simq submit --borrow --on-preempt freeze -- ./long-build.sh           # freeze | requeue (default) | cancel
simq submit -p powerflow --exclusive-timing -- ./benchmark.sh         # owner job whose cores are never lent
simq submit --test-only --borrow -c 2 -- /bin/true                    # where it would run now, and how
```

- A borrow job is restartable work. It may run on the cores of its own
  partition (`-p`, default `spheral`) and on lendable cores of every partition
  that lends (or of the partitions named by `--borrow-from`), in that order of
  preference, as long as no owner job of that partition is waiting for those
  cores. One job may span several partitions.
- `-c N` is the number of cores it needs (default: the partition default, as
  for owner jobs). `--max-borrow-cores M` lets it take up to `M` cores when that
  many are idle; with `--max-borrow-cores` and no `-c`, the minimum is 1. The
  job sees its cores in `SIMQ_CPU_LIST`/`SIMQ_NCPUS` and `nproc` (affinity),
  plus `SIMQ_BORROW=1` and `SIMQ_ON_PREEMPT`.
- Memory and time defaults and limits are those of its own partition
  (spheral's default memory is 18G: give builds and tests an explicit
  `--mem`, since a borrow job starts only when its memory fits the node's free
  schedulable memory). Borrow jobs cannot request GPUs (no GPU borrowing) or
  name explicit CPUs.
- Owners never wait: when an owner job is ready to start on cores a borrow job
  uses, the scheduler preempts the borrow job and starts the owner in the same
  pass. `--on-preempt` decides what happens to the borrow job:
  - `requeue` (default): the unit is stopped (SIGTERM, then SIGKILL after
    `kill_grace_seconds`); the job returns to the queue under the same ID and
    position and starts again from the beginning later.
  - `freeze`: every process of the job is frozen with the cgroup freezer
    (`systemctl --user freeze`); it keeps its memory and is thawed when its
    cores are lendable again. A job frozen longer than `freeze_max_s`
    (30 min) or within 60 s of its time limit is requeued instead; so is a
    frozen job whose memory an owner needs (freezing does not release memory).
  - `cancel`: the job ends `canceled` with result `preempted`.
- Every start, preemption (mode, reason, owner job), freeze, resume, requeue
  and end of a borrow job is recorded in the job history (`simq show ID`,
  table `job_events`); the owner job's history records which borrow jobs it
  reclaimed cores from. `simq` lists borrow jobs with `~` before the
  partition and shows `frozen` for frozen ones.

### Where borrow jobs may run (lendable cores)

| Class | Core | Lent when | Borrower preempted when |
|---|---|---|---|
| `free` | no queue job, no lease holder, and the CPU itself idle (work pinned there without a lease file keeps it) | partition `lend` (always for the job's own partition) and no owner job of that partition is waiting | an owner job needs it; anyone takes the core's lease file |
| `owner-idle` | allocated to a running owner job of the partition that leaves it unused | partition `lend`, job not `--exclusive-timing` | the owner job's process tree uses more than the cores it kept |
| `leased-idle` | lease file held outside the queue (any-core.sh, both-cores.sh, pf-launch), holder idle | partition `lend` and `lend_leased_idle` | the holder's process tree becomes busy, or another process takes the lease |

Idle and busy are measured by the worker (only while borrow jobs exist):

- a *holder tree* is every process that has the core's lease file open (from
  `/proc/<pid>/fd`; the pid in `/proc/locks` is only the process that took
  the lock and may have exited) plus all descendants; for an owner job it is
  the unit's wrapper and its descendants (including processes that
  `bounded-job.py` moves into its own scope);
- a tree is idle on `k` of its `n` cores when its CPU use over `idle_window_s`
  (30 s) leaves `k` cores unused, rounding usage up only beyond
  `idle_threshold` (0.1 core), and the lent CPU itself was at most
  `cpu_idle_threshold` (15 %) busy over the same window;
- a tree is busy when, over `busy_window_s` (2 s), it uses more than the cores
  not lent plus `busy_threshold` (0.25 core): with one leased core, any use
  above a quarter core. Detection takes one or two worker passes (about 1-3 s);
  meanwhile CPU priority already gives the holder the CPU.

Holding `borrow-core<K>.lock` exclusively from outside the queue (for example
`flock -x .../pic-profile/borrow-core6.lock sleep 3600`) marks core K as not
lendable; borrowers there are preempted.

### Partition policy (node configuration)

| Key | Default | Meaning |
|---|---|---|
| `lend` | `true` | other partitions' borrow jobs may use this partition's free cores and the unused cores of its running jobs |
| `lend_leased_idle` | `false` | also lend cores whose lease is held outside the queue while the holder tree is idle |
| `timing_sensitive` | `false` | never place any borrow job (not even the partition's own) on its cores or on their SMT siblings |

Per job, `--exclusive-timing` makes an owner job's cores and their SMT
siblings unlendable while it runs (and preempts a borrower on a sibling when
it starts). SimQueue never allocates SMT siblings (8-15) to any job, borrow
jobs included, so on zenobia the sibling clauses hold by construction; they
are enforced explicitly for configurations whose partitions list sibling
CPUs (`allow_smt_siblings`).

Node-wide settings live in the `borrow` section (`enabled`, `slice`,
`slice_cpu_weight`, `unit_cpu_weight`, `idle_window_s`, `busy_window_s`,
`idle_threshold`, `cpu_idle_threshold`, `busy_threshold`, `freeze_max_s`,
`default_on_preempt`) and `lock_files.borrow_path`.

**zenobia choices** (`config/zenobia-node.json`):

- `spheral`, `raikou-mfem`, `powerflow`: `lend: true`. A partition's free cores,
  and the unused cores of a job that only uses part of its allocation (tests,
  builds), go to borrow jobs; a full simulation uses its cores and lends
  nothing.
- `lend_leased_idle: false` everywhere. Work started with any-core.sh,
  both-cores.sh or pf-launch has no time limit and no way to mark a timing
  run, and an idle phase of such a tool (between two test runs) is not a
  promise that the core stays idle. Those cores are never lent.
- `powerflow` is not `timing_sensitive` as a whole: Power-Flow timing and
  benchmark runs submit with `--exclusive-timing`, which keeps their cores
  (and siblings) unlent while they run, while the partition's other core stays
  usable. Set `timing_sensitive: true` to keep borrowers off the whole
  partition instead.
- No GPU is lent.

### CPU priority

Owners keep the CPU even before a preemption takes effect:

- borrow units run in `simqborrow.slice`, a slice directly under the user
  manager (a sibling of `app.slice`, `session.slice` and `simq.slice`) with
  `CPUWeight=idle` (cgroup `cpu.idle=1`): its processes run only when no
  process of any sibling slice wants that CPU. The slice must not contain
  `-`: `simq-borrow.slice` would sit inside `simq.slice` and its idle weight
  would only apply against owner job units;
- each borrow unit has `CPUWeight=1`, `CPUSchedulingPolicy=idle`
  (`SCHED_IDLE`, checked and set again by `jobwrap.py` before the command
  starts, inherited by every process) and `IOSchedulingClass=idle`;
- the cpu controller is delegated to the user manager
  (`cgroup.controllers`: `cpu memory pids`); cpuset is not, so affinity
  remains `CPUAffinity` per process, as for owners.

Measured on zenobia (tests): a busy process outside the queue pinned to a
borrower's core left the borrower under 5 % of that core.
`simq doctor` reports the slice weight, the cpu controller and the freezer.

### Lease hygiene: `simq partitions`

```sh
simq partitions               # samples lease holders for 1 s
simq partitions --sample 3    # longer window; --sample 0 skips sampling; --json for scripts
```

Every core is shown as `free`, `job N` (owner job), `external lease busy` /
`external lease idle` (CPU use of the holder's process tree over the sample
window, its number of processes, and each holder's pid and command), `external
lease (not sampled)`, or `stale lease` (held according to `/proc/locks` but no
holder process is visible); borrow jobs on a core and `not lendable` markers
follow. The `POLICY` column shows `lend`, `lend,leased-idle`, `no-lend` or
`timing-sensitive`. In `--json`, each core keeps the 0.6.0 `state` text and
adds `class`, `holders`, `usage_cores`, `borrowers` and `no_borrow`.

## Running the worker as the node scheduler

`~/.config/systemd/user/simqueue.service` (installed by `switchover.py`):

```ini
[Service]
ExecStart=/home/treefor/.local/bin/simq worker --wait-for-lock
Environment=SIMQUEUE_PROJECT_ROOT=/home/treefor/projects
Environment=SIMQUEUE_STATE_DIR=/home/treefor/.local/state/simqueue
Environment=SIMQUEUE_CONFIG=/home/treefor/.local/state/simqueue/node.json
Nice=5
Restart=on-failure
KillMode=mixed
UMask=0077
```

There is no `CPUAffinity`: the worker only polls (about one `systemctl show`
per second) and places each job unit explicitly. `~/.local/bin/simq` exports the
same three variables and runs the 0.6 runtime with the pinned interpreter.

Each owner job is `simq-job-<id>-<token12>.service` in `simq.slice` (a borrow
job: in `simqborrow.slice`, with the properties above) with
`CPUAffinity=` (the allocated CPUs), `MemoryMax=`, `MemorySwapMax=`,
`OOMPolicy=kill`, `RuntimeMaxSec=` (time limit), `WorkingDirectory=`,
`KillMode=control-group`, `TimeoutStopSec=10s`, `--collect`, its log as
stdout/stderr, the `SIMQ_*` placement variables in `Environment=`, and
`ExecStopPost=` recording `$SERVICE_RESULT`/`$EXIT_CODE`/`$EXIT_STATUS`. The
main process is `source/simqueue/jobwrap.py` (stdlib only, imports nothing from
the scheduler); it checks affinity and working directory, takes the leases,
runs the command with the captured environment, waits until every process of
the unit has exited (as 0.5 did), and writes `status.json`. `AllowedCPUs=` is
not used because the user manager has no cpuset controller delegated
(`cgroup.controllers`: cpu memory pids); affinity is per process, so a job can
still widen its own mask: correct stale masks with `-e` as in 0.5.

The worker: on start (and every pass) it reconciles every active job with its
unit and status files: an active unit with the job's launch token is adopted;
a finished one is recorded with the wrapper's exit status, timestamps and
systemd's result (timeout, oom-kill, signal); a unit that never started is
requeued; a unit that vanished without a status record (reboot) is marked
`interrupted` and never restarted automatically. Resources of adopted jobs stay
reserved from their database allocation. SIGTERM, SIGKILL or a crash of the
worker leaves every job running.

### Upgrade procedure (stop worker, swap code, start worker)

```sh
systemctl --user stop simqueue.service      # jobs keep running in their own units
# install the new runtime directory next to the old one; edit ~/.local/bin/simq to point at it
systemctl --user start simqueue.service     # adopts running jobs, records jobs that ended meanwhile
simq; tail ~/.local/state/simqueue/worker.log   # look for "adopted job N"
```

### Upgrading 0.6.0 to 0.6.1: `upgrade.py`

```sh
cd /home/treefor/projects/runtimes/simqueue-0.6.1-20261031
./upgrade.py install                         # dry run: preflight and the exact actions
./upgrade.py install --execute               # upgrade; running job units keep running
./upgrade.py status
```

`install --execute`:

1. Preflight: 0.6.1 runtime sealed; its `config/zenobia-node.json` valid and
   identical to the live `node.json` apart from the 0.6.1 keys (`borrow`,
   `lend`, `lend_leased_idle`, `timing_sensitive`, `borrow_path`, notes);
   the wrapper runs 0.6.0; database schema 7; lists the active jobs.
2. Stops `simqueue.service` (job units keep running).
3. Backs up the database (sqlite3 backup API, source opened read-only), the
   wrapper, the `simqueue` link, the three simqueue unit files, `node.json` and
   `worker.log` to `~/.local/state/simqueue-backups/0.6.0-to-0.6.1-<UTC>/`
   with a manifest of sha256 sums.
4. Installs the 0.6.1 wrapper, link, `simqueue.service` and
   `simqborrow.slice`. From here on every `simq` call runs 0.6.1.
5. Migrates the database to schema 8 (new columns and the `job_events`
   table; nothing else changes) and checks that every schema-7 column of every
   existing row is unchanged; installs the 0.6.1 `node.json`.
6. daemon-reload, enable, start.
7. Verifies that the 0.6.1 worker started and adopted every job unit that was
   active before (or recorded its end) within `--adoption-timeout` (60 s);
   records `simq doctor`, `simq partitions` and the borrow slice weight in the
   receipt (`--require-doctor` also fails the install on a doctor problem).
   If the verification fails, it rolls back by itself. Receipt:
   `<backup>/install-receipt.json`.

Rollback (0.6.0 again; jobs keep running and are adopted by 0.6.0):

```sh
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.6.0-to-0.6.1-<UTC>                 # dry run
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.6.0-to-0.6.1-<UTC> --execute
./upgrade.py rollback --backup ... --execute --cancel-borrow-jobs   # when borrow jobs are queued or running
```

It refuses while borrow jobs are queued or active (0.6.0 would run them as
owner jobs) unless `--cancel-borrow-jobs`; stops the worker; keeps the
schema-8 database as `<backup>/queue-0.6.1-at-rollback.sqlite3`; restores the
wrapper, link, unit files and `node.json` byte for byte and removes
`simqborrow.slice`; then sets the database back to schema 7 in place (all
rows stay, including owner jobs started under 0.6.1 that still run; the 0.6.1
columns stay and are ignored by 0.6.0) and checks that the 0.6.0 code opens
it; starts the 0.6.0 worker and verifies that it adopted every active job.
Owner job units started by 0.6.1 run the 0.6.1 `jobwrap.py`, so the 0.6.1
runtime directory must stay until they end (as the 0.6.0 one must stay after
the upgrade).

- Keep the old runtime directory until its last job has finished: running units
  execute that runtime's `jobwrap.py`.
- A newer build must keep reading schema 7 (or migrate it in place); the worker
  refuses a database whose schema is newer than itself.
- `simq cancel` works while the worker is down (it stops the unit directly; the
  next worker records `canceled`).
- Restarting the user manager or rebooting ends job units: such jobs are marked
  `interrupted`. Linger keeps the user manager alive across logout.

### Recovery notes

- `simq doctor` checks the service, user manager, wrapper interpreter, lease
  directory and lists job units.
- `systemctl --user list-units 'simq-job-*'` shows running jobs; `systemctl --user status <unit>`
  shows the processes.
- A job that ends while the worker is down is recorded on the next start with
  its real end time.
- A job started by a 0.5 worker (only possible if a 0.5 job is still alive at
  switchover) cannot be adopted: 0.6 keeps its cores reserved, never signals it,
  and marks it `interrupted` once its processes are gone.

## Switchover from 0.5

(Historical: the 0.5 to 0.6.0 switchover ran at 2026-10-31T00:00Z. In the
0.6.1 runtime, `switchover.py` would install 0.6.1 directly over 0.5.)

State at time of writing: the 0.5 queue is idle (Spheral job 17 was canceled at
23:10Z on 2026-10-30; nothing queued). A 0.5 job runs as a child of the 0.5
worker, so stopping the 0.5 service ends it; the switchover therefore runs only
while 0.5 is idle (the script refuses otherwise). Run it only after agreement
with the Spheral session.

```sh
cd /home/treefor/projects/runtimes/simqueue-0.6.0-20261030
./switchover.py install                     # dry run: preflight and the exact actions
./switchover.py install --execute           # performs the switchover and the canaries
```

What `install --execute` does:

1. Preflight: sealed runtime files unchanged, node config valid, database
   schema 6, no active 0.5 job and no live 0.5 job process, nothing queued
   (unless `--allow-queued`), 0.6 not already installed.
2. Stops the idle 0.5 `simqueue.service` and checks that the queue did not
   change meanwhile.
3. Backs up the database (sqlite backup API), `~/.local/bin/simq`, the
   `~/.local/bin/simqueue` link target, the three simqueue unit files and
   `worker.log` to `~/.local/state/simqueue-backups/0.5-to-0.6-<UTC>/` with a
   manifest of sha256 sums.
4. Migrates the database to schema 7 with the 0.6 code and verifies that all
   0.5 columns of every row are unchanged.
5. Installs `node.json`, the 0.6 `~/.local/bin/simq`, `~/.local/bin/simqueue`
   (link to the 0.6 runtime launcher) and the 0.6 `simqueue.service`;
   daemon-reload, enable, start; runs `simq doctor` and `simq partitions`.
6. Canaries: `powerflow` runs a real one-core job that reports its affinity and
   runs `nvidia-smi -L` while holding both GPU lease files; `spheral` runs a
   four-core no-op that reports its affinity; `raikou-mfem` is validated with
   `simq submit --test-only` (no launch on cores 4-5). If cores 6-7 or the GPU
   are leased elsewhere, the powerflow canary waits up to `--canary-timeout`
   (600 s) and is then canceled and reported. Results:
   `<backup>/install-receipt.json`.

`simqueue-spheral-dispatch.path` keeps working unchanged: its service runs the
read-only `spheral-dispatch.py` from the 0.5 runtime, which calls
`~/.local/bin/simq submit --cores 4 --idempotency-key ...`; under 0.6 that is
the default partition `spheral` with cores 0-3. Its status file already records
job 7, so it does not submit again unless that file is reset.

Rollback (refuses while 0.6 job units are active):

```sh
./switchover.py rollback --backup ~/.local/state/simqueue-backups/0.5-to-0.6-<UTC>            # dry run
./switchover.py rollback --backup ~/.local/state/simqueue-backups/0.5-to-0.6-<UTC> --execute
```

It stops the 0.6 worker, saves the 0.6 database as
`<backup>/queue-0.6-at-rollback.sqlite3`, restores the schema-6 database and
the wrapper, link and unit files byte for byte, removes `node.json`, and starts
the 0.5 worker. Jobs submitted under 0.6 are not carried back (their history
stays in the saved 0.6 database).

Risks and how they are handled:

- A 0.5 job still running: stopping the 0.5 service would kill it. The script
  refuses while any 0.5 job is active or any process carries a 0.5 launch
  token.
- Database migration is one-way for the 0.5 code (0.5 refuses schema 7 and
  fails closed). Rollback restores the backup taken after the 0.5 worker
  stopped.
- A submission between preflight and stop: the script compares the queue
  before and after stopping and aborts before changing anything.
- Lease files held by Power-Flow or raikou-mfem work at switchover time: the
  0.6 scheduler simply waits for them; the powerflow canary may time out and is
  then reported, not failed.
- The 0.5 runtime directory must stay: 0.6 uses its interpreter, and
  `spheral-dispatch.py` lives there.

### Spheral resubmission

The 0.5 form works unchanged under 0.6 (default partition `spheral`, `-c 4` =
cores 0-3, memory 18G, time limit 2 days):

```sh
simq submit -c 4 -C /home/treefor/projects/runs/rz-spheral-setup-20261029 -- \
  /home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen \
  /home/treefor/projects/runtimes/simqueue-resource-20261030/bounded-job.py -- <command ...>
```

Recommended 0.6 form (explicit and Slurm-like; the time limit sits above the
command's own `--timeout-s 86400` so that the command's timeout fires first):

```sh
simq submit -p spheral -c 4 --mem 18G -t 1-01:00:00 -n zapdos-raikou-spheral \
  -C /home/treefor/projects/runs/rz-spheral-setup-20261029 -- \
  /home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen \
  /home/treefor/projects/runtimes/simqueue-resource-20261030/bounded-job.py -- <command ...>
```

The interpreter path does not change (0.6 uses the same `python-frozen`).
`bounded-job.py` sees `SIMQ_CPU_LIST=0-3`, `SIMQ_CPU_COUNT=4`, `SIMQ_JOB_ID`
and `SIMQ_LAUNCH_TOKEN` with their 0.5 meaning; the 0.6 suite runs the
unmodified `bounded-job.py` and its `scope-entry.py` gate inside 0.6 job
units. `bounded-job.py` creates its own 18G/1G scope; `--mem 18G` charges the
same amount to the node's schedulable memory. `simq cancel` stops the job unit;
`bounded-job.py`'s SIGTERM handler then stops its scope.

## Moving to Slurm later

| SimQueue 0.6 | Slurm |
|---|---|
| `simq submit -p P -c N -- cmd` | `sbatch -p P -c N --wrap 'cmd'` or `srun -p P -c N cmd` (one task with N CPUs) |
| `--cpus LIST` | no direct equivalent; `srun --cpu-bind=map_cpu:LIST` inside an allocation |
| `--gpus N` | `--gpus N` (`--gres=gpu:N`) |
| `--gpu-mem SIZE` (shared GPU) | no admission budget in Slurm; use MPS/shard gres (`--gres=shard:N`) or exclusive GPUs |
| `--mem SIZE` | `--mem SIZE` (same units: bare number = MiB) |
| `-t TIME` | `-t TIME` (same formats) |
| `--priority N` | `--nice` (inverse) or `scontrol update jobid=ID priority=N` |
| `-n/-J NAME`, `-C DIR`, `-e K=V` | `-J NAME`, `--chdir DIR`, `--export=ALL,K=V` |
| `--idempotency-key` | none; keep the controller ledger |
| `simq`, `simq partitions`, `simq show`, `simq cancel`, `simq log` | `squeue`, `sinfo`, `scontrol show job`, `scancel`, the `--output` file |
| `simq submit --test-only` | `sbatch --test-only` |
| `simq ensemble submit` | job array (`--array`) or one allocation running `srun` steps |
| priority/FIFO + EASY backfill | `SchedulerType=sched/backfill` |
| lease files | not needed inside Slurm (cgroup task/affinity plugins confine jobs) |

| Variable | Slurm equivalent |
|---|---|
| `SIMQ_JOB_ID` | `SLURM_JOB_ID` |
| `SIMQ_PARTITION` | `SLURM_JOB_PARTITION` |
| `SIMQ_NODE` | `SLURMD_NODENAME` |
| `SIMQ_CPU_LIST` | the step's affinity (`taskset -pc $$`); no variable |
| `SIMQ_NCPUS` / `SIMQ_CPU_COUNT` | `SLURM_CPUS_PER_TASK` (one task) / `SLURM_CPUS_ON_NODE` |
| `SIMQ_GPUS` | `SLURM_JOB_GPUS` (indices) and `CUDA_VISIBLE_DEVICES` |
| `SIMQ_MEM` | `SLURM_MEM_PER_NODE` (MiB) |
| `SIMQ_TIME_LIMIT` | `squeue -h -j $SLURM_JOB_ID -o %l` |
| `SIMQ_LOG` | `--output` path |

A script that runs under both can use
`JOB=${SLURM_JOB_ID:-$SIMQ_JOB_ID}`, `NCPUS=${SLURM_CPUS_PER_TASK:-$SIMQ_NCPUS}`.
pf-launch 1.1 already handles both allocations (Slurm first).

Partition configuration as `slurm.conf` (zenobia as one node; with
`SelectTypeParameters=CR_Core_Memory` and `--hint=nomultithread` each task gets
one thread per core and the SMT siblings stay idle, as in SimQueue):

```text
NodeName=zenobia Sockets=1 CoresPerSocket=8 ThreadsPerCore=2 RealMemory=30000 Gres=gpu:1
PartitionName=spheral     Nodes=zenobia Default=YES MaxTime=14-00:00:00 DefaultTime=2-00:00:00
PartitionName=raikou-mfem Nodes=zenobia MaxTime=7-00:00:00 DefaultTime=12:00:00
PartitionName=powerflow   Nodes=zenobia MaxTime=7-00:00:00 DefaultTime=12:00:00
```

Slurm partitions on one node share its cores; reproducing the fixed core sets
needs one node definition per core group (several `slurmd -N` instances or
`NodeName` with `CPUSpecList`), or accepting shared cores with per-partition
`MaxCPUsPerNode`.

## Validation of 0.6.1 (2026-10-31)

See `release.json` and
`/home/treefor/projects/.workspace/outputs/simqueue-0.6.1-20261031/`: test
receipts per suite, the observation of cores 6-7 through the scheduler while
another session held their leases (`receipts/powerflow-lease-observation.json`:
both cores classified busy with holder pids and commands, the borrow job never
placed, the real leases never taken), the CPU-priority measurement, the
owner-idle reclaim latency, the schema 7 to 8 migration of a read-only copy of
the live database, and the `upgrade.py` dry run against the live installation.
Real job units ran only on cores 0-3, in chunks started when two of those
cores were free of other work (`core-usage-notes.md`).

## Validation of 0.6.0 (2026-10-30)

See `release.json` and
`/home/treefor/sim-workspace/powerflow/outputs/simqueue-0.6.0-20261030/`:
test receipts, the migration check on a copy of the live database, and the
real-configuration validation (`validation-real-config-*`: powerflow jobs on
core 6 while core 7 was leased by another session, worker SIGKILL, adoption by
a relocated copy of the code, a spheral job that finished while the worker was
down, the GPU canary).

## Known limitations

0.6.1 borrowing:

- A frozen borrow job keeps its memory (and page cache, open files and file
  locks). Memory is reclaimed only by requeueing (done automatically when an
  owner needs it, after `freeze_max_s`, or near the time limit); a frozen
  build that holds a lock in a tree an owner also uses (for example a git or
  build-directory lock) blocks that owner until the borrower is requeued. Use
  `--on-preempt requeue` (default) for such work.
- Time spent frozen counts against the job's time limit.
- `requeue` and `cancel` give the borrower `kill_grace_seconds` (10 s) to
  exit; the owner starts at once, so for those seconds the dying borrower
  (at idle priority) and its memory overlap the owner.
- Idle/busy detection follows process trees. Work that a lease holder or an
  owner job starts outside its own process tree (for example through another
  systemd service) is not seen by the busy check; placement still requires
  the CPU itself to be idle, and CPU priority still favors the owner.
- CPU priority (idle slice, `SCHED_IDLE`) does not cover shared caches,
  memory bandwidth or SMT; a borrower on a neighboring core can still slow a
  memory-bound owner. Use `--exclusive-timing` or `timing_sensitive` for
  timing work.
- Processes outside the user manager (system services, login-session scopes)
  meet the borrow slice only through `user@1000.service`'s weight.
- Borrow jobs are placed greedily in queue order: a small borrow job can take
  cores a larger one is waiting for.
- `simq submit --test-only` for a borrow job sees only free cores (it has no
  usage history for leased or allocated cores).

0.6.0:

- Job units do not survive a reboot or a restart of the user manager; such jobs
  are marked `interrupted`.
- The time-limit shadow of a resource held outside the queue is unknown, so
  backfill cannot protect a waiting job from a later job that grabs the
  resource first once the outside holder releases it.
- `AllowedCPUs` is unavailable without cpuset delegation; containment is
  `CPUAffinity` plus the lease files.
- GPU memory claims are admission numbers, not enforced limits (as for 0.5
  ensembles).
- Priority has no aging or fair share.
