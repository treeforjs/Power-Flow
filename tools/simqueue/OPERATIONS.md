# SimQueue 0.6.3 on zenobia: operations

SimQueue 0.6 is a single-node, Slurm-shaped scheduler. One worker
(`simqueue.service`) schedules the whole node; every job runs in its own
transient systemd user service, independent of the worker.

0.6.3 removes the **partitions** on zenobia: one pool, `zenobia`, of every
physical core plus the GPU, from which every job requests. The GPU stays a
scheduled resource (and the GPU line stays), idle-core borrowing stays for
`--borrow` jobs, and the old partition names `spheral`, `raikou-mfem` and
`powerflow` are **aliases** of the pool: `-p powerflow` keeps working, with
the defaults that partition had and its old cores tried first, and jobs
recorded under the old names are adopted and scheduled in the pool. 0.6.3
also adds **`simq run`** (Slurm `srun`-shaped: submit, stream the job's log,
wait, exit with the job's exit code; a signal cancels the job), for GitHub
Actions runner jobs and terminals (see "The pool" and "`simq run`" below).
A configuration for the 16-core Ryzen 9 9950X3D2 is included
(`config/zenobia-16core-node.json`).

0.6.1 adds **opportunistic borrowing**: restartable work submitted with
`simq submit --borrow` (builds, tests) may run on idle cores of any partition
that lends them, at the lowest CPU priority, and is frozen, requeued or
canceled the moment an owner needs those cores. Owners never wait for a
borrower. With no `--borrow` job in the queue, 0.6.1 schedules exactly as
0.6.0 (see "Borrowing idle cores" below).

0.6.2 adds the **GPU line**: a queued GPU job that is ready except for the
GPU now actually wins the GPU lease files when they are released, instead of
losing every release to processes blocked in a raw `flock testing-gpu.lock`.
A small claim unit waits for the lease files in the kernel's queue on behalf
of the queue and hands the held files to the job unit, raw `flock` waiters
are moved behind the line, `gpu-lease` serves direct GPU users fairly, and
the line is published in `<pic-profile>/simq-gpu-line.json` (see "The GPU
line" below). Without queued GPU jobs, 0.6.2 schedules exactly as 0.6.1.

| Item | Path |
|---|---|
| Runtime (sealed, read-only once installed) | `/home/treefor/projects/runtimes/simqueue-0.6.3-20261105` |
| Previous runtime (kept: job units started by 0.6.2 execute its `jobwrap.py`) | `/home/treefor/projects/runtimes/simqueue-0.6.2-20261031` |
| Interpreter (unchanged from 0.5) | `/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen` |
| Code | `source/simqueue/` (`nodeworker.py`, `scheduler.py`, `borrow.py`, `telemetry.py`, `gpuqueue.py`, `gpuline.py`, `units.py`, `jobwrap.py`, `locks.py`, `nodeconfig.py`, plus the 0.5 modules) |
| Direct GPU use (0.6.2) | `bin/gpu-lease` (runs `source/simqueue/gpuline.py lease`) |
| GPU line state, tickets, log (0.6.2) | `<pic-profile>/simq-gpu-line.json`, `<pic-profile>/gpu-line.d/`, `<pic-profile>/simq-gpu-line.log` |
| GPU claims (0.6.2) | units `simq-job-gpuclaim-<token12>.service`; files `~/.local/state/simqueue/claims/<token>/{spec,control,status}.json`, `claim.log` |
| Node configuration (template) | `config/zenobia-node.json` (8 cores, the CPU until 2026-11-05); `config/zenobia-16core-node.json` (Ryzen 9 9950X3D2) |
| Node configuration (live) | `~/.local/state/simqueue/node.json` (`SIMQUEUE_CONFIG`) |
| Database | `~/.local/state/simqueue/queue.sqlite3` (schema 8, unchanged by 0.6.2 and 0.6.3; 0.6.0 used schema 7) |
| Borrow slice | `~/.config/systemd/user/simqborrow.slice` (`CPUWeight=idle`) |
| Per-job files | `~/.local/state/simqueue/jobs/NNNNNN-<token>/{spec,status,stop-post}.json` |
| Logs | `~/.local/state/simqueue/logs/NNNNNN.log`, `~/.local/state/simqueue/worker.log` |
| Upgrade 0.6.2 to 0.6.3 / rollback | `upgrade.py` (dry run unless `--execute`) |
| Switchover from 0.5 (historical) | `switchover.py` |
| Tests | `source/run-tests.py` (takes the cores' lease files first; stops when other work needs those cores; `--shared` on cores 6-7) |

## The pool on zenobia (0.6.3: no partitions)

| Name | What it is | Default cores | Default memory | Default time | Cores tried first |
|---|---|---|---|---|---|
| `zenobia` (default) | the pool: physical cores 0-7 (one logical CPU each) and the RTX 3080 Ti | 4 | 18G | 2-00:00:00 | lowest free |
| `spheral` | alias of `zenobia` (the spheral partition before 0.6.3) | 4 | 18G | 2-00:00:00 | 0-3 |
| `raikou-mfem` | alias of `zenobia` (the raikou-mfem partition before 0.6.3) | 2 | 4G | 12:00:00 | 4, 5 |
| `powerflow` | alias of `zenobia` (the powerflow partition before 0.6.3) | 2 | 4G | 12:00:00 | 6, 7 |

Limits are the pool's for every name: up to all 8 cores, 20G memory and
14-00:00:00 per job (the largest limits of the old partitions), 1G swap.
Every name may request the GPU. The defaults without `-p` are those of the old
default partition (`spheral`), so the 0.5 form `simq submit -c 4 -- ./run.sh`
behaves as before apart from which cores it gets. SMT siblings 8-15 are never
allocated. Memory is a counted resource: the sum of `--mem` over running jobs
stays within 28G (`memory_schedulable`). The GPU has an 11G admission budget:
a GPU job is exclusive unless it passes `--gpu-mem`, in which case jobs share
the GPU while their claims fit the budget. Time limits are enforced by systemd
(`RuntimeMaxSec`); an expired job ends `failed` with result `timeout`.

Aliases (node configuration, per partition):

```json
"aliases": [
  {"name": "powerflow", "default_cpus": 2, "default_mem": "4G", "default_time": "12:00:00",
   "prefer_cpus": [6, 7]},
  "plain-name"
]
```

An alias is a plain name or an object that may set `default_cpus` (`all` or
a number), `default_mem`, `default_time` and `prefer_cpus` (an ordered array,
or a CPU list); the cores, GPUs, limits and lending policy are the
partition's. Every lookup of a partition by name resolves aliases: `submit`,
`run`, `--test-only`, `--borrow-from`, `ensemble submit`, the scheduler, job
adoption, `simq partitions`, `simq doctor`. A job is recorded under the name it
was submitted with: `simq` lists it under that name, `simq show` shows
`powerflow (pool zenobia)`, and the job sees `SIMQ_PARTITION=powerflow`
(sealed tools that check `SIMQ_PARTITION` keep working). Jobs recorded under
the old names before 0.6.3 (queued, running, borrow jobs, GPU line members)
are scheduled, adopted and lent under the pool. The scheduler tries the
truly free cores first, then the alias's `prefer_cpus` in order, then the
lowest-numbered free cores; a job without an alias takes the lowest free
cores (borrow jobs take the highest free cores, as before). An alias is a
name, not a reservation: a `-p powerflow` job runs on cores 0-3 when 6-7 are
taken.

After the CPU change (Ryzen 9 9950X3D2, 16 cores, 32 threads; logical CPUs
0-15 are the first threads of cores 0-15, 16-31 their SMT siblings) install
`config/zenobia-16core-node.json` (the same pool, aliases, limits and GPU on
cores 0-15; core K's lease file is `testing-core<K>.lock`, core 6's
`testing-core.lock`, the borrow leases `borrow-core<K>.lock`; its
`smt_siblings` map is checked against sysfs). The topology check reads the
running CPU: the 16-core file is refused on the 8-core CPU, and the 8-core
file would leave cores 8-15 of the new CPU unused (still valid). See
"Upgrading 0.6.2 to 0.6.3" for the commands. `memory_schedulable` (28G)
assumes the same 31 GiB of RAM.

## Using the queue

```sh
simq submit -c 4 -- ./run.sh                          # 0.5 syntax: the pool, 4 cores (lowest free)
simq submit -c 8 --mem 16G -t 1-00:00:00 -- ./pilot.sh   # every core of the pool
simq submit -p powerflow -c 1 -t 2:00:00 --mem 2G -- ./case.sh   # old name: alias of the pool
simq submit -p powerflow --cpus 7 -- ./pinned.sh      # explicit core (any core of the pool)
simq submit -p raikou-mfem --gpus 1 -- ./gpu.sh       # exclusive GPU
simq submit -p powerflow --gpus 1 --gpu-mem 3G -- ./x # shared GPU, 3G claim
simq submit --test-only -p raikou-mfem --gpus 1 -- /bin/true   # where/when it would run; queues nothing
simq run -c 2 --mem 4G -t 30 -- ./test.sh             # 0.6.3: submit, stream the log, wait, exit with its code
simq                     # active and recent jobs (PARTITION, CPU, GPU, reason for waiting)
simq partitions          # free cores, aliases, GPU claims, lease files held outside the queue (alias: sinfo)
simq show ID             # request, allocation, unit, result, job files
simq log -f ID
simq cancel ID           # works with or without a running worker
simq priority ID N       # higher runs first; promote/demote still move to the ends
simq ensemble submit -p powerflow sweep.json          # an ensemble runs inside its own allocation
simq submit --borrow -p powerflow --max-borrow-cores 4 -- make -j"$(nproc)"   # 0.6.1: borrow idle cores
simq submit -p powerflow --exclusive-timing -- ./timing-run.sh               # 0.6.1: never lend its cores
```

Options: `-p/--partition` (a partition or an alias), `-c/--cores/--cpus-per-task N`, `--cpus LIST`,
`--gpus N`, `--gpu-mem SIZE`, `--mem SIZE` (`4G`; bare number = MiB, as
`sbatch --mem`), `-t/--time` (`MM`, `MM:SS`, `HH:MM:SS`, `D-HH[:MM[:SS]]`),
`--priority N`, `-n/-J/--name`, plus every 0.5 option (`-C`, `-e`, `--shell`,
`--idempotency-key`, `--json`). Without `-p` a job goes to the default
partition (the pool); without `-c` it gets the default of the name it used
(see the table above). The 0.5 options are accepted unchanged and an
idempotency key replayed with 0.5 syntax keeps its 0.5 fingerprint; a key
replayed with the same `-p` name keeps its fingerprint across 0.6.2 and 0.6.3
(the name given enters the fingerprint, not the pool it resolves to).

Scheduling: waiting jobs are ordered by priority, then queue order. Each worker
pass starts every job that fits, except that once a job is blocked, a later
job starts only if (given all time limits) it cannot delay that blocked job
(EASY backfill). With one pool, a job that needs all cores (an 8-core pilot)
reserves the pool when it is first blocked: later jobs start only if they end
before it can start, or fit beside it. A paused ensemble holds its position
and blocks later jobs of its own partition, which on zenobia is now every
later owner job (as in 0.5); `simq demote ID` moves it out of the way. Borrow
jobs (0.6.1) always sort after every owner job, whatever their priority, and
never take part in backfill reservations.

### `simq run` (0.6.3)

```sh
simq run -c 2 --mem 4G -t 30 -n pr-tests --borrow -- ./ci/test.sh
simq run -c 8 --mem 16G -t 1-00:00:00 -n pilot --idempotency-key "gh:$GITHUB_RUN_ID:$GITHUB_RUN_ATTEMPT:pilot" -- ./pilot.sh
```

`simq run` (alias `srun`) takes every option of `simq submit`, queues the job
exactly as `submit` does, prints `simq: job N (NAME) queued in partition ...`
on stderr, then copies the job's combined log (the header written by the job
wrapper, the command's stdout and stderr, the final `state:` line) to its
stdout as it grows, waits until the job has ended and exits with the job's own
exit code when one was recorded (a command killed by signal N: 128+N), else 0
for `succeeded` and 1 for any other end (`failed` by time limit or memory,
`canceled`, `interrupted`); a job that did not succeed never exits 0. Its own
messages go to stderr only (with `--json`, the submit result is one JSON line
on stderr). SIGINT, SIGTERM or SIGHUP cancel the job (as `simq cancel`: a
queued job ends at once, a running unit is stopped) and `simq run` exits with
130, 143 or 129 after the job has ended (at most 20 s; a second signal exits
at once; the cancellation is recorded either way). A borrow job that is
preempted and requeued keeps the client waiting; the restart's output follows
in the same log. Stdin is not forwarded: the command runs in its own job unit.

With `--idempotency-key KEY`, running the same `simq run` again (same name,
directory, command, options and environment, as for `submit`) reattaches to
the existing job instead of queuing a second one: it prints
`attached to existing job N`, streams the whole log from the beginning and
returns that job's exit code (at once if it has already ended). A different
request under the same key is refused, as for `submit`. Because the captured
environment is part of the request, a GitHub Actions step reattaches when it
is retried with the same environment (a retry loop inside the step); a re-run
of the workflow has a new `GITHUB_RUN_ATTEMPT` and step file paths, so give it
a key that includes `$GITHUB_RUN_ATTEMPT` and it submits a fresh job.

Exit statuses of `simq run` itself: the job's as above; 2 for an invalid
request (nothing queued); 1 when the protected runner of the 0.5 exclusive
mode cannot start (the job stays queued).


### Environment of every job

| Variable | Meaning |
|---|---|
| `SIMQ_JOB_ID` | job ID (as in 0.5) |
| `SIMQ_LAUNCH_TOKEN` | random per launch (as in 0.5; used by `bounded-job.py`) |
| `SIMQ_CPU_LIST` | allocated logical CPUs, one per physical core (as in 0.5), e.g. `0-3`, `6` |
| `SIMQ_CPU_COUNT` | their number (as in 0.5) |
| `SIMQ_NCPUS` | same as `SIMQ_CPU_COUNT` |
| `SIMQ_PARTITION`, `SIMQ_NODE` | the name the job was submitted under (a partition or an alias, 0.6.3: `zenobia`, `spheral`, `raikou-mfem`, `powerflow`), node name |
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

- Lock order on this machine (0.6.2): the GPU files first
  (`testing-gpu.lock`, then `testing-gpu-<UUID>.lock`), then the core leases.
  The reverse order deadlocked once; every tool that blocks on more than one
  lease must use this order.
- Each job unit's wrapper takes the lease files of its GPU (exclusive, or
  shared for a `--gpu-mem` job) and of its cores (exclusive) before the
  command starts, in that order and without blocking, and passes them to the
  job's processes, so the leases stay held while any process of the job
  lives, including across worker restarts. When the GPU was contended, the
  GPU files arrive already held from the queue's GPU claim (see "The GPU
  line").
- The worker never takes a lease itself. It reads `/proc/locks`; a core or
  GPU whose lease is held while no queue job owns it is busy, and
  `simq partitions` names the holder. (0.6.2: the GPU claim units, which are
  separate units like job units, wait for and hold the GPU files.)
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

## The GPU line (0.6.2)

`<pic-profile>` is `/home/treefor/projects/.workspace/outputs/pic-profile`
(the same directory that `any-core.sh` reaches through its own, older path).

**Problem fixed.** Before 0.6.2 the worker only observed the GPU lease files.
A process blocked in `flock testing-gpu.lock` gets a released flock from the
kernel at once, while the worker looked again up to a second later, so a
queued GPU job lost every release to workers that ran raw `flock` commands
back to back (job 169 waited more than 20 minutes this way).

**How it works.**

1. Each pass the worker finds the queued GPU jobs that are *ready except for
   the GPU*: an exclusive job that the planner would start now if everybody
   holding the GPU released it (outside holders, or the queue's own jobs on
   it, whose cores free at the same moment); a shared (`--gpu-mem`) job that
   it would start if only the outside holders released it, the queue's own
   GPU jobs keeping their cores and GPU budget. A shared job that waits for
   cores or GPU budget held by queue jobs is not waiting for the lease files
   (it joins their shared lock when it starts), so it gets no claim.
2. For such jobs it starts one *claim* per GPU:
   `simq-job-gpuclaim-<token12>.service`, a small unit (one Python process,
   `MemoryMax=256M`) that is independent of the worker like a job unit. The
   claim enters the kernel's wait queue on `testing-gpu.lock`; raw `flock`
   waiters that are ahead of it there are moved behind it (below). When it
   wins `testing-gpu.lock` it takes `testing-gpu-<UUID>.lock` without
   blocking.
3. The worker then starts the job unit with a hand-off: the job's wrapper
   receives the claim's open lock files over an abstract unix socket
   (SCM_RIGHTS, checked with a per-claim secret and the peer uid), so the lock
   is never released between claim and job. The job's time limit starts with
   its unit, after the GPU was won. A `--gpu-mem` (shared) job turns the
   received locks into shared ones once all its leases are held. If the job
   cannot take its core leases, it gives the GPU back to the claim and is
   requeued as before; if the hand-off itself fails, the wrapper takes the GPU
   files without blocking, as 0.6.1 did.
4. A claim is stopped when no queued job can use it (a held claim after 5 s,
   a waiting one at once), when the queue is paused, or when it held the GPU
   for `claim_hold_max_s` (300 s) without a hand-off. Claims survive worker
   restarts and upgrades; the next worker adopts them from
   `~/.local/state/simqueue/claims/`.

**Order of service** (the line). Members are served by their *since* time:
a claim's is the submit time of the oldest ready queued job (exact to the
microsecond: the creation time of the job's log file, which `simq submit`
creates in the same transaction); a `gpu-lease` user's is the moment it
started. Ties go to the queue. A member enters the kernel queue only when
every earlier member is in it or holds the GPU, and leaves it again when an
earlier member arrives, so the kernel's FIFO order equals the line order.
Raw `flock` users are served after every member that is in the line when the
lock is released.

**Raw waiters are moved behind the line.** A process blocked in a raw
`flock` on the GPU file ahead of a member in the kernel queue gets SIGSTOP and,
about 30 ms later, SIGCONT: the stop interrupts its blocked flock(2), the
continue restarts it, and the restarted request queues at the tail. The
process notices nothing else (`flock` neither sees the signals nor loses its
command or timeout). Only processes of this user whose command is listed in
`requeue_commands` (default `flock`) and that are blocked are touched; holders
are never touched, nor jobs of an interactive job-control shell (a terminal
would report them stopped). A member that dies between the two signals is
cleaned up by the next reader of its ticket (SIGCONT). Every move is logged in
`<pic-profile>/simq-gpu-line.log`. `gpu-lease --no-requeue-raw` and
`"requeue_raw_waiters": false` turn the moves off.

**Lock order and deadlocks.** Only the first GPU file is waited for in the
kernel; the second file and core leases are polled without blocking. A
member that finds the holder of a file it needs blocked on a file the member
holds (someone using the reversed order, or holding a core lease and waiting
for the GPU) releases everything and re-enters the line with its original
since time. Raw users of two different orders can still deadlock with each
other, whatever the queue does; use `gpu-lease`.

### Direct GPU use: `gpu-lease`

```sh
G=/home/treefor/projects/runtimes/simqueue-0.6.2-20261031/bin/gpu-lease
$G -- ./gpu-case.sh                       # GPU only
$G --core any -- ./gpu-case.sh            # GPU, then core 6 or 7 (as any-core.sh --wait)
$G --core both -- ./build-and-run.sh      # GPU, then cores 6 and 7 (as both-cores.sh)
$G --core 4 --label "raikou probe" -- ./probe.sh
$G --timeout 600 -- ./gpu-case.sh         # give up after 10 min of waiting (exit 75)
$G --status                               # holders, line, raw waiters, queue jobs (--json)
```

It waits for queued SimQueue GPU jobs that were waiting before it started
and for earlier `gpu-lease` users, then takes `testing-gpu.lock`,
`testing-gpu-<UUID>.lock` and the core lease (in that order), runs the command
(with `--core`: pinned to the core, nice 15, `PF_CORE`, `OMP_NUM_THREADS=1` as
`any-core.sh`), and releases everything when the command and its children
have closed the files. It exports `GPU_LEASE_GPU`, `GPU_LEASE_LOCKS` and
(with `--core`) `GPU_LEASE_CORES`. Exit status is the command's; 75 after
`--timeout`; SIGTERM/SIGHUP are passed to the command. It keeps working while
the worker is stopped (it then waits for nobody in the queue) and after a
rollback to 0.6.1. Raw `flock testing-gpu.lock ...` still works but is served
last; use `gpu-lease` instead.

### Published state: `<pic-profile>/simq-gpu-line.json`

Written by the worker (atomic replace) whenever it changes and at least every
10 s (`heartbeat_unix`); schema `simqueue-gpu-line-v1`. Per GPU (`gpus.<UUID>`):
`queue` (the queued GPU jobs that are ready except for the GPU, that is the
queue's part of the line, in dispatch order, with `job_id`, `name`,
`partition`, `mode`, `since_unix`/`since`, `ready_since_unix`, `waited_s`,
`state`; a GPU job still waiting for cores is not listed), `jobs` (queue
jobs on the GPU), `claim`, `holders`, `helpers`
(`gpu-lease` users), `raw_waiters`, `blocked` and `warnings`. A tool that
wants to yield should wait while `queue` is non-empty and the file is fresh
(`now - heartbeat_unix <= stale_after_s`, 120 s). When the worker stops it
sets `heartbeat_unix` to 0 (withdrawn).

### Display

- `simq partitions` / `simq sinfo`: per GPU, the holder (queue job, `gpu-lease`
  user, the queue's claim, or raw flock holder with pid and command), the
  line (queue jobs and `gpu-lease` users in order, with since times), raw
  waiters, and `WARNING: queue job N has waited ... for the GPU` once a ready
  job has waited longer than `starvation_warn_s` (600 s). `--json` adds
  `gpus[].line` and `gpu_line`.
- `simq show ID` (queued GPU job): `GPU line` with its position, waiting time,
  claim phase and holder; the job history records `gpu-wait`, `gpu-claim` and
  `gpu-handoff` (with `waited_s`).
- `simq` (status): the reason column says `waiting for the GPU (line position
  N of the queue, ...; claim queued; GPU held by ...)`.
- `simq doctor`: a `GPU line` line (state age, ticket directory, limits).
- `worker.log`: claims started, adopted, stopped and ended, and the starvation
  warning.

### Configuration (`gpu_line` in node.json)

| Key | zenobia | Meaning |
|---|---|---|
| `enabled` | true | default true when lease files and a GPU are configured |
| `state_path`, `ticket_dir`, `log_path` | see above | default: next to the GPU lease files |
| `starvation_warn_s` | 600 | warning threshold for a ready queued GPU job |
| `claim_hold_max_s` | 300 | a held claim that no job takes releases the GPU |
| `stale_after_s` | 120 | readers ignore an older published state |
| `requeue_raw_waiters`, `requeue_commands` | true, `["flock"]` | move raw waiters behind the line |
| `poll_s` | 0.1 | polling interval of claims and helpers |
| `helper_cores` | `[6, 7]` | cores of `gpu-lease --core any/both` |
| `claim_cpus` | (unset) | CPUAffinity of claim units (tests only) |

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
  partition (`-p`) and on lendable cores of every partition that lends (or of
  the partitions named by `--borrow-from`), in that order of preference, as
  long as no owner job of that partition is waiting for those cores (0.6.3,
  `backfill_waiting_owners`: unless the borrow job's time limit ends before
  that owner job can start). One job may span several partitions. On zenobia
  (0.6.3) there is one partition: every name, including `--borrow-from`
  names, resolves to the pool, so a borrow job may use any lendable core.
- `-c N` is the number of cores it needs (default: the partition default, as
  for owner jobs). `--max-borrow-cores M` lets it take up to `M` cores when that
  many are idle; with `--max-borrow-cores` and no `-c`, the minimum is 1. The
  job sees its cores in `SIMQ_CPU_LIST`/`SIMQ_NCPUS` and `nproc` (affinity),
  plus `SIMQ_BORROW=1` and `SIMQ_ON_PREEMPT`.
- Memory and time defaults and limits are those of the name it used (the
  pool's and spheral's default memory is 18G: give builds and tests an explicit
  `--mem`, since a borrow job starts only when its memory fits the node's free
  schedulable memory, and an explicit `-t`, since a short time limit lets it
  use cores a waiting owner job needs later). Borrow jobs cannot request GPUs (no GPU borrowing) or
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
| `free` | no queue job, no lease holder, and the CPU itself idle (work pinned there without a lease file keeps it) | partition `lend` (always for the job's own partition) and no owner job of that partition is waiting for it; with `backfill_waiting_owners` (0.6.3, on in zenobia's configuration) also while one waits, if the borrow job's time limit ends before that owner job's earliest start | an owner job needs it; anyone takes the core's lease file |
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
`default_on_preempt`, 0.6.3: `backfill_waiting_owners`) and
`lock_files.borrow_path`.

`backfill_waiting_owners` (0.6.3, default `false`, `true` on zenobia): EASY
backfill for borrow jobs. Without it (the 0.6.1 rule) a free core is never
lent while an owner job that needs it waits; in a single pool every waiting
owner job needs every core, so one waiting 8-core pilot would keep every
borrow job (PR tests) off all free cores for as long as it waits. With it, a
borrow job may take such a core when its time limit ends no later than the
owner job's earliest start (computed from the time limits of the running
jobs, as the reservation is; a job held for the reservation counts from the
reserved start). The owner still never waits: if it can start earlier than
expected (a job ended early), the borrow job is preempted as always. A borrow
job without a time limit, or with one past the owner's start, still waits;
an owner job blocked by a job without a time limit or a lease held outside the
queue has no start to protect, and its cores are lent.

**zenobia choices** (`config/zenobia-node.json`):

- The pool `zenobia`: `lend: true`. Free cores, and the unused cores of a job
  that only uses part of its allocation (tests, builds), go to borrow jobs; a
  full simulation uses its cores and lends nothing.
- `backfill_waiting_owners: true` (see above): PR tests with short time limits
  keep running while a pilot waits for every core.
- `lend_leased_idle: false` everywhere. Work started with any-core.sh,
  both-cores.sh or pf-launch has no time limit and no way to mark a timing
  run, and an idle phase of such a tool (between two test runs) is not a
  promise that the core stays idle. Those cores are never lent.
- The pool is not `timing_sensitive`: timing and benchmark runs submit with
  `--exclusive-timing`, which keeps their cores (and siblings) unlent while
  they run, while the other cores stay usable. `timing_sensitive: true` on the
  pool would keep borrowers off every core.
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

### Upgrading 0.6.2 to 0.6.3: `upgrade.py`

```sh
cd /home/treefor/projects/runtimes/simqueue-0.6.3-20261105
./upgrade.py install                         # dry run: preflight, the layout change and the exact actions
./upgrade.py install --execute               # upgrade; running job units keep running
./upgrade.py status
# after the CPU change to the 9950X3D2 (once 0.6.3 runs, or instead of the 8-core file):
./upgrade.py install --config config/zenobia-16core-node.json            # dry run
```

`install --execute`:

1. Preflight: 0.6.3 runtime sealed; the new node configuration
   (`config/zenobia-node.json`, or `--config FILE`) valid for 0.6.3 on this
   machine's CPU and identical to the live `node.json` apart from the
   partition layout (`partitions`, `default_partition`, `smt_siblings`), the
   0.6.3 key `borrow.backfill_waiting_owners` and the notes. The layout change
   must keep every core and GPU of the live configuration, and every
   partition name of the live configuration, and of every queued or active
   job in the database, must resolve in the new one (as a partition or an
   alias); the dry run prints both layouts, how each old name resolves, the
   queued and active jobs by recorded partition and the GPU claim units. The
   wrapper must run 0.6.2; database schema 8.
2. Waits until no job is in state `starting` (`--starting-timeout`, 60 s),
   stops `simqueue.service` (job and GPU claim units keep running) and checks
   again; a job that became `starting` meanwhile makes it start the 0.6.2
   worker again and stop without changing anything.
3. Backs up the database (sqlite3 backup API, source opened read-only), the
   wrapper, the `simqueue` link, the unit files (`simqueue.service`, the two
   spheral dispatch units, `simqborrow.slice`), `node.json` and `worker.log` to
   `~/.local/state/simqueue-backups/0.6.2-to-0.6.3-<UTC>/` with a manifest of
   sha256 sums (and the layout comparison and the active claim units).
4. Installs the 0.6.3 wrapper, link, `simqueue.service` and the new
   `node.json`. From here on every `simq` call runs 0.6.3. `simqborrow.slice`
   stays as it is.
5. Opens the database with the 0.6.3 code and checks that the schema (8),
   the columns and every column of every row are unchanged: 0.6.3 migrates
   nothing and does not rewrite the old partition names (the aliases resolve
   them).
6. daemon-reload, enable, start.
7. Verifies that the 0.6.3 worker started (its start line lists the pool and
   its aliases), adopted every job unit that was active before (or recorded
   its end), adopted every GPU claim unit that was active (or saw it end) and
   published `<pic-profile>/simq-gpu-line.json` within `--adoption-timeout`
   (60 s), and that `simq partitions` runs; records `simq doctor` and
   `simq partitions` in the receipt (`--require-doctor` also fails the
   install on a doctor problem). If the verification fails, it rolls back by
   itself. Receipt: `<backup>/install-receipt.json`.

Rollback (0.6.2 again; jobs keep running and are adopted by 0.6.2):

```sh
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.6.2-to-0.6.3-<UTC>             # dry run
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.6.2-to-0.6.3-<UTC> --execute
```

It checks the backup against its manifest, stops the worker (the 0.6.3
worker withdraws the published GPU line state), stops the GPU claim units
recorded in `~/.local/state/simqueue/claims/` that have not finished (the
0.6.2 worker plans with the restored partitions and sets up its own; job
units are never stopped), keeps the database as
`<backup>/queue-0.6.3-at-rollback.sqlite3`, restores the wrapper, link, unit
files and `node.json` byte for byte, checks that the 0.6.2 code opens the
database (schema 8, nothing to undo), starts the 0.6.2 worker and verifies
that it adopted every active job. A running job keeps its cores whatever
name it has. Queued jobs recorded under an old name run under 0.6.2 in that
partition. Queued jobs recorded under the pool's name (`zenobia`: submitted
without `-p` or with `-p zenobia` under 0.6.3) are unknown to the 0.6.2
configuration and stay queued ("unknown partition zenobia"); the rollback
lists them (`NOTE: job(s) [...]`, `stranded_jobs` in the receipt): cancel and
resubmit them with an old name, or install 0.6.3 again. Job units started by
0.6.3 run the 0.6.3 `jobwrap.py`, so the 0.6.3 runtime directory must stay
until they end (as the 0.6.2 one must stay after the upgrade). The 0.6.1 to
0.6.2 procedure is documented in the 0.6.2 runtime's `OPERATIONS.md` and
`upgrade.py`.

Moving to the 16-core configuration after the CPU change: with 0.6.3
installed, `node.json` is the only file to change. Stop the worker
(`systemctl --user stop simqueue.service`; job units keep running), copy
`config/zenobia-16core-node.json` to `~/.local/state/simqueue/node.json`
(mode 600), start the worker and check `simq partitions` (16 cores) and
`simq doctor`. If 0.6.2 is still installed at that point, pass
`--config config/zenobia-16core-node.json` to `upgrade.py install` instead.
The 0.6.2 configuration's partitions name cores 0-7 only, which remain valid
physical cores on the new CPU; the 8-core pool file keeps working there too
(cores 8-15 unused).

- Keep the old runtime directory until its last job has finished: running units
  execute that runtime's `jobwrap.py`.
- A newer build must keep reading schema 8 (or migrate it in place); the worker
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

The 0.5 form works unchanged under 0.6.3 (default partition: the pool, with
the old spheral defaults: memory 18G, time limit 2 days; `-c 4` = the four
lowest free cores, 0-3 on an idle node):

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
With `-p spheral` (an alias since 0.6.3) cores 0-3 are tried first.
`bounded-job.py` sees `SIMQ_CPU_LIST=0-3` (or the cores it got), `SIMQ_CPU_COUNT=4`, `SIMQ_JOB_ID`
and `SIMQ_LAUNCH_TOKEN` with their 0.5 meaning; the 0.6 suite runs the
unmodified `bounded-job.py` and its `scope-entry.py` gate inside 0.6 job
units. `bounded-job.py` creates its own 18G/1G scope; `--mem 18G` charges the
same amount to the node's schedulable memory. `simq cancel` stops the job unit;
`bounded-job.py`'s SIGTERM handler then stops its scope.

## Moving to Slurm later

| SimQueue 0.6 | Slurm |
|---|---|
| `simq submit -p P -c N -- cmd` | `sbatch -p P -c N --wrap 'cmd'` (one task with N CPUs) |
| `simq run -p P -c N -- cmd` (0.6.3) | `srun -p P -c N cmd` (foreground; exit code of the job; Ctrl-C cancels) |
| `--cpus LIST` | no direct equivalent; `srun --cpu-bind=map_cpu:LIST` inside an allocation |
| `--gpus N` | `--gpus N` (`--gres=gpu:N`) |
| `--gpu-mem SIZE` (shared GPU) | no admission budget in Slurm; use MPS/shard gres (`--gres=shard:N`) or exclusive GPUs |
| `--mem SIZE` | `--mem SIZE` (same units: bare number = MiB) |
| `-t TIME` | `-t TIME` (same formats) |
| `--priority N` | `--nice` (inverse) or `scontrol update jobid=ID priority=N` |
| `-n/-J NAME`, `-C DIR`, `-e K=V` | `-J NAME`, `--chdir DIR`, `--export=ALL,K=V` |
| `--idempotency-key` | none; keep the controller ledger (for `simq run`: reattach on retry) |
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

The pool as `slurm.conf` (zenobia as one node; with
`SelectTypeParameters=CR_Core_Memory` and `--hint=nomultithread` each task gets
one thread per core and the SMT siblings stay idle, as in SimQueue). Slurm
partitions over the same node share its cores, which is what the 0.6.3
aliases are: the old names keep their defaults over one pool (Slurm has no
per-partition core preference; `prefer_cpus` has no equivalent):

```text
NodeName=zenobia Sockets=1 CoresPerSocket=8 ThreadsPerCore=2 RealMemory=30000 Gres=gpu:1
PartitionName=zenobia     Nodes=zenobia Default=YES MaxTime=14-00:00:00 DefaultTime=2-00:00:00
PartitionName=spheral     Nodes=zenobia MaxTime=14-00:00:00 DefaultTime=2-00:00:00
PartitionName=raikou-mfem Nodes=zenobia MaxTime=14-00:00:00 DefaultTime=12:00:00
PartitionName=powerflow   Nodes=zenobia MaxTime=14-00:00:00 DefaultTime=12:00:00
```

After the CPU change: `CoresPerSocket=16`. Default memory per name maps to
`DefMemPerNode` per partition; borrow jobs map to a preemptible partition
(`PreemptType=preempt/partition_prio`, `PreemptMode=REQUEUE`).

## Validation of 0.6.2 (2026-10-31)

See `release.json` and
`/home/treefor/projects/.workspace/outputs/simqueue-0.6.2-20261031/`
(`receipts/`, `experiments/`):

- Kernel behavior the GPU line relies on, checked on private lock files
  (`experiments/kernel_flock.out.json`): blocked flock(2) waiters are served
  in FIFO order; SIGSTOP then SIGCONT of a blocked raw `flock` moves it to
  the tail of the queue; an exclusive lock turns shared atomically on the
  same open file description; a lock received over SCM_RIGHTS stays held
  after the sender exits.
- `test_gpuline.py` (24, processes on private lock files): which queued
  jobs are ready except for the GPU (pure), `/proc/locks`
  chain parsing, a queue claim and a `gpu-lease` user each win the first
  release over earlier back-to-back raw waiters, helpers are FIFO by start
  time, a helper yields to an earlier queued job, back-off from reversed
  lock orders and from core holders waiting for the GPU, cleanup after a
  member dies between SIGSTOP and SIGCONT, interactive jobs are not
  requeued, hand-off without a gap, wrong secret refused, atomic shared
  downgrade.
- `test_gpu_line_units.py` (8, real claim and job units on the leased test
  cores with private lock files and a private `<pic-profile>`): a queued GPU job starts
  at the first release despite raw `flock` and `gpu-lease` users; a helper
  that started before the job goes first; the time limit starts with the
  job; shared jobs end with a shared lease; cancel of a waiting job stops the
  claim; worker restart while the job waits and while it runs; no deadlock
  between GPU users and core leases; `simq partitions`, `simq show` and the
  published state.
- `test_upgrade.py` (7): 0.6.1 job units keep running through the install
  and are adopted by 0.6.2, 0.6.2 job units through the rollback (claims
  stopped), refusal while a job is starting, automatic rollback after a
  failed verification, dry run changes nothing.
- The 0.6.1 suites, unchanged, through `run-tests.py` (leases held, live
  guard on, private state and lock files; the real GPU lease files never
  touched): first on cores 0-1 (lent by Spheral) until Spheral asked for
  them back (`RETURN-CORES-0-3`, which stopped a run as designed), then on
  cores 6-7 with `--wait-any --shared`, between Power-Flow jobs under the
  same lease files. `test_exclusive_mode.py` runs as a
  pure suite (mocks only; its fixture needs four CPUs in the affinity
  mask). `test_shield_helper.py` cannot be imported without
  `source/privileged/simqueue-shield`, as for 0.6.0 and 0.6.1.
- The test temporary directory is used by its resolved path (the outputs
  directory lies behind a symbolic link): with the unresolved path the
  wrapper's working-directory check and the 0.5 submit fingerprint
  comparison fail in the tests, because 0.6.1 ran its tests in `/tmp`.
- `test_switchover.py` (the historical 0.5 to 0.6 switchover) runs only on
  a sealed runtime and with a temporary directory whose path does not
  contain `simqueue-0.6.` (`switchover.py` takes that text in the sandbox
  unit for an installed 0.6); it ran with
  `SIMQ_TEST_TMPDIR=<outputs>/simq062-test-tmp`.
- `test_borrow_units.py` places real borrow jobs by a CPU-idle sample of
  the leased cores; with only two cores, load of the test harness itself can
  make one of them look busy (one run placed the borrower on the other
  core). The test and the borrow code are unchanged since 0.6.1.
- `upgrade.py install` dry run against the live installation: no problem.
- Post-install check on the live queue (2026-11-01): a shared check job
  waiting for a core held by a shared queue GPU job was counted as ready
  except for the GPU and got an exclusive claim, which held back every
  queued GPU job until the GPU drained. The readiness rule was corrected
  (`gpuqueue.ready_except_gpu`, five pure tests in `test_gpuline.py`), the
  GPU suites were run again, and the worker was restarted on the corrected
  code (job and claim units are separate units and kept running).

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

0.6.3 pool and `simq run`:

- An alias is a name, not a reservation: under load a `-p powerflow` job can
  land on cores 0-3. Sealed Power-Flow tools that check
  `<pic-profile>/RETURN-CORES-0-3` refuse to run on cores 0-3 while that file
  exists (it does on 2026-11-04); remove the file once the pool is live, or
  expect such jobs to fail when they spill there. Tools that check
  `SIMQ_PARTITION == "powerflow"` keep working (the job keeps the name).
- Placement is lowest-free-core first (or the alias's preferred cores); it
  knows nothing about the 9950X3D2's two CCDs, so a job can be split across
  them when the pool is fragmented. Use `--cpus 0-7` / `--cpus 8-15` for a
  job that must stay on one CCD.
- A paused ensemble now holds every later owner job of the pool (it held
  only its partition's); demote or resume it.
- `--borrow-from NAME` no longer narrows anything on zenobia: every name is
  the pool.
- `simq run` polls the database and the log (0.2 s); output appears in
  chunks as the job writes its log. It does not forward stdin or
  propagate signals other than SIGINT/SIGTERM/SIGHUP (which cancel the job).
  If the client is killed with SIGKILL, the job keeps running; rerun with the
  same `--idempotency-key` to reattach.
- An idempotency key reattaches only for an identical request, captured
  environment included (as for `submit`); see "`simq run`".
- Rolling back to 0.6.2 leaves jobs queued under the pool's name waiting (see
  the upgrade section).

0.6.2 GPU line:

- The line never preempts. A process that already holds the GPU lease files
  (raw `flock`, `gpu-lease`, a pf-launch run) keeps them until it ends; the
  line only decides who gets them next. A long outside holder still delays a
  queued GPU job, which `simq partitions` then flags after
  `starvation_warn_s`.
- Only blocked waiters whose command (`/proc/PID/comm`) is listed in
  `requeue_commands` (`flock`) are moved behind the line. A program that
  blocks in flock(2) on a GPU lease file by itself (for example Python
  `fcntl.flock` without `LOCK_NB`), and any job-control job of an interactive
  shell, keeps its place in the kernel queue and can still take a release
  ahead of the line. Such waiters are listed as raw waiters in the published
  state and logged as "not requeueing".
- A queued job joins the line only when it is ready except for the GPU
  (cores, memory and the planner's order allow it). While it waits for
  cores or behind other jobs, `gpu-lease` users that start meanwhile go
  first; once it is ready its since time is its submit time, so helpers
  still waiting that started after the submit move behind it (a helper that
  already holds the GPU keeps it).
- A shared (`--gpu-mem`) job that waits for the cores of a queue job on the
  same GPU (or for GPU budget) gets no claim; when the last queue job on the
  GPU ends, a raw `flock` or `gpu-lease` user blocked at that moment can take
  the release first, and the shared job joins the line after that (one
  release cycle at most, since raw waiters are then moved behind its claim).
- While a claim waits for the GPU, queued GPU jobs do not start on it
  directly, even shared jobs that would fit next to the queue's shared jobs
  already on it: an exclusive job in the line waits until the GPU drains,
  and later shared jobs do not overtake it.
- Between the claim winning the GPU and the hand-off (normally well under a
  second) and after a failed start (the job's core lease taken outside the
  queue in the meantime), the claim holds the GPU without using it: up to
  5 s when no queued job can use it, up to `claim_hold_max_s` (300 s) while
  a job is still expected.
- The line order uses the creation time (statx birth time) of the job's log
  file; on a filesystem without birth times it falls back to the
  modification time within the `created_at` second, else the end of that
  second.
- Raw users of two different lock orders (for example core lease then GPU,
  and GPU then core lease) can deadlock with each other; claims and
  `gpu-lease` back off from such chains but cannot break them. Use
  `gpu-lease`, which always takes the GPU files first.
- The published state is advisory: a tool that neither reads it nor uses
  `gpu-lease` is affected only through the raw-waiter moves.
- After a rollback to 0.6.1 queued GPU jobs again lose releases to raw
  `flock` waiters (`gpu-lease` keeps serving its own users fairly).

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
