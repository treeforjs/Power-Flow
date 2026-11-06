# SimQueue 0.8.0 on zenobia: operations

SimQueue is a single-node, Slurm-shaped scheduler. One worker
(`simqueue.service`) schedules the whole node; every job runs in its own
transient systemd user service, independent of the worker.

0.8.0 (a minor release: schema 10, new scheduling inputs) is the version for
the 16-core node (Ryzen 9 9950X3D2 since 2026-11-06: cores 0-15, SMT siblings
16-31, two L3 domains CCD0 = cores 0-7 and CCD1 = cores 8-15, 96 MiB each),
see "0.8.0" below: a job's cores stay inside one CCD when they fit (heavy jobs
spread over the CCDs, light jobs pack, `--l3-domains`/`--split-domains`),
SMT siblings are lent to borrow jobs, `-t auto` derives time limits from
history, a partition's `cpus` may be `"auto"`, and `simq partitions` shows
the memory admission ceiling.

0.7.2 (a patch release: no schema change, no new scheduling input) contains
MPI jobs, see "0.7.2" below: every job gets its own temporary and MPI
rendezvous directory, so Open MPI's `prun`/`--dvm` can no longer attach to a
PRRTE DVM outside the job and run the ranks outside the job's cgroup; the
worker reports processes carrying a job's launch token outside its cgroup
and MPI launcher daemons outside every job.

0.7.1 (a patch release: no new scheduling input, no schema change, no new
side effect) changes defaults and fixes, see "0.7.1" below: a job submitted
without `--mem` gets `--mem auto` when its name has a recorded peak; the
`spheral` alias's default memory is auto instead of 18G; the load-sensitive
tests are robust; the secret values in the older upgrade backups are
redacted.

0.7.0 changes four things (see "0.7.0: event-driven worker, memory,
GitHub reporting" below):

- **The worker is event-driven and runs at low priority.** It sleeps until
  a submission, a cancellation or a job's end wakes it (a FIFO in the state
  directory), reads files between passes instead of starting `systemctl`,
  and runs with `Nice=10`, `CPUWeight=20`, idle I/O class and
  `MemoryMax=128M`. Measured: about 0.2 CPU-seconds per hour with an empty
  queue and about 0.5 with two running and ten queued jobs, against about 7
  and 20 for 0.6.3.
- **It tracks every job's memory and CPU time** from the job unit's cgroup
  (samples every 3 s while it runs; the unit's own `memory.peak` and CPU time
  at its end), shown by `simq show`, `simq list` and `simq partitions`.
- **It schedules with it, safely:** `--mem auto` derives a request from the
  peaks of earlier jobs of the same name; an admission guard holds a start
  that actual free memory (MemAvailable minus a margin) cannot take, even
  when the reservations would; memory pressure (MemAvailable below a
  high-water mark) stops dispatch and preempts borrow jobs, records which
  job grew, and never stops an owner job. Reservations (`memory_schedulable`)
  stay the primary limit and are never exceeded because of current use.
- **It reports job state to GitHub** for jobs submitted with
  `--github-status` (commit statuses on every transition) or
  `--github-comment` (one pull-request comment at the end, edited in place).

- **Security: the job environment is an allowlist.** Up to 0.6.3 the whole
  environment of `simq submit` was stored in the database and given to the
  job (a workflow's `GH_TOKEN` ended up in `queue.sqlite3`; jobs that
  inherited `RUNNER_TRACKING_ID` were killed by the Actions runner). 0.7.0
  passes only HOME, USER, LOGNAME, LANG, PATH, TZ and SimQueue's own
  `SIMQUEUE_*` settings, plus `-e NAME=VALUE` and `--export NAME,...`
  (`--export ALL`: everything); secret-looking names are never taken from
  the environment and are stored redacted when given with `-e` (see "Job
  environment" below). The migration redacts what older builds stored.

The database moves to schema 9 (columns added, stored secrets redacted;
`upgrade.py` migrates it and its rollback sets the schema back in place).

**Versioning rule.** A minor version (0.6 -> 0.7) is for a change of the
scheduling model, of state or schema, or of side effects outside the
machine; 0.7.0 is all three (measured memory as a scheduling input, schema
9, GitHub posts) and also changes the worker's architecture. A patch
version (0.7.0 -> 0.7.1) is for compatible additions only: no new
scheduling input, no schema change, no new external side effect.

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
| Runtime (sealed, read-only once installed) | `/home/treefor/projects/runtimes/simqueue-0.8.0-20261106` |
| Previous runtime (kept: job units started by 0.7.2 execute its `jobwrap.py`) | `/home/treefor/projects/runtimes/simqueue-0.7.2-20261105` |
| Interpreter (unchanged from 0.5) | `/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen` |
| Code | `source/simqueue/` (`nodeworker.py`, `memory.py` and `github.py` (0.7.0), `scheduler.py`, `borrow.py`, `telemetry.py`, `gpuqueue.py`, `gpuline.py`, `units.py`, `jobwrap.py`, `locks.py`, `nodeconfig.py`, plus the 0.5 modules) |
| Worker wake FIFO (0.7.0) | `~/.local/state/simqueue/worker.wake` (clients and job units write a byte; the worker wakes) |
| Direct GPU use (0.6.2) | `bin/gpu-lease` (runs `source/simqueue/gpuline.py lease`) |
| GPU line state, tickets, log (0.6.2) | `<pic-profile>/simq-gpu-line.json`, `<pic-profile>/gpu-line.d/`, `<pic-profile>/simq-gpu-line.log` |
| GPU claims (0.6.2) | units `simq-job-gpuclaim-<token12>.service`; files `~/.local/state/simqueue/claims/<token>/{spec,control,status}.json`, `claim.log` |
| Node configuration (template) | `config/zenobia-16core-node.json` (Ryzen 9 9950X3D2, live since 2026-11-06; `cpus: "auto"`); `config/zenobia-node.json` (8 cores, the CPU until 2026-11-05) |
| Node configuration (live) | `~/.local/state/simqueue/node.json` (`SIMQUEUE_CONFIG`) |
| Database | `~/.local/state/simqueue/queue.sqlite3` (schema 10 from 0.8.0: `l3_domains`, `split_domains`, `cpu_utilization`; schema 9 from 0.7.0: memory, CPU time and GitHub columns; schema 8 for 0.6.1 to 0.6.3; 0.6.0 used schema 7) |
| Borrow slice | `~/.config/systemd/user/simqborrow.slice` (`CPUWeight=idle`) |
| Per-job files | `~/.local/state/simqueue/jobs/NNNNNN-<token>/{spec,status,stop-post,escaped}.json` (`escaped.json`: 0.7.2, processes seen outside the job's cgroup) |
| Per-job temporary and MPI rendezvous directory (0.7.2) | `/tmp/simq.<job>.<token12>` (`TMPDIR`, `PMIX_SYSTEM_TMPDIR`, `PMIX_SERVER_TMPDIR`, `PRTE_MCA_prte_tmpdir_base`; removed at the job's end) |
| Outside MPI daemons (0.7.2) | `~/.local/state/simqueue/containment.json` (written by the worker's escape scan) |
| Admission ceiling (0.8.0) | `~/.local/state/simqueue/admission.json` (written by the worker; shown by `simq partitions`) |
| Logs | `~/.local/state/simqueue/logs/NNNNNN.log`, `~/.local/state/simqueue/worker.log` |
| Upgrade 0.7.x to 0.8.0 / rollback | `upgrade.py` (dry run unless `--execute`); 0.7.1 to 0.7.2: `../simqueue-0.7.2-20261105/upgrade.py`; 0.7.0 to 0.7.1: `../simqueue-0.7.1-20261105/upgrade.py`; 0.6.3 to 0.7.0: `../simqueue-0.7.0-20261105/upgrade.py` |
| Overhead benchmark (0.7.0) | `source/tools/bench-overhead.py` (0.6.3 against 0.7.0, each worker in its own unit) |
| Switchover from 0.5 (historical) | `switchover.py` |
| Tests | `source/run-tests.py`: `--in-job` inside a queue job (`simq run -c 2 ... -- run-tests.py --in-job`, 0.7.0), `--pure` without units; or holding the cores' lease files as before (`--shared` on cores 6-7) |

## 0.8.0: topology-aware placement on the 16-core node

Built on 0.7.2 (its MPI containment included) for the Ryzen 9 9950X3D2
installed on 2026-11-06 (BIOS 2503): cores 0-15 with SMT siblings 16-31 (c
and c+16), two L3 domains of 96 MiB -- CCD0 = cores 0-7 (logical 0-7,16-23),
CCD1 = cores 8-15 (8-15,24-31) -- and one NUMA node, as read from this
machine's sysfs (`test_topology.py` checks exactly this layout when it runs
on a 9950X3D2). The live 0.7.2 was switched to
`config/zenobia-16core-node.json` the same night (backup of the 8-core
`node.json`: `~/.local/state/simqueue-backups/node-8core-20261106T002055Z.json`).
`upgrade.py` installs 0.8.0 over the installed 0.7.x (0.7.2, 0.7.1 or 0.7.0,
read from the wrapper) with `config/zenobia-16core-node.json`, migrates the
database to schema 10 (`l3_domains`, `split_domains`, `cpu_utilization`;
every existing row kept) and rolls back by setting the schema back to 9 in
place; tested in a sandbox (install, adoption, GPU claim, rollback,
refusals).

What changes on this node (computed with the 0.8.0 planner on the real
sysfs and config; `scheduling-demo.py` in the validation outputs):

| Situation | 0.7.2 (lowest cores first) | 0.8.0 |
|---|---|---|
| empty node, two 4-rank jobs, then a 6-rank job | 0-3, 4-7, 8-13 | 0-3 (CCD0), 8-11 (CCD1: a heavy job goes to the emptiest CCD), 4-7+12-13 (fits no CCD any more) |
| CCD0 has 5 busy cores; a 1-, a 2- and an 8-core job | 5, 6-7, 8-15 | the same: light jobs pack into the fuller CCD, the 8-rank job gets CCD1 whole |
| core 0 busy; an 8-rank job | 1-8 (both CCDs) | 8-15 (CCD1) |
| coupled pilot, 16 cores | 0-15 | `--split-domains 8,8`: 0-15 with `SIMQ_L3_DOMAINS=0-7;8-15` (FE ranks on one CCD, Spheral on the other) |
| 8 ranks spread on purpose | - | `--l3-domains 2`: 0-3,8-11 |
| all 16 cores held by two 8-rank sims; a 4-core `--borrow` build | waits | SMT siblings 28-31 under `cpu.idle`; with the CCD1 sim on `--exclusive-timing`, siblings 20-23 (24-31 stay idle) |
| `-t auto` | - | 1.5 x the p95 runtime of the last 20 succeeded runs, at least 10 min: the full test suite -> 15 min instead of the 2-day default |

**Admission ceiling in `simq partitions`.** The worker publishes the memory
admission ceiling it parks against (`admission.json`, written when it moves
by 128 MiB or every 5 minutes): "Admission ceiling: ~X (a request above it is
parked, not reserved): best MemAvailable with the queue idle over the last
hour ~Y - margin Z, at most memory_schedulable - margin = W; this moment's
memory alone gives ~V". `simq partitions --json` has it under
`node_memory.admission`.

- **Topology from sysfs** (`topology.py`), no hard-coded counts: L3 domains
  from `cpu*/cache/index*/shared_cpu_list` (level 3), SMT siblings from
  `thread_siblings_list`, NUMA nodes from `/sys/devices/system/node/node*/cpulist`.
  The 9850X3D is one domain (96 MiB); the 9950X3D2 two 8-core domains, one
  per chiplet. A partition's `cpus` may be `"auto"`: every physical core
  (lowest logical CPU of each SMT group), so `config/zenobia-16core-node.json`
  serves both CPUs (`memory_schedulable` stays 28G).
- **Placement**: a job stays inside one L3 domain when one has room; a heavy
  job (at least `topology.heavy_cores`, default half a domain) goes to the
  domain with the most free cores, so concurrent heavy jobs spread over the
  caches; a light job to the fullest domain that still fits, so whole
  domains stay free; a job larger than any free domain takes the fewest
  domains. `--l3-domains N` spreads a job evenly over N domains;
  `--split-domains A,B` gives each part its own domain (the coupled pilot: FE
  ranks on one chiplet, Spheral ranks on the other). The job finds its cores
  per domain in `SIMQ_L3_DOMAINS` ("0-7;8-15", in the order of the parts),
  the domains in `SIMQ_L3_DOMAIN_IDS` and its NUMA nodes in
  `SIMQ_NUMA_NODES`. `simq partitions` shows each domain's cores, L3 size
  and occupancy. `topology.l3_placement: false` restores lowest-core-first.
- **Physical cores for sims, siblings for borrowers**: owner jobs get whole
  physical cores (siblings are never allocated to owner jobs). With
  `borrow.lend_siblings` (on in the 16-core config) the SMT siblings go to
  borrow jobs under `cpu.idle`, after every free, leased-idle and owner-idle
  core, unless the owner runs with `--exclusive-timing` (an exclusive start
  preempts a borrower on its siblings).
- **Utilization**: each finished job records `cpu_utilization` = CPU time /
  (runtime x allocated cores); `simq show` flags a job below 50% over at
  least a minute as UNDER-USED, `simq list` has a USE column. `simq
  partitions` measures every core and its sibling over 0.5 s and prints the
  node's used against idle cores; the worker keeps a per-CPU snapshot in
  `queue_meta.cpu_utilization` while jobs run.
- **Backfill starvation** (seen live: an 8-core job held the queue head for
  77 minutes while 6 cores idled, because every other job's time limit was
  too generous to backfill): `-t auto` takes 1.5 x the 95th-percentile
  runtime of the last 20 succeeded jobs of the same name (else name family),
  rounded up to 5 minutes, at least 10, at most the partition's limit, so
  limits are accurate enough to backfill; submit warns when an owner job
  below priority 10 asks for 75% or more of the pool's cores; with
  `borrow.lend_reserved` (on in the 16-core config) borrow jobs use the free
  cores held for the waiting head job whatever their time limit, and are
  preempted by their policy when it starts; the head job's reason reads
  "reserved to start by ...; head of the queue for N min" (`simq queue`).
- **Memory limit by priority** (prepared, not set): a partition may carry
  `priority_max_mem: {"min_priority": 10, "max_mem": "26G"}` so owner jobs
  of that priority may ask for more than `max_mem`. It stays unset until the
  full-machine pilot's memory (16 FE + 16 Spheral ranks) is measured.

Future work (kept in mind, not built): NUMA memory locality
(`cpuset.mems` -- the user manager has no cpuset controller delegated -- or
a memory policy set by `jobwrap.py`), per-NUMA-node memory accounting for
admission, several sockets and EPYC-sized topologies (12-16 CCDs, NPS1/2/4).
Discovery and placement already handle N domains and N nodes generically;
they are tested on one-, two- and four-domain single-node topologies.

### Upgrading 0.7.2 to 0.8.0: `upgrade.py`

```sh
cd /home/treefor/projects/runtimes/simqueue-0.8.0-20261106
./upgrade.py install                 # dry run
./upgrade.py install --execute
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.7.2-to-0.8.0-<UTC> --execute
```

The preflight checks the seal; that `config/zenobia-16core-node.json` is
valid on this machine and equal to the live `node.json` apart from the 0.8.0
keys (`topology`, `borrow.lend_siblings`, `borrow.lend_reserved`,
`partitions[].cpus` "auto", `priority_max_mem`, `smt_siblings` (now read from
sysfs), `containment`, descriptions, notes) without losing a core, a GPU or
a partition name; that the wrapper runs 0.7.x; and that a copy of the
database migrated to schema 10 keeps every row. The install stops the worker
(job and claim units keep running), backs up, switches the wrapper, link,
unit and `node.json`, migrates to schema 10 and verifies the rows, starts
0.8.0 and verifies that every running job unit and GPU claim unit is
adopted. Rollback restores the files byte for byte, sets the schema back to
9 in place (the 0.8.0 columns stay; 0.7.x ignores them) and 0.7.2 adopts
everything back. Job units started by 0.8.0 run its `jobwrap.py`: after a
rollback keep the 0.8.0 runtime until they have ended.

## 0.7.2: MPI containment

**What happened.** The live Raikou work reported that job 993's MPI ranks
appeared in the user slice instead of the job's cgroup. Escaped ranks bypass
the job's MemoryMax, its CPU set, the OOM handling and the kill at its end.
Reproduced on 2026-11-05 with a two-rank job on the live 0.7.1 queue (job
1175): with a PRRTE DVM running outside the queue (Open MPI 5.0.10, `prte
--daemonize` in its own scope), `prun -n 2` inside the job found the DVM
through its PMIx rendezvous files in `/tmp` and the DVM forked both ranks in
*its* scope (`app.slice/mpi-repro-outside-dvm.scope`); `mpiexec -n 2` in the
same job started its own runtime and kept both ranks in the job unit.
`prterun --dvm search|system|pid:N` behaves like `prun`. A job's own
`mpiexec` (prterun) is a PMIx server too, and before 0.7.2 its rendezvous
files sat in `/tmp`, where an outside tool could find it. (Job 993's own
scripts call plain `mpiexec`; which launcher or tool reached an outside
runtime there is not recorded, so 0.7.2 closes the discovery path for every
case and reports what it cannot prevent.)

**Per-job rendezvous.** Every job started by 0.7.2 gets its own directory
`<containment.tmp_base>/simq.<job>.<token12>` (default `/tmp`, the same
tmpfs jobs used for temporary files before; mode 700; made by the job
wrapper before the command starts, removed when every process of the job has
exited, and by the worker if the wrapper was killed). `TMPDIR`,
`PMIX_SYSTEM_TMPDIR`, `PMIX_SERVER_TMPDIR` and `PRTE_MCA_prte_tmpdir_base`
point at it, and `SIMQ_TMPDIR` names it. A DVM search from the job now finds
only runtimes of the same job: `prun` with an outside DVM fails with "prun
failed to initialize, likely due to no DVM being available" instead of
escaping, and `mpiexec` works as before with its session directory inside the
job's directory. A variable given with `-e` wins (`-e TMPDIR=/scratch/x`).
Only an explicit address (`--dvm file:PATH` or `uri:...`) still reaches an
outside DVM. `mktemp`, Python's `tempfile`, compilers and MPI now put their
temporary files in the job's directory, which goes with the job (files a
job leaves in its `TMPDIR` no longer pile up in `/tmp`; a job that must leave
files behind writes them elsewhere).

**Launcher identity variables.** Variables an MPI launcher sets for the
processes it starts (`PMIX_*` but not `PMIX_MCA_*`, `PRTE_*` but not
`PRTE_MCA_*`, `PMI_*`, `OMPI_COMM_WORLD_*`, `OMPI_UNIVERSE_SIZE`,
`OMPI_MCA_orte_hnp_uri`, ...) are dropped from a job's captured environment,
as `SLURM_*` are: inherited (for example with `--export ALL` from a shell
started by `mpiexec`), they make the job's launcher a client of a server
outside the job. A variable given with `-e` is kept; MCA parameters
(`OMPI_MCA_pml`, `PMIX_MCA_*`, `PRTE_MCA_*`) are kept.

**Detection (report only; nothing is killed).** Every
`containment.scan_seconds` (60) while job units run, and at each job's end,
the worker reads `/proc` (no process is spawned; a process's environment is
read once, since it cannot change):

- a process of this user whose environment carries a running job's
  `SIMQ_LAUNCH_TOKEN` but whose cgroup is not inside the job's unit (moved
  out with `systemd-run --scope` or a cgroup write) is recorded in the job's
  `escaped.json`, as an `escaped-process` event, in the worker log ("job N:
  process P (comm) carries the job's launch token but runs outside its
  cgroup ...") and in `simq show` (`ESCAPED`); one still alive after the job
  ended is marked "after the job ended";
- an MPI launcher daemon (`prte`, `prted`, `prterun`, `prun`, `orted`,
  `hydra_pmi_proxy`, `mpiexec.hydra`, ...) of this user outside every job's
  cgroup is recorded in `containment.json` of the state directory, in the
  worker log, in `simq partitions` ("OUTSIDE MPI DAEMON: prte pid N in ...")
  and in `simq doctor` ("Containment").

Ranks that an outside DVM forks do **not** carry the job's environment (PRRTE
forwards none of it, checked with `SIMQ_LAUNCH_TOKEN`, `OMPI_MCA_*`,
`PMIX_MCA_*` and `PRTE_MCA_*`), so they cannot be attributed to a job by its
token; they show up as children of the outside daemon. A scan of about 550 processes
takes 0.4 to 2 ms of CPU (3 to 5 ms the first time, when every environment is
read once), so at most about 0.13 CPU-s per hour at the default 60 s
(`containment-scan-cost.json` in the test receipts).

**Configuration** (`containment` in node.json, all optional):
`enabled` (true), `tmp_base` ("/tmp"), `scan_seconds` (60; 0 = only at each
job's end). Running jobs started by 0.7.1 keep their environment (`/tmp`)
until they end.

**Tests** (`test_containment.py`): the identity filter and the rendezvous
variables, the private directory (mode, reuse, refusal of a foreign or
symlinked one, removal that never follows symlinks), the scan against real
processes (token outside and inside its cgroup, left behind, an outside
`prte` with its children, cost), a job unit's environment and directory
(removed at the end), a process moved out of a job with `systemd-run --scope`
(reported while it runs and after the job ended), and the two-rank
reproduction: an outside DVM whose rendezvous directory the job shares (as
0.7.1 jobs shared `/tmp`) gets both `prun` ranks of a job without its own
rendezvous variables, in the DVM's cgroup; with them (and a second stray
DVM in the real `/tmp`) `prun` finds no DVM and nothing runs outside;
`mpiexec -n 2` keeps both ranks in the job's unit; the worker lists both DVMs
in `containment.json` and `simq partitions`. The shared directory is not the
real `/tmp` because there other jobs' servers may be found too: in the first
official 0.7.2 run (job 1198) `prun` from the test refused with "PMIx has
found multiple possible servers" -- the test's DVM and another server in
`/tmp`, most likely the two-rank `pore_plasma` run of live job 1063
(17:11-17:17 UTC): a 0.7.1 job's `mpiexec` puts its rendezvous files in
`/tmp`, where any `prun` on the node can find it. The reproduction needs an Open MPI 5 prefix with `prte`, `prun`
and `mpiexec` (`SIMQ_TEST_OMPI_PREFIX`, else the MOOSE conda environment) and
is skipped without one. Two older borrow tests were made independent of
where they run: the placement check names the partition of the free core the
borrower got, and the preemption, owner-idle and leased-idle tests turn off
the per-CPU idle guard (`cpu_idle_threshold` 1.0; their idleness still comes
from the owner's or holder's own cgroup): inside a queue job the two test
cores also carry the harness and unpinned work the busy node schedules
there, so they were often above 15% for long stretches and a borrower waited
past the timeouts (jobs 1230 and 1248).

### Upgrading 0.7.1 to 0.7.2: `upgrade.py`

```sh
cd /home/treefor/projects/runtimes/simqueue-0.7.2-20261105
./upgrade.py install                 # dry run
./upgrade.py install --execute
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.7.1-to-0.7.2-<UTC> --execute
```

No schema change: the preflight checks the seal, that the new
`config/zenobia-node.json` equals the live `node.json` apart from the new
`containment` section and the notes, that the wrapper runs 0.7.1, and that a
copy of the database opened by 0.7.2 keeps every row; it lists MPI launcher
daemons running outside the queue. The install stops the worker (job and
claim units keep running), backs up, switches the wrapper, link, unit and
`node.json`, starts 0.7.2 and verifies that every running job unit and GPU
claim unit is adopted; rollback restores the files byte for byte and 0.7.1
adopts everything back. Job units started by 0.7.2 run its `jobwrap.py`:
after a rollback keep the 0.7.2 runtime until they have ended.

## 0.7.1

- **CRITICAL: a job unit can no longer take down the user manager.** On
  2026-11-05 at 09:24 UTC a test job unit (OOMPolicy=kill, so
  `memory.oom.group` 1) hit its MemoryMax while the user manager had just
  spawned a process into its cgroup with vfork semantics (posix_spawn,
  CLONE_VM: until exec the child shares the manager's mm). The kernel's
  group kill (`mm/oom_kill.c`: `oom_kill_memcg_member` ->
  `__oom_kill_process`) SIGKILLs every process sharing the victim's mm, so
  it killed `user@1000.service`'s manager, and with it the desktop, SimQueue,
  the runners, docker and every agent. Live 0.7.0 job units had the same
  exposure (OOMPolicy=kill and `ExecStopPost=jobwrap.py stop-post` running in
  the job's cgroup after the job, possibly while MPI shared memory left in
  `/dev/shm` still held the cgroup at its limit). 0.7.1:
  - job units have **no ExecStopPost, ExecStartPost or ExecStop**: the
    manager spawns into a job cgroup only for ExecStart, while it is new and
    empty. The worker records a job's end from the wrapper's own
    `status.json` and, for a killed job, from the unit's properties
    (`Result`, `ExecMainStatus`, `MemoryPeak`, `CPUUsageNSec`): units are no
    longer `--collect`ed, so a failed unit stays loaded until the worker has
    read it and runs `systemctl --user reset-failed`. A unit that ended
    normally is unloaded at once, so its result (`success`, `exit-code`,
    `signal`) follows from the wrapper's exit status. The wrapper wakes the
    worker itself when the command ends;
  - **OOMPolicy=stop** (`memory.oom.group` 0): `oom_kill_process` kills the
    one victim `oom_badness` picks -- which never picks an `in_vfork` task --
    and only the processes sharing that victim's mm (job processes); systemd
    then stops the rest of the unit and records `Result=oom-kill`;
  - after a job ends the worker **removes its leftover MPI shared-memory
    segments** from `/dev/shm` (this user's, MPI names such as
    `sm_segment.*` and `vader_segment.*`, made after the job started, mapped
    or open in no process), so their shmem no longer stays charged;
  - tests never drive a unit into a real OOM: the OOM path is tested with a
    mocked unit state (`test_oom_safety.py`); the real-OOM tests run only with
    `SIMQ_TEST_REAL_OOM=1`, which `run-tests.py` refuses on zenobia; a test
    checks the job units' properties (no Exec*Post/Stop, OOMPolicy=stop,
    `memory.oom.group` 0, a failed unit recorded and unloaded);
  - the install removes the interim drop-in
    `~/.config/systemd/user/simq-job-.service.d/50-oom-policy-stop.conf`
    (backed up; a rollback restores it). Job units started by 0.7.0 keep
    their ExecStopPost until they end (the drop-in gives them
    OOMPolicy=stop meanwhile).
- **Orphans are reaped while the job runs.** The job wrapper is a
  subreaper, so processes a job orphans (daemonized helpers, MPI launchers)
  become its children; 0.7.0 and earlier waited only for the leader, so
  those that exited stayed zombies until the leader ended. The 0.7.1 wrapper
  reaps any child as it exits (`test_oom_safety.py` checks it with a
  double-forked process).
- **Memory without `--mem`.** With `memory.auto_default` (on for zenobia) a
  job submitted without `--mem` gets `--mem auto` when its name (or name
  family) has a recorded peak, else the partition's or alias's default (4G).
  `simq show` says "Memory request auto (default without --mem): ...". An
  explicit `--mem` always wins, and `--mem auto` works as in 0.7.0. The
  idempotency fingerprint of a submission without `--mem` does not depend on
  the history.
- **`spheral` alias:** `default_mem: "auto"` (it was 18G): auto from the
  history even without `memory.auto_default`, else the pool's 4G. `simq
  partitions` shows "auto (else 4G)". Any alias may say so.
- **Memory head-of-line fix** (seen live on 0.7.0: an 8-core 20G job held
  the queue head for hours because the admission guard could never admit
  it, while about 15.8G was admissible): the planner now accounts *actual*
  memory as a resource next to the reservations -- the starts of a pass may
  take MemAvailable minus the margin minus what young jobs have not
  allocated, and a running job's end gives back what it uses now -- so a head
  job that does not fit gets a reservation from the memory running jobs will
  free, and backfill works against it ("memory busy: needs 12G, 9.8G actually
  free now; reserved to start by ..."). A job whose request exceeds what the
  node can ever admit (the best MemAvailable seen with the queue's jobs
  taken away over the last hour, or `memory_schedulable`, minus the margin)
  is **parked**: left out of planning, so it holds no reservation, with the
  reason "parked: mem 20G exceeds admissible ~15.8G on zenobia (...)", and
  looked at again on every pass. `simq submit` warns when the request
  exceeds what the node can admit now; `simq list` shows the reason.
- **Tests:** the borrow placement check accepts any free test core; the CPU
  priority check takes the best of three 2 s windows; the ensemble launch-gate
  members have tails long enough for a loaded node; owner-idle lending gets
  60 s. `run-tests.py` removes every inherited `SIMQUEUE_*` variable (inside
  a queue job they name the live installation) and points `XDG_STATE_HOME`
  and `SIMQUEUE_STATE_DIR` into the test directory: the first official 0.7.1
  run (job 1057) failed 43 tests because the 0.5-mode and ensemble tests ran
  against the live `node.json`, saw every live lease held and never started
  their jobs (nothing reached the live queue). The second (job 1114) failed
  only where the tests met the job's own surroundings: the exclusive-mode
  fixture built overlapping P/E sets from a 2-CPU affinity (now one CPU each),
  two restart tests counted the killed orphan's zombie (its parent, the live
  0.7.0 wrapper, reaps only at the end) as alive (a zombie now counts as
  exited), and `test_shield_helper.py` cannot be imported (as before; the
  official runs list the suites without it).
- **Older backups:** the secret-pattern environment values in the database
  copies of the four older upgrade backups were redacted in place on
  2026-11-05 with the migration's redaction (7, 60, 205 and 705 job rows; 0
  left); their schemas and every other column are unchanged, so they stay
  usable for the rollbacks they belong to.
- **Idempotency:** as documented for 0.7.0, a key first submitted under 0.6.x
  is refused on replay ("already bound to a different submission request")
  instead of running twice.

### Upgrading 0.7.0 to 0.7.1: `upgrade.py`

```sh
cd /home/treefor/projects/runtimes/simqueue-0.7.1-20261105
./upgrade.py install                 # dry run
./upgrade.py install --execute
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.7.0-to-0.7.1-<UTC> --execute
```

No schema change: the preflight checks the seal, that the new
`config/zenobia-node.json` equals the live `node.json` apart from
`memory.auto_default`, the aliases' `default_mem` and the notes, that the
wrapper runs 0.7.0, and that a copy of the database opened by 0.7.1 keeps
every row. The install stops the worker (job and claim units keep running),
backs up, switches the wrapper, link, unit and `node.json`, starts 0.7.1 and
verifies that every running job unit and GPU claim unit is adopted; rollback
restores the files byte for byte and 0.7.0 adopts everything back (it also
makes GitHub posts that are still due).

### 0.7.0 install history (2026-11-05)

At 06:16:55 UTC a 0.7.0 test run migrated the live database to schema 9 (the
ensemble tests' supervisors briefly ran without `SIMQUEUE_STATE_DIR` and
fell back to the live state directory): new 0.6.3 clients were refused
until 0.7.0 was installed at 06:34 UTC (`--allow-unsealed`; backup
`0.6.3-to-0.7.0-20261105T063438Z`). The cause is fixed (`SIMQUEUE_*` is
allowlisted; the test runner isolates the state directory).

## 0.7.0: event-driven worker, memory, GitHub reporting

### The event-driven worker

Up to 0.6.3 the worker woke ten times a second and made a full pass every
`poll_seconds` (1 s), each starting one `systemctl --user show` for the
running jobs. 0.7.0 sleeps in `select()` on the FIFO
`~/.local/state/simqueue/worker.wake` and wakes for:

- **clients**: `simq submit`, `run`, `cancel`, `retry`, `promote`, `demote`,
  `priority`, `pause`, `resume` and the ensemble commands write a byte to the
  FIFO after their database change (never blocking; without a worker nothing
  happens);
- **job ends**: the 0.7.0 job wrapper's `ExecStopPost` hook writes a byte;
  the worker then looks for the unit's cgroup to disappear every 0.2 s for
  5 s and reconciles the job at once;
- **other writers**: every 3 s while jobs run (5 s otherwise) it compares the
  database's `PRAGMA data_version` on its own persistent connection, which
  changes only when another connection committed (an ensemble controller, an
  older client);
- **timers**: a pass every `launch.poll_seconds` (1 s) only while something
  must be watched -- a job starting or being stopped, borrow jobs (the idle
  telemetry), the GPU line, queued jobs held by leases outside the queue or
  by the memory guard, memory pressure -- else at most every
  `launch.idle_seconds` (30 s); every 10 idle intervals (5 min) every job
  unit is queried from systemd as a safety net.

Between passes it only reads files: each running job's cgroup directory
(`memory.current`, `memory.peak`; its disappearance is how the end of a job
started by an older wrapper is noticed) and `/proc/meminfo`. A pass queries
systemd only for the jobs that need it (starting, being stopped, cgroup
gone, not verified yet, borrow job whose freezer state differs from the
record). Dispatch latency: a submission is planned within milliseconds of
the client's write; a job's end is reconciled about 0.2 s after its unit is
gone (adopted 0.6.3 jobs: within 3 s).

The worker log's last line on stop counts the passes and the systemd
queries, for example `node worker stopped ...; 13 pass(es), 2 systemd unit
quer(ies)`.

Measured overhead (`source/tools/bench-overhead.py`: each worker in its own
transient unit with the properties it is installed with, cgroup `cpu.stat`
over 10 minutes after a 30 s warm-up, all four at once; busy = a two-core
pool with two running and ten queued `sleep` jobs; sandbox config like
zenobia's with lease files, a GPU with the GPU line and borrowing on):

| Worker | Queue | CPU-s per hour | Share of one core | RSS | cgroup memory (peak) |
|---|---|---|---|---|---|
| 0.6.3 | idle | 7.08 | 0.197% | 31.1 MiB | 18.7 MiB (19.2) |
| 0.6.3 | busy | 19.46 | 0.540% | 31.1 MiB | 18.9 MiB (20.5) |
| 0.7.0 | idle | 0.18 | 0.005% | 31.3 MiB | 19.4 MiB (20.3) |
| 0.7.0 | busy | 0.69 | 0.019% | 31.5 MiB | 19.8 MiB (21.1) |

The live 0.6.3 worker measured 11.9 CPU-seconds over 77 minutes (9.3 per hour, MemoryCurrent 21 MB) with real jobs; its `systemctl show` calls also cost the user manager CPU, which these figures do not include. Receipt: `<outputs>/simqueue-0.7.0-20261105/bench/bench-overhead.json`. One sample of eight job cgroups and MemAvailable costs about 0.14 ms (`test_memory.py`), 0.16 CPU-s per hour at the 3 s interval.

### Memory tracking

| What | Where |
|---|---|
| current memory of a running job | read from its cgroup every `memory.sample_seconds` (3 s); written to the database every `memory.record_seconds` (60 s) as `mem_current_bytes`; `simq show` and `simq partitions` read the cgroup directly |
| history | min/max per minute (`mem_series_json`), at most `series_buckets` (240) entries: older buckets merge in pairs when a job runs longer, so a job's whole history stays bounded |
| peak | `mem_peak_bytes`: the unit's own `memory.peak` recorded by the 0.7.0 job wrapper when the job ends (source `unit`), else the highest `memory.peak` the worker sampled (`sampled`, for jobs adopted from 0.6.3); over all runs of a requeued job |
| CPU time | `cpu_usage_usec`: the unit's `cpu.stat` usage at its end (all processes of the job) |
| node | `/proc/meminfo` MemTotal and MemAvailable at every sample |

`simq show ID` adds "Memory now" (running: live cgroup reading and its share
of the limit), "Memory peak" (and its source), "CPU time" (and average
cores), "Memory history" (a sparkline of the per-minute maxima) and, for
`--mem auto`, how the request was derived. `simq list` has a MEM column:
now/request while running, peak/request after. `simq partitions` prints the
node's available memory, the high- and low-water marks, requested against
actual memory of every running job and the memory used outside the queue,
and MEMORY PRESSURE while it lasts (`--json`: `node_memory`).

A job using at least `near_limit_fraction` (90%) of its limit is logged once
and gets a `memory-near-limit` event (`simq show` History): an
out-of-memory kill is likely.

### `--mem auto`

`simq submit|run --mem auto` (also `ensemble submit`) takes `auto_factor`
(1.25) times the largest peak of the last `auto_history` (5) finished jobs
with the same name -- else the same *name family*, the name without trailing
index tokens (`fe-audit-a0-n4` -> `fe-audit`,
`lv-interp_300K_xe1e-6` -> `lv-interp`, `pore1p5-cpu-a` -> `pore1p5-cpu`) --
preferring jobs with the same core count, rounded up to `auto_round` (256M)
and bounded by the partition's `min_mem` (256M) and `max_mem` (20G). A job
killed at its memory limit counts with its limit (it needed more), so the
next auto request grows. Without any history: the partition's (or alias's)
default. The request is fixed at submission (`simq show`: "Memory request
auto: 1.25 x largest peak 2.9G of 3 earlier job(s) with the same name");
an idempotent replay with `--mem auto` finds the same job however the
history has grown since.

Over-request warning: `simq submit --mem SIZE` prints `simq: warning: --mem
18G is 6.2x the largest recent peak ...` (and `warnings` in `--json`) when
the request is at least `warn_factor` (2x) and `warn_min_excess` (2G) above
the largest recent peak of the name (family). The job is queued as asked.

### The admission guard

Reservations stay the primary constraint: the planner never lets the sum of
`--mem` of running jobs exceed `memory_schedulable`, whatever the jobs use
now (usage grows). On top of that, a start the planner chose is held when

    request > MemAvailable - margin - (unallocated requests of jobs started
                                       less than ramp_seconds ago)

`margin` is `memory.admission_margin`, by default the high-water mark (a job
that then grows to its full request does not by itself cause memory
pressure); jobs started within `ramp_seconds` (60 s) count their unallocated
request as used, so a burst of starts cannot all see the same MemAvailable.
Processes outside the queue (the desktop, other sessions, builds run
directly) are thereby accounted. A held job keeps its place: smaller jobs
behind it may still start, and since the reservations already include it,
they cannot push it out; it starts when memory outside the queue is freed.
Its reason reads `pending: memory: needs 18G but MemAvailable 12.0G - margin
2.2G leaves 9.8G`. Starts are re-checked every `poll_seconds` while one is
held.

### Memory pressure

Every sample compares MemAvailable with the high-water mark
`max(high_water_available, high_water_percent of MemTotal)` (2G or 7%: 2.2G
on zenobia). Below it:

1. dispatch stops (queued jobs read `pending: memory pressure since ...`);
2. every running borrow job is preempted by its `--on-preempt` policy
   (freeze stops its growth at once; requeue and cancel free its memory); a
   job frozen for pressure is requeued when the pressure lasts
   `pressure_escalate_seconds` (120 s), because a frozen job keeps its
   memory;
3. the worker log, the jobs' History (`memory-pressure` event: current
   memory, request, growth over the last 60 s, culprit) and `simq partitions`
   record which job grew and by how much against its request, and how much
   memory is used outside the queue;
4. when MemAvailable is back at or above the low-water mark
   `max(low_water_available, low_water_percent)` (4G or 12%: 4G), dispatch
   resumes and frozen borrow jobs are thawed when their cores are lendable.

Owner jobs are never stopped for memory pressure: their own `MemoryMax`
bounds each of them; they are recorded and warned about.

### Job environment (security fix)

What a job sees from its submitter: HOME, USER, LOGNAME, LANG, PATH, TZ and
the `SIMQUEUE_*` settings (so a job running `simq` stays on the same queue),
the names given with `--export NAME,...` (`--export ALL`: the whole
environment), and `-e NAME=VALUE`. The worker adds the user manager's
`XDG_RUNTIME_DIR` and `DBUS_SESSION_BUS_ADDRESS` (so `systemctl --user` works
in a job) and the `SIMQ_*` placement variables.

A name matching `*TOKEN*`, `*SECRET*`, `*PASSWORD*`, `*KEY*`, `GITHUB_*`,
`ACTIONS_*` or `RUNNER_*` (case-insensitive; SimQueue's own `SIMQ_*`
excepted) is never taken from the environment, whatever `--export` says;
`simq submit` lists such names it left out (names only). Given with `-e`, it
reaches the job, but the database stores `<redacted>`: the value waits in
`~/.local/state/simqueue/secrets/<job>.json` (mode 600) until the job ends,
when the worker deletes that file and redacts the job's `spec.json`. A
`simq retry` of such a job runs without it (pass it again with a new
submission). `simq show` and `simq list` print names, never values.

The schema-9 migration (`upgrade.py install`) redacts secret-pattern values
in every stored job environment and in finished ensemble members; a queued
job's secret that had been given with `-e` moves to its secrets file. The
install also redacts the `spec.json` of every finished job and its own
backup copy of the database. Dry-run counts on 2026-11-05 (never values;
`<outputs>/simqueue-0.7.0-20261105/env-scrub-dry-run.json`): the live
database had no secret-pattern values left; 859 of 864 job spec files had
some (857 `*TOKEN*`, 156 `GITHUB_*`, 32 `RUNNER_*`, 4 `ACTIONS_*` values);
the database copies in the four older upgrade backups hold 7, 60, 205 and
705 such rows (not touched by the install).

CI scripts that submit through `env -i` keep working unchanged.
Idempotency fingerprints cover the passed environment only, so a replay of
a key that was submitted under 0.6.x is refused ("already bound to a
different submission request") rather than run twice: 0.6.x fingerprinted
the whole environment, including `PWD`, `SHLVL` and `LC_CTYPE`.

### Configuration (`memory` in node.json)

| Key | Default | Meaning |
|---|---|---|
| `enabled` | true | false turns tracking, the guard and pressure handling off |
| `sample_seconds` | 3 | cgroup and MemAvailable sampling while jobs run (0.5-60) |
| `record_seconds` | 60 | database writes of the samples |
| `series_buckets` | 240 | bound of each job's min/max history |
| `admission_margin` | null | MemAvailable kept free by a start; null = the high-water mark |
| `ramp_seconds` | 60 | young jobs count their unallocated request |
| `high_water_available`, `high_water_percent` | 2G, 7 | pressure below the larger |
| `low_water_available`, `low_water_percent` | 4G, 12 | pressure over at the larger |
| `pressure_escalate_seconds` | 120 | frozen-for-pressure borrow job requeued after |
| `near_limit_fraction` | 0.9 | near-limit warning |
| `auto_factor`, `auto_history`, `auto_round` | 1.25, 5, 256M | `--mem auto` |
| `warn_factor`, `warn_min_excess` | 2, 2G | over-request warning |
| `test_node_total` | null | tests only: a simulated node of this size used only by queue jobs |

Partitions take `min_mem` (default 256M); `launch` takes `idle_seconds`
(default 30).

### GitHub reporting

The worker, not the job, posts:

```sh
simq run -c 2 --mem auto -t 1:00:00 -n rz-pr-tests \
    --github-status treeforjs/plasma-solid-shocks@$SHA:simqueue/rz-pr-tests \
    -- ./run-pr-tests.sh
simq submit -c 8 -n pore-pilot --github-status treeforjs/plasma-solid-shocks@$SHA:simqueue/pilot \
    --github-comment treeforjs/plasma-solid-shocks#123 --github-summary pilot-summary.md -- ./pilot.sh
```

- `--github-status OWNER/REPO@SHA:CONTEXT` (repeatable; full commit SHA):
  a commit status for every transition -- `pending` "SimQueue job N queued
  (C cores, M)", `pending` "SimQueue job N running on CPUs X", then
  `success` (exit 0) / `failure` (non-zero exit) / `error` (canceled, time
  limit exceeded, out of memory at its limit, preempted and canceled, lost
  unit), with runtime and peak memory, at most 140 characters. A job that
  starts before the worker's first post shows running directly; a borrow
  job requeued after preemption is pending again.
- `--github-comment OWNER/REPO#PR` (opt-in): one comment at the end with the
  job ID, state, runtime, peak memory, cores, the last `log_lines` (20) log
  lines and the job's summary.
- `--github-summary PATH`: a markdown file the job may write; exported to
  the job as `SIMQ_SUMMARY` (relative paths are resolved against the job's
  working directory); appended to the comment (at most `summary_max_bytes`).

**Notifications.** Commit statuses never notify anyone. A comment can only
notify when it is first created; it carries the hidden marker
`<!-- simq-job:N -->` and every later post for the job finds it and edits it
(PATCH), which does not notify. Nothing SimQueue writes can mention anyone:
in descriptions, comments and summaries every `@` and every `#` before a
digit gets a zero-width space, so a log line or a summary cannot ping a
person or team or cross-reference an issue (`#123`, `owner/repo#123`).
GitHub does not email the token's owner about their own comments unless
"Include your own updates" is enabled in that account's notification
settings, so with the default settings the token's owner gets no email from
SimQueue's reporting at all.

**Token.** `github.token_file` (a file holding a token, mode 600) if set,
else `gh auth token` run by the worker at post time; kept only in the
posting thread's memory for `token_cache_seconds` (600) and sent only in the
Authorization header. It is never written to the database, the logs, the
job's environment or any unit file (tested). `gh` is found on the user
manager's PATH (`/usr/bin/gh`; `gh auth status`: logged in as treeforjs,
keyring).

**Never blocking.** Submissions record the targets in the job row
(`github_json`); the worker compares each job's phase (queued, running,
final) with the last one posted after every pass and hands due posts to one
daemon thread of the worker process (no extra process), which makes the
HTTP requests with a `timeout_seconds` (10 s) timeout and hands results back
through a queue. A failure (GitHub unreachable, HTTP 5xx, rate limit,
missing token) is recorded (`simq show`: "GitHub error", History
`github-failed`) and retried with exponential backoff from
`backoff_seconds` (15 s) up to `backoff_max_seconds` (15 min), at most
`max_attempts` (12) times per phase; a permanent error (HTTP 4xx such as a
wrong SHA or repository) is recorded and not retried. Scheduling never waits
for any of it. Posts still due when the worker stops are made by the next
worker.

Configuration (`github` in node.json): `enabled` (true), `api_base`
(`https://api.github.com`), `token_file` (null), `gh` (`gh`),
`timeout_seconds`, `token_cache_seconds`, `backoff_seconds`,
`backoff_max_seconds`, `max_attempts`, `log_lines`, `summary_max_bytes`.
`enabled: false` makes the submit options an error.

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
| `SIMQ_SUMMARY` (0.7.0) | with `--github-summary`: the markdown file the job may write for its PR comment |
| `SIMQ_TMPDIR` (0.7.2) | the job's own temporary and MPI rendezvous directory `/tmp/simq.<job>.<token12>` (mode 700, removed at the job's end) |
| `TMPDIR`, `PMIX_SYSTEM_TMPDIR`, `PMIX_SERVER_TMPDIR`, `PRTE_MCA_prte_tmpdir_base` (0.7.2) | `SIMQ_TMPDIR` (a value given with `-e` wins); see "0.7.2: MPI containment" |
| (0.7.0) everything else | only HOME, USER, LOGNAME, LANG, PATH, TZ, `SIMQUEUE_*`, `--export` names and `-e` values from the submitter, plus XDG_RUNTIME_DIR and DBUS_SESSION_BUS_ADDRESS; see "Job environment" |
| `CUDA_VISIBLE_DEVICES` | the allocated GPU UUIDs (empty when none); `CUDA_DEVICE_ORDER=PCI_BUS_ID` |

The rest of the environment is the submitting shell's, captured at submit time
(as in 0.5), except that inherited `SLURM_*` variables are removed unless given
with `-e`, and (0.7.2) so are inherited MPI launcher identity variables
(`PMIX_*` but `PMIX_MCA_*`, `PRTE_*` but `PRTE_MCA_*`, `PMI_*`,
`OMPI_COMM_WORLD_*`, ...).

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

`~/.config/systemd/user/simqueue.service` (as `upgrade.py` installs it for 0.7.0):

```ini
[Service]
ExecStart=/home/treefor/.local/bin/simq worker --wait-for-lock
Environment=SIMQUEUE_PROJECT_ROOT=/home/treefor/projects
Environment=SIMQUEUE_STATE_DIR=/home/treefor/.local/state/simqueue
Environment=SIMQUEUE_CONFIG=/home/treefor/.local/state/simqueue/node.json
Nice=10
CPUWeight=20
IOSchedulingClass=idle
MemoryHigh=96M
MemoryMax=128M
Restart=on-failure
KillMode=mixed
UMask=0077
```

There is no `CPUAffinity`: the worker places each job unit explicitly. Up to
0.6.3 it polled (one `systemctl show` per second while jobs ran); 0.7.0 is
event-driven (see below). Job units and GPU claim units are started by the
user manager, not by the worker, so they inherit none of the worker's
priority or memory limit. `~/.local/bin/simq` exports the same three
variables and runs the runtime with the pinned interpreter.

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

### Upgrading 0.6.3 to 0.7.0: `upgrade.py`

```sh
cd /home/treefor/projects/runtimes/simqueue-0.7.0-20261105
./upgrade.py install                 # dry run: preflight, migration rehearsed on a copy, nothing changed
./upgrade.py install --execute       # stop worker, back up, switch files, migrate, start, verify
./upgrade.py rollback --backup ~/.local/state/simqueue-backups/0.6.3-to-0.7.0-<UTC> --execute
./upgrade.py status
```

Install (see the docstring of `upgrade.py` for every step): the preflight
checks the seal, that `config/zenobia-node.json` equals the live
`node.json` apart from the 0.7.0 keys (`memory`, `github`,
`launch.idle_seconds`, `partitions[].min_mem`, notes; `default_mem` stays
4G), that the layout and every recorded partition name resolve, that the
wrapper runs 0.6.3 and the database is schema 8, and migrates a copy of the
database with the 0.7.0 code to prove that every schema-8 column of every
row stays as it is. Then it waits for no job to be starting, stops
`simqueue.service` (job and claim units keep running), backs up the
database (sqlite3 backup API), wrapper, link, unit files, `node.json` and
`worker.log`, installs the 0.7.0 wrapper, link, unit (`Nice=10`,
`CPUWeight=20`, `IOSchedulingClass=idle`, `MemoryHigh=96M`,
`MemoryMax=128M`) and `node.json`, migrates the database to schema 9 and
checks the rows again, starts the worker and verifies that it adopted every
running job unit and GPU claim unit, published the GPU line state and
created the wake FIFO, and that `simq doctor` and `simq partitions` run. A
failed verification rolls back by itself.

Rollback restores the files byte for byte, keeps a copy of the schema-9
database in the backup directory, sets the schema back to 8 **in place**
(rows, including jobs that 0.7.0 queued or ran, stay; the new columns stay
and 0.6.3 ignores them), starts 0.6.3 and verifies that it adopted every
running job unit and GPU claim unit. Jobs whose GitHub posts are still due
are listed: 0.6.3 does not post; installing 0.7.0 again posts them.

During the switch: a `simq` started before the wrapper changed that opens
the database after the migration (a few seconds) fails with "schema 9 is
newer"; retry it. Supervisors of running ensemble jobs (0.6.3 code) opened
the database before and keep working. Job units started by 0.6.3 keep
running its `jobwrap.py`: keep the 0.6.3 runtime until they have ended
(their memory is sampled by the 0.7.0 worker; their peak is the highest
sample, source `sampled`). The first 0.7.0 worker clears any stale memory
pressure record.

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

## Validation of 0.7.2 (2026-11-05)

- Full suite inside live queue job 1254 (`simq run -c 2 --mem 1G -t 0:40:00 --priority 20 -- run-tests.py --in-job`, the test units on the job's cores [0, 7]): 391 tests, 0 failure(s), 0 error(s), 12 skipped, 541 s.
- Real-OOM tests not run (2 skipped: `SIMQ_TEST_REAL_OOM` unset, and refused on zenobia); the OOM path ran against mocked unit states (`test_oom_safety.py`), and the unit-property test checked live-started test units (OOMPolicy=stop, no Exec*Post/Stop, `memory.oom.group` 0).
- Every suite without units, with the live-database copy: 375 tests, 0 failures (64 skipped: the unit-level tests).
- Historical migration suites: 7 tests, 0 failures.
- Borrower CPU share under contention (three windows): [0.07773934788254602, 0.07814236119953473, 0.07864425399262283].
- Escape scan: 556 processes, 3.4 ms CPU the first time, 0.4 ms after.
- Two-rank MPI reproduction (`test_containment.py`): ran (Open MPI 5.0.10 prefix of the MOOSE environment).
- `test_shield_helper.py` cannot be imported, as before.
- Receipts: `/home/treefor/projects/.workspace/outputs/simqueue-0.7.2-20261105/receipts/`.

## Validation of 0.7.1 (2026-11-05)

- Full suite inside live queue job 1123 (`simq run -c 2 --mem 1G -t 0:40:00 --priority 20 -- run-tests.py --in-job`, the test units on the job's cores [0, 1]): 379 tests, 0 failure(s), 0 error(s), 12 skipped, 594 s.
- Real-OOM tests not run (2 skipped: `SIMQ_TEST_REAL_OOM` unset, and refused on zenobia); the OOM path ran against mocked unit states (`test_oom_safety.py`), and the unit-property test checked live-started test units (OOMPolicy=stop, no Exec*Post/Stop, `memory.oom.group` 0).
- Every suite without units, with the live-database copy: 363 tests, 0 failures (61 skipped: the unit-level tests).
- Historical migration suites: 7 tests, 0 failures.
- Borrower CPU share under contention (three windows): [0.1179767492535054, 0.11677841506968292, 0.117988022014696].
- `test_shield_helper.py` cannot be imported, as before.
- Receipts: `/home/treefor/projects/.workspace/outputs/simqueue-0.7.1-20261105/receipts/`.

## Validation of 0.6.3 (2026-11-04)

Built from 0.6.2 on the Power-Flow branch `simqueue-0.6.3`; not installed
(`installed: false` in `release.json`). Receipts and logs in
`/home/treefor/projects/.workspace/outputs/simqueue-0.6.3-20261105/`.

- Full suite through `run-tests.py` on cores 6-7 (`--wait-any --min-cores 2
  --shared`; it waited while a Spheral run held all eight core leases):
  299 run, 0 failures, 1 error, 6 skipped (`receipts/full-suite-1.json`).
  The error was a test bug in `test_upgrade.py` (a helper shadowed by a
  fixture attribute) after the install, the adoption of the running 0.6.2
  job and the pool run of the job queued under the old name had passed.
  After the fix `test_upgrade.py` ran 10/10 (`receipts/suite-upgrade-3.json`;
  the rerun before it, `suite-upgrade-2.json`, showed that 0.6.2 sets up its
  own GPU claim right after a rollback, which the test now expects).
- Pure: `test_exclusive_mode.py` (9, needs four CPUs in the affinity mask)
  and `test_pool_aliases.LiveCopyTests` on a read-only copy of the live
  database (`live-queue-copy-20261104.sqlite3`, 715 jobs recorded as NULL,
  spheral, raikou-mfem or powerflow: every one resolves to the pool, plans
  without "unknown partition", and opening it with 0.6.3 changes no row):
  10/10 (`receipts/pure-exclusive-and-live-copy.json`).
- Per file: 308 tests, 303 ok, 5 skipped (the two historical live-copy
  migration tests, which expect a schema 6/7 live database; three
  `test_node_worker` tests that need two-core partitions or pf-launch 1.0),
  0 failed; `test_shield_helper.py` (16) cannot be imported, as before.
- New: `test_pool_aliases.py` (24: alias parsing and validation, the zenobia
  pool and its 16-core variant, lease paths, planning in the pool with
  preferred cores, an 8-core owner job reserving the pool, borrow jobs next to
  an 8-core owner job (owner-idle lending, preemption, backfill against a
  waiting owner with and without `backfill_waiting_owners`), `--borrow-from`
  through aliases, jobs recorded under old names, CLI submissions under
  aliases, `simq run` exit statuses) and `test_run_units.py` (6 with real
  units: jobs recorded under the old partition names before the switch are
  adopted and scheduled in the pool; a borrow job under an alias and an owner
  job of every pool core that preempts it at once; `simq run` streaming and
  exit codes, reattach by idempotency key, SIGINT/SIGTERM while running and
  SIGTERM while queued). `test_upgrade.py` rewritten for 0.6.2 to 0.6.3.
- `./upgrade.py install` (dry run) against the live state: no PROBLEM line;
  layout `spheral=0-3, raikou-mfem=4-5+1gpu, powerflow=6-7+1gpu` to
  `zenobia=0-7+1gpu (aliases spheral, raikou-mfem, powerflow)`. With
  `--config config/zenobia-16core-node.json` it is refused on the 8-core CPU
  ("CPU 8 is an SMT sibling of CPU 0"), as intended.
- Not validated: the 16-core file on the real 9950X3D2 (parsed in the tests
  without the topology check; `test_pool_aliases` loads it with the check
  when the machine has 32 logical CPUs), and a live install.

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

0.7.0 worker, memory, GitHub:

- The admission guard sees MemAvailable, which counts reclaimable page cache
  as available; a job whose request fits may still push other processes'
  cache out. It does not count the unused part of running jobs' requests
  (that would block almost every start on zenobia, where reservations far
  exceed use); the reservations (`memory_schedulable`) bound that.
- Memory pressure stops dispatch and preempts borrow jobs only. If an owner
  job (or something outside the queue) keeps growing, the pressure lasts
  until the kernel's OOM handling or the job's own `MemoryMax` ends it; the
  record names the job that grew.
- `--mem auto` learns only from jobs that ended under 0.7.0 (earlier jobs
  have no measured peak) and only by name / name family; two different jobs
  with the same name share a history. A request is fixed when the job is
  queued.
- A job adopted from 0.6.3 has no unit record of its peak: the highest
  sampled `memory.peak` is kept (source `sampled`); it ends at most 3 s
  before the worker notices (no stop-post wake-up from the old wrapper).
- `IOSchedulingClass=idle` on the worker has no effect while the state
  directory is on NVMe with the `none` I/O scheduler (as on zenobia); it
  takes effect on a disk with BFQ or mq-deadline.
- GitHub: only the latest phase is posted (a job that starts before the
  first post shows running directly). Comments are edited in place but the
  first creation can notify watchers of the pull request depending on their
  settings; statuses never notify. A job submitted with GitHub options
  under 0.7.0 is not reported by 0.6.3 after a rollback. The real API was
  exercised once (statuses on a throwaway branch); the comment path was
  tested against the local fake API only.

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
