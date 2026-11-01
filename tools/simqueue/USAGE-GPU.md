# Using the GPU on zenobia (SimQueue 0.6.2)

For every session and worker that runs GPU commands on zenobia.

**Stop taking the GPU lease files with raw `flock`.** Run GPU work either
through the queue (`simq submit --gpus 1 ...`) or with `gpu-lease`:

```sh
G=/home/treefor/projects/runtimes/simqueue-0.6.2-20261031/bin/gpu-lease
```

| Before | Now |
|---|---|
| `flock $L/testing-gpu.lock flock $L/testing-gpu-<UUID>.lock any-core.sh --wait CMD` | `$G --core any -- CMD` |
| `flock $L/testing-gpu.lock flock $L/testing-gpu-<UUID>.lock both-cores.sh CMD` | `$G --core both -- CMD` |
| `flock $L/testing-gpu.lock flock $L/testing-gpu-<UUID>.lock CMD` (GPU only) | `$G -- CMD` |
| GPU plus one specific core K | `$G --core K -- CMD` |

`$L` is `/home/treefor/projects/.workspace/outputs/pic-profile`.

## What `gpu-lease` does

1. Waits its turn: queued SimQueue GPU jobs that were waiting before it
   started go first, other `gpu-lease` users are served in start order.
2. Takes `testing-gpu.lock`, then `testing-gpu-<UUID>.lock`, then the core
   lease (`--core`), always in that order, so it cannot deadlock with the
   queue or with other `gpu-lease` users.
3. Runs the command. With `--core` it is pinned to the core(s), at nice 15,
   with `PF_CORE` and `OMP_NUM_THREADS=1`, as `any-core.sh` and
   `both-cores.sh` do. It also exports `GPU_LEASE_GPU`, `GPU_LEASE_LOCKS`
   and (with `--core`) `GPU_LEASE_CORES`.
4. Releases everything when the command and its children have closed the
   lease files. The exit status is the command's.

Useful options: `--timeout S` gives up after S seconds of waiting (exit 75,
as `any-core.sh`); `--label TEXT` names the request in `simq partitions`;
`-q` silences the progress messages.

## Why

The kernel gives a released flock to the first blocked waiter at once. Raw
`flock` chains run back to back therefore took every release, and queued GPU
jobs starved (job 169 waited more than 20 minutes). Since 0.6.2 a queued GPU
job that is ready joins the GPU line and wins the next release; a process
blocked in raw `flock` on the GPU file is moved behind the line (it is
stopped and continued, which restarts its wait at the tail; the command
itself is not affected). Raw `flock` still works, but it is served last.

## Seeing the line

```sh
$G --status            # holder, line, raw waiters, queued jobs (--json for scripts)
simq partitions        # the same per GPU, with a warning when a queued job starves
```

Tools that decide by themselves whether to start GPU work can read
`$L/simq-gpu-line.json` (schema `simqueue-gpu-line-v1`): yield while
`gpus.<UUID>.queue` is non-empty and the file is fresh
(`now - heartbeat_unix <= stale_after_s`).

## Do not

- Do not call `gpu-lease` (or `flock` on the GPU files) from inside a
  SimQueue job submitted with `--gpus`/`--gpu-mem`, or from inside another
  `gpu-lease` command: the outer one already holds the GPU files, and the
  inner request would wait for its own parent.
- Do not take a core lease first and then the GPU files: the GPU files come
  first, always.

Details: `OPERATIONS.md`, section "The GPU line (0.6.2)".
