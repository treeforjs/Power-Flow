"""Print one JSON line describing this rank's placement (affinity, nice, thread and GPU variables).

Optional arguments: directories; for each, report whether it is the current root (frozen-root check).
"""
import json
import os
import socket
import sys


def inside(root):
    try:
        return os.path.samestat(os.stat("/"), os.stat(root))
    except OSError:
        return False


def ancestors(depth=5):
    """Name, allowed CPUs and nice of the parent chain up to the pf-launch process
    (shows where the rank wrapper, mpirun and the sandbox run)."""
    chain, pid = [], os.getppid()
    for _ in range(depth):
        if pid <= 0:
            break
        try:
            with open(f"/proc/{pid}/status") as stream:
                status = dict(line.split(":", 1) for line in stream if ":" in line)
            with open(f"/proc/{pid}/stat") as stream:
                nice = int(stream.read().rsplit(")", 1)[1].split()[16])
        except OSError:
            break
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as stream:
                argv = stream.read().split(b"\0")
        except OSError:
            argv = []
        outer = any(a.endswith(b"pf_launch.py") for a in argv) and b"--rank-exec" not in argv
        chain.append({"pid": pid, "name": status["Name"].strip(), "cpus_allowed": status["Cpus_allowed_list"].strip(),
                      "nice": nice, "role": "pf-launch" if outer else None})
        if outer:
            break  # stop at the launcher; its parents are outside the run
        pid = int(status["PPid"])
    return chain


keys = ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "SLURM_PROCID", "PF_LAUNCH_RANK", "PF_LAUNCH_NRANKS", "PF_LAUNCH_CPU",
        "PF_CORE", "CUDA_VISIBLE_DEVICES", "CUDA_DEVICE_ORDER", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS")
print(json.dumps({
    "probe": "rank_probe.py",
    "host": socket.gethostname(),
    "pid": os.getpid(),
    "affinity": sorted(os.sched_getaffinity(0)),
    "nice": os.getpriority(os.PRIO_PROCESS, 0),
    "env": {k: os.environ.get(k) for k in keys},
    "python": sys.version.split()[0],
    "inside_root": {root: inside(root) for root in sys.argv[1:]},
    "ancestors": ancestors(),
}, sort_keys=True), flush=True)
