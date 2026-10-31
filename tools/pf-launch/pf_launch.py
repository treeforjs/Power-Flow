#!/usr/bin/env python3
"""pf-launch: portable launcher with per-run core and GPU selection.

    pf-launch [--site SITE.json] [--cores N|auto | --cpus LIST] [--gpus none|auto|LIST|count:N]
              [--ranks-per-gpu K] [--cores-available] [--wait] [--launcher auto|exec|mpirun|srun]
              [--execution site|host] [--prefix CMD] [--dry-run] [--json] [--receipt PATH]
              -- PROGRAM ARGS...
    pf-launch --detect [--site PATH] [--force]

One rank per selected physical core (never two SMT siblings), pinned to one logical CPU;
one thread per rank; nice level, reservations, leases, execution prefix and identity
admission from a per-site JSON file; an execution envelope receipt per run.
Python 3.9+ standard library only.  See README.md.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import datetime as _dt
import fcntl
import hashlib
import io
import json
import math
import os
import re
import shlex
import signal
import socket
import subprocess
import sys
import time

VERSION = "1.1.0"
SITE_SCHEMA = "pf-launch-site/1"
RANK_PLAN_SCHEMA = "pf-launch-rank-plan/1"

EXIT_USAGE = 64        # malformed request (bad list, SMT siblings, too many GPUs for the ranks)
EXIT_ABSENT = 69       # requested CPU or GPU does not exist on this host
EXIT_SOFTWARE = 70     # launcher could not start the program
EXIT_BUSY = 75         # requested resources reserved or leased
EXIT_CONFIG = 78       # site file problem or identity admission failure

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
RANK_ENV_KEYS = ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "SLURM_PROCID")
LOCAL_RANK_ENV_KEYS = ("OMPI_COMM_WORLD_LOCAL_RANK", "SLURM_LOCALID", "MPI_LOCALRANKID")
SLURM_GPU_ENV_KEYS = ("CUDA_VISIBLE_DEVICES", "SLURM_STEP_GPUS", "SLURM_JOB_GPUS", "GPU_DEVICE_ORDINAL")
PER_RANK_ENV_KEYS = ("CUDA_VISIBLE_DEVICES", "PF_LAUNCH_RANK", "PF_LAUNCH_NRANKS", "PF_LAUNCH_CPU", "PF_CORE")
GPU_QUERY_FIELDS = "index,uuid,name,compute_cap,memory.total,driver_version"
IDENTITY_KEYS = ("kernel_release", "nvidia_driver_version", "nvidia_srcversion", "gpu_uuids", "cpu_model",
                 "frozen_root_manifest_sha256")
RECEIPT_FIELDS = ("host", "site_file_digest", "selected_cores", "selected_gpus", "rank_count", "rank_to_core",
                  "rank_to_gpu", "gpus", "launcher_argv", "environment_digest", "root_or_image_digest",
                  "start_time", "end_time", "exit_codes")
WAIT_POLL_SECONDS = 3.0
SELF_PATH = os.path.abspath(__file__)


class PfError(Exception):
    """A refusal or configuration error with a process exit code."""

    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code
        self.message = message


def utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


# --------------------------------------------------------------------------- CPU lists


def parse_cpu_list(text: str) -> list:
    """Parse a Linux CPU list ("4-7", "0,2,5", "0-3,8-11") preserving order; reject duplicates."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty CPU list")
    cpus = []
    for part in text.split(","):
        part = part.strip()
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not match:
            raise ValueError(f"invalid CPU list element {part!r} in {text!r}")
        low = int(match.group(1))
        high = int(match.group(2)) if match.group(2) is not None else low
        if high < low:
            raise ValueError(f"descending CPU range {part!r} in {text!r}")
        cpus.extend(range(low, high + 1))
    seen, duplicates = set(), set()
    for cpu in cpus:
        (duplicates if cpu in seen else seen).add(cpu)
    if duplicates:
        raise ValueError(f"CPU(s) listed more than once: {format_cpu_list(duplicates)}")
    return cpus


def format_cpu_list(cpus) -> str:
    values = sorted(set(int(c) for c in cpus))
    parts, i = [], 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[j + 1] == values[j] + 1:
            j += 1
        parts.append(str(values[i]) if i == j else f"{values[i]}-{values[j]}")
        i = j + 1
    return ",".join(parts)


# --------------------------------------------------------------------------- topology


class Topology:
    """Physical cores keyed by their primary (lowest) logical CPU."""

    def __init__(self, cores, cpu_model=None):
        self.cpu_model = cpu_model
        self.cores = {}
        self.cpu_core = {}
        for entry in cores:
            cpus = sorted(int(c) for c in entry["cpus"])
            if not cpus:
                continue
            primary = cpus[0]
            for cpu in cpus:
                if cpu in self.cpu_core:
                    raise ValueError(f"CPU {cpu} appears in two SMT sibling groups")
                self.cpu_core[cpu] = primary
            self.cores[primary] = {
                "core": primary,
                "cpus": cpus,
                "package": int(entry.get("package", 0)),
                "core_id": int(entry.get("core_id", primary)),
                "numa": int(entry.get("numa", 0)),
            }

    @classmethod
    def from_siblings(cls, groups, numa_of=None, cpu_model=None):
        return cls([{"cpus": list(group), "core_id": i, "numa": numa_of(min(group)) if numa_of else 0}
                    for i, group in enumerate(groups)], cpu_model)

    def primaries(self) -> list:
        return sorted(self.cores)

    def logical_cpus(self) -> list:
        return sorted(self.cpu_core)

    def siblings(self, cpu: int) -> list:
        return self.cores[self.cpu_core[cpu]]["cpus"]

    def to_json(self) -> dict:
        numa = {}
        for rec in self.cores.values():
            numa.setdefault(str(rec["numa"]), []).extend(rec["cpus"])
        return {
            "cpu_model": self.cpu_model,
            "logical_cpus": format_cpu_list(self.logical_cpus()),
            "physical_core_count": len(self.cores),
            "smt_threads_per_core": max((len(r["cpus"]) for r in self.cores.values()), default=0),
            "numa_nodes": {k: format_cpu_list(v) for k, v in sorted(numa.items(), key=lambda kv: int(kv[0]))},
            "physical_cores": [self.cores[c] for c in self.primaries()],
        }


def _read_text(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as stream:
            return stream.read().strip()
    except OSError:
        return None


def read_cpu_model(procfs="/proc"):
    text = _read_text(os.path.join(procfs, "cpuinfo")) or ""
    for line in text.splitlines():
        key, _, value = line.partition(":")
        if key.strip() in ("model name", "Model name", "cpu model"):
            return value.strip()
    return None


def read_topology(sysfs="/sys", procfs="/proc") -> Topology:
    """Physical cores from /sys/devices/system/cpu/cpu*/topology/thread_siblings_list."""
    base = os.path.join(sysfs, "devices/system/cpu")
    online_text = _read_text(os.path.join(base, "online"))
    if online_text:
        online = parse_cpu_list(online_text)
    else:
        online = sorted(int(n[3:]) for n in os.listdir(base) if re.fullmatch(r"cpu\d+", n))
    online_set = set(online)
    numa = {}
    node_base = os.path.join(sysfs, "devices/system/node")
    if os.path.isdir(node_base):
        for name in os.listdir(node_base):
            match = re.fullmatch(r"node(\d+)", name)
            text = _read_text(os.path.join(node_base, name, "cpulist")) if match else None
            if text:
                for cpu in parse_cpu_list(text):
                    numa[cpu] = int(match.group(1))
    assigned, groups = set(), []
    for cpu in online:
        if cpu in assigned:
            continue
        topo = os.path.join(base, f"cpu{cpu}", "topology")
        text = _read_text(os.path.join(topo, "thread_siblings_list")) or _read_text(os.path.join(topo, "core_cpus_list"))
        siblings = [c for c in parse_cpu_list(text) if c in online_set and c not in assigned] if text else []
        if cpu not in siblings:
            siblings.append(cpu)
        assigned.update(siblings)
        package = _read_text(os.path.join(topo, "physical_package_id"))
        core_id = _read_text(os.path.join(topo, "core_id"))
        groups.append({"cpus": sorted(siblings), "package": int(package) if package not in (None, "") else 0,
                       "core_id": int(core_id) if core_id not in (None, "") else cpu,
                       "numa": numa.get(min(siblings), 0)})
    return Topology(groups, cpu_model=read_cpu_model(procfs))


def core_order(topology: Topology, preference=None) -> list:
    """Site preference first (mapped to physical cores), then remaining cores by (NUMA node, core)."""
    ordered = []
    for cpu in preference or []:
        core = topology.cpu_core.get(int(cpu))
        if core is not None and core not in ordered:
            ordered.append(core)
    rest = sorted((c for c in topology.primaries() if c not in ordered), key=lambda c: (topology.cores[c]["numa"], c))
    return ordered + rest


# --------------------------------------------------------------------------- GPUs


def parse_gpu_csv(text: str) -> list:
    """Parse `nvidia-smi --query-gpu=index,uuid,name,compute_cap,memory.total,driver_version --format=csv,noheader`."""
    gpus = []
    for row in csv.reader(io.StringIO(text), skipinitialspace=True):
        row = [cell.strip() for cell in row]
        if not row or row == [""]:
            continue
        if len(row) != 6:
            raise ValueError(f"unexpected nvidia-smi row {row!r}")
        index, uuid, name, capability, memory, driver = row
        gpus.append({"index": int(index), "uuid": uuid, "name": name, "compute_cap": capability,
                     "memory_total": memory, "driver_version": driver, "pci_bus_id": None, "numa": None})
    return sorted(gpus, key=lambda g: g["index"])


def parse_gpu_request(text: str):
    """none | auto | count:N | comma list of indices and/or GPU-UUIDs."""
    value = (text or "").strip()
    lowered = value.lower()
    if lowered == "none":
        return ("none", None)
    if lowered == "auto":
        return ("auto", None)
    match = re.fullmatch(r"count:(\d+)", lowered)
    if match:
        return ("count", int(match.group(1)))
    tokens = [t.strip() for t in value.split(",")]
    if not value or any(not t for t in tokens):
        raise ValueError(f"--gpus: empty element in {text!r}")
    out = []
    for token in tokens:
        if re.fullmatch(r"\d+", token):
            out.append(int(token))
        elif re.fullmatch(r"GPU-[0-9A-Fa-f-]+", token):
            out.append(token)
        else:
            raise ValueError(f"--gpus: invalid element {token!r}; expected none, auto, count:N, "
                             "GPU indices (0,1) or UUIDs (GPU-...)")
    if len(set(str(t).lower() for t in out)) != len(out):
        raise ValueError(f"--gpus: GPU listed more than once in {text!r}")
    return ("list", out)


def find_gpu(inventory, token):
    for gpu in inventory:
        if isinstance(token, int) and gpu["index"] == token:
            return gpu
        if isinstance(token, str) and gpu["uuid"].lower() == token.lower():
            return gpu
    return None


def describe_inventory(inventory) -> str:
    if not inventory:
        return "no GPUs detected"
    return ", ".join(f"{g['index']}={g['uuid']}" for g in inventory)


def rank_gpu_map(nranks: int, gpus, ranks_per_gpu: int) -> list:
    """Rank r uses gpus[r // K]; ranks beyond len(gpus) * K run CPU-only (None)."""
    out = []
    for rank in range(nranks):
        slot = rank // ranks_per_gpu
        out.append(gpus[slot]["uuid"] if slot < len(gpus) else None)
    return out


# --------------------------------------------------------------------------- system access


class System:
    """Live host access.  Tests substitute a subclass with fake topology, GPUs, probes and leases."""

    def __init__(self, sysfs="/sys", procfs="/proc"):
        self.sysfs = sysfs
        self.procfs = procfs
        self.environ = dict(os.environ)
        self._sha_cache = {}

    def hostname(self):
        return socket.gethostname()

    def topology(self):
        return read_topology(self.sysfs, self.procfs)

    def cpu_model(self):
        return read_cpu_model(self.procfs)

    def kernel_release(self):
        return _read_text(os.path.join(self.procfs, "sys/kernel/osrelease"))

    def read_text(self, path):
        return _read_text(path)

    def exists(self, path):
        return os.path.exists(path)

    def sha256_file(self, path):
        """sha256 of a file, cached by (path, mtime, size) so --wait re-planning does not re-hash."""
        try:
            st = os.stat(path)
            key = (os.path.abspath(path), st.st_mtime_ns, st.st_size)
            if key not in self._sha_cache:
                self._sha_cache[key] = sha256_file(path)
            return self._sha_cache[key]
        except OSError:
            return None

    def run(self, argv, timeout):
        done = subprocess.run(argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, timeout=timeout, check=False)
        return done.returncode, done.stdout, done.stderr

    def gpu_inventory(self, nvidia_smi="nvidia-smi"):
        """Return (gpus, note).  An absent or failing nvidia-smi yields an empty inventory and a note."""
        try:
            rc, out, err = self.run([nvidia_smi, f"--query-gpu={GPU_QUERY_FIELDS}", "--format=csv,noheader"], 30)
        except FileNotFoundError:
            return [], f"{nvidia_smi} not found"
        except (OSError, subprocess.SubprocessError) as error:
            return [], f"{nvidia_smi} failed: {error}"
        if rc != 0:
            return [], f"{nvidia_smi} exited {rc}: {(err or out).strip()[:200]}"
        try:
            gpus = parse_gpu_csv(out)
        except ValueError as error:
            return [], f"{nvidia_smi} output not understood: {error}"
        try:
            rc, out, _ = self.run([nvidia_smi, "--query-gpu=uuid,pci.bus_id", "--format=csv,noheader"], 30)
            if rc == 0:
                for row in csv.reader(io.StringIO(out), skipinitialspace=True):
                    if len(row) == 2:
                        gpu = find_gpu(gpus, row[0].strip())
                        if gpu is not None:
                            gpu["pci_bus_id"] = row[1].strip()
                            gpu["numa"] = self._pci_numa(row[1].strip())
        except (OSError, subprocess.SubprocessError):
            pass
        return gpus, None

    def _pci_numa(self, bus_id):
        match = re.fullmatch(r"([0-9A-Fa-f]+):([0-9A-Fa-f]{2}):([0-9A-Fa-f]{2})\.([0-7])", bus_id)
        if not match:
            return None
        name = f"{int(match.group(1), 16):04x}:{match.group(2).lower()}:{match.group(3).lower()}.{match.group(4)}"
        text = _read_text(os.path.join(self.sysfs, "bus/pci/devices", name, "numa_node"))
        try:
            return int(text)
        except (TypeError, ValueError):
            return None

    def lease_busy(self, path):
        return lease_probe_busy(path)

    def current_nice(self):
        return os.getpriority(os.PRIO_PROCESS, 0)

    def cgroup_cpus(self):
        """CPUs permitted by the cgroup v2 cpuset of this process, or None when unknown."""
        text = _read_text(os.path.join(self.procfs, "self/cgroup")) or ""
        rel = None
        for line in text.splitlines():
            if line.startswith("0::"):
                rel = line[3:].strip().lstrip("/")
        if rel is None:
            return None
        path = os.path.join(self.sysfs, "fs/cgroup", rel)
        root = os.path.join(self.sysfs, "fs/cgroup")
        while True:
            value = _read_text(os.path.join(path, "cpuset.cpus.effective"))
            if value:
                try:
                    return set(parse_cpu_list(value))
                except ValueError:
                    return None
            if os.path.normpath(path) == os.path.normpath(root):
                return None
            path = os.path.dirname(path)

    def same_dir(self, a, b):
        try:
            return os.path.samestat(os.stat(a), os.stat(b))
        except OSError:
            return False

    def now(self):
        return utc_now()


# --------------------------------------------------------------------------- leases


def lease_probe_busy(path: str) -> bool:
    """True when another open file description holds an flock on PATH.  A missing file is free."""
    try:
        fd = os.open(path, os.O_RDONLY | os.O_CLOEXEC)
    except FileNotFoundError:
        return False
    except OSError:
        return True
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(fd, fcntl.LOCK_UN)
        return False
    finally:
        os.close(fd)


class LeaseSet:
    """flock(2) leases (compatible with flock(1)), held until release() or process exit."""

    def __init__(self):
        self.held = []

    def acquire(self, path: str, block: bool) -> bool:
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        fd = os.open(path, os.O_RDONLY | os.O_CREAT, 0o644)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | (0 if block else fcntl.LOCK_NB))
        except BlockingIOError:
            os.close(fd)
            return False
        except BaseException:
            os.close(fd)
            raise
        self.held.append((path, fd))
        return True

    def paths(self):
        return [p for p, _ in self.held]

    def fds(self):
        return [fd for _, fd in self.held]

    def release(self):
        while self.held:
            _, fd = self.held.pop()
            try:
                fcntl.flock(fd, fcntl.LOCK_UN)
            finally:
                os.close(fd)


def core_lease_path(site: dict, core: int):
    cfg = site.get("leases") or {}
    overrides = cfg.get("core_overrides") or {}
    if str(core) in overrides:
        return os.path.expanduser(overrides[str(core)])
    if cfg.get("core_path"):
        return os.path.expanduser(cfg["core_path"].format(core=core))
    return None


def gpu_lease_paths(site: dict, gpu: dict) -> list:
    cfg = site.get("leases") or {}
    paths = []
    if cfg.get("gpu_path"):
        paths.append(os.path.expanduser(cfg["gpu_path"].format(uuid=gpu["uuid"], index=gpu["index"])))
    for extra in (cfg.get("gpu_extra") or {}).get(gpu["uuid"], []) or []:
        paths.append(os.path.expanduser(extra))
    return paths


def lease_order(site: dict, cores, gpus) -> list:
    """Acquisition order shared by every pf-launch instance: cores ascending, then GPUs by index, then extras."""
    order = []
    for core in sorted(set(cores)):
        path = core_lease_path(site, core)
        if path and path not in order:
            order.append(path)
    extras = []
    for gpu in sorted(gpus, key=lambda g: g["index"]):
        paths = gpu_lease_paths(site, gpu)
        cfg = site.get("leases") or {}
        primary = paths[:1] if cfg.get("gpu_path") else []
        for path in primary:
            if path not in order:
                order.append(path)
        for path in paths[len(primary):]:
            if path not in order and path not in extras:
                extras.append(path)
    return order + extras


# --------------------------------------------------------------------------- reservations


def run_reservation_probe(spec: dict, system: System):
    """Return (reserved, detail).  Missing paths, errors, timeouts, non-zero exits and
    unrecognized output all count as reserved."""
    probe = spec.get("probe") or {}
    for path in probe.get("require_paths", []) or []:
        if not system.exists(os.path.expanduser(path)):
            return True, f"required path missing: {path} (treated as reserved)"
    argv = [os.path.expanduser(a) for a in probe.get("argv") or []]
    if not argv:
        return True, "probe has no argv (treated as reserved)"
    timeout = float(probe.get("timeout_s", 20))
    try:
        rc, out, _ = system.run(argv, timeout)
    except subprocess.TimeoutExpired:
        return True, f"probe timed out after {timeout:g} s (treated as reserved)"
    except (OSError, subprocess.SubprocessError) as error:
        return True, f"probe failed: {error} (treated as reserved)"
    if rc != 0:
        return True, f"probe exited {rc} (treated as reserved)"
    lines = [line.rstrip() for line in (out or "").splitlines() if line.strip()]
    try:
        busy_re = probe.get("busy_regex")
        if busy_re:
            busy = [line for line in lines if re.search(busy_re, line)]
            if busy:
                return True, f"active ({len(busy)} matching line(s); first: {busy[0].strip()!r})"
        free_re = probe.get("free_regex")
        if free_re:
            if any(re.fullmatch(free_re, line.strip()) for line in lines):
                return False, "idle"
            return True, f"unrecognized probe output (treated as reserved): {(out or '').strip()[:120]!r}"
    except re.error as error:
        return True, f"invalid probe regex: {error} (treated as reserved)"
    return False, "probe exited 0"


def evaluate_reservations(site: dict, topology: Topology, system: System, inventory) -> dict:
    """Return {"cores": {core: [reasons]}, "gpus": {uuid: [reasons]}, "probes": [...]}."""
    cores, gpus, probes = {}, {}, []

    def add_core(cpu, reason):
        core = topology.cpu_core.get(int(cpu), int(cpu))
        cores.setdefault(core, [])
        if reason not in cores[core]:
            cores[core].append(reason)

    def add_gpu(token, reason):
        gpu = find_gpu(inventory, token)
        if gpu is not None:
            gpus.setdefault(gpu["uuid"], []).append(reason)

    static = site.get("reserved_cores") or {}
    if isinstance(static, list):
        static = {"cores": static}
    for cpu in static.get("cores", []) or []:
        add_core(cpu, f"static: {static.get('reason') or 'reserved in site file'}")
    gstatic = site.get("reserved_gpus") or {}
    if isinstance(gstatic, list):
        gstatic = {"gpus": gstatic}
    for token in gstatic.get("gpus", []) or []:
        add_gpu(token, f"static: {gstatic.get('reason') or 'reserved in site file'}")
    for spec in site.get("reservations", []) or []:
        name = spec.get("name", "probe")
        reserved, detail = run_reservation_probe(spec, system)
        probes.append({"name": name, "cores": spec.get("cores", []), "gpus": spec.get("gpus", []),
                       "reserved": reserved, "detail": detail})
        if reserved:
            for cpu in spec.get("cores", []) or []:
                add_core(cpu, f"{name}: {detail}")
            for token in spec.get("gpus", []) or []:
                add_gpu(token, f"{name}: {detail}")
    return {"cores": cores, "gpus": gpus, "probes": probes}


def describe_reasons(reason_map: dict) -> str:
    groups = {}
    for core in sorted(reason_map):
        for reason in reason_map[core]:
            groups.setdefault(reason, []).append(core)
    return "; ".join(f"cores {format_cpu_list(cs)} ({reason})" for reason, cs in groups.items())


def describe_unavailable(reserved: dict, busy: dict) -> str:
    parts = []
    if reserved:
        parts.append("reserved: " + describe_reasons(reserved))
    if busy:
        parts.append("leased: " + "; ".join(f"core {c} ({busy[c]})" for c in sorted(busy)))
    return "; ".join(parts) if parts else "no reservations or leases"


# --------------------------------------------------------------------------- selection


def _refusal(code, message, retryable=False):
    return {"code": code, "message": message, "retryable": retryable}


def select_cores(topology: Topology, request, reserved=None, busy=None, preference=None,
                 allow_fewer=False, wait=False) -> dict:
    """Resolve a core request into ranks.

    request: ("auto", None) | ("count", N) | ("list", [cpu, ...]).
    reserved: {core: [reasons]}; busy: {core: lease path held elsewhere}.
    Returns {"ranks": [{"rank", "cpu", "core", "siblings"}], "refusals": [...], "notes": [...]}.
    """
    reserved = reserved or {}
    busy = busy or {}
    out = {"ranks": [], "refusals": [], "notes": []}
    refuse = lambda *a, **k: out["refusals"].append(_refusal(*a, **k))  # noqa: E731
    kind, value = request
    if kind == "list":
        cpus = list(value)
        unknown = [c for c in cpus if c not in topology.cpu_core]
        if unknown:
            refuse(EXIT_ABSENT, f"--cpus: CPU(s) {format_cpu_list(unknown)} not present or offline on this host "
                                f"(online CPUs: {format_cpu_list(topology.logical_cpus())})")
        by_core = {}
        for cpu in cpus:
            if cpu in topology.cpu_core:
                by_core.setdefault(topology.cpu_core[cpu], []).append(cpu)
        for core, members in by_core.items():
            if len(members) > 1:
                refuse(EXIT_USAGE, f"--cpus: CPUs {' and '.join(map(str, members))} are SMT siblings of physical "
                                   f"core {core} (thread_siblings_list "
                                   f"{format_cpu_list(topology.cores[core]['cpus'])}); one rank per physical core, "
                                   f"so name at most one of them")
        checked = set()
        for cpu in cpus:
            core = topology.cpu_core.get(cpu)
            if core is None or core in checked:
                continue
            checked.add(core)
            if core in reserved:
                refuse(EXIT_BUSY, f"--cpus: CPU {cpu} (physical core {core}) is reserved: {'; '.join(reserved[core])}")
            elif core in busy:
                if wait:
                    out["notes"].append(f"CPU {cpu} (physical core {core}) is leased ({busy[core]}); --wait queues for it")
                else:
                    refuse(EXIT_BUSY, f"--cpus: CPU {cpu} (physical core {core}) is leased: lock held on "
                                      f"{busy[core]}; use --wait to queue for it", retryable=True)
        if not out["refusals"]:
            out["ranks"] = [{"rank": i, "cpu": cpu, "core": topology.cpu_core[cpu],
                             "siblings": topology.siblings(cpu)} for i, cpu in enumerate(cpus)]
        return out

    order = core_order(topology, preference)
    total = len(order)
    free = [c for c in order if c not in reserved and c not in busy]
    unreserved = [c for c in order if c not in reserved]
    why = describe_unavailable({c: reserved[c] for c in order if c in reserved},
                               {c: busy[c] for c in order if c in busy and c not in reserved})
    if kind == "auto":
        if not free:
            refuse(EXIT_BUSY, f"--cores auto: no free physical core of {total}; {why}", retryable=bool(unreserved))
            return out
        chosen = free
        if len(free) < total:
            out["notes"].append(f"--cores auto: {len(free)} of {total} physical cores free ({format_cpu_list(free)}); {why}")
    elif kind == "count":
        n = int(value)
        if n < 1:
            refuse(EXIT_USAGE, "--cores must be at least 1")
            return out
        if n > total and not allow_fewer:
            refuse(EXIT_ABSENT, f"--cores {n}: this host has {total} physical cores "
                                f"({format_cpu_list(order)}); use --cores-available to accept fewer")
            return out
        if len(free) >= n:
            chosen = free[:n]
        elif allow_fewer and free:
            chosen = free
            out["notes"].append(f"--cores {n}: {len(free)} physical core(s) available "
                                f"({format_cpu_list(free)}); running {len(free)} rank(s) (--cores-available); {why}")
        else:
            enough_unreserved = len(unreserved) >= min(n, total)
            if enough_unreserved:
                hint = "--wait queues until the leased cores are released"
                if free:
                    hint += f", or use --cores {len(free)} / --cores-available to accept fewer"
            else:
                hint = (f"at most {len(unreserved)} physical core(s) are unreserved on this host "
                        f"({format_cpu_list(unreserved) or 'none'})")
                if unreserved:
                    hint += f"; use --cores {len(free) or len(unreserved)} or --cores-available to accept fewer"
            refuse(EXIT_BUSY, f"--cores {n}: only {len(free)} of {total} physical cores available "
                              f"({format_cpu_list(free) or 'none'}); {why}; {hint}", retryable=enough_unreserved)
            return out
    else:
        raise ValueError(f"unknown core request {kind!r}")
    out["ranks"] = [{"rank": i, "cpu": core, "core": core, "siblings": topology.cores[core]["cpus"]}
                    for i, core in enumerate(chosen)]
    return out


# --------------------------------------------------------------------------- SimQueue allocations

SIMQ_ENV_KEYS = ("SIMQ_JOB_ID", "SIMQ_PARTITION", "SIMQ_NODE", "SIMQ_CPU_LIST", "SIMQ_CPU_COUNT", "SIMQ_NCPUS",
                 "SIMQ_GPUS", "SIMQ_GPU_MODE", "SIMQ_GPU_MEM", "SIMQ_MEM", "SIMQ_TIME_LIMIT", "SIMQ_UNIT")


def simqueue_allocation(env: dict):
    """The SimQueue allocation this process runs in, or None.

    SimQueue 0.6 exports SIMQ_JOB_ID with SIMQ_CPU_LIST and SIMQ_GPUS (empty when no GPU is allocated).
    A SimQueue 0.5 job has SIMQ_CPU_LIST but no SIMQ_GPUS: its cores are the allocation and GPUs are
    still selected (and leased) locally.
    """
    if not env.get("SIMQ_JOB_ID") or ("SIMQ_CPU_LIST" not in env and "SIMQ_GPUS" not in env):
        return None
    allocation = {"env": {k: env.get(k) for k in SIMQ_ENV_KEYS if k in env}, "error": None}
    try:
        allocation["cpus"] = parse_cpu_list(env.get("SIMQ_CPU_LIST", ""))
    except ValueError as error:
        allocation["cpus"] = []
        allocation["error"] = f"invalid SIMQ_CPU_LIST {env.get('SIMQ_CPU_LIST')!r}: {error}"
    allocation["gpus_managed"] = "SIMQ_GPUS" in env
    allocation["gpus"] = [g.strip() for g in env.get("SIMQ_GPUS", "").split(",") if g.strip()]
    return allocation


def select_simq_cores(topology: Topology, request, allocation: dict, allow_fewer=False) -> dict:
    """Ranks from the SimQueue allocation: pf-launch picks no cores of its own."""
    out = {"ranks": [], "refusals": [], "notes": []}
    refuse = lambda *a, **k: out["refusals"].append(_refusal(*a, **k))  # noqa: E731
    cpus = list(allocation["cpus"])
    job = allocation["env"].get("SIMQ_JOB_ID")
    if allocation.get("error"):
        refuse(EXIT_CONFIG, allocation["error"])
        return out
    if not cpus:
        refuse(EXIT_CONFIG, f"SimQueue job {job}: SIMQ_CPU_LIST is empty")
        return out
    unknown = [c for c in cpus if c not in topology.cpu_core]
    if unknown:
        refuse(EXIT_ABSENT, f"SimQueue allocation CPU(s) {format_cpu_list(unknown)} not present on this host")
        return out
    kind, value = request
    if kind == "list":
        outside = [c for c in value if c not in cpus]
        if outside:
            refuse(EXIT_BUSY, f"--cpus: CPU(s) {format_cpu_list(outside)} are not in SimQueue job {job}'s "
                              f"allocation (SIMQ_CPU_LIST={format_cpu_list(cpus)})")
            return out
        chosen = list(value)
    elif kind == "auto":
        chosen = cpus
    else:
        n = int(value)
        if n < 1:
            refuse(EXIT_USAGE, "--cores must be at least 1")
            return out
        if n <= len(cpus):
            chosen = cpus[:n]
        elif allow_fewer:
            chosen = cpus
            out["notes"].append(f"--cores {n}: SimQueue allocation has {len(cpus)} core(s); running {len(cpus)} "
                                f"rank(s) (--cores-available)")
        else:
            refuse(EXIT_BUSY, f"--cores {n} exceeds SimQueue job {job}'s allocation of {len(cpus)} core(s) "
                              f"(SIMQ_CPU_LIST={format_cpu_list(cpus)}); use --cores {len(cpus)} or --cores-available")
            return out
    by_core = {}
    for cpu in chosen:
        by_core.setdefault(topology.cpu_core[cpu], []).append(cpu)
    for core, members in by_core.items():
        if len(members) > 1:
            refuse(EXIT_USAGE, f"CPUs {' and '.join(map(str, members))} are SMT siblings of physical core {core}; "
                               f"one rank per physical core")
    if out["refusals"]:
        return out
    out["notes"].append(f"SimQueue allocation (job {job}, partition {allocation['env'].get('SIMQ_PARTITION', '?')}): "
                        f"CPUs {format_cpu_list(cpus)}; the queue's job unit holds the core leases")
    out["ranks"] = [{"rank": i, "cpu": cpu, "core": topology.cpu_core[cpu], "siblings": topology.siblings(cpu)}
                    for i, cpu in enumerate(chosen)]
    return out


def select_simq_gpus(inventory, request, nranks: int, ranks_per_gpu: int, allocation: dict) -> dict:
    """GPUs from the SimQueue allocation (SIMQ_GPUS); no GPU leases are taken."""
    out = {"gpus": [], "refusals": [], "notes": []}
    refuse = lambda *a, **k: out["refusals"].append(_refusal(*a, **k))  # noqa: E731
    allocated = []
    for uuid in allocation["gpus"]:
        gpu = find_gpu(inventory, uuid)
        allocated.append(gpu if gpu is not None else {"index": None, "uuid": uuid, "name": None})
    usable = math.ceil(nranks / ranks_per_gpu) if nranks > 0 else 0
    job = allocation["env"].get("SIMQ_JOB_ID")
    held = ",".join(allocation["gpus"]) or "none"
    kind, value = request
    if kind == "none":
        out["notes"].append("--gpus none: every rank runs with CUDA_VISIBLE_DEVICES empty")
    elif kind == "auto":
        out["gpus"] = allocated[:usable]
        if not allocated:
            out["notes"].append(f"--gpus auto: SimQueue job {job} has no GPU allocated; ranks run CPU-only")
    elif kind == "count":
        n = int(value)
        if n > len(allocated):
            refuse(EXIT_BUSY, f"--gpus count:{n}: SimQueue job {job} has {len(allocated)} GPU(s) allocated "
                              f"(SIMQ_GPUS={held}); submit with --gpus {n}")
        elif n > usable:
            refuse(EXIT_USAGE, f"--gpus count:{n}: {nranks} rank(s) at --ranks-per-gpu {ranks_per_gpu} use at most "
                               f"{usable} GPU(s)")
        else:
            out["gpus"] = allocated[:n]
    elif kind == "list":
        selected = []
        for token in value:
            gpu = find_gpu(inventory, token)
            uuid = gpu["uuid"] if gpu is not None else (token if str(token).startswith("GPU-") else None)
            match = next((g for g in allocated if g["uuid"] == uuid), None)
            if match is None:
                refuse(EXIT_BUSY, f"--gpus: GPU {token} is not in SimQueue job {job}'s allocation (SIMQ_GPUS={held})")
                continue
            if match in selected:
                refuse(EXIT_USAGE, f"--gpus: GPU {token} named twice")
                continue
            selected.append(match)
        if len(selected) > usable:
            refuse(EXIT_USAGE, f"--gpus: {len(selected)} GPU(s) selected but {nranks} rank(s) at --ranks-per-gpu "
                               f"{ranks_per_gpu} use at most {usable}")
        if not out["refusals"]:
            out["gpus"] = selected
    else:
        raise ValueError(f"unknown GPU request {kind!r}")
    return out


def select_gpus(inventory, request, nranks: int, ranks_per_gpu: int = 1, reserved=None, busy=None,
                wait=False, inventory_note=None) -> dict:
    """Resolve a GPU request.  Returns {"gpus": [...], "refusals": [...], "notes": [...]}.

    At most ceil(nranks / ranks_per_gpu) GPUs can be used; `auto` trims to that count,
    explicit lists and counts beyond it are refused.
    """
    reserved = reserved or {}
    busy = busy or {}
    out = {"gpus": [], "refusals": [], "notes": []}
    refuse = lambda *a, **k: out["refusals"].append(_refusal(*a, **k))  # noqa: E731
    kind, value = request
    usable = math.ceil(nranks / ranks_per_gpu) if nranks > 0 else 0
    free = [g for g in inventory if g["uuid"] not in reserved and g["uuid"] not in busy]
    unreserved = [g for g in inventory if g["uuid"] not in reserved]
    suffix = f"; {inventory_note}" if inventory_note else ""

    def unavailable_text():
        parts = []
        for gpu in inventory:
            if gpu["uuid"] in reserved:
                parts.append(f"GPU {gpu['index']} reserved ({'; '.join(reserved[gpu['uuid']])})")
            elif gpu["uuid"] in busy:
                parts.append(f"GPU {gpu['index']} leased ({busy[gpu['uuid']]})")
        return "; ".join(parts) or "none reserved or leased"

    if kind == "none":
        out["notes"].append("--gpus none: every rank runs with CUDA_VISIBLE_DEVICES empty")
    elif kind == "auto":
        if not inventory:
            out["notes"].append(f"--gpus auto: no GPU detected{suffix}; ranks run CPU-only")
        elif not free:
            out["notes"].append(f"--gpus auto: no free GPU ({unavailable_text()}); ranks run CPU-only")
        else:
            out["gpus"] = free[:usable]
            if len(free) > usable:
                out["notes"].append(f"--gpus auto: {len(free)} free GPU(s); {usable} used by {nranks} rank(s) at "
                                    f"--ranks-per-gpu {ranks_per_gpu}")
    elif kind == "count":
        n = int(value)
        if n == 0:
            out["notes"].append("--gpus count:0: CPU-only")
        elif n > len(inventory):
            refuse(EXIT_ABSENT, f"--gpus count:{n}: this host has {len(inventory)} GPU(s){suffix}")
        elif n > usable:
            refuse(EXIT_USAGE, f"--gpus count:{n}: {nranks} rank(s) at --ranks-per-gpu {ranks_per_gpu} use at most "
                               f"{usable} GPU(s); select fewer GPUs or more cores")
        elif len(free) < n:
            refuse(EXIT_BUSY, f"--gpus count:{n}: only {len(free)} GPU(s) free ({unavailable_text()})",
                   retryable=len(unreserved) >= n)
        else:
            out["gpus"] = free[:n]
    elif kind == "list":
        selected = []
        for token in value:
            gpu = find_gpu(inventory, token)
            if gpu is None:
                refuse(EXIT_ABSENT, f"--gpus: GPU {token} not present (inventory: {describe_inventory(inventory)}){suffix}")
                continue
            if any(g["uuid"] == gpu["uuid"] for g in selected):
                refuse(EXIT_USAGE, f"--gpus: GPU {gpu['index']} ({gpu['uuid']}) named twice")
                continue
            if gpu["uuid"] in reserved:
                refuse(EXIT_BUSY, f"--gpus: GPU {gpu['index']} ({gpu['uuid']}) is reserved: {'; '.join(reserved[gpu['uuid']])}")
            elif gpu["uuid"] in busy:
                if wait:
                    out["notes"].append(f"GPU {gpu['index']} is leased ({busy[gpu['uuid']]}); --wait queues for it")
                else:
                    refuse(EXIT_BUSY, f"--gpus: GPU {gpu['index']} ({gpu['uuid']}) is leased: lock held on "
                                      f"{busy[gpu['uuid']]}; use --wait to queue for it", retryable=True)
            selected.append(gpu)
        if len(selected) > usable:
            refuse(EXIT_USAGE, f"--gpus: {len(selected)} GPU(s) selected but {nranks} rank(s) at --ranks-per-gpu "
                               f"{ranks_per_gpu} use at most {usable}; select fewer GPUs or more cores")
        if not out["refusals"]:
            out["gpus"] = selected
    else:
        raise ValueError(f"unknown GPU request {kind!r}")
    return out


# --------------------------------------------------------------------------- identity admission


def observe_identity(system: System, site: dict, inventory) -> dict:
    execution = site.get("execution") or {}
    manifest = execution.get("manifest")
    return {
        "kernel_release": system.kernel_release(),
        "nvidia_driver_version": system.read_text(os.path.join(system.sysfs, "module/nvidia/version")),
        "nvidia_srcversion": system.read_text(os.path.join(system.sysfs, "module/nvidia/srcversion")),
        "gpu_uuids": sorted(g["uuid"] for g in inventory),
        "cpu_model": system.cpu_model(),
        "frozen_root_manifest_sha256": system.sha256_file(os.path.expanduser(manifest)) if manifest else None,
    }


def admit_identity(identity_cfg: dict, observed: dict) -> dict:
    """strict: every recorded (non-null) value must match the observed one; record-only: never refuses."""
    identity_cfg = identity_cfg or {}
    policy = identity_cfg.get("policy", "record-only")
    if policy not in ("strict", "record-only"):
        raise PfError(EXIT_CONFIG, f"identity policy must be 'strict' or 'record-only', not {policy!r}")
    checks = []
    for key, expected in sorted((identity_cfg.get("record") or {}).items()):
        if expected is None:
            continue
        got = observed.get(key)
        if key == "gpu_uuids":
            match = sorted(expected) == sorted(got or [])
        else:
            match = expected == got
        checks.append({"key": key, "recorded": expected, "observed": got, "match": match})
    mismatched = [c for c in checks if not c["match"]]
    return {"policy": policy, "admitted": policy != "strict" or not mismatched, "checks": checks,
            "observed": observed}


# --------------------------------------------------------------------------- execution environment


def resolve_execution(site: dict, opts, system: System) -> dict:
    execution = dict(site.get("execution") or {})
    notes = []
    if getattr(opts, "prefix", None) is not None:
        prefix = shlex.split(opts.prefix)
        execution = {"kind": "container" if prefix else "host", "prefix": prefix, "mpi_inside_prefix": False}
        notes.append(f"--prefix: per-rank command prefix {prefix!r}")
    elif getattr(opts, "execution", "site") == "host":
        if execution.get("prefix"):
            notes.append("--execution host: site execution prefix not applied")
        execution = {"kind": "host", "prefix": [], "mpi_inside_prefix": False}
    kind = execution.get("kind", "host")
    prefix = [os.path.expanduser(str(p)) for p in execution.get("prefix") or []]
    inside = False
    root = execution.get("root")
    if kind == "frozen-root" and root and system.same_dir("/", os.path.expanduser(root)):
        inside = True
        prefix = []
        notes.append("already running inside the frozen root; execution prefix not applied")
    return {
        "kind": kind,
        "prefix": prefix,
        "mpi_inside_prefix": bool(execution.get("mpi_inside_prefix", False)),
        "wrapper_python": execution.get("wrapper_python"),
        "root": root,
        "manifest": execution.get("manifest"),
        "image": execution.get("image"),
        "image_digest": execution.get("image_digest"),
        "hash_image": bool(execution.get("hash_image", False)),
        "inside_root": inside,
        "notes": notes,
    }


def root_or_image_digest(execution: dict, observed: dict, system: System) -> dict:
    kind = execution["kind"]
    if kind == "frozen-root":
        launcher = execution["prefix"][0] if execution["prefix"] else None
        return {"kind": kind, "root": execution["root"], "manifest": execution["manifest"],
                "manifest_sha256": observed.get("frozen_root_manifest_sha256"),
                "launcher": launcher, "launcher_sha256": system.sha256_file(launcher) if launcher else None,
                "inside_root": execution["inside_root"]}
    if kind == "container":
        digest = execution.get("image_digest")
        image = execution.get("image")
        if execution.get("hash_image") and image:
            digest = system.sha256_file(os.path.expanduser(image)) or digest
        return {"kind": kind, "prefix": execution["prefix"], "image": image, "image_digest": digest}
    return {"kind": "host", "manifest_sha256": observed.get("frozen_root_manifest_sha256")}


def build_base_env(base: dict, site: dict) -> dict:
    env = dict(base)
    cfg = site.get("environment") or {}
    for key in cfg.get("unset", []) or []:
        env.pop(key, None)
    for key, value in (cfg.get("set") or {}).items():
        env[str(key)] = str(value)
    for key in THREAD_VARS:
        env[key] = "1"
    env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    return env


def env_digest(env: dict) -> str:
    return sha256_bytes(canonical_json(sorted(env.items())).encode())


# --------------------------------------------------------------------------- plan


def load_site(path: str):
    try:
        with open(path, "rb") as stream:
            raw = stream.read()
    except OSError as error:
        raise PfError(EXIT_CONFIG, f"cannot read site file {path}: {error.strerror}; "
                                   "generate one with: pf-launch --detect") from None
    try:
        site = json.loads(raw)
    except ValueError as error:
        raise PfError(EXIT_CONFIG, f"site file {path} is not valid JSON: {error}") from None
    if not isinstance(site, dict) or site.get("schema") != SITE_SCHEMA:
        raise PfError(EXIT_CONFIG, f"site file {path} does not declare schema {SITE_SCHEMA!r}")
    nice = site.get("nice", 0)
    if not isinstance(nice, int) or not 0 <= nice <= 19:
        raise PfError(EXIT_CONFIG, f"site file {path}: nice must be an integer 0..19")
    return site, {"path": os.path.abspath(path), "sha256": sha256_bytes(raw), "schema": SITE_SCHEMA,
                  "site": site.get("site")}


def default_site_path(system: System) -> str:
    env_path = system.environ.get("PF_LAUNCH_SITE")
    if env_path:
        return env_path
    return os.path.join(os.path.dirname(SELF_PATH), "sites", f"{system.hostname()}.json")


def parse_core_request(opts):
    if opts.cpus is not None:
        return ("list", parse_cpu_list(opts.cpus))
    if opts.cores is None or str(opts.cores).strip().lower() == "auto":
        return ("auto", None)
    text = str(opts.cores).strip()
    if not re.fullmatch(r"\d+", text) or int(text) < 1:
        raise ValueError(f"--cores expects a positive integer or 'auto', not {opts.cores!r}")
    return ("count", int(text))


def request_text(opts) -> str:
    if opts.cpus is not None:
        return f"--cpus {opts.cpus}"
    return f"--cores {opts.cores if opts.cores is not None else 'auto'}"


def build_plan(opts, site: dict, site_info: dict, system: System) -> dict:
    """Resolve every decision of a launch without side effects (probes and lease tests only)."""
    env = system.environ
    plan = {"refusals": [], "notes": [], "site": site_info, "pf_launch_argv": list(getattr(opts, "argv", []) or [])}
    refuse = lambda *a, **k: plan["refusals"].append(_refusal(*a, **k))  # noqa: E731
    launcher_cfg = site.get("launcher") or {}

    try:
        core_request = parse_core_request(opts)
    except ValueError as error:
        refuse(EXIT_USAGE, str(error))
        core_request = None
    if opts.gpus is None:
        opts.gpus = site.get("default_gpus") or "auto"
    try:
        gpu_request = parse_gpu_request(opts.gpus)
    except ValueError as error:
        refuse(EXIT_USAGE, str(error))
        gpu_request = None
    ranks_per_gpu = int(opts.ranks_per_gpu)
    if ranks_per_gpu < 1:
        refuse(EXIT_USAGE, "--ranks-per-gpu must be at least 1")
        ranks_per_gpu = 1
    plan["core_request"] = request_text(opts)
    plan["gpu_request"] = opts.gpus
    plan["ranks_per_gpu"] = ranks_per_gpu

    topology = system.topology()
    inventory, inventory_note = system.gpu_inventory(launcher_cfg.get("nvidia_smi", "nvidia-smi"))
    plan["inventory"] = inventory
    if inventory_note:
        plan["notes"].append(f"GPU inventory: {inventory_note}")
    plan["topology"] = {"logical_cpus": format_cpu_list(topology.logical_cpus()),
                        "physical_cores": len(topology.cores), "cpu_model": topology.cpu_model}

    observed = observe_identity(system, site, inventory)
    try:
        identity = admit_identity(site.get("identity") or {}, observed)
    except PfError as error:
        refuse(error.code, error.message)
        identity = {"policy": None, "admitted": False, "checks": [], "observed": observed}
    plan["identity"] = identity
    if not identity["admitted"]:
        bad = [f"{c['key']}: recorded {c['recorded']!r}, observed {c['observed']!r}"
               for c in identity["checks"] if not c["match"]]
        refuse(EXIT_CONFIG, "identity admission (strict) failed: " + "; ".join(bad))

    execution = resolve_execution(site, opts, system)
    plan["execution"] = execution
    plan["notes"].extend(execution["notes"])
    if execution["prefix"] and "/" in execution["prefix"][0] and not system.exists(execution["prefix"][0]):
        refuse(EXIT_CONFIG, f"execution prefix {execution['prefix'][0]} not found")

    in_slurm = bool(env.get("SLURM_JOB_ID"))
    launcher = opts.launcher
    if launcher == "auto" and in_slurm and launcher_cfg.get("srun_in_allocation", True):
        launcher = "srun"
    slurm_mode = launcher == "srun"
    plan["slurm"] = None
    use_leases = not slurm_mode  # inside an allocation the scheduler owns cores and GPUs
    simq = None if slurm_mode or in_slurm else simqueue_allocation(env)
    plan["simqueue"] = None if simq is None else dict(simq["env"], gpus_managed=simq["gpus_managed"])
    use_core_leases = use_leases and simq is None
    use_gpu_leases = use_leases and (simq is None or not simq["gpus_managed"])

    # reservations and leases (local placement only)
    if slurm_mode:
        reservations = {"cores": {}, "gpus": {}, "probes": []}
        plan["notes"].append("Slurm allocation: site reservations not evaluated; the allocation owns the cores")
    elif simq is not None:
        reservations = {"cores": {}, "gpus": {}, "probes": []}
        plan["notes"].append(f"SimQueue allocation (SIMQ_JOB_ID={simq['env'].get('SIMQ_JOB_ID')}): site reservations "
                             "not evaluated and no core leases taken; the queue's job unit holds them")
    else:
        reservations = evaluate_reservations(site, topology, system, inventory)
        cg = system.cgroup_cpus()
        if cg is not None:
            for core, rec in topology.cores.items():
                if core not in cg:
                    reservations["cores"].setdefault(core, []).append("outside this process's cgroup cpuset")
    plan["reservations"] = reservations
    busy_cores, busy_gpus = {}, {}
    if use_core_leases:
        for core in topology.primaries():
            path = core_lease_path(site, core)
            if path and system.lease_busy(path):
                busy_cores[core] = path
    if use_gpu_leases:
        for gpu in inventory:
            for path in gpu_lease_paths(site, gpu):
                if system.lease_busy(path):
                    busy_gpus[gpu["uuid"]] = path
                    break
    plan["busy"] = {"cores": busy_cores, "gpus": busy_gpus}

    # cores
    ranks = []
    if core_request is not None and slurm_mode:
        ntasks_text = env.get("SLURM_NTASKS") or env.get("SLURM_NPROCS")
        ntasks = int(ntasks_text) if ntasks_text and ntasks_text.isdigit() else None
        if core_request[0] == "list":
            refuse(EXIT_USAGE, f"--cpus is not available inside a Slurm allocation (SLURM_JOB_ID="
                               f"{env.get('SLURM_JOB_ID')}): the allocation and srun binding choose the cores; "
                               f"use --cores N")
            n = 0
        elif core_request[0] == "count":
            n = core_request[1]
            if ntasks is not None and n > ntasks:
                if opts.cores_available:
                    plan["notes"].append(f"--cores {n}: allocation has {ntasks} task(s); running {ntasks} "
                                         f"(--cores-available)")
                    n = ntasks
                else:
                    refuse(EXIT_BUSY, f"--cores {n} exceeds the allocation's {ntasks} task(s) (SLURM_NTASKS); "
                                      f"use --cores {ntasks} or --cores-available")
        else:
            if ntasks is not None:
                n = ntasks
            else:
                allowed = set(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else set()
                n = max(1, len({topology.cpu_core[c] for c in allowed if c in topology.cpu_core}))
                plan["notes"].append(f"--cores auto: SLURM_NTASKS unset; {n} rank(s) from this node's "
                                     f"allowed physical cores")
        ranks = [{"rank": r, "cpu": None, "core": None, "siblings": None} for r in range(n)]
        plan["slurm"] = {k: env.get(k) for k in ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NNODES",
                                                  "SLURM_NTASKS", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE",
                                                  "SLURM_JOB_GPUS", "SLURM_STEP_GPUS")}
    elif core_request is not None and simq is not None:
        selection = select_simq_cores(topology, core_request, simq, allow_fewer=opts.cores_available)
        plan["refusals"].extend(selection["refusals"])
        plan["notes"].extend(selection["notes"])
        ranks = selection["ranks"]
    elif core_request is not None:
        selection = select_cores(topology, core_request, reservations["cores"], busy_cores,
                                 preference=site.get("core_preference"), allow_fewer=opts.cores_available,
                                 wait=opts.wait)
        plan["refusals"].extend(selection["refusals"])
        plan["notes"].extend(selection["notes"])
        ranks = selection["ranks"]
    ncores = len(ranks)

    # launcher and process layout
    if launcher == "auto":
        launcher = "exec" if ncores <= 1 else "mpirun"
    if launcher == "srun" and not in_slurm:
        plan["notes"].append("--launcher srun outside an allocation: srun requests one")
    plan["launcher"] = launcher
    group = launcher == "exec" and ncores > 1
    plan["group"] = group
    if group:
        plan["notes"].append(f"--launcher exec with {ncores} cores: one process pinned to CPUs "
                             f"{','.join(str(r['cpu']) for r in ranks)} (not MPI ranks)")

    # GPUs
    if gpu_request is None:
        gpu_sel = {"gpus": [], "refusals": [], "notes": []}
    elif slurm_mode:
        gpu_sel = {"gpus": [], "refusals": [], "notes": [
            f"Slurm: GPUs assigned per node from the step's visible GPUs, rank SLURM_LOCALID // {ranks_per_gpu}"]}
    elif simq is not None and simq["gpus_managed"]:
        gpu_sel = select_simq_gpus(inventory, gpu_request, ncores, ranks_per_gpu, simq)
    else:
        n_for_gpus = ncores
        if not n_for_gpus and core_request is not None:
            n_for_gpus = core_request[1] if core_request[0] == "count" else len(core_request[1] or [None])
        gpu_sel = select_gpus(inventory, gpu_request, n_for_gpus, ranks_per_gpu, reservations["gpus"], busy_gpus,
                              wait=opts.wait, inventory_note=inventory_note)
    plan["refusals"].extend(gpu_sel["refusals"])
    plan["notes"].extend(gpu_sel["notes"])
    gpus = gpu_sel["gpus"]
    plan["gpus"] = gpus

    # processes (one per rank, or one pinned group)
    processes = []
    if group:
        uuids = [g["uuid"] for g in gpus]
        processes.append({"rank": 0, "cpus": [r["cpu"] for r in ranks], "cores": [r["core"] for r in ranks],
                          "siblings": [r["siblings"] for r in ranks], "gpus": uuids,
                          "cuda_visible_devices": ",".join(uuids)})
    else:
        gpu_map = rank_gpu_map(len(ranks), gpus, ranks_per_gpu)
        for r, uuid in zip(ranks, gpu_map):
            processes.append({"rank": r["rank"], "cpus": [r["cpu"]] if r["cpu"] is not None else None,
                              "cores": [r["core"]] if r["core"] is not None else None,
                              "siblings": [r["siblings"]] if r["siblings"] is not None else None,
                              "gpus": [uuid] if uuid else [], "cuda_visible_devices": uuid or ""})
    plan["processes"] = processes
    plan["ranks"] = ranks
    if slurm_mode and gpu_request is not None:
        plan["slurm_gpu_request"] = {"kind": gpu_request[0], "value": gpu_request[1]}

    # leases to take, nice, environment
    order = (lease_order(site, [r["core"] for r in ranks] if use_core_leases else [],
                         gpus if use_gpu_leases else []) if use_leases else [])
    core_paths = {core_lease_path(site, r["core"]) for r in ranks if r["core"] is not None}
    plan["core_leases"] = [p for p in order if p in core_paths]
    plan["gpu_leases"] = [p for p in order if p not in core_paths]
    plan["leases"] = order
    plan["lease_block"] = bool(opts.wait)
    plan["waits"] = [p for p in order if system.lease_busy(p)] if opts.wait else []
    target_nice = int(site.get("nice", 0))
    current_nice = system.current_nice()
    plan["nice"] = {"site": target_nice, "current": current_nice, "effective": max(target_nice, current_nice)}
    if current_nice > target_nice:
        plan["notes"].append(f"current nice {current_nice} is above the site's {target_nice}; it cannot be lowered")
    base_env = build_base_env(env, site)
    nprocs = len(processes)
    base_env["PF_LAUNCH_NRANKS"] = str(nprocs)
    if group:
        base_env["MAKEFLAGS"] = f"-j{ncores}"
        base_env["CMAKE_BUILD_PARALLEL_LEVEL"] = str(ncores)
    plan["env"] = base_env
    env_cfg = site.get("environment") or {}
    plan["env_record"] = {"site_set": env_cfg.get("set") or {}, "site_unset": env_cfg.get("unset") or [],
                          "forced": {**{k: "1" for k in THREAD_VARS}, "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
                                     **({"MAKEFLAGS": f"-j{ncores}", "CMAKE_BUILD_PARALLEL_LEVEL": str(ncores)}
                                        if group else {})},
                          "per_rank": list(PER_RANK_ENV_KEYS)}
    plan["launcher_cfg"] = launcher_cfg
    return plan


def process_env(plan: dict, process: dict) -> dict:
    env = dict(plan["env"])
    env["CUDA_VISIBLE_DEVICES"] = process["cuda_visible_devices"]
    env["PF_LAUNCH_RANK"] = str(process["rank"])
    if process["cpus"]:
        env["PF_LAUNCH_CPU"] = ",".join(str(c) for c in process["cpus"])
        env["PF_CORE"] = env["PF_LAUNCH_CPU"]
    return env


def rank_plan_payload(plan: dict, report_dir) -> dict:
    payload = {"schema": RANK_PLAN_SCHEMA, "mode": plan["launcher"], "nranks": len(plan["processes"]),
               "nice": plan["nice"]["effective"], "report_dir": report_dir, "ranks_per_gpu": plan["ranks_per_gpu"]}
    if plan["launcher"] == "srun":
        payload["gpu_request"] = plan.get("slurm_gpu_request") or {"kind": "none", "value": None}
    else:
        payload["cpus"] = [p["cpus"][0] for p in plan["processes"]]
        payload["gpus"] = [p["cuda_visible_devices"] for p in plan["processes"]]
    return payload


def build_child_argv(plan: dict, program: list, report_dir=None, self_path=SELF_PATH) -> list:
    execution = plan["execution"]
    prefix = list(execution["prefix"])
    launcher = plan["launcher"]
    cfg = plan.get("launcher_cfg") or {}
    program = list(program) if program else ["<PROGRAM>"]
    if launcher == "exec":
        return prefix + program
    inside = bool(prefix) and execution["mpi_inside_prefix"] and launcher == "mpirun"
    if inside:
        wrapper_python = execution.get("wrapper_python") or "python3"
    else:
        wrapper_python = cfg.get("wrapper_python") or sys.executable or "python3"
    wrapper = [wrapper_python, self_path, "--rank-exec", canonical_json(rank_plan_payload(plan, report_dir)), "--"]
    n = str(len(plan["processes"]))
    if launcher == "mpirun":
        mpi = [cfg.get("mpirun") or "mpirun", "--bind-to", "none", "-np", n] + list(cfg.get("mpirun_args") or [])
        if inside:
            return prefix + mpi + wrapper + program
        return mpi + wrapper + prefix + program
    if launcher == "srun":
        srun = [cfg.get("srun") or "srun", f"--ntasks={n}", "--cpus-per-task=1"] + list(cfg.get("srun_args") or [])
        return srun + wrapper + prefix + program
    raise PfError(EXIT_USAGE, f"unknown launcher {launcher!r}")


def plan_exit_code(plan: dict) -> int:
    """First refusal's code; EXIT_BUSY when the plan is valid but --wait would block on held leases."""
    if plan["refusals"]:
        return plan["refusals"][0]["code"]
    return EXIT_BUSY if plan.get("waits") else 0


def dry_run_status(plan: dict) -> str:
    if plan["refusals"]:
        return "dry-run-refused"
    return "dry-run-wait" if plan.get("waits") else "dry-run"


# --------------------------------------------------------------------------- receipt


def build_receipt(plan: dict, *, status: str, pf_code: int, child_argv, start_time, end_time,
                  launcher_rc=None, rank_reports=None, messages=None, leases_held=None, system=None,
                  observed_child=None) -> dict:
    system = system or System()
    identity = plan.get("identity") or {}
    observed = identity.get("observed") or {}
    rank_reports = rank_reports or {}
    processes = plan.get("processes") or []
    rank_to_core, rank_to_gpu = [], []
    for proc in processes:
        report = rank_reports.get(proc["rank"])
        if report is None and observed_child is not None and proc["rank"] == 0:
            report = observed_child
        rank_to_core.append({"rank": proc["rank"], "cpus": proc["cpus"], "physical_cores": proc["cores"],
                             "smt_siblings": proc["siblings"],
                             "binding": ("srun --cpu-bind (allocation)" if plan.get("launcher") == "srun"
                                         else "sched_setaffinity"),
                             "observed": report})
        rank_to_gpu.append({"rank": proc["rank"], "gpus": proc["gpus"],
                            "cuda_visible_devices": (None if plan.get("launcher") == "srun"
                                                     else proc["cuda_visible_devices"]),
                            "observed_cuda_visible_devices": (report or {}).get("cuda_visible_devices")})
    rank_codes = None
    if rank_reports:
        rank_codes = [{"rank": r, "host": rep.get("host"), "exit_code": rep.get("exit_code")}
                      for r, rep in sorted(rank_reports.items())]
    elif plan.get("launcher") == "exec" and launcher_rc is not None:
        rank_codes = [{"rank": 0, "host": system.hostname(), "exit_code": launcher_rc}]
    held = leases_held if leases_held is not None else []
    execution = plan.get("execution") or {"kind": "host", "prefix": [], "inside_root": False, "root": None,
                                          "manifest": None}
    env = plan.get("env") or {}
    return {
        "host": {"name": system.hostname(), "site": (plan.get("site") or {}).get("site"),
                 "kernel_release": observed.get("kernel_release"), "cpu_model": observed.get("cpu_model"),
                 "topology": plan.get("topology"),
                 "identity": {"policy": identity.get("policy"), "admitted": identity.get("admitted"),
                              "checks": identity.get("checks", [])}},
        "site_file_digest": plan.get("site"),
        "selected_cores": {"request": plan.get("core_request"), "cpus": [r["cpu"] for r in plan.get("ranks", [])],
                           "physical_cores": [r["core"] for r in plan.get("ranks", [])],
                           "leases": plan.get("core_leases", []),
                           "leases_held": [p for p in held if p in plan.get("core_leases", [])],
                           "reserved": {str(k): v for k, v in (plan.get("reservations") or {}).get("cores", {}).items()},
                           "reservation_probes": (plan.get("reservations") or {}).get("probes", []),
                           "leased_elsewhere": {str(k): v for k, v in (plan.get("busy") or {}).get("cores", {}).items()},
                           "slurm": plan.get("slurm"),
                           "simqueue": plan.get("simqueue"),
                           "nice": plan.get("nice")},
        "selected_gpus": {"request": plan.get("gpu_request"), "uuids": [g["uuid"] for g in plan.get("gpus", [])],
                          "ranks_per_gpu": plan.get("ranks_per_gpu"),
                          "leases": plan.get("gpu_leases", []),
                          "leases_held": [p for p in held if p in plan.get("gpu_leases", [])],
                          "leased_elsewhere": (plan.get("busy") or {}).get("gpus", {})},
        "rank_count": len(processes),
        "rank_to_core": rank_to_core,
        "rank_to_gpu": rank_to_gpu,
        "gpus": [{k: g.get(k) for k in ("index", "uuid", "name", "compute_cap", "memory_total", "driver_version",
                                          "pci_bus_id", "numa")} for g in plan.get("gpus", [])],
        "launcher_argv": {"pf_launch": plan.get("pf_launch_argv"), "pf_launch_version": VERSION,
                          "pf_launch_sha256": system.sha256_file(SELF_PATH), "launcher": plan.get("launcher"),
                          "argv": child_argv, "cwd": os.getcwd(),
                          "dry_run": status.startswith("dry-run")},
        "environment_digest": {"sha256": env_digest(env) if env else None,
                               "variables": plan.get("env_record"),
                               "notes": plan.get("notes", [])},
        "root_or_image_digest": root_or_image_digest(execution, observed, system) if plan.get("execution") else None,
        "start_time": start_time,
        "end_time": end_time,
        "exit_codes": {"status": status, "pf_launch": pf_code, "launcher": launcher_rc, "ranks": rank_codes,
                       "messages": messages or [r["message"] for r in plan.get("refusals", [])]},
    }


def plan_cores(plan):
    return [r["core"] for r in plan.get("ranks", []) if r.get("core") is not None]


def write_json_atomic(path: str, value) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    tmp = f"{path}.tmp-{os.getpid()}"
    with open(tmp, "w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=False)
        stream.write("\n")
    os.replace(tmp, path)


# --------------------------------------------------------------------------- output


def format_plan(plan: dict, program, dry_run: bool) -> str:
    lines = [f"pf-launch {VERSION} plan{' (dry run)' if dry_run else ''}"]
    site = plan.get("site") or {}
    lines.append(f"  site       {site.get('site')}  {site.get('path')}  sha256 {str(site.get('sha256'))[:16]}")
    identity = plan.get("identity") or {}
    checks = identity.get("checks", [])
    ok = [c["key"] for c in checks if c["match"]]
    bad = [c["key"] for c in checks if not c["match"]]
    lines.append(f"  identity   {identity.get('policy')}: {'admitted' if identity.get('admitted') else 'REFUSED'}"
                 f" (match: {', '.join(ok) or '-'}{'; mismatch: ' + ', '.join(bad) if bad else ''})")
    execution = plan.get("execution") or {}
    where = execution.get("kind")
    if execution.get("prefix"):
        where += f" via {' '.join(execution['prefix'])}"
    if plan.get("processes"):
        lines.append(f"  launcher   {plan.get('launcher')}, {len(plan['processes'])} process(es); execution {where}")
    else:
        lines.append(f"  launcher   not resolved (no ranks); execution {where}")
    lines.append(f"  cores      {plan.get('core_request')} -> CPUs "
                 f"{','.join(str(r['cpu']) for r in plan.get('ranks', []) if r['cpu'] is not None) or '-'}")
    simq = plan.get("simqueue")
    if simq:
        lines.append(f"             SimQueue job {simq.get('SIMQ_JOB_ID')} partition {simq.get('SIMQ_PARTITION')}: "
                     f"CPUs {simq.get('SIMQ_CPU_LIST')}, GPUs {simq.get('SIMQ_GPUS') or 'none'}"
                     f"{'' if simq.get('gpus_managed') else ' (0.5 queue: GPUs leased locally)'}")
    reservations = plan.get("reservations") or {}
    if reservations.get("cores"):
        lines.append(f"             reserved: {describe_reasons(reservations['cores'])}")
    for probe in reservations.get("probes", []):
        lines.append(f"             probe {probe['name']} (cores {format_cpu_list(probe['cores']) or '-'}): "
                     f"{'RESERVED' if probe['reserved'] else 'free'}: {probe['detail']}")
    busy = plan.get("busy") or {}
    if busy.get("cores"):
        lines.append("             leased elsewhere: " + "; ".join(f"core {c} ({p})" for c, p in sorted(busy["cores"].items())))
    lines.append(f"  gpus       {plan.get('gpu_request')} -> "
                 + (", ".join(f"{g.get('index')} {g['uuid']} ({g.get('name')}, cc {g.get('compute_cap')}, {g.get('memory_total')}, "
                              f"driver {g.get('driver_version')})" for g in plan.get("gpus", [])) or "none"))
    for proc in plan.get("processes") or []:
        cpus = ",".join(str(c) for c in proc["cpus"]) if proc["cpus"] else "srun binding"
        sib = ""
        if proc.get("siblings"):
            sib = f" (physical core {','.join(str(c) for c in proc['cores'])}; siblings " + \
                  " ".join(format_cpu_list(s) for s in proc["siblings"]) + ")"
        cvd = proc["cuda_visible_devices"] if plan.get("launcher") != "srun" else "per node"
        lines.append(f"  rank {proc['rank']:<4}  cpu {cpus}{sib}  CUDA_VISIBLE_DEVICES={cvd}")
    nice = plan.get("nice") or {}
    lines.append(f"  nice       {nice.get('effective')}; threads per rank 1 ({'='.join(THREAD_VARS)}=1)")
    if plan.get("leases") and plan.get("processes"):
        lines.append(f"  leases     {'blocking (--wait)' if plan.get('lease_block') else 'non-blocking'}, in order: "
                     + ", ".join(plan["leases"]))
    if plan.get("processes"):
        lines.append("  argv       " + shlex.join(build_child_argv(plan, program)))
    for note in plan.get("notes", []):
        lines.append(f"  note       {note}")
    if plan.get("refusals"):
        for refusal in plan["refusals"]:
            lines.append(f"  REFUSED    (exit {refusal['code']}) {refusal['message']}")
        if plan.get("lease_block") and all(r["retryable"] for r in plan["refusals"]):
            lines.append(f"  wait       with --wait a real run re-plans every {WAIT_POLL_SECONDS:g} s until "
                         f"enough cores/GPUs are free")
    elif plan.get("waits"):
        lines.append(f"  result     WAIT (exit {EXIT_BUSY}): with --wait a real run blocks on "
                     + ", ".join(plan["waits"]))
    else:
        lines.append("  result     OK")
    return "\n".join(lines)


# --------------------------------------------------------------------------- run


def _set_pdeathsig():
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        libc.prctl(1, int(signal.SIGKILL), 0, 0, 0)  # PR_SET_PDEATHSIG
    except Exception:  # noqa: BLE001 - best effort
        pass


def _normalize_rc(rc: int) -> int:
    return 128 + (-rc) if rc < 0 else rc


def _wait_forwarding(proc) -> int:
    def forward(signum, _frame):
        if proc.poll() is None:
            try:
                proc.send_signal(signum)
            except OSError:
                pass

    previous = {}
    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP, signal.SIGUSR1, signal.SIGUSR2):
        try:
            previous[sig] = signal.signal(sig, forward)
        except (ValueError, OSError):
            pass
    try:
        return proc.wait()
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def default_receipt_path(site: dict, system: System) -> str:
    directory = os.path.expanduser(site.get("receipt_dir") or os.path.join(os.getcwd(), "pf-launch-receipts"))
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return os.path.join(directory, f"pf-launch-{system.hostname()}-{stamp}-{os.getpid()}.json")


def collect_rank_reports(report_dir) -> dict:
    reports = {}
    if not report_dir or not os.path.isdir(report_dir):
        return reports
    for name in sorted(os.listdir(report_dir)):
        match = re.fullmatch(r"rank-(\d+)\.json", name)
        if not match:
            continue
        path = os.path.join(report_dir, name)
        try:
            with open(path, encoding="utf-8") as stream:
                reports[int(match.group(1))] = json.load(stream)
            os.unlink(path)
        except (OSError, ValueError):
            pass
    try:
        os.rmdir(report_dir)
    except OSError:
        pass
    return reports


def check_rank_overlap(reports: dict, topology: Topology, hostname: str) -> list:
    """Warn when two ranks on this host share a physical core (possible under srun misconfiguration)."""
    owners, warnings = {}, []
    for rank, report in sorted(reports.items()):
        if report.get("host") != hostname:
            continue
        for cpu in report.get("affinity") or []:
            core = topology.cpu_core.get(cpu)
            if core is None:
                continue
            other = owners.setdefault(core, rank)
            if other != rank:
                warnings.append(f"ranks {other} and {rank} share physical core {core}")
    return sorted(set(warnings))


def acquire_plan_leases(plan: dict, block: bool, err):
    """Take the plan's leases in its fixed order.  Returns (LeaseSet, None) or (None, busy_path)."""
    leases = LeaseSet()
    try:
        for path in plan["leases"]:
            if leases.acquire(path, block=False):
                continue
            if not block:
                leases.release()
                return None, path
            print(f"pf-launch: waiting for lease {path}", file=err)
            leases.acquire(path, block=True)
    except BaseException:
        leases.release()
        raise
    return leases, None


def launch(opts, program, site, site_info, system, out, err) -> int:
    announced = False
    while True:
        plan = build_plan(opts, site, site_info, system)
        if plan["refusals"] and opts.wait and all(r["retryable"] for r in plan["refusals"]):
            if not announced:
                print("pf-launch: waiting: " + "; ".join(r["message"] for r in plan["refusals"]), file=err)
                announced = True
            time.sleep(WAIT_POLL_SECONDS)
            continue
        if plan["refusals"]:
            print(format_plan(plan, program, dry_run=False), file=err)
            if opts.receipt:
                write_json_atomic(opts.receipt, build_receipt(plan, status="refused", pf_code=plan_exit_code(plan),
                                                              child_argv=None, start_time=system.now(),
                                                              end_time=None, system=system))
            return plan_exit_code(plan)
        try:
            leases, busy_path = acquire_plan_leases(plan, opts.wait, err)
        except OSError as error:
            print(f"pf-launch: refused: cannot take lease: {error}", file=err)
            return EXIT_CONFIG
        if leases is None:
            print(f"pf-launch: refused: lease {busy_path} was taken by another process after planning; "
                  f"use --wait to queue", file=err)
            return EXIT_BUSY
        break

    if opts.wait and not plan.get("simqueue"):
        recheck = evaluate_reservations(site, system.topology(), system, plan["inventory"])
        now_reserved = [c for c in plan_cores(plan) if c in recheck["cores"]]
        if now_reserved:
            leases.release()
            print(f"pf-launch: refused: cores {format_cpu_list(now_reserved)} became reserved while waiting: "
                  f"{describe_reasons({c: recheck['cores'][c] for c in now_reserved})}", file=err)
            return EXIT_BUSY

    receipt_path = opts.receipt or default_receipt_path(site, system)
    report_dir = None
    if plan["launcher"] in ("mpirun", "srun"):
        report_dir = receipt_path + ".ranks"
        os.makedirs(report_dir, exist_ok=True)
    argv = build_child_argv(plan, program, report_dir)
    if plan["launcher"] == "exec":
        child_env = process_env(plan, plan["processes"][0])
    else:
        child_env = dict(plan["env"])
        child_env["CUDA_VISIBLE_DEVICES"] = ""
    union = sorted({c for p in plan["processes"] for c in (p["cpus"] or [])})
    nice = plan["nice"]["effective"]
    procs_text = ",".join(map(str, union)) if union else "srun binding"
    print(f"pf-launch: {plan['launcher']} {len(plan['processes'])} process(es) on CPUs {procs_text}, GPUs "
          f"{','.join(g['uuid'] for g in plan['gpus']) or 'none'}, nice {nice}; receipt {receipt_path}", file=err)

    def preexec():
        os.setpriority(os.PRIO_PROCESS, 0, nice)
        if union:
            os.sched_setaffinity(0, union)

    start = system.now()
    observed_child = None
    try:
        proc = subprocess.Popen(argv, env=child_env, preexec_fn=preexec, pass_fds=leases.fds())
    except (OSError, subprocess.SubprocessError) as error:
        end = system.now()
        leases_held = leases.paths()
        leases.release()
        message = f"cannot start {argv[0]!r}: {error}"
        print(f"pf-launch: {message}", file=err)
        write_json_atomic(receipt_path, build_receipt(plan, status="failed-to-start", pf_code=EXIT_SOFTWARE,
                                                      child_argv=argv, start_time=start, end_time=end,
                                                      messages=[message], leases_held=leases_held, system=system))
        return EXIT_SOFTWARE
    if plan["launcher"] == "exec":
        try:
            observed_child = {"pid": proc.pid, "host": system.hostname(),
                              "affinity": sorted(os.sched_getaffinity(proc.pid)),
                              "nice": os.getpriority(os.PRIO_PROCESS, proc.pid),
                              "cuda_visible_devices": child_env.get("CUDA_VISIBLE_DEVICES"),
                              "source": "launcher view of the child process"}
        except OSError:
            observed_child = None
    rc = _wait_forwarding(proc)
    end = system.now()
    reports = collect_rank_reports(report_dir)
    messages = check_rank_overlap(reports, system.topology(), system.hostname()) if reports else []
    if report_dir and not reports:
        messages.append(f"no rank reports in {report_dir} (directory not visible to the ranks?)")
    leases_held = leases.paths()
    receipt = build_receipt(plan, status="completed", pf_code=_normalize_rc(rc), child_argv=argv, start_time=start,
                            end_time=end, launcher_rc=rc, rank_reports=reports, messages=messages,
                            leases_held=leases_held, system=system, observed_child=observed_child)
    write_json_atomic(receipt_path, receipt)
    leases.release()
    print(f"pf-launch: exit {rc}; receipt {receipt_path}", file=err)
    return _normalize_rc(rc)


# --------------------------------------------------------------------------- per-rank wrapper


def first_int(env: dict, keys):
    for key in keys:
        value = env.get(key)
        if value is not None and value.strip().isdigit():
            return int(value)
    return None


def slurm_visible_gpus(env: dict) -> list:
    for key in SLURM_GPU_ENV_KEYS:
        value = env.get(key)
        if value is not None and value.strip() and value.strip() != "NoDevFiles":
            return [v.strip() for v in value.split(",") if v.strip()]
    return []


def slurm_gpu_for_rank(local_rank: int, ranks_per_gpu: int, request: dict, env: dict) -> str:
    kind, value = request.get("kind"), request.get("value")
    if kind == "none":
        return ""
    visible = slurm_visible_gpus(env)
    if kind == "auto":
        pool = visible
    elif kind == "count":
        pool = visible[: int(value)]
    else:
        pool = []
        for token in value or []:
            if isinstance(token, int):
                if token < len(visible):
                    pool.append(visible[token])
            else:
                pool.append(token)
    slot = local_rank // max(1, ranks_per_gpu)
    return pool[slot] if slot < len(pool) else ""


def rank_assignment(rank_plan: dict, env: dict) -> dict:
    if rank_plan.get("schema") != RANK_PLAN_SCHEMA:
        raise PfError(EXIT_USAGE, "rank plan schema mismatch")
    nranks = int(rank_plan["nranks"])
    rank = first_int(env, RANK_ENV_KEYS)
    if rank is None:
        if nranks == 1:
            rank = 0
        else:
            raise PfError(EXIT_USAGE, f"no rank variable ({', '.join(RANK_ENV_KEYS)}) in the environment")
    local = first_int(env, LOCAL_RANK_ENV_KEYS)
    local = rank if local is None else local
    if rank >= nranks:
        raise PfError(EXIT_USAGE, f"rank {rank} outside the plan's {nranks} rank(s)")
    if rank_plan["mode"] == "srun":
        cpu = None
        cuda = slurm_gpu_for_rank(local, int(rank_plan.get("ranks_per_gpu", 1)), rank_plan.get("gpu_request") or {}, env)
    else:
        cpu = int(rank_plan["cpus"][rank])
        cuda = rank_plan["gpus"][rank] or ""
    updates = {k: "1" for k in THREAD_VARS}
    updates.update({"CUDA_VISIBLE_DEVICES": cuda, "PF_LAUNCH_RANK": str(rank), "PF_LAUNCH_NRANKS": str(nranks)})
    if cpu is not None:
        updates["PF_LAUNCH_CPU"] = str(cpu)
        updates["PF_CORE"] = str(cpu)
    return {"rank": rank, "local_rank": local, "cpu": cpu, "cuda_visible_devices": cuda, "env": updates}


def rank_exec_main(args) -> int:
    if len(args) < 3 or args[1] != "--":
        print("pf-launch --rank-exec: internal use: --rank-exec PLAN_JSON -- PROGRAM ARGS", file=sys.stderr)
        return EXIT_USAGE
    try:
        rank_plan = json.loads(args[0])
        assignment = rank_assignment(rank_plan, dict(os.environ))
    except (ValueError, KeyError, IndexError) as error:
        print(f"pf-launch rank wrapper: bad rank plan: {error}", file=sys.stderr)
        return EXIT_USAGE
    except PfError as error:
        print(f"pf-launch rank wrapper: {error.message}", file=sys.stderr)
        return error.code
    program = args[2:]
    try:
        if assignment["cpu"] is not None:
            os.sched_setaffinity(0, {assignment["cpu"]})
        target = int(rank_plan.get("nice", 0))
        if os.getpriority(os.PRIO_PROCESS, 0) < target:
            os.setpriority(os.PRIO_PROCESS, 0, target)
    except OSError as error:
        print(f"pf-launch rank wrapper: rank {assignment['rank']}: cannot pin to CPU {assignment['cpu']} "
              f"or set nice: {error}", file=sys.stderr)
        return EXIT_SOFTWARE
    env = dict(os.environ)
    env.update(assignment["env"])
    report_dir = rank_plan.get("report_dir")
    report = {"rank": assignment["rank"], "local_rank": assignment["local_rank"], "host": socket.gethostname(),
              "wrapper_pid": os.getpid(), "affinity": sorted(os.sched_getaffinity(0)),
              "nice": os.getpriority(os.PRIO_PROCESS, 0), "cuda_visible_devices": assignment["cuda_visible_devices"],
              "start_time": utc_now(), "exit_code": None, "end_time": None}
    report_path = os.path.join(report_dir, f"rank-{assignment['rank']}.json") if report_dir else None

    def save():
        if report_path:
            try:
                write_json_atomic(report_path, report)
            except OSError as error:
                print(f"pf-launch rank wrapper: cannot write {report_path}: {error}", file=sys.stderr)

    save()
    parent = os.getpid()

    def preexec():
        _set_pdeathsig()
        if os.getppid() != parent:
            os._exit(EXIT_SOFTWARE)

    try:
        proc = subprocess.Popen(program, env=env, preexec_fn=preexec)
    except OSError as error:
        print(f"pf-launch rank wrapper: cannot start {program[0]!r}: {error}", file=sys.stderr)
        report["exit_code"] = 127
        report["end_time"] = utc_now()
        save()
        return 127
    report["pid"] = proc.pid
    rc = _wait_forwarding(proc)
    report["exit_code"] = rc
    report["end_time"] = utc_now()
    save()
    return _normalize_rc(rc)


# --------------------------------------------------------------------------- detect


def detect_site(system: System, nvidia_smi="nvidia-smi") -> dict:
    topology = system.topology()
    inventory, note = system.gpu_inventory(nvidia_smi)
    host = system.hostname()
    observed = observe_identity(system, {}, inventory)
    record = {k: observed[k] for k in IDENTITY_KEYS if k != "frozen_root_manifest_sha256"}
    return {
        "schema": SITE_SCHEMA,
        "site": host,
        "hostname": host,
        "detected": {"at": system.now(), "by": f"pf-launch {VERSION} --detect", "gpu_inventory": note or "ok"},
        "topology": topology.to_json(),
        "gpus": inventory,
        "nice": 0,
        "default_gpus": "auto",
        "core_preference": [],
        "reserved_cores": {"cores": [], "reason": ""},
        "reserved_gpus": {"gpus": [], "reason": ""},
        "reservations": [],
        "leases": {"core_path": "~/.local/state/pf-launch/leases/core{core}.lock", "core_overrides": {},
                   "gpu_path": "~/.local/state/pf-launch/leases/gpu-{uuid}.lock", "gpu_extra": {}},
        "launcher": {"mpirun": "mpirun", "mpirun_args": [], "srun": "srun",
                     "srun_args": ["--cpu-bind=cores", "--ntasks-per-core=1", "--kill-on-bad-exit=1"],
                     "srun_in_allocation": True, "nvidia_smi": nvidia_smi, "wrapper_python": None},
        "execution": {"kind": "host", "prefix": [], "mpi_inside_prefix": False, "wrapper_python": None,
                      "root": None, "manifest": None, "image": None, "image_digest": None, "hash_image": False},
        "environment": {"set": {}, "unset": []},
        "identity": {"policy": "record-only", "record": record},
        "receipt_dir": None,
    }


# --------------------------------------------------------------------------- CLI


def make_parser():
    parser = argparse.ArgumentParser(
        prog="pf-launch", usage="pf-launch [options] -- PROGRAM ARGS...\n       pf-launch --detect [--site PATH] [--force]",
        description="Launch PROGRAM with one pinned single-threaded rank per selected physical core.")
    parser.add_argument("--site", help="site file (default: $PF_LAUNCH_SITE or sites/<hostname>.json next to pf_launch.py)")
    parser.add_argument("--detect", action="store_true", help="write a site file from /sys, NUMA and nvidia-smi")
    parser.add_argument("--force", action="store_true", help="--detect: overwrite an existing site file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--cores", help="N physical cores, or auto (default: every free physical core)")
    group.add_argument("--cpus", help="exact logical CPUs, one per rank, e.g. 6-7 or 0,2,5 (no SMT siblings)")
    parser.add_argument("--cores-available", action="store_true", help="with --cores N: accept fewer free cores")
    parser.add_argument("--gpus", default=None,
                        help="none | auto | indices 0,1 | GPU-UUID,... | count:N (default: the site's default_gpus, "
                             "else auto); a bare integer is an index, not a count")
    parser.add_argument("--ranks-per-gpu", type=int, default=1, help="ranks sharing each GPU (default 1)")
    parser.add_argument("--wait", action="store_true", help="queue for leased cores/GPUs instead of refusing")
    parser.add_argument("--launcher", choices=("auto", "exec", "mpirun", "srun"), default="auto",
                        help="auto: srun inside Slurm, exec for one core, else mpirun")
    parser.add_argument("--execution", choices=("site", "host"), default="site",
                        help="site: use the site's frozen launcher or container prefix; host: run directly")
    parser.add_argument("--prefix", help="command prefix per rank (e.g. 'apptainer exec --nv img.sif')")
    parser.add_argument("--dry-run", action="store_true", help="print the resolved plan; refuse if not launchable")
    parser.add_argument("--json", action="store_true", help="--dry-run: print the plan receipt as JSON")
    parser.add_argument("--receipt", help="execution envelope receipt path")
    parser.add_argument("--version", action="version", version=f"pf-launch {VERSION}")
    return parser


def main(argv=None, system=None, out=None, err=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    out = out or sys.stdout
    err = err or sys.stderr
    if argv[:1] == ["--rank-exec"]:
        return rank_exec_main(argv[1:])
    if "--" in argv:
        split = argv.index("--")
        options, program = argv[:split], argv[split + 1:]
    else:
        options, program = argv, []
    parser = make_parser()
    opts = parser.parse_args(options)
    opts.argv = ["pf-launch"] + argv
    system = system or System()

    if opts.detect:
        path = opts.site or default_site_path(system)
        if os.path.exists(path) and not opts.force:
            print(f"pf-launch: {path} exists; edit it or pass --force to overwrite", file=err)
            return EXIT_CONFIG
        site = detect_site(system)
        write_json_atomic(path, site)
        topo = site["topology"]
        print(f"pf-launch: wrote {path}: {topo['physical_core_count']} physical cores "
              f"(SMT {topo['smt_threads_per_core']}), NUMA {topo['numa_nodes']}, "
              f"{len(site['gpus'])} GPU(s); review reservations, leases, nice, execution and identity", file=out)
        return 0

    try:
        site, site_info = load_site(opts.site or default_site_path(system))
    except PfError as error:
        print(f"pf-launch: {error.message}", file=err)
        return error.code

    if opts.dry_run:
        plan = build_plan(opts, site, site_info, system)
        code = plan_exit_code(plan)
        receipt = build_receipt(plan, status=dry_run_status(plan), pf_code=code,
                                child_argv=build_child_argv(plan, program), start_time=system.now(), end_time=None,
                                system=system)
        if opts.json:
            print(json.dumps(receipt, indent=2), file=out)
        else:
            print(format_plan(plan, program, dry_run=True), file=out)
        if opts.receipt:
            write_json_atomic(opts.receipt, receipt)
        return code

    if not program:
        parser.error("PROGRAM required after -- (or use --dry-run / --detect)")
    if opts.receipt and os.path.exists(opts.receipt):
        print(f"pf-launch: receipt {opts.receipt} exists; choose a new path", file=err)
        return EXIT_USAGE
    try:
        return launch(opts, program, site, site_info, system, out, err)
    except KeyboardInterrupt:
        print("pf-launch: interrupted before launch; no leases held", file=err)
        return 130


if __name__ == "__main__":
    sys.exit(main())
