"""Job containment for MPI launchers (0.7.2).

A job's processes are contained by its unit's cgroup: MemoryMax, the CPU
set, the OOM handling and the kill at its end all apply to that cgroup only.
Open MPI 5 launchers can put a job's ranks somewhere else: ``prun`` and
``prterun --dvm search|system|pid:N`` look for a running PRRTE DVM through the
PMIx rendezvous files under ``$TMPDIR`` (``/tmp`` when unset) and
``$PMIX_SYSTEM_TMPDIR``, connect to it as a tool and ask it to spawn the
ranks -- the DVM forks them in its own cgroup (reproduced on zenobia: two
``prun`` ranks from a job ran in the outside DVM's scope).  A job's own
``mpiexec`` (prterun) is such a server too: its rendezvous files in ``/tmp``
let an outside tool attach to it.

0.7.2 therefore gives every job its own rendezvous and temporary directory
``<tmp_base>/simq.<job>.<token12>`` (mode 700, made by the job wrapper,
removed when the job ends) and points ``TMPDIR``, ``PMIX_SYSTEM_TMPDIR``,
``PMIX_SERVER_TMPDIR`` and ``PRTE_MCA_prte_tmpdir_base`` at it: a search finds
only DVMs of the same job, an outside one makes ``prun`` fail ("no DVM being
available") instead of escaping, and ``mpiexec`` works as before.  Launcher
identity variables captured at submission (``PMIX_*`` but ``PMIX_MCA_*``,
``PRTE_*`` but ``PRTE_MCA_*``, ``PMI_*``, ``OMPI_COMM_WORLD_*`` and the like)
are dropped, as ``SLURM_*`` are: they would make a launcher act as a client
of a server outside the job.  A variable given with ``-e`` always wins.

Only an explicit address (``--dvm file:PATH``/``uri:``) can still reach an
outside DVM.  The worker reports what it can see (``EscapeScanner``):

* processes carrying a running job's ``SIMQ_LAUNCH_TOKEN`` outside the job's
  cgroup (moved out with ``systemd-run --scope`` or a cgroup write, or left
  behind after the job ended), per job in ``escaped.json`` of the job
  directory, in ``simq show`` and in the worker log;
* MPI launcher daemons (``prte``, ``prted``, ``prterun``, ``orted``,
  ``hydra_pmi_proxy``, ...) of this user running outside every job's cgroup,
  node-wide in ``containment.json`` of the state directory, in ``simq
  partitions`` and ``simq doctor``.  Ranks such a daemon spawns for a job do
  not carry the job's environment (PRRTE forwards none of it), so they are
  reported with their daemon, not per job.

The scan reads ``/proc`` only (no process is spawned): a process's
environment is read once (it cannot change), its cgroup on every scan while
it carries a token.
"""

from __future__ import annotations

import json
import os
import shutil
import stat
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

TOKEN_KEY = b"SIMQ_LAUNCH_TOKEN="
# the rendezvous and temporary directory variables set for every job
RENDEZVOUS_VARIABLES = ("TMPDIR", "PMIX_SYSTEM_TMPDIR", "PMIX_SERVER_TMPDIR", "PRTE_MCA_prte_tmpdir_base")
# processes of MPI launchers that may spawn ranks for others
MPI_DAEMONS = frozenset({
    "prte", "prted", "prterun", "prun", "orted", "orterun", "ompi-server",
    "hydra_pmi_proxy", "mpiexec.hydra", "hydra_nameserver", "pmix_server", "psrvr",
})
_IDENTITY_PREFIXES = ("PMIX_", "PRTE_", "PMI_", "OMPI_COMM_WORLD_", "OMPI_MCA_ess_base_", "OMPI_MCA_orte_")
_IDENTITY_KEEP_PREFIXES = ("PMIX_MCA_", "PRTE_MCA_")
_IDENTITY_NAMES = frozenset({
    "OMPI_UNIVERSE_SIZE", "OMPI_APP_CTX_NUM_PROCS", "OMPI_FILE_LOCATION", "OMPI_FIRST_RANKS",
    "OMPI_NUM_APP_CTX", "OMPI_ARGV", "OMPI_COMMAND", "OMPI_MCA_initial_wdir", "OMPI_MCA_num_procs",
    "OMPI_MCA_orte_hnp_uri", "OMPI_MCA_pmix", "OMPI_MCA_ess",
})


@dataclass(frozen=True)
class ContainmentConfig:
    """``containment`` in node.json (all optional)."""

    enabled: bool = True
    # parent of the per-job directories (same filesystem jobs used for /tmp)
    tmp_base: str = "/tmp"
    # escape scan while job units run; 0: only when a job ends
    scan_seconds: float = 60.0


def parse_containment_config(value: Any) -> ContainmentConfig:
    if value is None:
        return ContainmentConfig()
    if not isinstance(value, dict):
        raise ValueError("containment must be an object")
    unknown = set(value) - {"enabled", "tmp_base", "scan_seconds", "notes"}
    if unknown:
        raise ValueError(f"containment: unknown key(s) {', '.join(sorted(unknown))}")
    enabled = value.get("enabled", True)
    if not isinstance(enabled, bool):
        raise ValueError("containment.enabled must be true or false")
    base = value.get("tmp_base", "/tmp")
    if not isinstance(base, str) or not base.startswith("/") or len(base) > 60:
        raise ValueError("containment.tmp_base must be an absolute path of at most 60 characters")
    try:
        scan = float(value.get("scan_seconds", 60.0))
    except (TypeError, ValueError) as error:
        raise ValueError("containment.scan_seconds must be a number") from error
    if scan != 0 and not 5 <= scan <= 3600:
        raise ValueError("containment.scan_seconds must be 0 or 5 to 3600")
    return ContainmentConfig(enabled=enabled, tmp_base=base.rstrip("/") or "/", scan_seconds=scan)


def is_launcher_identity(name: str) -> bool:
    """A variable an MPI launcher sets for the processes it starts: inherited
    by a job, it makes the job's launcher a client of that outside server."""
    if name.startswith(_IDENTITY_KEEP_PREFIXES):
        return False
    return name in _IDENTITY_NAMES or name.startswith(_IDENTITY_PREFIXES)


def job_tmpdir(config: ContainmentConfig, job_id: int, token: str) -> str:
    """Short on purpose: PMIx and MPI build socket and file names below it."""
    clean = "".join(ch for ch in token if ch.isalnum())[:12]
    return f"{config.tmp_base.rstrip('/')}/simq.{int(job_id)}.{clean}"


def rendezvous_environment(path: str) -> dict[str, str]:
    environment = {name: path for name in RENDEZVOUS_VARIABLES}
    environment["SIMQ_TMPDIR"] = path
    return environment


def make_job_tmpdir(path: str) -> None:
    """Create the job's directory (mode 700).  An existing one is reused only
    if it is a real directory of this user with mode 700 (a restarted
    wrapper); anything else under that name is refused."""
    try:
        os.mkdir(path, 0o700)
        os.chmod(path, 0o700)
        return
    except FileExistsError:
        pass
    info = os.lstat(path)
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) != 0o700:
        raise PermissionError(f"{path} exists and is not this user's private directory")


def remove_job_tmpdir(path: str | None) -> int:
    """Remove a job's directory (never following symlinks); returns the
    bytes freed (best effort).  Only a directory of this user named like a
    job directory is touched."""
    if not path or not Path(path).name.startswith("simq."):
        return 0
    try:
        info = os.lstat(path)
    except OSError:
        return 0
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid():
        return 0
    size = 0
    for root, dirs, files in os.walk(path, followlinks=False):
        for name in files:
            try:
                size += os.lstat(os.path.join(root, name)).st_blocks * 512
            except OSError:
                pass
    shutil.rmtree(path, ignore_errors=True)
    return size


# ----- escape scan ----------------------------------------------------------------


def _read(path: str, limit: int = 1 << 20) -> bytes | None:
    try:
        with open(path, "rb") as handle:
            return handle.read(limit)
    except OSError:
        return None


def _cgroup_of(pid: int) -> str | None:
    data = _read(f"/proc/{pid}/cgroup", 4096)
    if data is None:
        return None
    for line in data.decode("utf-8", "replace").splitlines():
        if line.startswith("0::"):
            return line[3:]
    return None


def _cmdline(pid: int) -> str:
    data = _read(f"/proc/{pid}/cmdline", 4096) or b""
    return data.replace(b"\0", b" ").decode("utf-8", "replace").strip()[:200]


def _comm(pid: int) -> str:
    data = _read(f"/proc/{pid}/comm", 64) or b""
    return data.decode("utf-8", "replace").strip()


def _ppid(pid: int) -> int | None:
    data = _read(f"/proc/{pid}/stat", 1024)
    if not data:
        return None
    try:
        return int(data.rsplit(b")", 1)[1].split()[1])
    except (IndexError, ValueError):
        return None


def _within(cgroup: str | None, root: str) -> bool:
    return cgroup is not None and (cgroup == root or cgroup.startswith(root.rstrip("/") + "/"))


@dataclass
class JobScope:
    job_id: int
    token: str
    cgroup: str  # as in /proc/<pid>/cgroup, e.g. /user.slice/.../simq-job-1-abc.service
    name: str = ""


@dataclass
class Escaped:
    pid: int
    job_id: int
    comm: str
    cmdline: str
    cgroup: str
    first_seen: float
    last_seen: float
    after_end: bool = False

    def to_json(self) -> dict[str, Any]:
        return {"pid": self.pid, "comm": self.comm, "cmdline": self.cmdline, "cgroup": self.cgroup,
                "first_seen": _iso(self.first_seen), "last_seen": _iso(self.last_seen),
                "after_end": self.after_end}


def _iso(unix: float) -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(unix))


@dataclass
class ScanResult:
    escaped: dict[int, list[Escaped]] = field(default_factory=dict)  # job -> processes
    daemons: list[dict[str, Any]] = field(default_factory=list)
    new: list[Escaped] = field(default_factory=list)  # first seen in this scan
    seconds: float = 0.0
    processes: int = 0


class EscapeScanner:
    """Finds processes of running (or just ended) jobs outside their cgroup,
    and MPI launcher daemons outside every job's cgroup."""

    def __init__(self, uid: int | None = None) -> None:
        self.uid = os.getuid() if uid is None else uid
        # pid -> (token or None, comm); a process's environment is fixed
        self._seen: dict[int, tuple[str | None, str]] = {}
        self._first: dict[tuple[int, int], float] = {}
        self.scans = 0
        self.cpu_seconds = 0.0

    def _classify(self, pid: int) -> tuple[str | None, str] | None:
        known = self._seen.get(pid)
        if known is not None:
            return known
        try:
            if os.stat(f"/proc/{pid}").st_uid != self.uid:
                self._seen[pid] = (None, "")
                return self._seen[pid]
        except OSError:
            return None
        data = _read(f"/proc/{pid}/environ")
        token = None
        if data:
            at = data.find(TOKEN_KEY)
            if at == 0 or (at > 0 and data[at - 1:at] == b"\0"):
                end = data.find(b"\0", at)
                token = data[at + len(TOKEN_KEY): end if end >= 0 else len(data)].decode("ascii", "replace")
        comm = _comm(pid)
        self._seen[pid] = (token, comm)
        return self._seen[pid]

    def scan(self, jobs: Iterable[JobScope], ended: Iterable[JobScope] = ()) -> ScanResult:
        """``jobs``: running job units; ``ended``: jobs whose unit just ended
        (any process of theirs still alive is reported as left behind)."""
        start = time.process_time()
        result = ScanResult()
        by_token = {scope.token: scope for scope in jobs}
        ended_by_token = {scope.token: scope for scope in ended}
        roots = [scope.cgroup for scope in by_token.values()]
        try:
            pids = [int(entry.name) for entry in os.scandir("/proc") if entry.name.isdigit()]
        except OSError:
            return result
        alive = set(pids)
        for pid in list(self._seen):
            if pid not in alive:
                self._seen.pop(pid, None)
        now = time.time()
        own = os.getpid()
        for pid in pids:
            if pid == own:
                continue
            known = self._classify(pid)
            if known is None:
                continue
            token, comm = known
            if token is None and comm not in MPI_DAEMONS:
                continue
            cgroup = _cgroup_of(pid)
            if cgroup is None:
                continue
            if token is not None:
                scope = by_token.get(token)
                after_end = False
                if scope is None:
                    scope = ended_by_token.get(token)
                    after_end = scope is not None
                if scope is not None and (after_end or not _within(cgroup, scope.cgroup)):
                    key = (scope.job_id, pid)
                    first = self._first.setdefault(key, now)
                    item = Escaped(pid=pid, job_id=scope.job_id, comm=comm, cmdline=_cmdline(pid), cgroup=cgroup,
                                   first_seen=first, last_seen=now, after_end=after_end)
                    result.escaped.setdefault(scope.job_id, []).append(item)
                    if first == now:
                        result.new.append(item)
                    continue
            if comm in MPI_DAEMONS and not any(_within(cgroup, root) for root in roots):
                result.daemons.append({"pid": pid, "comm": comm, "cmdline": _cmdline(pid), "cgroup": cgroup,
                                       "children": self._children(pid, pids)})
        result.processes = len(pids)
        result.seconds = time.process_time() - start
        self.scans += 1
        self.cpu_seconds += result.seconds
        return result

    def _children(self, parent: int, pids: list[int]) -> list[dict[str, Any]]:
        children = []
        for pid in pids:
            if _ppid(pid) == parent:
                children.append({"pid": pid, "comm": _comm(pid), "cgroup": _cgroup_of(pid) or ""})
                if len(children) >= 16:
                    break
        return children

    def forget_job(self, job_id: int) -> None:
        for key in [key for key in self._first if key[0] == job_id]:
            self._first.pop(key, None)


def write_escaped(job_dir: str | Path, items: list[Escaped]) -> None:
    """``escaped.json`` of the job directory: every process seen, merged."""
    path = Path(job_dir) / "escaped.json"
    previous: dict[int, dict[str, Any]] = {}
    try:
        for entry in json.loads(path.read_text(encoding="utf-8")).get("processes", []):
            previous[int(entry["pid"])] = entry
    except (OSError, ValueError, KeyError, TypeError):
        pass
    for item in items:
        entry = item.to_json()
        if item.pid in previous:
            entry["first_seen"] = previous[item.pid].get("first_seen", entry["first_seen"])
            entry["after_end"] = bool(previous[item.pid].get("after_end")) or item.after_end
        previous[item.pid] = entry
    document = {"schema": "simqueue-escaped-v1", "processes": sorted(previous.values(), key=lambda e: e["pid"])}
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    os.replace(temporary, path)


def read_escaped(job_dir: str | Path | None) -> list[dict[str, Any]]:
    if not job_dir:
        return []
    try:
        document = json.loads((Path(job_dir) / "escaped.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    processes = document.get("processes") if isinstance(document, dict) else None
    return [p for p in processes or [] if isinstance(p, dict)]


def write_daemons(path: str | Path, daemons: list[dict[str, Any]]) -> None:
    document = {"schema": "simqueue-containment-v1", "scanned_at": _iso(time.time()), "outside_mpi_daemons": daemons}
    target = Path(path)
    temporary = target.with_name(target.name + ".tmp")
    temporary.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    os.replace(temporary, target)


def read_daemons(path: str | Path) -> dict[str, Any] | None:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return document if isinstance(document, dict) else None


def describe_escaped(processes: list[dict[str, Any]], limit: int = 3) -> str:
    parts = []
    for entry in processes[:limit]:
        where = str(entry.get("cgroup", "")).rsplit("/", 1)[-1] or "?"
        parts.append(f"pid {entry.get('pid')} {entry.get('comm') or '?'} in {where}"
                     + (" (after the job ended)" if entry.get("after_end") else ""))
    more = len(processes) - limit
    return "; ".join(parts) + (f"; and {more} more" if more > 0 else "")


def environment_with_containment(environment: Mapping[str, str], explicit: set[str],
                                 path: str | None) -> dict[str, str]:
    """``environment`` without inherited launcher identity variables, with the
    job's rendezvous directory (an explicit ``-e`` value wins)."""
    result = {k: v for k, v in environment.items() if k in explicit or not is_launcher_identity(k)}
    if path:
        for key, value in rendezvous_environment(path).items():
            if key not in explicit or key == "SIMQ_TMPDIR":
                result[key] = value
    return result
