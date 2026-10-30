#!/usr/bin/env python3
"""Fixed-pin component dispatcher. Default CLI operation verifies only."""
from __future__ import annotations

import argparse
from contextlib import ExitStack, contextmanager
import fcntl
import hashlib
import json
import os
import posixpath
from pathlib import Path, PurePosixPath
import re
import signal
import stat
import subprocess
import sys
import time

SCHEMA = "pic_generic_test_dispatch_v1"
CPU_LOCK = "<workspace>/outputs/pic-profile/testing-core7.lock"
GPU_LOCK = "<workspace>/outputs/pic-profile/testing-gpu.lock"
MAX_MANIFEST = 1024 * 1024
MAX_FILES = 10000
MAX_BYTES = 2 * 1024**3
HEX = re.compile(r"[0-9a-f]{64}\Z")
BASE_ENV = {
    "HOME": "/home/treefor", "PATH": "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
    "LANG": "C", "LC_ALL": "C", "TZ": "UTC", "PYTHONNOUSERSITE": "1",
    "PYTHONHASHSEED": "0", "PYTHONOPTIMIZE": "0", "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "BLIS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "MAKEFLAGS": "-j1",
    "CMAKE_BUILD_PARALLEL_LEVEL": "1", "MAX_JOBS": "1",
}
EXTRA_ENV = {"RAIKOU_PIC_FIELD_BACKEND", "LD_LIBRARY_PATH", "TMPDIR"}
HOST_PATHS = {"kernel_release": "/proc/sys/kernel/osrelease", "kernel_notes": "/sys/kernel/notes"}
for _module in ("nvidia", "nvidia_uvm", "nvidia_modeset", "nvidia_drm"):
    HOST_PATHS[_module + "_version"] = "/sys/module/" + _module + "/version"
    HOST_PATHS[_module + "_build_id"] = "/sys/module/" + _module + "/notes/.note.gnu.build-id"
HOST_PATHS["nvidia_srcversion"] = "/sys/module/nvidia/srcversion"
GPU_HOST = set(HOST_PATHS) - {"kernel_release", "kernel_notes"}
GPU_HOST |= {"gpu_information", "gpu_device"}


class Refusal(Exception):
    def __init__(self, code, detail):
        self.code = code
        super().__init__(detail)


def need(ok, code, detail):
    if not ok:
        raise Refusal(code, detail)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def exact_keys(value, keys, where):
    need(type(value) is dict and set(value) == set(keys), "schema", where)


def integer(value, lo, hi, where):
    need(type(value) is int and lo <= value <= hi, "integer", where)


def text(value, limit, where):
    need(type(value) is str and 0 < len(value) <= limit and "\0" not in value, "text", where)


def absolute(value):
    text(value, 4096, "path")
    need(value.startswith("/") and ".." not in PurePosixPath(value).parts and
         str(PurePosixPath(value)) == value, "path", value)


def strict_json(raw):
    def pairs(rows):
        out = {}
        for key, value in rows:
            need(key not in out, "duplicate_key", key)
            out[key] = value
        return out
    try:
        return json.loads(raw, object_pairs_hook=pairs,
                          parse_constant=lambda x: (_ for _ in ()).throw(Refusal("number", x)))
    except (ValueError, UnicodeError, RecursionError) as exc:
        raise Refusal("json", str(exc)) from exc


def record_schema(row, host=False):
    need(type(row) is dict, "schema", "identity object")
    kind = row.get("kind")
    prefix = {"label"} if host else {"role"}
    shapes = {"file": {"path", "kind", "bytes", "sha256"},
              "symlink": {"path", "kind", "target"},
              "char_device": {"path", "kind", "major", "minor"}}
    need(kind in shapes, "schema", "identity kind")
    exact_keys(row, prefix | shapes[kind], "identity fields")
    absolute(row["path"])
    if host:
        text(row["label"], 64, "host label")
        need(kind in {"file", "char_device"}, "schema", "host identity kind")
    else:
        need(row["role"] in {"source", "build", "runtime", "input", "tool"}, "schema", "identity role")
        need(kind != "char_device", "schema", "device only belongs to host admission")
        need(not row["path"].startswith(("/proc/", "/sys/", "/dev/")), "path", "host path in file closure")
    if kind == "file":
        integer(row["bytes"], 0, 65536 if host else 1024**3, "file bytes")
        need(type(row["sha256"]) is str and HEX.fullmatch(row["sha256"]), "digest", "file SHA256")
    elif kind == "symlink":
        text(row["target"], 4096, "symlink target")
    else:
        integer(row["major"], 0, 2**32 - 1, "device major")
        integer(row["minor"], 0, 2**32 - 1, "device minor")


def validate_manifest(m):
    exact_keys(m, {"schema", "mode", "stage", "command", "cwd", "env", "identities", "host",
                   "limits", "closure", "sealed_roots"}, "manifest fields")
    need(m["schema"] == SCHEMA, "schema", "schema version")
    need(m["mode"] in {"cpu", "gpu"}, "schema", "explicit cpu/gpu mode")
    need(m["stage"] in {"build", "primitive_test", "component_test"}, "schema", "bounded component scope")
    need(m["closure"] == "reviewed_complete_source_build_runtime_and_transitive_inputs_v1",
         "closure", "explicit reviewed closure declaration required")
    absolute(m["cwd"])
    need(type(m["command"]) is list and 1 <= len(m["command"]) <= 256, "schema", "command argv")
    for value in m["command"]:
        text(value, 4096, "argument")
        need(value not in {"-O", "-OO"} and not value.startswith("PYTHONOPTIMIZE="),
             "optimization", "optimization override in command refused")
    absolute(m["command"][0])
    need(sum(map(len, m["command"])) <= 65536, "resource", "argv bytes")
    need(type(m["identities"]) is list and 4 <= len(m["identities"]) <= MAX_FILES,
         "resource", "identity count")
    paths, roles, total = set(), set(), 0
    for row in m["identities"]:
        record_schema(row)
        need(row["path"] not in paths, "duplicate_identity", row["path"])
        paths.add(row["path"])
        roles.add(row["role"])
        total += row.get("bytes", 0)
    need(total <= MAX_BYTES, "resource", "identity byte budget")
    need({"source", "build", "runtime", "tool"} <= roles, "closure", "missing closure role")
    need(any(r["path"] == m["command"][0] and r["role"] == "tool" for r in m["identities"]),
         "closure", "command executable is not pinned as tool")
    index = {r["path"]: r for r in m["identities"]}
    for row in m["identities"]:
        seen, current = set(), row
        while current["kind"] == "symlink":
            need(current["path"] not in seen, "closure", "symlink cycle")
            seen.add(current["path"])
            target = current["target"]
            resolved = posixpath.normpath(target if target.startswith("/") else
                                         posixpath.join(posixpath.dirname(current["path"]), target))
            need(resolved in index, "closure", "unlisted symlink target: " + resolved)
            current = index[resolved]
    need(type(m["sealed_roots"]) is list and 1 <= len(m["sealed_roots"]) <= 32,
         "closure", "bounded sealed source roots required")
    need(len(set(m["sealed_roots"])) == len(m["sealed_roots"]), "closure", "duplicate sealed roots")
    for root in m["sealed_roots"]:
        absolute(root)
        need(any(r["path"].startswith(root + "/") for r in m["identities"]),
             "closure", "sealed root has no pins")
    need(all(any(r["path"].startswith(root + "/") for root in m["sealed_roots"])
             for r in m["identities"] if r["role"] == "source"), "closure", "source outside sealed roots")
    need(type(m["host"]) is list and len(m["host"]) <= 16, "resource", "host identity count")
    labels = set()
    for row in m["host"]:
        record_schema(row, True)
        label = row["label"]
        need(label not in labels, "duplicate_identity", label)
        labels.add(label)
        if label in HOST_PATHS:
            need(row["path"] == HOST_PATHS[label] and row["kind"] == "file", "host_path", label)
        elif label == "gpu_information":
            need(re.fullmatch(r"/proc/driver/nvidia/gpus/[0-9a-f]{4}:[0-9a-f]{2}:[0-9a-f]{2}\.[0-7]/information", row["path"])
                 and row["kind"] == "file", "host_path", label)
        elif label == "gpu_device":
            need(re.fullmatch(r"/dev/nvidia[0-9]+", row["path"]) and row["kind"] == "char_device", "host_path", label)
        else:
            raise Refusal("host_path", label)
    required = {"kernel_release", "kernel_notes"} | (GPU_HOST if m["mode"] == "gpu" else set())
    need(labels == required, "host_admission", "missing or unexpected host identities")
    exact_keys(m["limits"], {"wall_s", "output_bytes"}, "limits")
    integer(m["limits"]["wall_s"], 1, 600, "wall seconds")
    integer(m["limits"]["output_bytes"], 1, 64 * 1024**2, "output bytes")
    need(type(m["env"]) is dict and set(m["env"]) <= EXTRA_ENV, "environment", "extra environment keys")
    for key, value in m["env"].items():
        text(value, 16384, "environment value")
        if key == "LD_LIBRARY_PATH":
            for path in value.split(":"):
                absolute(path)
                need(path in m["sealed_roots"], "closure", "library search directory is not sealed")
                need(any(r["role"] == "runtime" and r["path"].startswith(path + "/") for r in m["identities"]),
                     "closure", "library directory has no declared runtime pin")
        if key == "TMPDIR":
            absolute(value)
    return m


class RealFS:
    def members(self, root, cap):
        pending, out, visited = [root], set(), 0
        while pending:
            directory = pending.pop()
            need(self.info(directory)["kind"] == "directory", "closure", "sealed root is not a directory")
            with os.scandir(directory) as entries:
                for entry in entries:
                    visited += 1
                    need(visited <= cap, "resource", "sealed tree entry cap")
                    if entry.is_dir(follow_symlinks=False):
                        pending.append(entry.path)
                    else:
                        out.add(entry.path)
        return out

    def read(self, path, cap):
        fd = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        with os.fdopen(fd, "rb") as stream:
            data = stream.read(cap + 1)
        need(len(data) <= cap, "resource", "bounded read: " + path)
        return data

    def info(self, path):
        s = os.lstat(path)
        if stat.S_ISREG(s.st_mode):
            return {"kind": "file", "bytes": s.st_size}
        if stat.S_ISLNK(s.st_mode):
            return {"kind": "symlink", "target": os.readlink(path)}
        if stat.S_ISCHR(s.st_mode):
            return {"kind": "char_device", "major": os.major(s.st_rdev), "minor": os.minor(s.st_rdev)}
        if stat.S_ISDIR(s.st_mode):
            return {"kind": "directory"}
        return {"kind": "unsupported"}

    def digest(self, path, size, host):
        fd = os.open(path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        with os.fdopen(fd, "rb") as stream:
            before = os.fstat(stream.fileno())
            need(stat.S_ISREG(before.st_mode), "identity", "not regular: " + path)
            h, count = hashlib.sha256(), 0
            while True:
                block = stream.read(min(1024**2, size + 1 - count))
                if not block:
                    break
                count += len(block)
                need(count <= size, "identity", "larger than expected: " + path)
                h.update(block)
            after = os.fstat(stream.fileno())
        need(count == size, "identity", "byte count: " + path)
        need((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
             (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), "identity", "changed during hash: " + path)
        return h.hexdigest()


class RealLeases:
    @contextmanager
    def acquire(self, path):
        fd = os.open(path, os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW)
        try:
            need(stat.S_ISREG(os.fstat(fd).st_mode), "lock", "not regular lock file")
            deadline = time.monotonic() + 10
            while True:
                try:
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    need(time.monotonic() < deadline, "lock_busy", path)
                    time.sleep(.05)
            yield
        finally:
            os.close(fd)


class ProcessRunner:
    """Session adapter; arbitrary daemonizing commands require a separate owned-process guard."""
    def __init__(self, log_path):
        self.log_path = log_path

    def run(self, argv, cwd, env, limits):
        import resource
        def child_limits():
            resource.setrlimit(resource.RLIMIT_FSIZE, (limits["output_bytes"], limits["output_bytes"]))
            resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        fd = os.open(self.log_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, "wb") as log:
            proc = subprocess.Popen(argv, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True, close_fds=True,
                                    preexec_fn=child_limits)
            try:
                code = proc.wait(timeout=limits["wall_s"])
                return {"returncode": code, "timed_out": False}
            except subprocess.TimeoutExpired:
                # No poll/wait between timeout and these signals: the unreaped
                # session leader retains this PID/PGID and prevents reuse.
                for sig in (signal.SIGTERM, signal.SIGKILL):
                    try:
                        os.killpg(proc.pid, sig)
                    except ProcessLookupError:
                        pass
                    if sig == signal.SIGTERM:
                        time.sleep(.25)
                proc.wait(timeout=2)
                return {"returncode": proc.returncode, "timed_out": True}
            except BaseException:
                if proc.returncode is None:
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    proc.wait(timeout=2)
                raise


def real_context():
    return {"optimize": sys.flags.optimize, "affinity": sorted(os.sched_getaffinity(0)),
            "nice": os.getpriority(os.PRIO_PROCESS, 0), "env": dict(os.environ)}


class Dispatcher:
    def __init__(self, fs, context, leases, runner):
        self.fs, self.context, self.leases, self.runner = fs, context, leases, runner

    def verify(self, rows, host):
        for row in rows:
            actual = self.fs.info(row["path"])
            need(actual["kind"] == row["kind"], "identity", "kind: " + row["path"])
            if row["kind"] == "file":
                # proc/sys file st_size is commonly zero; the bounded stream
                # length/hash is authoritative for those explicit host pins.
                if not host:
                    need(actual.get("bytes") == row["bytes"], "identity", "size: " + row["path"])
                need(self.fs.digest(row["path"], row["bytes"], host) == row["sha256"],
                     "identity", "SHA256: " + row["path"])
            elif row["kind"] == "symlink":
                need(actual.get("target") == row["target"], "identity", "symlink: " + row["path"])
            else:
                need(actual.get("major") == row["major"] and actual.get("minor") == row["minor"],
                     "identity", "device: " + row["path"])

    def verify_roots(self, m):
        for root in m["sealed_roots"]:
            expected = {r["path"] for r in m["identities"] if r["path"].startswith(root + "/")}
            need(self.fs.members(root, MAX_FILES) == expected, "closure", "changed sealed tree members: " + root)

    def gpu_environment(self, host):
        index = {r["label"]: r for r in host}
        row = index["gpu_information"]
        value = self.fs.read(row["path"], 65536).decode("ascii")
        need(sha(value.encode("ascii")) == row["sha256"], "identity", "GPU information changed")
        fields = {}
        for line in value.splitlines():
            if ":" in line:
                key, item = line.split(":", 1)
                need(key.strip() not in fields, "host_admission", "duplicate GPU information field")
                fields[key.strip()] = item.strip()
        uuid = fields.get("GPU UUID", "")
        need(re.fullmatch(r"GPU-[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", uuid),
             "host_admission", "GPU UUID format")
        need(fields.get("Bus Location") == row["path"].split("/")[-2], "host_admission", "GPU PCI mapping")
        device = index["gpu_device"]
        need(fields.get("Device Minor") == str(device["minor"]) == device["path"].removeprefix("/dev/nvidia"),
             "host_admission", "GPU minor mapping")
        return {"CUDA_VISIBLE_DEVICES": uuid, "CUDA_DEVICE_ORDER": "PCI_BUS_ID"}

    def dispatch(self, manifest_path, expected_sha256, execute=False):
        receipt = {"schema": SCHEMA + "_receipt", "status": "refused", "child_started": False,
                   "expected_manifest_sha256": expected_sha256, "lock_order": [], "pins_before": False,
                   "pins_after": False, "integration_qualified": False, "automatic_deployment": False}
        phase = "bootstrap"
        try:
            c = self.context
            need(type(c["optimize"]) is int and c["optimize"] == 0 and c["env"].get("PYTHONOPTIMIZE", "") in {"", "0"},
                 "optimization", "optimized dispatcher execution refused")
            need(c["affinity"] == [7] and c["nice"] == 15, "placement", "CPU7/nice15 required")
            need(type(expected_sha256) is str and HEX.fullmatch(expected_sha256), "digest", "trusted manifest digest required")
            raw = self.fs.read(manifest_path, MAX_MANIFEST)
            need(sha(raw) == expected_sha256, "manifest_identity", "manifest differs from reviewed bytes")
            m = validate_manifest(strict_json(raw))
            receipt.update(mode=m["mode"], stage=m["stage"], command=m["command"], cpu=[7], nice=15,
                           identity_records=len(m["identities"]), host_records=len(m["host"]))
            need(self.fs.info(m["cwd"])["kind"] == "directory", "cwd", "cwd is not a real directory")
            env = dict(BASE_ENV)
            env.update(m["env"])
            receipt["child_environment"] = env
            phase = "locks"
            order = ([GPU_LOCK] if m["mode"] == "gpu" else []) + [CPU_LOCK]
            with ExitStack() as stack:
                for path in order:
                    stack.enter_context(self.leases.acquire(path))
                    receipt["lock_order"].append(path)
                phase = "pre_admission"
                self.verify(m["host"], True)
                self.verify(m["identities"], False)
                self.verify_roots(m)
                receipt["pins_before"] = True
                env.update(self.gpu_environment(m["host"]) if m["mode"] == "gpu" else {"CUDA_VISIBLE_DEVICES": "-1"})
                if not execute:
                    receipt["status"] = "verified_only_no_child"
                else:
                    phase = "child"
                    receipt["child_dispatch_attempted"] = True
                    receipt["child_started"] = None
                    try:
                        result = self.runner.run(m["command"], m["cwd"], env, m["limits"])
                        receipt["child_started"] = True
                        exact_keys(result, {"returncode", "timed_out"}, "runner result")
                        integer(result["returncode"], -255, 255, "child exit")
                        need(type(result["timed_out"]) is bool, "schema", "runner timeout flag")
                        receipt["child"] = result
                    finally:
                        phase = "post_admission"
                        self.verify(m["host"], True)
                        self.verify(m["identities"], False)
                        self.verify_roots(m)
                        receipt["pins_after"] = True
                    need(not result["timed_out"] and result["returncode"] == 0, "child_failed", "nonzero or timed-out child")
                    receipt["status"] = "component_command_completed_pins_unchanged"
            return receipt
        except BaseException as exc:
            receipt.update(status="refused", failed_phase=phase,
                           error_code=exc.code if isinstance(exc, Refusal) else
                           "interrupted" if isinstance(exc, (KeyboardInterrupt, SystemExit, InterruptedError)) else "adapter_error",
                           error=str(exc))
            return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--execute", action="store_true", help="explicit reviewed component command only")
    args = parser.parse_args()
    absolute(args.receipt)
    absolute(args.manifest)
    def interrupted(signum, frame):
        raise InterruptedError("Dispatcher signal " + str(signum))
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGHUP, interrupted)
    # Reserve the fresh receipt before a child is even considered.
    fd = os.open(args.receipt, os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w") as stream:
        receipt = Dispatcher(RealFS(), real_context(), RealLeases(),
                             ProcessRunner(args.receipt + ".child.log")).dispatch(args.manifest, args.expected_sha256, args.execute)
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    return 0 if receipt["status"] != "refused" else 78


if __name__ == "__main__":
    raise SystemExit(main())
