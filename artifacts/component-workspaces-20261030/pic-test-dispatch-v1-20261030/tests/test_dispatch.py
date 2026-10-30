"""Injected primitive contract tests; no real locks, child processes or host pins."""
import ast
from contextlib import contextmanager
import copy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "source"))
import test_dispatch_v1 as d

checks = []


def check(ok, name):
    if not ok:
        raise RuntimeError(name)
    checks.append(name)


class MemoryFS:
    def __init__(self):
        self.files, self.links, self.devices = {}, {}, {}
        self.directories = {"/owned/work", "/owned/source", "/sealed/runtime"}

    def read(self, path, cap):
        value = self.files[path]
        d.need(len(value) <= cap, "resource", "bounded injected read")
        return value

    def info(self, path):
        if path in self.links:
            return {"kind": "symlink", "target": self.links[path]}
        if path in self.files:
            return {"kind": "file", "bytes": len(self.files[path])}
        if path in self.devices:
            return {"kind": "char_device", "major": self.devices[path][0], "minor": self.devices[path][1]}
        if path in self.directories:
            return {"kind": "directory"}
        raise FileNotFoundError(path)

    def digest(self, path, size, host):
        value = self.files[path]
        d.need(len(value) == size, "identity", "injected stream length")
        return d.sha(value)

    def members(self, root, cap):
        rows = {p for p in set(self.files) | set(self.links) | set(self.devices) if p.startswith(root + "/")}
        d.need(len(rows) <= cap, "resource", "injected tree cap")
        return rows


class Leases:
    def __init__(self, fail=None):
        self.trace, self.held, self.fail = [], [], fail

    @contextmanager
    def acquire(self, path):
        self.trace.append(["attempt", path])
        if path == self.fail:
            raise d.Refusal("lock_busy", path)
        self.held.append(path)
        self.trace.append(["acquire", path])
        try:
            yield
        finally:
            if self.held[-1] != path:
                raise RuntimeError("non-LIFO lease release")
            self.held.pop()
            self.trace.append(["release", path])


class Runner:
    def __init__(self, callback=None, result=None):
        self.calls, self.callback = [], callback
        self.result = result or {"returncode": 0, "timed_out": False}

    def run(self, command, cwd, env, limits):
        self.calls.append({"command": list(command), "cwd": cwd, "env": dict(env), "limits": dict(limits)})
        if self.callback:
            self.callback()
        return self.result


def fixture(gpu=False):
    fs = MemoryFS()
    identities = []
    for path, role, value in [("/owned/source/input.cpp", "source", b"source"),
                              ("/owned/source/build.py", "build", b"recipe"),
                              ("/sealed/runtime/lib.so", "runtime", b"library"),
                              ("/sealed/tool", "tool", b"executable")]:
        fs.files[path] = value
        identities.append({"path": path, "role": role, "kind": "file", "bytes": len(value), "sha256": d.sha(value)})
    labels = {"kernel_release", "kernel_notes"} | (d.GPU_HOST if gpu else set())
    host = []
    for label in sorted(labels):
        if label == "gpu_device":
            fs.devices["/dev/nvidia0"] = (195, 0)
            host.append({"path": "/dev/nvidia0", "label": label, "kind": "char_device", "major": 195, "minor": 0})
            continue
        if label == "gpu_information":
            path = "/proc/driver/nvidia/gpus/0000:01:00.0/information"
            value = b"GPU UUID: GPU-12345678-1234-5678-1234-567812345678\nBus Location: 0000:01:00.0\nDevice Minor: 0\n"
        else:
            path = d.HOST_PATHS[label]
            value = (label + ":fixed-manufactured-pin\n").encode()
        fs.files[path] = value
        host.append({"path": path, "label": label, "kind": "file", "bytes": len(value), "sha256": d.sha(value)})
    manifest = {"schema": d.SCHEMA, "mode": "gpu" if gpu else "cpu", "stage": "primitive_test",
                "command": ["/sealed/tool", "--primitive"], "cwd": "/owned/work", "env": {},
                "identities": identities, "host": host, "limits": {"wall_s": 10, "output_bytes": 4096},
                "closure": "reviewed_complete_source_build_runtime_and_transitive_inputs_v1",
                "sealed_roots": ["/owned/source"]}
    context = {"optimize": 0, "affinity": [7], "nice": 15,
               "env": {"SECRET": "must-not-leak", "OMPI_CXX": "/foreign/compiler", "LD_PRELOAD": "/foreign.so",
                       "PYTHONPATH": "/foreign/python", "OPENBLAS_NUM_THREADS": "99"}}
    return fs, manifest, context


def invoke(fs, m, c, execute=True, leases=None, runner=None, raw=None, digest=None):
    data = raw if raw is not None else json.dumps(m, sort_keys=True).encode()
    fs.files["/manifest.json"] = data
    locks, child = leases or Leases(), runner or Runner()
    result = d.Dispatcher(fs, c, locks, child).dispatch("/manifest.json", digest or d.sha(data), execute)
    check(not locks.held, "all injected leases released")
    return result, locks, child


def refusal(mutate, name, gpu=False, code=None):
    fs, m, c = fixture(gpu)
    mutate(fs, m, c)
    r, _, child = invoke(fs, m, c)
    check(r["status"] == "refused" and not child.calls, name)
    if code:
        check(r["error_code"] == code, name + " reason")


def main(output):
    fs, m, c = fixture()
    r, locks, child = invoke(fs, m, c)
    check(r["status"] == "component_command_completed_pins_unchanged", "CPU injected dispatch succeeds")
    check(r["pins_before"] and r["pins_after"], "both pin barriers pass")
    check(r["lock_order"] == [d.CPU_LOCK], "CPU-only lease excludes GPU lease")
    env = child.calls[0]["env"]
    check(all(env.get(k) == "1" for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                         "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "MAX_JOBS"]), "one-thread/job environment")
    check(env["CMAKE_BUILD_PARALLEL_LEVEL"] == "1" and env["MAKEFLAGS"] == "-j1", "build parallelism defaults")
    check(not ({"SECRET", "OMPI_CXX", "LD_PRELOAD", "PYTHONPATH"} & set(env)), "ambient overrides/secrets do not reach child")
    check(env["PYTHONOPTIMIZE"] == "0" and env["CUDA_VISIBLE_DEVICES"] == "-1", "unoptimized CPU child and CUDA hidden")
    check("must-not-leak" not in json.dumps(r), "ambient secret absent from receipt")
    fs, m, c = fixture(True)
    r, locks, child = invoke(fs, m, c)
    check(r["status"] == "component_command_completed_pins_unchanged", "explicit injected GPU mode succeeds without GPU work")
    check(r["lock_order"] == [d.GPU_LOCK, d.CPU_LOCK], "GPU acquired before CPU")
    check([p for action, p in locks.trace if action == "release"] == [d.CPU_LOCK, d.GPU_LOCK], "reverse-order lease release")
    check(child.calls[0]["env"]["CUDA_VISIBLE_DEVICES"] == "GPU-12345678-1234-5678-1234-567812345678", "GPU UUID selects admitted device")
    fs, m, c = fixture()
    r, _, child = invoke(fs, m, c, execute=False)
    check(r["status"] == "verified_only_no_child" and not child.calls and not r["child_started"], "default verification launches nothing")
    for optimize in [1, 2]:
        refusal(lambda f, m, c, v=optimize: c.update(optimize=v), "optimized interpreter refuses " + str(optimize), code="optimization")
    for optimize in ["1", "2", "garbage"]:
        refusal(lambda f, m, c, v=optimize: c["env"].update(PYTHONOPTIMIZE=v), "inherited optimization refuses " + optimize, code="optimization")
    for affinity in [[6], [6, 7], []]:
        refusal(lambda f, m, c, v=affinity: c.update(affinity=v), "wrong affinity refuses " + str(affinity), code="placement")
    for nice in [0, 14, 16]:
        refusal(lambda f, m, c, v=nice: c.update(nice=v), "wrong nice refuses " + str(nice), code="placement")
    refusal(lambda f, m, c: m.update(extra="unknown"), "unknown manifest field refuses", code="schema")
    refusal(lambda f, m, c: m.update(mode="automatic"), "mode must be explicit", code="schema")
    refusal(lambda f, m, c: m.update(stage="full_simulation"), "full simulation scope refuses", code="schema")
    refusal(lambda f, m, c: m.update(closure="assumed"), "undeclared closure refuses", code="closure")
    refusal(lambda f, m, c: m["limits"].update(wall_s=True), "boolean timeout refuses", code="integer")
    refusal(lambda f, m, c: m["limits"].update(wall_s=1.5), "fractional timeout refuses", code="integer")
    refusal(lambda f, m, c: m["limits"].update(wall_s=601), "unbounded timeout refuses", code="integer")
    refusal(lambda f, m, c: m["limits"].update(output_bytes=64 * 1024**2 + 1), "oversize output refuses", code="integer")
    refusal(lambda f, m, c: m["command"].__setitem__(0, "/unlisted/tool"), "unlisted executable refuses", code="closure")
    for option in ["-O", "-OO", "PYTHONOPTIMIZE=1"]:
        refusal(lambda f, m, c, v=option: m["command"].append(v), "command optimization refuses " + option, code="optimization")
    refusal(lambda f, m, c: m["env"].update(PYTHONPATH="/foreign"), "explicit forbidden environment refuses", code="environment")
    refusal(lambda f, m, c: m["env"].update(LD_LIBRARY_PATH="/unlisted"), "unlisted library directory refuses", code="closure")
    refusal(lambda f, m, c: m["env"].update(LD_LIBRARY_PATH="/sealed/runtime"), "unsealed library directory refuses", code="closure")
    refusal(lambda f, m, c: m["env"].update(LD_LIBRARY_PATH="/sealed/runtime:"), "empty library search component refuses")
    refusal(lambda f, m, c: m["host"].pop(), "missing kernel host identity refuses", code="host_admission")
    refusal(lambda f, m, c: m["host"].pop(), "missing GPU host identity refuses", gpu=True, code="host_admission")
    refusal(lambda f, m, c: f.files.update({d.HOST_PATHS["kernel_release"]: b"changed kernel"}), "kernel drift refuses", code="identity")
    refusal(lambda f, m, c: f.files.update({d.HOST_PATHS["nvidia_version"]: b"changed driver"}), "driver drift refuses", gpu=True, code="identity")
    refusal(lambda f, m, c: f.files.update({d.HOST_PATHS["nvidia_build_id"]: b"changed module"}), "module build drift refuses", gpu=True, code="identity")
    refusal(lambda f, m, c: f.devices.update({"/dev/nvidia0": (195, 1)}), "GPU device drift refuses", gpu=True, code="identity")
    for role in ["source", "build", "runtime", "tool"]:
        def change(f, m, c, role=role):
            path = next(r["path"] for r in m["identities"] if r["role"] == role)
            f.files[path] = b"changed"
        refusal(change, role + " drift refuses", code="identity")
    refusal(lambda f, m, c: f.files.update({"/owned/source/extra.hpp": b"injected"}), "extra source-tree file refuses", code="closure")
    refusal(lambda f, m, c: m.update(sealed_roots=[]), "missing source membership barrier refuses", code="closure")
    refusal(lambda f, m, c: m["identities"].append(copy.deepcopy(m["identities"][0])), "duplicate identity refuses", code="duplicate_identity")
    fs, m, c = fixture()
    m["env"]["LD_LIBRARY_PATH"] = "/sealed/runtime"
    m["sealed_roots"].append("/sealed/runtime")
    r, _, child = invoke(fs, m, c)
    check(r["status"] != "refused" and child.calls[0]["env"]["LD_LIBRARY_PATH"] == "/sealed/runtime",
          "explicit sealed library search path accepted")
    fs.files["/sealed/runtime/unlisted-plugin.so"] = b"unlisted"
    r, _, child = invoke(fs, m, c)
    check(r["error_code"] == "closure" and not child.calls, "added runtime search member refuses")
    for field, replacement in [("Bus Location", "0000:02:00.0"), ("Device Minor", "1"), ("GPU UUID", "not-a-uuid")]:
        fs, m, c = fixture(True)
        row = next(r for r in m["host"] if r["label"] == "gpu_information")
        lines = fs.files[row["path"]].decode().splitlines()
        value = ("\n".join(field + ": " + replacement if line.startswith(field + ":") else line for line in lines) + "\n").encode()
        fs.files[row["path"]] = value
        row["bytes"], row["sha256"] = len(value), d.sha(value)
        r, _, child = invoke(fs, m, c)
        check(r["error_code"] == "host_admission" and not child.calls, "internally inconsistent GPU identity refuses " + field)
    fs, m, c = fixture()
    raw = json.dumps(m).replace('"schema":', '"schema": "duplicate", "schema":', 1).encode()
    r, _, child = invoke(fs, m, c, raw=raw)
    check(r["error_code"] == "duplicate_key" and not child.calls, "duplicate JSON key refuses")
    fs, m, c = fixture()
    r, _, child = invoke(fs, m, c, digest="0" * 64)
    check(r["error_code"] == "manifest_identity" and not child.calls, "manifest byte drift refuses")
    fs, m, c = fixture()
    r, _, child = invoke(fs, m, c, raw=b" " * (d.MAX_MANIFEST + 1))
    check(r["error_code"] == "resource" and not child.calls, "manifest byte cap before parsing")
    fs, m, c = fixture(True)
    r, leases, child = invoke(fs, m, c, leases=Leases(d.CPU_LOCK))
    check(r["error_code"] == "lock_busy" and not child.calls, "busy CPU lease refuses after GPU lease")
    check([p for action, p in leases.trace if action == "release"] == [d.GPU_LOCK], "GPU lease released on CPU contention")
    for gpu, path in [(False, "/owned/source/input.cpp"), (True, d.HOST_PATHS["nvidia_build_id"])]:
        fs, m, c = fixture(gpu)
        child = Runner(lambda p=path: fs.files.update({p: b"post-run drift"}))
        r, _, child = invoke(fs, m, c, runner=child)
        check(r["status"] == "refused" and r["failed_phase"] == "post_admission" and len(child.calls) == 1,
              "post-command drift refuses rather than refreshing pin " + path)
    for result in [{"returncode": 7, "timed_out": False}, {"returncode": -9, "timed_out": True}]:
        fs, m, c = fixture()
        r, _, _ = invoke(fs, m, c, runner=Runner(result=result))
        check(r["error_code"] == "child_failed" and r["pins_after"], "child failure still checks final pins")
    fs, m, c = fixture()
    def interrupt():
        raise KeyboardInterrupt()
    r, _, _ = invoke(fs, m, c, runner=Runner(callback=interrupt))
    check(r["error_code"] == "interrupted" and r["pins_after"] and r["status"] == "refused",
          "injected interruption retains final pin checks and releases leases")
    fs, m, c = fixture()
    m["identities"].append({"role": "tool", "path": "/sealed/link", "kind": "symlink", "target": "tool"})
    fs.links["/sealed/link"] = "tool"
    m["command"][0] = "/sealed/link"
    r, _, _ = invoke(fs, m, c)
    check(r["status"] != "refused", "complete declared symlink chain works")
    fs.links["/sealed/link"] = "foreign"
    r, _, child = invoke(fs, m, c)
    check(r["error_code"] == "identity" and not child.calls, "symlink target drift refuses")
    m["identities"][-1]["target"] = "foreign"
    r, _, child = invoke(fs, m, c)
    check(r["error_code"] == "closure" and not child.calls, "unlisted symlink target refuses")
    source = Path(d.__file__).read_text()
    check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "dispatcher contains no removable assert statements")
    check("write_text" not in source and "write_bytes" not in source, "no pin refresh write path")
    result = {"schema": "test_dispatch_v1_primitive_contract_results", "status": "passed",
              "checks": len(checks), "check_names": checks, "real_child_processes": 0,
              "real_GPU_locks_or_queries": 0, "numerical_operators_executed": 0,
              "dependency_injected_adapters": ["filesystem", "process_environment", "leases", "runner"],
              "real_process_adapter_qualified": False, "full_integration_qualified": False}
    Path(output).write_text(json.dumps(result, indent=2) + "\n")
    print("Passed", len(checks), "injected primitive checks; zero real child dispatches.")


if __name__ == "__main__":
    main(sys.argv[1])
