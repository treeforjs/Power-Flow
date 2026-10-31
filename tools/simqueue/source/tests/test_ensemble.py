from __future__ import annotations

import ctypes
import json
import os
import shutil
import signal
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing, contextmanager
from pathlib import Path
from typing import Any
from unittest.mock import patch

from simqueue import ensemble as ensemble_module
from simqueue.config import format_cpu_list
from simqueue.ensemble import EnsembleSupervisor
from simqueue.store import SCHEMA_VERSION, QueueStore


REPOSITORY = Path(__file__).resolve().parents[1]


@contextmanager
def preserve_child_subreaper() -> Iterator[None]:
    """Keep direct supervisor tests from changing the test runner's reaper role."""
    libc = ctypes.CDLL(None, use_errno=True)
    previous = ctypes.c_int()
    if libc.prctl(37, ctypes.byref(previous), 0, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), "PR_GET_CHILD_SUBREAPER failed")
    try:
        yield
    finally:
        if libc.prctl(36, previous.value, 0, 0, 0) != 0:
            raise OSError(ctypes.get_errno(), "PR_SET_CHILD_SUBREAPER restore failed")


MEMBER_PROGRAM = r"""
import json
import os
import sys
import time
from pathlib import Path

name = sys.argv[1]
delay = float(sys.argv[2])
exit_code = int(sys.argv[3])
uvm_marker = sys.argv[4]
started = time.time()
Path("member.pid").write_text(str(os.getpid()), encoding="utf-8")
record = {
    "name": name,
    "pid": os.getpid(),
    "started": started,
    "affinity": sorted(os.sched_getaffinity(0)),
    "environment": {
        key: os.environ.get(key)
        for key in (
            "SIMQ_ENSEMBLE_JOB_ID",
            "SIMQ_ENSEMBLE_MEMBER_ID",
            "SIMQ_ENSEMBLE_MEMBER_NAME",
            "SIMQ_LANE",
            "SIMQ_LANE_CPU",
            "SIMQ_CPU_LIST",
            "SIMQ_CPU_COUNT",
            "CUDA_VISIBLE_DEVICES",
            "CUDA_MPS_PIPE_DIRECTORY",
            "OMP_NUM_THREADS",
            "OMP_PROC_BIND",
            "OMP_PLACES",
            "GOMP_CPU_AFFINITY",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "KOKKOS_NUM_THREADS",
            "ENSEMBLE_DEFAULT",
            "MEMBER_ONLY",
        )
    },
}
print(f"member-output name={name}", flush=True)
time.sleep(delay)
if uvm_marker:
    print(uvm_marker, flush=True)
record["finished"] = time.time()
Path("member-result.json").write_text(json.dumps(record), encoding="utf-8")
raise SystemExit(exit_code)
"""


UVM_FINAL_MARKER = (
    "MSIC_KOKKOS_CUDAUVM_BACKING_V1 event=finalize policy=managed-pool "
    "limit_bytes=536870912 live_allocations=0 live_bytes=0 "
    "peak_live_bytes=284464608 live_reserved_bytes=0 "
    "peak_reserved_bytes=284679424 allocation_requests=487 "
    "raw_allocations=1 pool_reclaims=3 teardown_ok=1 release=success"
)


def process_is_live(pid: int) -> bool:
    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
    except OSError:
        return False
    fields = stat[stat.rfind(")") + 2 :].split()
    return bool(fields) and fields[0] != "Z"


class EnsembleFixture:
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.projects = self.base / "projects"
        self.projects.mkdir()
        allowed = sorted(os.sched_getaffinity(0))
        self.cpus = tuple(allowed[: min(2, len(allowed))])
        self.environment = {
            "SIMQUEUE_STATE_DIR": str(self.base / "state"),
            "SIMQUEUE_PROJECT_ROOT": str(self.projects),
            "SIMQUEUE_PCORES": format_cpu_list(self.cpus),
            "SIMQUEUE_NO_AUTOSTART": "1",
            "SIMQUEUE_CANCEL_GRACE_SECONDS": "0.5",
        }
        self.environment_patch = patch.dict(
            os.environ, self.environment, clear=False
        )
        self.environment_patch.start()
        self.launcher = self.base / "simq-test-launcher"
        self.launcher.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "from pathlib import Path\n"
            f"sys.path.insert(0, {str(REPOSITORY)!r})\n"
            "from simqueue import exclusive\n"
            f"exclusive.CONFIG_PATH = Path({str(self.base / 'no-exclusive-config')!r})\n"
            "from simqueue.cli import main\n"
            "raise SystemExit(main())\n",
            encoding="utf-8",
        )
        self.launcher.chmod(0o700)
        self.store = QueueStore()

    def tearDown(self) -> None:
        self.environment_patch.stop()
        self.temporary.cleanup()

    def run_cli(
        self, *arguments: str, timeout: float = 15
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(self.launcher), *arguments],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )

    def run_worker(self, *, timeout: float = 20) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )

    def member(
        self,
        name: str,
        *,
        delay: float = 0.05,
        exit_code: int = 0,
        gpu_memory: str = "32MiB",
        uvm_marker: str = "",
        environment: dict[str, str] | None = None,
        cwd: Path | None = None,
    ) -> dict[str, Any]:
        directory = cwd or self.projects / name
        directory.mkdir(parents=True, exist_ok=True)
        result: dict[str, Any] = {
            "name": name,
            "cwd": str(directory),
            "argv": [
                sys.executable,
                "-c",
                MEMBER_PROGRAM,
                name,
                str(delay),
                str(exit_code),
                uvm_marker,
            ],
            "gpu_memory": gpu_memory,
        }
        if environment is not None:
            result["env"] = environment
        return result

    def write_manifest(
        self,
        members: list[dict[str, Any]],
        *,
        name: str = "test-ensemble",
        defaults: dict[str, Any] | None = None,
        max_parallel: int | None = None,
        allow_shared_cwd: bool | None = None,
        gpu: dict[str, Any] | None = None,
        launch_gate: dict[str, Any] | None = None,
    ) -> Path:
        document: dict[str, Any] = {
            "schema_version": "simqueue-ensemble-v1",
            "name": name,
            "members": members,
        }
        if defaults is not None:
            document["defaults"] = defaults
        if max_parallel is not None:
            document["max_parallel"] = max_parallel
        if allow_shared_cwd is not None:
            document["allow_shared_cwd"] = allow_shared_cwd
        if gpu is not None:
            document["gpu"] = gpu
        if launch_gate is not None:
            document["launch_gate"] = launch_gate
        path = self.base / f"manifest-{time.monotonic_ns()}.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def submit(
        self,
        manifest: Path,
        *extra: str,
        expect_success: bool = True,
    ) -> tuple[subprocess.CompletedProcess[str], int | None]:
        before = {int(job["id"]) for job in self.store.list_jobs(include_all=True)}
        result = self.run_cli("ensemble", "submit", str(manifest), *extra)
        after = self.store.list_jobs(include_all=True)
        new_ids = sorted(int(job["id"]) for job in after if int(job["id"]) not in before)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(new_ids), 1, result.stdout + result.stderr)
            return result, new_ids[0]
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertEqual(new_ids, [], "invalid manifest left a durable queue job")
        return result, None

    def status(self, job_id: int) -> dict[str, Any]:
        result = self.run_cli("ensemble", "status", str(job_id), "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        return dict(json.loads(result.stdout))


class EnsembleValidationTests(EnsembleFixture, unittest.TestCase):
    def test_binary_size_units_are_canonicalized_to_bytes(self) -> None:
        manifest = self.write_manifest(
            [self.member("one", gpu_memory="512MiB")]
        )
        result, job_id = self.submit(
            manifest,
            "--gpu-budget",
            "1GiB",
            "--gpu-device",
            "7",
        )
        assert job_id is not None
        self.assertIn(f"Queued ensemble {job_id}", result.stdout)
        payload = self.status(job_id)
        self.assertEqual(payload["ensemble"]["gpu_budget_bytes"], 1024**3)
        self.assertEqual(payload["ensemble"]["gpu_device"], "7")
        self.assertEqual(payload["members"][0]["gpu_admission_bytes"], 512 * 1024**2)

    def test_invalid_size_is_rejected_before_queue_mutation(self) -> None:
        manifest = self.write_manifest([self.member("one")])
        result, _job_id = self.submit(
            manifest, "--gpu-budget", "twelve-gib", expect_success=False
        )
        self.assertIn("size", result.stderr.lower())

    def test_nonfinite_diagnostic_intervals_fail_cleanly(self) -> None:
        for raw_value in ("nan", "inf", "-inf"):
            with self.subTest(value=raw_value):
                manifest = self.write_manifest([self.member(f"interval-{raw_value}")])
                result, _job_id = self.submit(
                    manifest,
                    "--gpu-budget",
                    "32MiB",
                    f"--diagnostic-interval={raw_value}",
                    expect_success=False,
                )
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertIn("diagnostics_interval", result.stderr.lower())
                self.assertIn("positive number", result.stderr.lower())

    def test_launch_gate_is_canonicalized_and_invalid_forms_are_rejected(self) -> None:
        gate = {
            "kind": "member-log-substring",
            "contains": "event=end phase=stk_build_objects rank=0",
            "timeout_seconds": 12,
        }
        manifest = self.write_manifest(
            [self.member("gated")],
            launch_gate=gate,
            gpu={"device": "0", "budget": "32MiB", "mps": "off"},
        )
        result = self.run_cli("ensemble", "validate", str(manifest), "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        canonical = json.loads(result.stdout)
        self.assertEqual(
            canonical["launch_gate"],
            {**gate, "timeout_seconds": 12.0},
        )

        invalid_gates = (
            {"kind": "unknown", "contains": "ready"},
            {"kind": "member-log-substring", "contains": ""},
            {"kind": "member-log-substring", "contains": "bad\nmarker"},
            {
                "kind": "member-log-substring",
                "contains": "ready",
                "timeout_seconds": float("nan"),
            },
            {
                "kind": "member-log-substring",
                "contains": "ready",
                "unexpected": True,
            },
        )
        for index, invalid in enumerate(invalid_gates):
            with self.subTest(index=index):
                manifest = self.write_manifest(
                    [self.member(f"invalid-gate-{index}")],
                    launch_gate=invalid,
                )
                result, _job_id = self.submit(
                    manifest,
                    "--gpu-budget",
                    "32MiB",
                    expect_success=False,
                )
                self.assertIn("launch_gate", result.stderr)

    def test_sizes_beyond_sqlite_int64_fail_cleanly(self) -> None:
        too_large = str(2**63)
        cases = (
            (
                "member-claim",
                self.member("huge-member-claim", gpu_memory=f"{too_large}B"),
                {"device": "0", "budget": f"{too_large}B", "mps": "off"},
            ),
            (
                "gpu-budget",
                self.member("huge-gpu-budget", gpu_memory="0B"),
                {"device": "0", "budget": f"{too_large}B", "mps": "off"},
            ),
            (
                "mapped-host-uvm",
                self.member(
                    "huge-mapped-host-uvm",
                    gpu_memory="0B",
                    environment={
                        "MSIC_KOKKOS_CUDAUVM_BACKING": "mapped-host",
                        "MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES": too_large,
                    },
                ),
                {"device": "0", "budget": "0B", "mps": "off"},
            ),
        )
        for label, member, gpu in cases:
            with self.subTest(case=label):
                manifest = self.write_manifest([member], gpu=gpu)
                result, _job_id = self.submit(manifest, expect_success=False)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertIn("maximum supported size", result.stderr.lower())
                self.assertIn(str(2**63 - 1), result.stderr)

    def test_schema_and_duplicate_member_names_are_rejected(self) -> None:
        first = self.member("duplicate")
        second = self.member("duplicate", cwd=self.projects / "duplicate-two")
        manifest = self.write_manifest([first, second])
        result, _job_id = self.submit(manifest, expect_success=False)
        self.assertIn("duplicate", result.stderr.lower())

        invalid_schema = self.write_manifest([self.member("schema")])
        document = json.loads(invalid_schema.read_text(encoding="utf-8"))
        document["schema_version"] = "simqueue-ensemble-v999"
        invalid_schema.write_text(json.dumps(document), encoding="utf-8")
        result, _job_id = self.submit(invalid_schema, expect_success=False)
        self.assertIn("schema", result.stderr.lower())

    def test_member_working_directory_cannot_escape_project_root(self) -> None:
        outside = self.base / "outside"
        outside.mkdir()
        manifest = self.write_manifest([self.member("escape", cwd=outside)])
        result, _job_id = self.submit(manifest, expect_success=False)
        self.assertIn("inside", result.stderr.lower())

        symlink = self.projects / "outside-link"
        symlink.symlink_to(outside, target_is_directory=True)
        member = self.member("symlink-member")
        member["cwd"] = str(symlink)
        manifest = self.write_manifest([member])
        result, _job_id = self.submit(manifest, expect_success=False)
        self.assertIn("inside", result.stderr.lower())

    def test_shared_working_directory_requires_explicit_opt_in(self) -> None:
        shared = self.projects / "shared"
        shared.mkdir()
        members = [
            self.member("first", cwd=shared),
            self.member("second", cwd=shared),
        ]
        manifest = self.write_manifest(members)
        result, _job_id = self.submit(manifest, expect_success=False)
        self.assertIn("shared", result.stderr.lower())

        manifest = self.write_manifest(members, allow_shared_cwd=True)
        _result, job_id = self.submit(manifest, "--gpu-budget", "64MiB")
        self.assertIsNotNone(job_id)

    def test_uvm_arena_must_be_positive_and_covered_by_gpu_claim(self) -> None:
        uvm_environment = {
            "MSIC_KOKKOS_CUDAUVM_BACKING": "managed-pool",
            "MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES": str(512 * 1024**2),
        }
        underdeclared = self.member(
            "underdeclared",
            gpu_memory="256MiB",
            environment=uvm_environment,
        )
        manifest = self.write_manifest([underdeclared])
        result, _job_id = self.submit(
            manifest, "--gpu-budget", "1GiB", expect_success=False
        )
        self.assertIn("uvm", result.stderr.lower())

        covered = self.member(
            "covered",
            gpu_memory="512MiB",
            environment=uvm_environment,
        )
        manifest = self.write_manifest([covered])
        _result, job_id = self.submit(manifest, "--gpu-budget", "1GiB")
        assert job_id is not None
        member = self.status(job_id)["members"][0]
        self.assertEqual(member["uvm_arena_bytes"], 512 * 1024**2)

        invalid_environment = dict(uvm_environment)
        invalid_environment["MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES"] = "0"
        invalid = self.member(
            "zero-arena",
            gpu_memory="512MiB",
            environment=invalid_environment,
        )
        manifest = self.write_manifest([invalid])
        result, _job_id = self.submit(manifest, expect_success=False)
        self.assertIn("positive", result.stderr.lower())


class EnsembleDurabilityTests(EnsembleFixture, unittest.TestCase):
    def test_idempotent_ensemble_json_replay_reuses_parent_and_artifacts(self) -> None:
        manifest = self.write_manifest(
            [self.member("idempotent-member")],
            gpu={"device": "-1", "budget": "32MiB", "mps": "off"},
        )
        arguments = (
            "ensemble",
            "submit",
            str(manifest),
            "--idempotency-key",
            "msc080-eq:v1:campaign:production:000",
            "--json",
        )
        first = self.run_cli(*arguments)
        replay = self.run_cli(*arguments)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(replay.returncode, 0, replay.stderr)
        first_payload = json.loads(first.stdout)
        replay_payload = json.loads(replay.stdout)
        self.assertTrue(first_payload["created"])
        self.assertFalse(replay_payload["created"])
        self.assertEqual(replay_payload["job_id"], first_payload["job_id"])
        self.assertEqual(
            replay_payload["artifact_dir"], first_payload["artifact_dir"]
        )
        job_id = int(first_payload["job_id"])
        self.assertEqual(len(self.store.list_jobs(include_all=True)), 1)
        self.assertEqual(len(self.store.list_ensemble_members(job_id)), 1)
        artifacts = Path(first_payload["artifact_dir"])
        self.assertEqual(len(list((artifacts / "members").iterdir())), 1)

    def test_idempotent_ensemble_key_rejects_changed_canonical_manifest(self) -> None:
        manifest = self.write_manifest(
            [self.member("collision-member")],
            gpu={"device": "-1", "budget": "32MiB", "mps": "off"},
        )
        arguments = (
            "ensemble",
            "submit",
            str(manifest),
            "--idempotency-key",
            "ensemble-collision",
            "--json",
        )
        first = self.run_cli(*arguments)
        self.assertEqual(first.returncode, 0, first.stderr)
        document = json.loads(manifest.read_text(encoding="utf-8"))
        document["members"][0]["argv"].append("changed")
        manifest.write_text(json.dumps(document), encoding="utf-8")
        collision = self.run_cli(*arguments)
        self.assertEqual(collision.returncode, 2, collision.stderr)
        self.assertIn("different submission request", collision.stderr)
        self.assertEqual(len(self.store.list_jobs(include_all=True)), 1)

    def test_concurrent_idempotent_ensemble_cli_submits_one_parent(self) -> None:
        manifest = self.write_manifest(
            [self.member("concurrent-member")],
            gpu={"device": "-1", "budget": "32MiB", "mps": "off"},
        )
        arguments = (
            "ensemble",
            "submit",
            str(manifest),
            "--idempotency-key",
            "concurrent-ensemble",
            "--json",
        )
        with ThreadPoolExecutor(max_workers=8) as pool:
            results = list(pool.map(lambda _index: self.run_cli(*arguments), range(16)))
        for result in results:
            self.assertEqual(result.returncode, 0, result.stderr)
        payloads = [json.loads(result.stdout) for result in results]
        self.assertEqual({int(payload["job_id"]) for payload in payloads}, {1})
        self.assertEqual(sum(bool(payload["created"]) for payload in payloads), 1)
        self.assertEqual(len(self.store.list_jobs(include_all=True)), 1)
        self.assertEqual(len(self.store.list_ensemble_members(1)), 1)

    def test_retry_clones_stored_snapshot_with_fresh_logs_and_tokens(self) -> None:
        original_member = self.member(
            "immutable-member",
            environment={"MEMBER_ONLY": "stored-value"},
        )
        manifest = self.write_manifest(
            [original_member],
            defaults={"env": {"ENSEMBLE_DEFAULT": "stored-default"}},
            gpu={"device": "4", "budget": "64MiB", "mps": "off"},
        )
        _result, original_id = self.submit(
            manifest, "--diagnostic-interval", "0.2"
        )
        assert original_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(original_id)["state"], "succeeded")

        original_job = self.store.require(original_id)
        original_ensemble = self.store.require_ensemble(original_id)
        original_rows = self.store.list_ensemble_members(original_id)
        original_row = original_rows[0]
        original_log = Path(str(original_row["log_path"]))
        original_log_contents = original_log.read_bytes()

        # Retry must not reread this now-invalid manifest or adopt its new data.
        manifest.write_text(
            json.dumps(
                {
                    "schema_version": "not-the-submitted-schema",
                    "name": "mutated",
                    "members": [],
                }
            ),
            encoding="utf-8",
        )
        before_ids = {
            int(job["id"]) for job in self.store.list_jobs(include_all=True)
        }
        retry_result = self.run_cli("retry", str(original_id))
        self.assertEqual(retry_result.returncode, 0, retry_result.stderr)
        retry_ids = {
            int(job["id"])
            for job in self.store.list_jobs(include_all=True)
            if int(job["id"]) not in before_ids
        }
        self.assertEqual(len(retry_ids), 1)
        retry_id = retry_ids.pop()

        retry_job = self.store.require(retry_id)
        retry_ensemble = self.store.require_ensemble(retry_id)
        retry_rows = self.store.list_ensemble_members(retry_id)
        retry_row = retry_rows[0]
        self.assertEqual(retry_job["retry_of"], original_id)
        self.assertEqual(retry_job["attempt"], 2)
        self.assertEqual(retry_ensemble["spec_json"], original_ensemble["spec_json"])
        for field in (
            "schema_version",
            "max_parallel",
            "gpu_device",
            "gpu_budget_bytes",
            "diagnostics_interval",
            "mps_mode",
            "mps_pipe_directory",
            "fail_fast",
            "allow_shared_cwd",
        ):
            self.assertEqual(retry_ensemble[field], original_ensemble[field], field)
        for field in (
            "name",
            "cwd",
            "argv_json",
            "env_json",
            "gpu_admission_bytes",
            "uvm_arena_bytes",
        ):
            self.assertEqual(retry_row[field], original_row[field], field)
        self.assertEqual(retry_row["state"], "queued")
        self.assertNotEqual(retry_row["member_token"], original_row["member_token"])
        self.assertNotEqual(retry_row["log_path"], original_row["log_path"])
        self.assertNotEqual(retry_ensemble["artifact_dir"], original_ensemble["artifact_dir"])
        self.assertEqual(Path(str(retry_row["log_path"])).read_bytes(), b"")
        self.assertEqual(original_log.read_bytes(), original_log_contents)

        result_path = self.projects / "immutable-member" / "member-result.json"
        result_path.unlink()
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(retry_id)["state"], "succeeded")
        self.assertTrue(result_path.is_file())
        retry_record = json.loads(result_path.read_text(encoding="utf-8"))
        self.assertEqual(retry_record["environment"]["MEMBER_ONLY"], "stored-value")
        self.assertEqual(
            retry_record["environment"]["ENSEMBLE_DEFAULT"], "stored-default"
        )
        self.assertIn(
            "member-output name=immutable-member",
            Path(str(retry_row["log_path"])).read_text(encoding="utf-8"),
        )
        self.assertEqual(original_log.read_bytes(), original_log_contents)

    def test_queued_cancel_atomically_reconciles_every_member(self) -> None:
        members = [self.member("queued-a"), self.member("queued-b")]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest, "--max-parallel", "1", "--gpu-budget", "32MiB"
        )
        assert job_id is not None

        cancel = self.run_cli("cancel", str(job_id))
        self.assertEqual(cancel.returncode, 0, cancel.stderr)
        parent = self.store.require(job_id)
        rows = self.store.list_ensemble_members(job_id)
        self.assertEqual(parent["state"], "canceled")
        self.assertIsNotNone(parent["finished_at"])
        self.assertEqual([row["state"] for row in rows], ["canceled", "canceled"])
        for row in rows:
            self.assertIsNotNone(row["finished_at"])
            self.assertIsNone(row["pid"])
            self.assertIn("before ensemble launch", str(row["message"]))
            self.assertEqual(Path(str(row["log_path"])).read_bytes(), b"")

        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertFalse((self.projects / "queued-a" / "member.pid").exists())
        self.assertFalse((self.projects / "queued-b" / "member.pid").exists())
        payload = self.status(job_id)
        self.assertEqual(payload["summary"]["canceled"], 2)

    def test_worker_restart_marks_stale_parent_and_all_members_interrupted(self) -> None:
        members = [self.member("stale-running"), self.member("stale-queued")]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest, "--max-parallel", "1", "--gpu-budget", "32MiB"
        )
        assert job_id is not None
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        self.assertEqual(int(claimed["id"]), job_id)
        self.store.mark_member_running(
            job_id,
            0,
            lane=0,
            cpu=self.cpus[0],
            pid=999_999_991,
            pgid=999_999_991,
        )

        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        parent = self.store.require(job_id)
        self.assertEqual(parent["state"], "interrupted")
        self.assertIn("worker restarted", str(parent["message"]))
        rows = self.store.list_ensemble_members(job_id)
        self.assertEqual(
            [row["state"] for row in rows], ["interrupted", "interrupted"]
        )
        self.assertTrue(
            all("parent job" in str(row["message"]) for row in rows), rows
        )
        self.assertFalse((self.projects / "stale-running" / "member.pid").exists())
        self.assertFalse((self.projects / "stale-queued" / "member.pid").exists())


class EnsemblePauseControlTests(EnsembleFixture, unittest.TestCase):
    def test_pause_control_is_durable_idempotent_and_generation_counted(self) -> None:
        manifest = self.write_manifest([self.member("pause-control", gpu_memory="0B")])
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None

        initial = self.store.ensemble_control_snapshot(job_id)
        self.assertFalse(initial["launch_paused"])
        self.assertEqual(initial["control_generation"], 0)
        self.assertIsNone(initial["control_changed_at"])
        self.assertEqual(initial["queued"], 1)

        paused = self.store.set_ensemble_paused(job_id, True)
        self.assertTrue(paused["launch_paused"])
        self.assertEqual(paused["control_generation"], 1)
        self.assertIsNotNone(paused["control_changed_at"])

        duplicate_pause = self.store.set_ensemble_paused(job_id, True)
        self.assertEqual(duplicate_pause, paused)
        self.assertEqual(
            QueueStore().ensemble_control_snapshot(job_id),
            paused,
            "a fresh store connection did not observe the durable pause",
        )

        resumed = self.store.set_ensemble_paused(job_id, False)
        self.assertFalse(resumed["launch_paused"])
        self.assertEqual(resumed["control_generation"], 2)
        self.assertIsNotNone(resumed["control_changed_at"])

        duplicate_resume = self.store.set_ensemble_paused(job_id, False)
        self.assertEqual(duplicate_resume, resumed)

    def test_pause_closes_launch_reservation_and_resume_reopens_it(self) -> None:
        manifest = self.write_manifest(
            [self.member("reservation-control", gpu_memory="0B")]
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None

        self.store.set_ensemble_paused(job_id, True)
        self.assertFalse(
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            ),
            "a pause committed before reservation did not close admission",
        )
        self.assertEqual(
            self.store.require_ensemble_member(job_id, 0)["state"], "queued"
        )

        self.store.set_ensemble_paused(job_id, False)
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.assertEqual(int(claimed["id"]), job_id)
        self.store.register_ensemble_supervisor(
            job_id,
            pid=os.getpid(),
            launch_token=str(claimed["launch_token"]),
        )

        self.store.set_ensemble_paused(job_id, True)
        self.assertFalse(
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            ),
            "a running supervisor admitted a member after pause",
        )
        self.store.set_ensemble_paused(job_id, False)
        self.assertTrue(
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            )
        )
        reserved = self.store.require_ensemble_member(job_id, 0)
        self.assertEqual(reserved["state"], "running")
        self.assertIsNone(reserved["pid"])
        self.assertEqual(reserved["lane"], 0)
        self.assertEqual(reserved["cpu"], self.cpus[0])

        self.store.record_member_process(job_id, 0, pid=987_654_321, pgid=987_654_321)
        tracked = self.store.require_ensemble_member(job_id, 0)
        self.assertEqual(tracked["pid"], 987_654_321)
        self.assertEqual(tracked["pgid"], 987_654_321)

        with self.assertRaisesRegex(RuntimeError, "not queued at reservation"):
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            )

    def test_drained_parent_requeues_same_id_and_preserves_member_progress(self) -> None:
        manifest = self.write_manifest(
            [
                self.member("drained-complete", gpu_memory="0B"),
                self.member("drained-pending", gpu_memory="0B"),
            ],
            max_parallel=1,
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.store.register_ensemble_supervisor(
            job_id,
            pid=os.getpid(),
            launch_token=str(claimed["launch_token"]),
        )
        self.assertTrue(
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            )
        )
        self.store.record_member_process(job_id, 0, pid=123_456_789, pgid=123_456_789)

        paused = self.store.set_ensemble_paused(job_id, True)
        generation = int(paused["control_generation"])
        self.store.finish_member(job_id, 0, state="succeeded", exit_code=0)
        self.assertTrue(self.store.mark_ensemble_drained(job_id, generation))
        self.assertTrue(self.store.requeue_drained_ensemble(job_id))
        self.assertFalse(
            self.store.requeue_drained_ensemble(job_id),
            "a drained parent was requeued more than once",
        )

        suspended_parent = self.store.require(job_id)
        self.assertEqual(suspended_parent["state"], "queued")
        self.assertIsNone(suspended_parent["pid"])
        self.assertIsNone(suspended_parent["pgid"])
        self.assertIsNone(suspended_parent["cpu_list"])
        self.assertIsNone(suspended_parent["finished_at"])
        self.assertIsNone(suspended_parent["exit_code"])
        self.assertIn("paused after draining", str(suspended_parent["message"]))
        suspended_members = self.store.list_ensemble_members(job_id)
        self.assertEqual(
            [member["state"] for member in suspended_members],
            ["succeeded", "queued"],
        )

        self.store.set_ensemble_paused(job_id, False)
        resumed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(resumed)
        assert resumed is not None
        self.assertEqual(int(resumed["id"]), job_id)
        self.assertEqual(len(self.store.list_jobs(include_all=True)), 1)

    def test_resume_then_repause_before_requeue_preserves_drain_marker(self) -> None:
        manifest = self.write_manifest(
            [
                self.member("rapid-control-complete", gpu_memory="0B"),
                self.member("rapid-control-pending", gpu_memory="0B"),
            ],
            max_parallel=1,
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.store.register_ensemble_supervisor(
            job_id,
            pid=os.getpid(),
            launch_token=str(claimed["launch_token"]),
        )
        self.assertTrue(
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            )
        )
        self.store.record_member_process(job_id, 0, pid=111_222_333, pgid=111_222_333)
        first_pause = self.store.set_ensemble_paused(job_id, True)
        self.assertEqual(first_pause["control_generation"], 1)
        self.store.finish_member(job_id, 0, state="succeeded", exit_code=0)
        self.assertTrue(self.store.mark_ensemble_drained(job_id, 1))

        resumed_before_requeue = self.store.set_ensemble_paused(job_id, False)
        self.assertEqual(resumed_before_requeue["control_generation"], 2)
        self.assertEqual(resumed_before_requeue["drained_generation"], 1)
        repaused_before_requeue = self.store.set_ensemble_paused(job_id, True)
        self.assertEqual(repaused_before_requeue["control_generation"], 3)
        self.assertEqual(
            repaused_before_requeue["drained_generation"],
            1,
            "rapid re-pause discarded the supervisor's committed drain marker",
        )
        self.assertTrue(self.store.requeue_drained_ensemble(job_id))
        self.assertEqual(self.store.require(job_id)["state"], "queued")

        resumed = self.store.set_ensemble_paused(job_id, False)
        self.assertEqual(resumed["control_generation"], 4)
        same_parent = self.store.claim_next(self.cpus)
        self.assertIsNotNone(same_parent)
        assert same_parent is not None
        self.assertEqual(int(same_parent["id"]), job_id)
        claimed_control = self.store.ensemble_control_snapshot(job_id)
        self.assertIsNone(claimed_control["drained_generation"])
        self.assertIsNone(claimed_control["drained_at"])
        self.assertEqual(claimed_control["supervisor_protocol"], 1)

    def test_worker_restart_recovers_a_committed_drain_without_losing_members(
        self,
    ) -> None:
        manifest = self.write_manifest(
            [
                self.member("drain-crash-complete", gpu_memory="0B"),
                self.member("drain-crash-pending", gpu_memory="0B"),
            ],
            max_parallel=1,
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.store.register_ensemble_supervisor(
            job_id,
            pid=os.getpid(),
            launch_token=str(claimed["launch_token"]),
        )
        self.assertTrue(
            self.store.reserve_member_launch(
                job_id, 0, lane=0, cpu=self.cpus[0]
            )
        )
        self.store.finish_member(job_id, 0, state="succeeded", exit_code=0)
        paused = self.store.set_ensemble_paused(job_id, True)
        self.assertTrue(
            self.store.mark_ensemble_drained(
                job_id, int(paused["control_generation"])
            )
        )

        recovered = self.run_worker()

        self.assertEqual(recovered.returncode, 0, recovered.stdout + recovered.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "queued")
        self.assertEqual(
            [row["state"] for row in self.store.list_ensemble_members(job_id)],
            ["succeeded", "queued"],
        )
        self.assertIn(
            "recovered drained ensemble",
            (self.base / "state" / "worker.log").read_text(encoding="utf-8"),
        )

    def test_new_supervisor_holds_safely_when_outer_worker_lacks_drain_protocol(
        self,
    ) -> None:
        manifest = self.write_manifest(
            [
                self.member("legacy-worker-active", delay=0.45, gpu_memory="0B"),
                self.member("legacy-worker-pending", delay=0.05, gpu_memory="0B"),
            ],
            max_parallel=1,
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        environment = dict(os.environ)
        environment.update(
            {
                "SIMQ_JOB_ID": str(job_id),
                "SIMQ_LAUNCH_TOKEN": str(claimed["launch_token"]),
                "SIMQ_CPU_LIST": str(self.cpus[0]),
            }
        )
        environment.pop("SIMQ_ENSEMBLE_DRAIN_PROTOCOL", None)
        supervisor = subprocess.Popen(
            [sys.executable, str(self.launcher), "_run-ensemble", str(job_id)],
            cwd=REPOSITORY,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            first_pid = self.projects / "legacy-worker-active" / "member.pid"
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline and not first_pid.exists():
                time.sleep(0.02)
            self.assertTrue(first_pid.exists(), "first compatibility member did not start")

            self.store.set_ensemble_paused(job_id, True)
            first_result = (
                self.projects / "legacy-worker-active" / "member-result.json"
            )
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline and not first_result.exists():
                time.sleep(0.02)
            self.assertTrue(first_result.exists(), "admitted compatibility member did not finish")
            time.sleep(0.25)
            self.assertIsNone(
                supervisor.poll(),
                "supervisor returned a drain exit to a legacy outer worker",
            )
            self.assertFalse(
                (self.projects / "legacy-worker-pending" / "member.pid").exists()
            )
            held = self.store.ensemble_control_snapshot(job_id)
            self.assertEqual((held["succeeded"], held["queued"]), (1, 1))
            self.assertIsNone(held["drained_generation"])

            self.store.set_ensemble_paused(job_id, False)
            stdout, stderr = supervisor.communicate(timeout=8)
            self.assertEqual(supervisor.returncode, 0, stdout + stderr)
            self.assertEqual(
                [row["state"] for row in self.store.list_ensemble_members(job_id)],
                ["succeeded", "succeeded"],
            )
            self.assertIn("pause-held-for-legacy-worker", stdout)
        finally:
            if supervisor.poll() is None:
                supervisor.kill()
                supervisor.wait()

    def test_prelaunch_pause_holds_fifo_head(self) -> None:
        first_manifest = self.write_manifest(
            [self.member("paused-fifo-head", gpu_memory="0B")],
            name="paused-fifo-head",
        )
        _result, first_id = self.submit(first_manifest, "--gpu-budget", "0B")
        assert first_id is not None
        second_manifest = self.write_manifest(
            [self.member("younger-fifo-member", gpu_memory="0B")],
            name="younger-fifo-ensemble",
        )
        _result, second_id = self.submit(second_manifest, "--gpu-budget", "0B")
        assert second_id is not None

        self.store.set_ensemble_paused(first_id, True)
        self.assertIsNone(
            self.store.claim_next(self.cpus),
            "dispatch skipped a paused FIFO head and launched younger work",
        )
        self.assertEqual(self.store.require(first_id)["state"], "queued")
        self.assertEqual(self.store.require(second_id)["state"], "queued")
        self.assertEqual(
            self.store.dispatch_blocker(),
            {
                "kind": "ensemble-paused",
                "job_id": first_id,
                "name": "paused-fifo-head",
            },
        )

        promoted_younger = self.store.promote(second_id)
        self.assertTrue(promoted_younger["changed"])
        self.assertEqual(promoted_younger["queue_position"], 1)
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [second_id, first_id],
        )
        self.assertIsNone(
            self.store.dispatch_blocker(),
            "a paused ensemble behind the dispatch head blocked promoted work",
        )

        self.store.demote(second_id)
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [first_id, second_id],
        )
        self.assertEqual(self.store.dispatch_blocker()["job_id"], first_id)

        demoted_paused = self.store.demote(first_id)
        self.assertTrue(demoted_paused["changed"])
        self.assertEqual(demoted_paused["queue_position"], 2)
        self.assertEqual(
            [job["id"] for job in self.store.queued_jobs()],
            [second_id, first_id],
        )
        self.assertIsNone(self.store.dispatch_blocker())

        self.store.promote(first_id)
        self.assertEqual(self.store.dispatch_blocker()["job_id"], first_id)

        self.store.set_ensemble_paused(first_id, False)
        claimed = self.store.claim_next(self.cpus)
        self.assertIsNotNone(claimed)
        assert claimed is not None
        self.assertEqual(int(claimed["id"]), first_id)
        self.assertEqual(
            self.store.ensemble_control_snapshot(first_id)["supervisor_protocol"],
            1,
            "claim did not advertise pause support before supervisor registration",
        )
        self.assertEqual(self.store.require(second_id)["state"], "queued")

    def test_pause_resume_cli_json_and_status_are_idempotent(self) -> None:
        manifest = self.write_manifest(
            [self.member("pause-cli-member", gpu_memory="0B")]
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None

        pause = self.run_cli("ensemble", "pause", str(job_id), "--json")
        self.assertEqual(pause.returncode, 0, pause.stderr)
        pause_payload = json.loads(pause.stdout)
        self.assertEqual(pause_payload["schema_version"], "simqueue-ensemble-control-v1")
        self.assertEqual(pause_payload["action"], "pause")
        self.assertTrue(pause_payload["changed"])
        self.assertTrue(pause_payload["launch_paused"])
        self.assertEqual(pause_payload["control_generation"], 1)

        duplicate_pause = self.run_cli(
            "ensemble", "pause", str(job_id), "--json"
        )
        self.assertEqual(duplicate_pause.returncode, 0, duplicate_pause.stderr)
        duplicate_payload = json.loads(duplicate_pause.stdout)
        self.assertFalse(duplicate_payload["changed"])
        self.assertEqual(duplicate_payload["control_generation"], 1)

        status = self.status(job_id)
        self.assertTrue(status["ensemble"]["launch_paused"])
        self.assertEqual(status["ensemble"]["control_state"], "paused")
        self.assertTrue(status["ensemble"]["allocation_released"])
        human_status = self.run_cli("ensemble", "status", str(job_id))
        self.assertEqual(human_status.returncode, 0, human_status.stderr)
        self.assertIn(
            "Launch control: paused; generation=1; allocation=released",
            human_status.stdout,
        )

        resume = self.run_cli("ensemble", "resume", str(job_id), "--json")
        self.assertEqual(resume.returncode, 0, resume.stderr)
        resume_payload = json.loads(resume.stdout)
        self.assertEqual(resume_payload["action"], "resume")
        self.assertTrue(resume_payload["changed"])
        self.assertFalse(resume_payload["launch_paused"])
        self.assertEqual(resume_payload["control_generation"], 2)
        self.assertEqual(resume_payload["worker_mode"], "disabled")

        duplicate_resume = self.run_cli(
            "ensemble", "resume", str(job_id), "--json"
        )
        self.assertEqual(duplicate_resume.returncode, 0, duplicate_resume.stderr)
        duplicate_resume_payload = json.loads(duplicate_resume.stdout)
        self.assertFalse(duplicate_resume_payload["changed"])
        self.assertEqual(duplicate_resume_payload["control_generation"], 2)

    def test_end_to_end_pause_drains_without_backfill_then_resumes_same_id(
        self,
    ) -> None:
        if len(self.cpus) < 2:
            self.skipTest("pause drain test requires two available CPUs")
        members = [
            self.member("pause-active-a", delay=0.8, gpu_memory="0B"),
            self.member("pause-active-b", delay=0.8, gpu_memory="0B"),
            self.member("pause-pending", delay=0.05, gpu_memory="0B"),
        ]
        manifest = self.write_manifest(members, max_parallel=2)
        _result, job_id = self.submit(manifest, "--gpu-budget", "0B")
        assert job_id is not None
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            active_pid_paths = [
                self.projects / "pause-active-a" / "member.pid",
                self.projects / "pause-active-b" / "member.pid",
            ]
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline:
                if all(path.is_file() for path in active_pid_paths):
                    break
                time.sleep(0.02)
            else:
                self.fail("both initial ensemble lanes did not start")

            paused = self.store.set_ensemble_paused(job_id, True)
            self.assertTrue(paused["launch_paused"])
            stdout, stderr = worker.communicate(timeout=10)
            self.assertEqual(worker.returncode, 0, stdout + stderr)

            parent = self.store.require(job_id)
            drained_rows = self.store.list_ensemble_members(job_id)
            self.assertEqual(parent["state"], "queued")
            self.assertIsNone(parent["pid"])
            self.assertIsNone(parent["pgid"])
            self.assertEqual(
                [row["state"] for row in drained_rows],
                ["succeeded", "succeeded", "queued"],
            )
            self.assertFalse(
                (self.projects / "pause-pending" / "member.pid").exists(),
                "a freed lane backfilled after the pause was committed",
            )
            completed_snapshots = [
                {
                    key: row[key]
                    for key in ("started_at", "finished_at", "pid", "pgid", "exit_code")
                }
                for row in drained_rows[:2]
            ]
            completed_results = [
                (self.projects / name / "member-result.json").read_bytes()
                for name in ("pause-active-a", "pause-active-b")
            ]

            self.store.set_ensemble_paused(job_id, False)
            resumed_worker = self.run_worker(timeout=15)
            self.assertEqual(
                resumed_worker.returncode,
                0,
                resumed_worker.stdout + resumed_worker.stderr,
            )
            self.assertEqual(self.store.require(job_id)["state"], "succeeded")
            final_rows = self.store.list_ensemble_members(job_id)
            self.assertEqual(
                [row["state"] for row in final_rows],
                ["succeeded", "succeeded", "succeeded"],
            )
            for index in range(2):
                self.assertEqual(
                    {
                        key: final_rows[index][key]
                        for key in (
                            "started_at",
                            "finished_at",
                            "pid",
                            "pgid",
                            "exit_code",
                        )
                    },
                    completed_snapshots[index],
                    f"completed member {index} was changed during resume",
                )
            self.assertEqual(
                [
                    (self.projects / name / "member-result.json").read_bytes()
                    for name in ("pause-active-a", "pause-active-b")
                ],
                completed_results,
            )
            events_path = Path(
                str(self.store.require_ensemble(job_id)["artifact_dir"])
            ) / "events.jsonl"
            events = [
                json.loads(line)
                for line in events_path.read_text(encoding="utf-8").splitlines()
                if line
            ]
            launch_counts = {
                index: sum(
                    event.get("event") == "member-launch"
                    and event.get("member_index") == index
                    for event in events
                )
                for index in range(3)
            }
            self.assertEqual(launch_counts, {0: 1, 1: 1, 2: 1})
            self.assertTrue(any(event.get("event") == "pause-drained" for event in events))
        finally:
            if worker.poll() is None:
                try:
                    self.store.cancel(job_id)
                except (KeyError, RuntimeError):
                    pass
                try:
                    worker.communicate(timeout=3)
                except subprocess.TimeoutExpired:
                    worker.kill()
                    worker.wait()


class EnsembleExecutionTests(EnsembleFixture, unittest.TestCase):
    def _read_result(self, name: str) -> dict[str, Any]:
        return dict(
            json.loads(
                (self.projects / name / "member-result.json").read_text(
                    encoding="utf-8"
                )
            )
        )

    def test_parallel_members_use_one_distinct_cpu_and_one_host_thread(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("parallel lane test requires two available CPUs")
        members = [
            self.member("lane-a", delay=0.45),
            self.member("lane-b", delay=0.45),
            self.member("lane-c", delay=0.05),
        ]
        manifest = self.write_manifest(
            members,
            defaults={"env": {"ENSEMBLE_DEFAULT": "from-defaults"}},
        )
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "64MiB",
            "--gpu-device",
            "3",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "succeeded")

        records = [self._read_result(name) for name in ("lane-a", "lane-b", "lane-c")]
        for record in records:
            self.assertEqual(len(record["affinity"]), 1, record)
            environment = record["environment"]
            self.assertEqual(environment["SIMQ_LANE_CPU"], str(record["affinity"][0]))
            self.assertEqual(environment["SIMQ_CPU_LIST"], str(record["affinity"][0]))
            self.assertEqual(environment["SIMQ_CPU_COUNT"], "1")
            self.assertEqual(environment["CUDA_VISIBLE_DEVICES"], "3")
            for variable in (
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "KOKKOS_NUM_THREADS",
            ):
                self.assertEqual(environment[variable], "1", (variable, record))
            self.assertEqual(environment["ENSEMBLE_DEFAULT"], "from-defaults")

        first, second = records[:2]
        overlap = min(first["finished"], second["finished"]) - max(
            first["started"], second["started"]
        )
        self.assertGreater(overlap, 0.25, records)
        self.assertNotEqual(first["affinity"], second["affinity"])
        payload = self.status(job_id)
        self.assertEqual(payload["summary"]["succeeded"], 3)
        self.assertTrue(all(member["state"] == "succeeded" for member in payload["members"]))

    def test_launch_gate_serializes_startup_but_overlaps_ready_member_tail(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("launch-gate overlap test requires two available CPUs")
        marker = "TEST_STARTUP_READY_MARKER"
        leader_directory = self.projects / "gate-leader"
        leader_directory.mkdir()
        leader_program = r"""
import json
import os
import time
from pathlib import Path

record = {
    "started": time.time(),
    "affinity": sorted(os.sched_getaffinity(0)),
}
time.sleep(0.2)
print("TEST_STARTUP_", end="", flush=True)
time.sleep(0.05)
record["ready"] = time.time()
print("READY_MARKER", flush=True)
time.sleep(0.35)
record["finished"] = time.time()
Path("member-result.json").write_text(json.dumps(record), encoding="utf-8")
"""
        leader = {
            "name": "gate-leader",
            "cwd": str(leader_directory),
            "argv": [sys.executable, "-c", leader_program],
            "gpu_memory": "32MiB",
        }
        follower = self.member("gate-follower", delay=0.05, uvm_marker=marker)
        manifest = self.write_manifest(
            [leader, follower],
            max_parallel=2,
            gpu={"device": "0", "budget": "64MiB", "mps": "off"},
            launch_gate={
                "kind": "member-log-substring",
                "contains": marker,
                "timeout_seconds": 2,
            },
        )
        _result, job_id = self.submit(
            manifest,
            "--diagnostic-interval",
            "0.05",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        leader_record = self._read_result("gate-leader")
        follower_record = self._read_result("gate-follower")
        self.assertGreaterEqual(follower_record["started"], leader_record["ready"])
        self.assertLess(follower_record["started"], leader_record["finished"])
        events_path = Path(
            str(self.store.require_ensemble(job_id)["artifact_dir"])
        ) / "events.jsonl"
        events = [
            json.loads(line)
            for line in events_path.read_text(encoding="utf-8").splitlines()
            if line
        ]
        self.assertTrue(any(row["event"] == "launch-gate-wait" for row in events))
        self.assertTrue(any(row["event"] == "launch-gate-ready" for row in events))
        payload = self.status(job_id)
        self.assertEqual(payload["ensemble"]["launch_gate"]["contains"], marker)
        self.assertEqual(
            payload["ensemble"]["launch_gate_runtime"]["state"], "complete"
        )
        self.assertEqual(
            payload["ensemble"]["launch_gate_runtime"]["passed_count"], 2
        )
        self.assertEqual(payload["summary"]["succeeded"], 2)

    def test_launch_gate_without_marker_fails_before_follower_launch(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("launch-gate fallback test requires two available CPUs")
        leader = self.member("gate-no-marker-leader", delay=0.25)
        follower = self.member("gate-no-marker-follower", delay=0.05)
        manifest = self.write_manifest(
            [leader, follower],
            max_parallel=2,
            gpu={"device": "0", "budget": "64MiB", "mps": "off"},
            launch_gate={
                "kind": "member-log-substring",
                "contains": "THIS_MARKER_IS_NEVER_EMITTED",
                "timeout_seconds": 2,
            },
        )
        _result, job_id = self.submit(manifest)
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        rows = self.store.list_ensemble_members(job_id)
        self.assertEqual([row["state"] for row in rows], ["failed", "canceled"])
        self.assertFalse(
            (self.projects / "gate-no-marker-follower" / "member-result.json").exists()
        )
        events_path = Path(
            str(self.store.require_ensemble(job_id)["artifact_dir"])
        ) / "events.jsonl"
        events = [
            json.loads(line)
            for line in events_path.read_text(encoding="utf-8").splitlines()
            if line
        ]
        self.assertTrue(
            any(row["event"] == "launch-gate-marker-missing" for row in events)
        )

    def test_launch_gate_drains_large_pre_marker_output_before_reaping(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("launch-gate log-drain test requires two available CPUs")
        marker = "LARGE_LOG_READY_MARKER"
        leader_directory = self.projects / "gate-large-log-leader"
        leader_directory.mkdir()
        leader_program = (
            "import sys; "
            f"sys.stdout.write('X' * (256 * 1024) + {marker!r} + '\\n'); "
            "sys.stdout.flush()"
        )
        leader = {
            "name": "gate-large-log-leader",
            "cwd": str(leader_directory),
            "argv": [sys.executable, "-c", leader_program],
            "gpu_memory": "32MiB",
        }
        follower = self.member(
            "gate-large-log-follower", delay=0.05, uvm_marker=marker
        )
        manifest = self.write_manifest(
            [leader, follower],
            max_parallel=2,
            gpu={"device": "0", "budget": "64MiB", "mps": "off"},
            launch_gate={
                "kind": "member-log-substring",
                "contains": marker,
                "timeout_seconds": 2,
            },
        )
        _result, job_id = self.submit(manifest)
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        rows = self.store.list_ensemble_members(job_id)
        self.assertEqual([row["state"] for row in rows], ["succeeded", "succeeded"])
        self.assertTrue(
            (self.projects / "gate-large-log-follower" / "member-result.json").exists()
        )

    def test_launch_gate_rolls_across_three_members(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("rolling launch-gate test requires two available CPUs")
        marker = "ROLLING_SETUP_READY_MARKER"
        program = r"""
import json
import sys
import time
from pathlib import Path

name = sys.argv[1]
marker = sys.argv[2]
setup_delay = float(sys.argv[3])
tail_delay = float(sys.argv[4])
record = {"started": time.time()}
time.sleep(setup_delay)
record["ready"] = time.time()
print(marker, flush=True)
time.sleep(tail_delay)
record["finished"] = time.time()
Path("member-result.json").write_text(json.dumps(record), encoding="utf-8")
"""
        members: list[dict[str, Any]] = []
        # The first member deliberately releases a lane well before the
        # second finishes, making third-vs-second overlap deterministic.
        for index, tail_delay in enumerate((0.25, 0.8, 0.05)):
            directory = self.projects / f"rolling-gate-{index}"
            directory.mkdir()
            members.append(
                {
                    "name": f"rolling-gate-{index}",
                    "cwd": str(directory),
                    "argv": [
                        sys.executable,
                        "-c",
                        program,
                        f"rolling-gate-{index}",
                        marker,
                        "0.15",
                        str(tail_delay),
                    ],
                    "gpu_memory": "32MiB",
                }
            )
        manifest = self.write_manifest(
            members,
            max_parallel=2,
            gpu={"device": "0", "budget": "64MiB", "mps": "off"},
            launch_gate={
                "kind": "member-log-substring",
                "contains": marker,
                "timeout_seconds": 2,
            },
        )
        _result, job_id = self.submit(
            manifest,
            "--diagnostic-interval",
            "0.05",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        records = [self._read_result(f"rolling-gate-{index}") for index in range(3)]
        self.assertGreaterEqual(records[1]["started"], records[0]["ready"])
        self.assertGreaterEqual(records[2]["started"], records[1]["ready"])
        self.assertLess(records[1]["started"], records[0]["finished"])
        self.assertLess(records[2]["started"], records[1]["finished"])
        payload = self.status(job_id)
        self.assertEqual(payload["summary"]["succeeded"], 3)
        self.assertEqual(
            payload["ensemble"]["launch_gate_runtime"]["passed_count"], 3
        )

    def test_launch_gate_timeout_fails_before_follower_launch(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("launch-gate timeout test requires two available CPUs")
        leader = self.member("gate-timeout-leader", delay=1.0)
        follower = self.member("gate-timeout-follower", delay=0.05)
        manifest = self.write_manifest(
            [leader, follower],
            max_parallel=2,
            gpu={"device": "0", "budget": "64MiB", "mps": "off"},
            launch_gate={
                "kind": "member-log-substring",
                "contains": "THIS_MARKER_IS_NEVER_EMITTED",
                "timeout_seconds": 0.15,
            },
        )
        _result, job_id = self.submit(manifest)
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "failed")
        rows = self.store.list_ensemble_members(job_id)
        self.assertEqual(rows[0]["state"], "interrupted")
        self.assertEqual(rows[1]["state"], "interrupted")
        self.assertFalse(
            (self.projects / "gate-timeout-follower" / "member-result.json").exists()
        )
        outer_log = Path(str(self.store.require(job_id)["log_path"]))
        self.assertIn("launch-gate-timeout", outer_log.read_text(encoding="utf-8"))

    def test_gpu_budget_prevents_oversubscribed_members_from_overlapping(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("GPU admission test requires two available CPUs")
        members = [
            self.member("budget-a", delay=0.35, gpu_memory="96MiB"),
            self.member("budget-b", delay=0.35, gpu_memory="96MiB"),
            self.member("budget-c", delay=0.35, gpu_memory="32MiB"),
        ]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "128MiB",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        records = {name: self._read_result(name) for name in ("budget-a", "budget-b", "budget-c")}
        claims = {"budget-a": 96, "budget-b": 96, "budget-c": 32}
        events: list[tuple[float, int, str]] = []
        for name, record in records.items():
            events.append((float(record["started"]), 1, name))
            events.append((float(record["finished"]), -1, name))
        # End events sort before start events at an identical timestamp.
        events.sort(key=lambda event: (event[0], event[1]))
        active: set[str] = set()
        peak_members = 0
        for _timestamp, change, name in events:
            if change < 0:
                active.discard(name)
            else:
                active.add(name)
            peak_members = max(peak_members, len(active))
            admitted = sum(claims[item] for item in active)
            self.assertLessEqual(admitted, 128, (active, records))
            self.assertLessEqual(len(active), 2, (active, records))
        self.assertGreaterEqual(peak_members, 1)

    def test_token_descendant_keeps_lane_and_gpu_claim_after_leader_exits(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("descendant admission test requires two available CPUs")
        leader_directory = self.projects / "descendant-leader"
        leader_directory.mkdir()
        child_record = leader_directory / "child-timing.json"
        child_program = r"""
import json
import os
import sys
import time
from pathlib import Path

path = Path(sys.argv[1])
record = {"pid": os.getpid(), "started": time.time()}
path.write_text(json.dumps(record), encoding="utf-8")
time.sleep(0.7)
record["finished"] = time.time()
path.write_text(json.dumps(record), encoding="utf-8")
"""
        leader_program = r"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

child = subprocess.Popen([sys.executable, "-c", sys.argv[1], sys.argv[2]])
Path("leader-timing.json").write_text(json.dumps({
    "pid": os.getpid(),
    "child_pid": child.pid,
    "exited_at": time.time(),
}), encoding="utf-8")
"""
        leader = {
            "name": "descendant-leader",
            "cwd": str(leader_directory),
            "argv": [
                sys.executable,
                "-c",
                leader_program,
                child_program,
                str(child_record),
            ],
            "gpu_memory": "96MiB",
        }
        follower = self.member(
            "descendant-follower", delay=0.02, gpu_memory="96MiB"
        )
        manifest = self.write_manifest([leader, follower])
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "96MiB",
            "--diagnostic-interval",
            "0.05",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "succeeded")

        leader_timing = json.loads(
            (leader_directory / "leader-timing.json").read_text(encoding="utf-8")
        )
        child_timing = json.loads(child_record.read_text(encoding="utf-8"))
        follower_timing = self._read_result("descendant-follower")
        self.assertEqual(leader_timing["child_pid"], child_timing["pid"])
        self.assertLess(leader_timing["exited_at"], child_timing["finished"])
        self.assertGreaterEqual(
            follower_timing["started"],
            child_timing["finished"],
            "the follower started while the prior member's token-inheriting "
            "child still held its lane and GPU admission claim",
        )
        events = Path(
            str(self.store.require_ensemble(job_id)["artifact_dir"])
        ) / "events.jsonl"
        event_rows = [
            json.loads(line)
            for line in events.read_text(encoding="utf-8").splitlines()
            if line
        ]
        leader_finish = next(
            row
            for row in event_rows
            if row.get("event") == "member-finish"
            and row.get("member_name") == "descendant-leader"
        )
        follower_launch = next(
            row
            for row in event_rows
            if row.get("event") == "member-launch"
            and row.get("member_name") == "descendant-follower"
        )
        self.assertLessEqual(
            event_rows.index(leader_finish), event_rows.index(follower_launch)
        )

    def test_environmentless_new_session_orphan_keeps_lane_and_gpu_claim(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("orphan admission test requires two available CPUs")
        leader_directory = self.projects / "orphan-leader"
        leader_directory.mkdir()
        child_record = leader_directory / "orphan-child-timing.json"
        child_program = r"""
import json
import os
import sys
import time
from pathlib import Path

path = Path(sys.argv[1])
record = {"pid": os.getpid(), "started": time.time()}
path.write_text(json.dumps(record), encoding="utf-8")
time.sleep(0.7)
record["finished"] = time.time()
path.write_text(json.dumps(record), encoding="utf-8")
"""
        leader_program = r"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

child = subprocess.Popen(
    [sys.executable, "-c", sys.argv[1], sys.argv[2]],
    env={},
    start_new_session=True,
)
Path("orphan-leader-timing.json").write_text(json.dumps({
    "pid": os.getpid(),
    "child_pid": child.pid,
    "exited_at": time.time(),
}), encoding="utf-8")
"""
        leader = {
            "name": "orphan-leader",
            "cwd": str(leader_directory),
            "argv": [
                sys.executable,
                "-c",
                leader_program,
                child_program,
                str(child_record),
            ],
            "gpu_memory": "96MiB",
        }
        follower = self.member(
            "orphan-follower", delay=0.02, gpu_memory="96MiB"
        )
        manifest = self.write_manifest([leader, follower])
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "96MiB",
            "--diagnostic-interval",
            "0.05",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "succeeded")

        leader_timing = json.loads(
            (leader_directory / "orphan-leader-timing.json").read_text(
                encoding="utf-8"
            )
        )
        child_timing = json.loads(child_record.read_text(encoding="utf-8"))
        follower_timing = self._read_result("orphan-follower")
        self.assertEqual(leader_timing["child_pid"], child_timing["pid"])
        self.assertLess(leader_timing["exited_at"], child_timing["finished"])
        self.assertGreaterEqual(
            follower_timing["started"],
            child_timing["finished"],
            "the follower started while an environment-cleared, new-session "
            "orphan still held the prior member's lane and GPU claim",
        )
        events = Path(
            str(self.store.require_ensemble(job_id)["artifact_dir"])
        ) / "events.jsonl"
        waiting_events = [
            json.loads(line)
            for line in events.read_text(encoding="utf-8").splitlines()
            if "member-leader-exited-waiting-descendants" in line
        ]
        self.assertTrue(waiting_events)
        self.assertIn(child_timing["pid"], waiting_events[0]["descendant_pids"])

    def test_fork_after_process_snapshot_does_not_release_lane_or_gpu_claim(
        self,
    ) -> None:
        leader_directory = self.projects / "snapshot-race-leader"
        follower_directory = self.projects / "snapshot-race-follower"
        leader_directory.mkdir()
        follower_directory.mkdir()
        ready = leader_directory / "ready"
        release = leader_directory / "release"
        child_record = leader_directory / "child-timing.json"
        follower_record = follower_directory / "follower-timing.json"
        child_program = r"""
import json
import os
import sys
import time
from pathlib import Path

path = Path(sys.argv[1])
record = {"pid": os.getpid(), "started": time.time()}
path.write_text(json.dumps(record), encoding="utf-8")
time.sleep(0.7)
record["finished"] = time.time()
path.write_text(json.dumps(record), encoding="utf-8")
"""
        leader_program = r"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ready = Path(sys.argv[1])
release = Path(sys.argv[2])
ready.write_text("ready", encoding="utf-8")
while not release.exists():
    time.sleep(0.001)
child = subprocess.Popen(
    [sys.executable, "-c", sys.argv[3], sys.argv[4]],
    env={},
    start_new_session=True,
)
Path("leader-timing.json").write_text(json.dumps({
    "pid": os.getpid(),
    "child_pid": child.pid,
    "exited_at": time.time(),
}), encoding="utf-8")
"""
        follower_program = r"""
import json
import os
import sys
import time
from pathlib import Path

path = Path(sys.argv[1])
record = {"pid": os.getpid(), "started": time.time()}
path.write_text(json.dumps(record), encoding="utf-8")
time.sleep(0.02)
record["finished"] = time.time()
path.write_text(json.dumps(record), encoding="utf-8")
"""
        leader = {
            "name": "snapshot-race-leader",
            "cwd": str(leader_directory),
            "argv": [
                sys.executable,
                "-c",
                leader_program,
                str(ready),
                str(release),
                child_program,
                str(child_record),
            ],
            "gpu_memory": "96MiB",
        }
        follower = {
            "name": "snapshot-race-follower",
            "cwd": str(follower_directory),
            "argv": [
                sys.executable,
                "-c",
                follower_program,
                str(follower_record),
            ],
            "gpu_memory": "96MiB",
        }
        manifest = self.write_manifest([leader, follower])
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "1",
            "--gpu-budget",
            "96MiB",
        )
        assert job_id is not None
        job = self.store.claim_next(self.cpus)
        self.assertIsNotNone(job)
        assert job is not None
        self.assertEqual(int(job["id"]), job_id)
        protected_environment = {
            "SIMQ_JOB_ID": str(job_id),
            "SIMQ_LAUNCH_TOKEN": str(job["launch_token"]),
            "SIMQ_CPU_LIST": str(self.cpus[0]),
        }
        with preserve_child_subreaper(), patch.dict(
            os.environ, protected_environment, clear=False
        ):
            supervisor = EnsembleSupervisor(job_id, store=self.store)
            original_cohort_processes = supervisor._cohort_processes
            released = False

            def snapshot_then_release() -> tuple[dict[int, set[int]], set[int]]:
                nonlocal released
                snapshot = original_cohort_processes()
                if not released and supervisor.active and ready.exists():
                    released = True
                    release.write_text("go", encoding="utf-8")
                    deadline = time.monotonic() + 3
                    while time.monotonic() < deadline:
                        active = next(iter(supervisor.active.values()))
                        try:
                            stat = Path(f"/proc/{active.process.pid}/stat").read_text(
                                encoding="utf-8"
                            )
                            state = stat[stat.rfind(")") + 2 :].split()[0]
                        except OSError:
                            state = "gone"
                        if child_record.exists() and state in {"Z", "gone"}:
                            break
                        time.sleep(0.005)
                    else:
                        self.fail("leader did not fork and exit inside the widened race")
                return snapshot

            supervisor._cohort_processes = snapshot_then_release
            supervisor._install_signal_handlers = lambda: None
            supervisor._sample = lambda: None
            return_code = supervisor.run()

        self.assertTrue(released)
        self.assertEqual(return_code, 0)
        child_timing = json.loads(child_record.read_text(encoding="utf-8"))
        follower_timing = json.loads(follower_record.read_text(encoding="utf-8"))
        self.assertGreaterEqual(
            follower_timing["started"],
            child_timing["finished"],
            "the lane and GPU claim were released from a stale process snapshot",
        )

    def test_adopted_descendant_zombies_are_reaped_while_supervisor_is_live(
        self,
    ) -> None:
        orphan_directory = self.projects / "zombie-orphans"
        follower_directory = self.projects / "zombie-follower"
        orphan_directory.mkdir()
        follower_directory.mkdir()
        pid_record = orphan_directory / "child-pids.json"
        follower_ready = follower_directory / "ready"
        child_program = "import time; time.sleep(0.08)"
        orphan_program = r"""
import json
import subprocess
import sys
from pathlib import Path

children = [
    subprocess.Popen(
        [sys.executable, "-c", sys.argv[2]],
        env={},
        start_new_session=True,
    )
    for _index in range(4)
]
Path(sys.argv[1]).write_text(
    json.dumps([child.pid for child in children]),
    encoding="utf-8",
)
"""
        follower_program = r"""
import sys
import time
from pathlib import Path

Path(sys.argv[1]).write_text("ready", encoding="utf-8")
time.sleep(0.7)
"""
        members = [
            {
                "name": "zombie-orphans",
                "cwd": str(orphan_directory),
                "argv": [
                    sys.executable,
                    "-c",
                    orphan_program,
                    str(pid_record),
                    child_program,
                ],
                "gpu_memory": "1MiB",
            },
            {
                "name": "zombie-follower",
                "cwd": str(follower_directory),
                "argv": [
                    sys.executable,
                    "-c",
                    follower_program,
                    str(follower_ready),
                ],
                "gpu_memory": "1MiB",
            },
        ]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "1",
            "--gpu-budget",
            "1MiB",
        )
        assert job_id is not None
        job = self.store.claim_next(self.cpus)
        self.assertIsNotNone(job)
        assert job is not None
        self.assertEqual(int(job["id"]), job_id)
        observed: dict[str, Any] = {}

        def observe_adopted_children() -> None:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline and not follower_ready.exists():
                time.sleep(0.01)
            if not follower_ready.exists():
                observed["error"] = "follower did not start"
                return
            child_pids = json.loads(pid_record.read_text(encoding="utf-8"))
            records: list[dict[str, Any]] = []
            for pid in child_pids:
                try:
                    stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
                    fields = stat[stat.rfind(")") + 2 :].split()
                    records.append(
                        {"pid": pid, "state": fields[0], "parent_pid": int(fields[1])}
                    )
                except OSError:
                    records.append({"pid": pid, "state": "gone"})
            observed["records"] = records

        observer = threading.Thread(target=observe_adopted_children)
        observer.start()
        protected_environment = {
            "SIMQ_JOB_ID": str(job_id),
            "SIMQ_LAUNCH_TOKEN": str(job["launch_token"]),
            "SIMQ_CPU_LIST": str(self.cpus[0]),
        }
        try:
            with preserve_child_subreaper(), patch.dict(
                os.environ, protected_environment, clear=False
            ):
                supervisor = EnsembleSupervisor(job_id, store=self.store)
                supervisor._install_signal_handlers = lambda: None
                supervisor._sample = lambda: None
                return_code = supervisor.run()
        finally:
            observer.join(timeout=6)
            if pid_record.exists():
                for pid in json.loads(pid_record.read_text(encoding="utf-8")):
                    try:
                        os.waitpid(pid, os.WNOHANG)
                    except (ChildProcessError, OSError):
                        pass

        self.assertFalse(observer.is_alive())
        self.assertNotIn("error", observed)
        self.assertEqual(return_code, 0)
        self.assertTrue(observed.get("records"))
        self.assertTrue(
            all(record["state"] == "gone" for record in observed["records"]),
            observed,
        )

    def test_parent_and_manifest_affinity_overrides_cannot_escape_lane_cpu(self) -> None:
        if len(self.cpus) < 2:
            self.skipTest("affinity override test requires two available CPUs")
        os.environ["GOMP_CPU_AFFINITY"] = str(self.cpus[1])
        parent_member = self.member("parent-affinity", delay=0.25)
        manifest_member = self.member(
            "manifest-affinity",
            delay=0.25,
            environment={"GOMP_CPU_AFFINITY": str(self.cpus[0])},
        )
        manifest = self.write_manifest([parent_member, manifest_member])
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "64MiB",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "succeeded")

        parent_record = self._read_result("parent-affinity")
        manifest_record = self._read_result("manifest-affinity")
        for record, assigned_cpu in (
            (parent_record, self.cpus[0]),
            (manifest_record, self.cpus[1]),
        ):
            self.assertEqual(record["affinity"], [assigned_cpu], record)
            self.assertEqual(
                record["environment"]["SIMQ_LANE_CPU"], str(assigned_cpu)
            )
            self.assertIsNone(record["environment"]["GOMP_CPU_AFFINITY"])
            self.assertIsNone(record["environment"]["OMP_PLACES"])
            self.assertEqual(record["environment"]["OMP_PROC_BIND"], "FALSE")

    def test_openmpi_rank_one_inherits_second_lane_singleton_affinity(self) -> None:
        mpiexec = shutil.which("mpiexec")
        if mpiexec is None:
            self.skipTest("mpiexec is not installed")
        if len(self.cpus) < 2:
            self.skipTest("Open MPI lane test requires two available CPUs")
        os.environ["OMPI_MCA_hwloc_base_binding_policy"] = "core"
        os.environ["PRTE_MCA_hwloc_default_binding_policy"] = "core"
        os.environ["PRTE_MCA_bindto"] = "core"
        mpi_directory = self.projects / "openmpi-rank-one"
        mpi_directory.mkdir()
        rank_program = r"""
import json
import os
from pathlib import Path

Path("mpi-rank-result.json").write_text(json.dumps({
    "rank": os.environ.get("OMPI_COMM_WORLD_RANK"),
    "affinity": sorted(os.sched_getaffinity(0)),
    "lane_cpu": os.environ.get("SIMQ_LANE_CPU"),
    "ompi_binding": os.environ.get("OMPI_MCA_hwloc_base_binding_policy"),
    "prte_bindto": os.environ.get("PRTE_MCA_bindto"),
    "legacy_prte_binding": os.environ.get("PRTE_MCA_hwloc_default_binding_policy"),
}), encoding="utf-8")
"""
        occupy_lane_zero = self.member("openmpi-lane-zero", delay=1.0)
        rank_one = {
            "name": "openmpi-rank-one",
            "cwd": str(mpi_directory),
            "argv": [mpiexec, "-np", "1", sys.executable, "-c", rank_program],
            "env": {
                "OMPI_MCA_hwloc_base_binding_policy": "hwthread",
                "PRTE_MCA_bindto": "hwthread",
                "PRTE_MCA_hwloc_default_binding_policy": "hwthread",
            },
            "gpu_memory": "32MiB",
        }
        manifest = self.write_manifest([occupy_lane_zero, rank_one])
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "64MiB",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        parent = self.store.require(job_id)
        self.assertEqual(
            parent["state"],
            "succeeded",
            Path(str(parent["log_path"])).read_text(encoding="utf-8"),
        )

        lane_zero_record = self._read_result("openmpi-lane-zero")
        self.assertEqual(lane_zero_record["affinity"], [self.cpus[0]])
        rank_record = json.loads(
            (mpi_directory / "mpi-rank-result.json").read_text(encoding="utf-8")
        )
        self.assertEqual(rank_record["rank"], "0")
        self.assertEqual(rank_record["lane_cpu"], str(self.cpus[1]))
        self.assertEqual(rank_record["affinity"], [self.cpus[1]])
        self.assertEqual(rank_record["ompi_binding"], "none")
        self.assertIn(rank_record["prte_bindto"], (None, "none"))
        self.assertIsNone(rank_record["legacy_prte_binding"])
        mpi_row = next(
            row
            for row in self.store.list_ensemble_members(job_id)
            if row["name"] == "openmpi-rank-one"
        )
        mpi_log = Path(str(mpi_row["log_path"])).read_text(encoding="utf-8")
        self.assertIn("environment: OMPI_MCA_hwloc_base_binding_policy=none", mpi_log)
        self.assertIn("environment: PRTE_MCA_bindto=none", mpi_log)
        self.assertNotIn("hwthread", mpi_log)

    def test_member_environment_overrides_are_merged_with_defaults(self) -> None:
        member = self.member(
            "environment",
            environment={"MEMBER_ONLY": "member-value"},
        )
        manifest = self.write_manifest(
            [member],
            defaults={"env": {"ENSEMBLE_DEFAULT": "default-value"}},
        )
        _result, job_id = self.submit(manifest, "--gpu-budget", "32MiB")
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        environment = self._read_result("environment")["environment"]
        self.assertEqual(environment["ENSEMBLE_DEFAULT"], "default-value")
        self.assertEqual(environment["MEMBER_ONLY"], "member-value")

    def test_launch_rechecks_cwd_and_rejects_symlink_escape(self) -> None:
        submitted_directory = self.projects / "runtime-cwd"
        member = self.member("runtime-cwd", cwd=submitted_directory)
        manifest = self.write_manifest([member])
        _result, job_id = self.submit(manifest, "--gpu-budget", "32MiB")
        assert job_id is not None

        retained_directory = self.projects / "runtime-cwd-original"
        submitted_directory.rename(retained_directory)
        outside = self.base / "runtime-cwd-outside"
        outside.mkdir()
        submitted_directory.symlink_to(outside, target_is_directory=True)
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        parent = self.store.require(job_id)
        row = self.store.list_ensemble_members(job_id)[0]
        self.assertEqual(parent["state"], "failed")
        self.assertEqual(row["state"], "failed")
        self.assertIn("resolves outside", str(row["message"]))
        self.assertFalse((outside / "member.pid").exists())
        self.assertIn(
            "resolves outside",
            Path(str(row["log_path"])).read_text(encoding="utf-8"),
        )

    def test_failure_keeps_going_by_default(self) -> None:
        members = [
            self.member("fails", exit_code=19),
            self.member("still-runs"),
        ]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest, "--max-parallel", "1", "--gpu-budget", "32MiB"
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertEqual(self.store.require(job_id)["state"], "failed")
        self.assertTrue((self.projects / "still-runs" / "member-result.json").is_file())
        payload = self.status(job_id)
        states = {member["name"]: member["state"] for member in payload["members"]}
        self.assertEqual(states, {"fails": "failed", "still-runs": "succeeded"})
        self.assertEqual(payload["summary"]["failed"], 1)
        self.assertEqual(payload["summary"]["succeeded"], 1)

    def test_fail_fast_cancels_members_that_have_not_started(self) -> None:
        members = [
            self.member("fails-fast", exit_code=23),
            self.member("must-not-run"),
        ]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "1",
            "--gpu-budget",
            "32MiB",
            "--fail-fast",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)
        self.assertFalse((self.projects / "must-not-run" / "member.pid").exists())
        payload = self.status(job_id)
        states = {member["name"]: member["state"] for member in payload["members"]}
        self.assertEqual(states["fails-fast"], "failed")
        self.assertEqual(states["must-not-run"], "canceled")

    def test_fail_fast_rescans_after_kill_snapshot_for_new_session_child(
        self,
    ) -> None:
        if len(self.cpus) < 2:
            self.skipTest("fail-fast kill race test requires two available CPUs")
        resistant_directory = self.projects / "kill-race-resistant"
        failure_directory = self.projects / "kill-race-failure"
        resistant_directory.mkdir()
        failure_directory.mkdir()
        resistant_ready = resistant_directory / "ready"
        spawn_release = resistant_directory / "spawn-release"
        child_record = resistant_directory / "child.json"
        child_finished = resistant_directory / "child-finished"
        parent_record = resistant_directory / "parent.json"
        child_program = r"""
import json
import os
import signal
import sys
import time
from pathlib import Path

signal.signal(signal.SIGTERM, signal.SIG_IGN)
Path(sys.argv[1]).write_text(json.dumps({
    "pid": os.getpid(),
    "started": time.time(),
}), encoding="utf-8")
time.sleep(1.5)
Path(sys.argv[2]).write_text("finished", encoding="utf-8")
"""
        resistant_program = r"""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

ready = Path(sys.argv[1])
release = Path(sys.argv[2])
signal.signal(signal.SIGTERM, signal.SIG_IGN)
ready.write_text("ready", encoding="utf-8")
while not release.exists():
    time.sleep(0.001)
child = subprocess.Popen(
    [sys.executable, "-c", sys.argv[3], sys.argv[4], sys.argv[5]],
    env={},
    start_new_session=True,
)
Path(sys.argv[6]).write_text(json.dumps({
    "pid": os.getpid(),
    "child_pid": child.pid,
}), encoding="utf-8")
while True:
    time.sleep(1)
"""
        failure_program = r"""
import sys
import time
from pathlib import Path

ready = Path(sys.argv[1])
while not ready.exists():
    time.sleep(0.001)
raise SystemExit(37)
"""
        members = [
            {
                "name": "kill-race-resistant",
                "cwd": str(resistant_directory),
                "argv": [
                    sys.executable,
                    "-c",
                    resistant_program,
                    str(resistant_ready),
                    str(spawn_release),
                    child_program,
                    str(child_record),
                    str(child_finished),
                    str(parent_record),
                ],
                "gpu_memory": "1MiB",
            },
            {
                "name": "kill-race-failure",
                "cwd": str(failure_directory),
                "argv": [
                    sys.executable,
                    "-c",
                    failure_program,
                    str(resistant_ready),
                ],
                "gpu_memory": "1MiB",
            },
        ]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "2MiB",
            "--fail-fast",
        )
        assert job_id is not None
        job = self.store.claim_next(self.cpus)
        self.assertIsNotNone(job)
        assert job is not None
        self.assertEqual(int(job["id"]), job_id)
        protected_environment = {
            "SIMQ_JOB_ID": str(job_id),
            "SIMQ_LAUNCH_TOKEN": str(job["launch_token"]),
            "SIMQ_CPU_LIST": format_cpu_list(self.cpus[:2]),
        }
        child_pid: int | None = None
        kill_calls = 0
        first_kill_snapshot_released = False
        with preserve_child_subreaper(), patch.dict(
            os.environ, protected_environment, clear=False
        ):
            supervisor = EnsembleSupervisor(job_id, store=self.store)
            original_signal_active = supervisor._signal_active
            original_process_records = ensemble_module._process_records
            current_signal: signal.Signals | None = None

            def signal_with_context(sig: signal.Signals) -> None:
                nonlocal current_signal, kill_calls
                current_signal = sig
                if sig == signal.SIGKILL:
                    kill_calls += 1
                try:
                    original_signal_active(sig)
                finally:
                    current_signal = None

            def release_after_first_kill_snapshot() -> dict[int, tuple[str, int, int]]:
                nonlocal first_kill_snapshot_released, child_pid
                records = original_process_records()
                if current_signal == signal.SIGKILL and not first_kill_snapshot_released:
                    first_kill_snapshot_released = True
                    spawn_release.write_text("spawn", encoding="utf-8")
                    deadline = time.monotonic() + 3
                    while time.monotonic() < deadline:
                        if child_record.exists() and parent_record.exists():
                            try:
                                child_pid = int(
                                    json.loads(
                                        child_record.read_text(encoding="utf-8")
                                    )["pid"]
                                )
                            except (OSError, json.JSONDecodeError, KeyError, ValueError):
                                pass
                            else:
                                break
                        time.sleep(0.005)
                    else:
                        self.fail("resistant member did not create its new-session child")
                return records

            supervisor._signal_active = signal_with_context
            supervisor._install_signal_handlers = lambda: None
            supervisor._sample = lambda: None
            try:
                with patch.object(
                    ensemble_module, "FAIL_FAST_GRACE_SECONDS", 0.1
                ), patch.object(
                    ensemble_module,
                    "_process_records",
                    side_effect=release_after_first_kill_snapshot,
                ):
                    return_code = supervisor.run()
            finally:
                if child_pid is not None and process_is_live(child_pid):
                    try:
                        os.killpg(child_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass

        self.assertEqual(return_code, 1)
        self.assertTrue(first_kill_snapshot_released)
        self.assertGreaterEqual(kill_calls, 2)
        assert child_pid is not None
        self.assertFalse(process_is_live(child_pid))
        self.assertFalse(
            child_finished.exists(),
            "the child survived until natural completion instead of a repeated SIGKILL scan",
        )
        states = {
            member["name"]: member["state"]
            for member in self.store.list_ensemble_members(job_id)
        }
        self.assertEqual(states["kill-race-failure"], "failed")
        self.assertEqual(states["kill-race-resistant"], "canceled")

    def test_exception_cleanup_rescans_after_kill_snapshot_escape(self) -> None:
        member_directory = self.projects / "exception-kill-race"
        member_directory.mkdir()
        spawn_release = member_directory / "spawn-release"
        child_record = member_directory / "child.json"
        child_finished = member_directory / "child-finished"
        parent_record = member_directory / "parent.json"
        child_program = r"""
import json
import os
import sys
import time
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "pid": os.getpid(),
    "started": time.time(),
}), encoding="utf-8")
time.sleep(1.5)
Path(sys.argv[2]).write_text("finished", encoding="utf-8")
"""
        member_program = r"""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

release = Path(sys.argv[1])
signal.signal(signal.SIGTERM, signal.SIG_IGN)
while not release.exists():
    time.sleep(0.001)
child = subprocess.Popen(
    [sys.executable, "-c", sys.argv[2], sys.argv[3], sys.argv[4]],
    env={},
    start_new_session=True,
)
Path(sys.argv[5]).write_text(json.dumps({
    "pid": os.getpid(),
    "child_pid": child.pid,
}), encoding="utf-8")
while True:
    time.sleep(1)
"""
        member = {
            "name": "exception-kill-race",
            "cwd": str(member_directory),
            "argv": [
                sys.executable,
                "-c",
                member_program,
                str(spawn_release),
                child_program,
                str(child_record),
                str(child_finished),
                str(parent_record),
            ],
            "gpu_memory": "1MiB",
        }
        manifest = self.write_manifest([member])
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "1",
            "--gpu-budget",
            "1MiB",
        )
        assert job_id is not None
        job = self.store.claim_next(self.cpus)
        self.assertIsNotNone(job)
        assert job is not None
        self.assertEqual(int(job["id"]), job_id)
        protected_environment = {
            "SIMQ_JOB_ID": str(job_id),
            "SIMQ_LAUNCH_TOKEN": str(job["launch_token"]),
            "SIMQ_CPU_LIST": str(self.cpus[0]),
        }
        child_pid: int | None = None
        kill_calls = 0
        first_kill_snapshot_released = False
        with preserve_child_subreaper(), patch.dict(
            os.environ, protected_environment, clear=False
        ):
            supervisor = EnsembleSupervisor(job_id, store=self.store)
            original_signal_active = supervisor._signal_active
            original_process_records = ensemble_module._process_records
            current_signal: signal.Signals | None = None

            def signal_with_context(sig: signal.Signals) -> None:
                nonlocal current_signal, kill_calls
                current_signal = sig
                if sig == signal.SIGKILL:
                    kill_calls += 1
                try:
                    original_signal_active(sig)
                finally:
                    current_signal = None

            def release_after_first_kill_snapshot() -> dict[int, tuple[str, int, int]]:
                nonlocal first_kill_snapshot_released, child_pid
                records = original_process_records()
                if current_signal == signal.SIGKILL and not first_kill_snapshot_released:
                    first_kill_snapshot_released = True
                    spawn_release.write_text("spawn", encoding="utf-8")
                    deadline = time.monotonic() + 3
                    while time.monotonic() < deadline:
                        if child_record.exists() and parent_record.exists():
                            try:
                                child_pid = int(
                                    json.loads(
                                        child_record.read_text(encoding="utf-8")
                                    )["pid"]
                                )
                            except (OSError, json.JSONDecodeError, KeyError, ValueError):
                                pass
                            else:
                                break
                        time.sleep(0.005)
                    else:
                        self.fail("member did not create its cleanup-race child")
                return records

            def fail_sample() -> None:
                raise RuntimeError("injected process-sample failure")

            supervisor._signal_active = signal_with_context
            supervisor._install_signal_handlers = lambda: None
            supervisor._sample = fail_sample
            try:
                with patch.object(
                    ensemble_module,
                    "_process_records",
                    side_effect=release_after_first_kill_snapshot,
                ):
                    with self.assertRaisesRegex(
                        RuntimeError, "injected process-sample failure"
                    ):
                        supervisor.run()
            finally:
                if child_pid is not None and process_is_live(child_pid):
                    try:
                        os.killpg(child_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass

        self.assertTrue(first_kill_snapshot_released)
        self.assertGreaterEqual(kill_calls, 2)
        assert child_pid is not None
        self.assertFalse(process_is_live(child_pid))
        self.assertFalse(
            child_finished.exists(),
            "exception cleanup let the escaped child finish naturally",
        )
        row = self.store.list_ensemble_members(job_id)[0]
        self.assertEqual(row["state"], "interrupted")

    def test_fail_fast_preserves_success_already_complete_in_same_reap_pass(
        self,
    ) -> None:
        if len(self.cpus) < 2:
            self.skipTest("simultaneous completion test requires two available CPUs")
        release = self.projects / "simultaneous-release"
        completion_program = r"""
import os
import sys
import time
from pathlib import Path

ready = Path(sys.argv[1])
release = Path(sys.argv[2])
exit_code = int(sys.argv[3])
Path("member.pid").write_text(str(os.getpid()), encoding="utf-8")
ready.write_text("ready", encoding="utf-8")
while not release.exists():
    time.sleep(0.005)
raise SystemExit(exit_code)
"""
        members: list[dict[str, Any]] = []
        for name, exit_code in (("simultaneous-fail", 31), ("simultaneous-pass", 0)):
            directory = self.projects / name
            directory.mkdir()
            members.append(
                {
                    "name": name,
                    "cwd": str(directory),
                    "argv": [
                        sys.executable,
                        "-c",
                        completion_program,
                        str(directory / "ready"),
                        str(release),
                        str(exit_code),
                    ],
                    "gpu_memory": "32MiB",
                }
            )
        members.append(self.member("simultaneous-pending", gpu_memory="32MiB"))
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "64MiB",
            "--fail-fast",
        )
        assert job_id is not None
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        supervisor_pid: int | None = None
        supervisor_stopped = False
        try:
            ready_paths = [
                self.projects / "simultaneous-fail" / "ready",
                self.projects / "simultaneous-pass" / "ready",
            ]
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline:
                parent = self.store.require(job_id)
                if all(path.is_file() for path in ready_paths) and parent["pid"]:
                    supervisor_pid = int(parent["pid"])
                    break
                time.sleep(0.02)
            else:
                self.fail("both active members did not reach the completion barrier")

            assert supervisor_pid is not None
            os.kill(supervisor_pid, signal.SIGSTOP)
            supervisor_stopped = True
            deadline = time.monotonic() + 3
            while time.monotonic() < deadline:
                try:
                    stat = Path(f"/proc/{supervisor_pid}/stat").read_text(
                        encoding="utf-8"
                    )
                except OSError:
                    self.fail("ensemble supervisor exited before SIGSTOP took effect")
                state = stat[stat.rfind(")") + 2 :].split()[0]
                if state in {"T", "t"}:
                    break
                time.sleep(0.01)
            else:
                self.fail("ensemble supervisor did not stop at the reap barrier")

            member_pids = [
                int(
                    (self.projects / name / "member.pid").read_text(
                        encoding="utf-8"
                    )
                )
                for name in ("simultaneous-fail", "simultaneous-pass")
            ]
            release.touch()
            deadline = time.monotonic() + 3
            while time.monotonic() < deadline:
                if all(not process_is_live(pid) for pid in member_pids):
                    break
                time.sleep(0.01)
            else:
                self.fail("both active members did not finish behind the reap barrier")

            os.kill(supervisor_pid, signal.SIGCONT)
            supervisor_stopped = False
            stdout, stderr = worker.communicate(timeout=10)
            self.assertEqual(worker.returncode, 0, stdout + stderr)
            rows = self.store.list_ensemble_members(job_id)
            states = {row["name"]: row["state"] for row in rows}
            self.assertEqual(
                states,
                {
                    "simultaneous-fail": "failed",
                    "simultaneous-pass": "succeeded",
                    "simultaneous-pending": "canceled",
                },
                {
                    "parent_state": self.store.require(job_id)["state"],
                    "parent_message": self.store.require(job_id)["message"],
                    "stdout": stdout,
                    "stderr": stderr,
                    "parent_log": Path(
                        str(self.store.require(job_id)["log_path"])
                    ).read_text(encoding="utf-8"),
                },
            )
            self.assertFalse(
                (self.projects / "simultaneous-pending" / "member.pid").exists()
            )
            events_path = Path(
                str(self.store.require_ensemble(job_id)["artifact_dir"])
            ) / "events.jsonl"
            events = [
                json.loads(line)
                for line in events_path.read_text(encoding="utf-8").splitlines()
                if line
            ]
            fail_fast_index = next(
                index
                for index, event in enumerate(events)
                if event.get("event") == "fail-fast-triggered"
            )
            completed_before_fail_fast = {
                event.get("member_name"): event.get("state")
                for event in events[:fail_fast_index]
                if event.get("event") == "member-finish"
            }
            self.assertEqual(
                completed_before_fail_fast,
                {
                    "simultaneous-fail": "failed",
                    "simultaneous-pass": "succeeded",
                },
            )
        finally:
            if supervisor_stopped and supervisor_pid is not None:
                try:
                    os.kill(supervisor_pid, signal.SIGCONT)
                except ProcessLookupError:
                    pass
            if worker.poll() is None:
                try:
                    self.store.cancel(job_id)
                except (KeyError, RuntimeError):
                    pass
                try:
                    worker.communicate(timeout=3)
                except subprocess.TimeoutExpired:
                    if supervisor_pid is not None:
                        try:
                            os.killpg(supervisor_pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    worker.kill()
                    worker.wait()


class EnsembleDiagnosticsTests(EnsembleFixture, unittest.TestCase):
    def _install_fake_nvidia_smi(self) -> None:
        binary_directory = self.base / "bin"
        binary_directory.mkdir()
        binary = binary_directory / "nvidia-smi"
        binary.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "mapping = {\n"
            " 'timestamp': '2026/08/24 12:34:56.000',\n"
            " 'index': '0', 'uuid': 'GPU-test-uuid',\n"
            " 'memory.total': '16384', 'memory.used': '2048',\n"
            " 'memory.free': '14336', 'utilization.gpu': '42',\n"
            " 'utilization.memory': '17', 'pstate': 'P2',\n"
            " 'temperature.gpu': '51', 'power.draw': '87.50',\n"
            "}\n"
            "query = next((a.split('=', 1)[1] for a in sys.argv "
            "if a.startswith('--query-gpu=')), '')\n"
            "if query:\n"
            " print(', '.join(mapping.get(field, '0') for field in query.split(',')))\n"
            "elif '-L' in sys.argv:\n"
            " print('GPU 0: Fake GPU (UUID: GPU-test-uuid)')\n",
            encoding="utf-8",
        )
        binary.chmod(0o700)
        os.environ["PATH"] = f"{binary_directory}{os.pathsep}{os.environ.get('PATH', '')}"

    def test_diagnostic_artifacts_and_uvm_final_marker_are_exposed(self) -> None:
        self._install_fake_nvidia_smi()
        uvm_environment = {
            "MSIC_KOKKOS_CUDAUVM_BACKING": "managed-pool",
            "MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES": str(512 * 1024**2),
        }
        member = self.member(
            "diagnostic",
            delay=0.45,
            gpu_memory="512MiB",
            uvm_marker=UVM_FINAL_MARKER,
            environment=uvm_environment,
        )
        manifest = self.write_manifest([member])
        _result, job_id = self.submit(
            manifest,
            "--gpu-budget",
            "768MiB",
            "--diagnostic-interval",
            "0.1",
        )
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        artifact_directory = self.base / "state" / "ensembles" / f"{job_id:06d}"
        events = artifact_directory / "events.jsonl"
        gpu_samples = artifact_directory / "gpu_samples.jsonl"
        process_samples = artifact_directory / "process_samples.jsonl"
        for path in (events, gpu_samples, process_samples):
            self.assertTrue(path.is_file(), path)
            lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
            self.assertTrue(lines, path)
            for line in lines:
                self.assertIsInstance(json.loads(line), dict)
        event_text = events.read_text(encoding="utf-8")
        self.assertIn("diagnostic", event_text)
        self.assertIn("launch", event_text.lower())
        self.assertIn("finish", event_text.lower())
        gpu_text = gpu_samples.read_text(encoding="utf-8")
        self.assertIn("GPU-test-uuid", gpu_text)
        self.assertIn("2048", gpu_text)

        payload = self.status(job_id)
        row = payload["members"][0]
        self.assertEqual(row["uvm_peak_live_bytes"], 284464608)
        self.assertEqual(row["uvm_peak_reserved_bytes"], 284679424)
        self.assertEqual(row["uvm_pool_reclaims"], 3)
        self.assertEqual(row["uvm_teardown_ok"], 1)
        self.assertEqual(row["uvm_release"], "success")
        member_log = Path(row["log_path"])
        self.assertTrue(member_log.is_file())
        log_text = member_log.read_text(encoding="utf-8")
        self.assertIn("member-output name=diagnostic", log_text)
        self.assertIn(UVM_FINAL_MARKER, log_text)
        self.assertIn("MSIC_KOKKOS_CUDAUVM_LIMIT_BYTES", log_text)
        self.assertEqual(member_log.stat().st_mode & 0o777, 0o600)

    def test_status_and_log_cli_support_names_and_zero_based_indices(self) -> None:
        manifest = self.write_manifest([self.member("visible-member")])
        _result, job_id = self.submit(manifest, "--gpu-budget", "32MiB")
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        human = self.run_cli("ensemble", "status", str(job_id))
        self.assertEqual(human.returncode, 0, human.stderr)
        self.assertIn("visible-member", human.stdout)
        self.assertIn("succeeded", human.stdout)

        by_name = self.run_cli(
            "ensemble", "log", str(job_id), "visible-member", "-n", "20"
        )
        self.assertEqual(by_name.returncode, 0, by_name.stderr)
        self.assertIn("member-output name=visible-member", by_name.stdout)
        by_index = self.run_cli("ensemble", "log", str(job_id), "0", "-n", "20")
        self.assertEqual(by_index.returncode, 0, by_index.stderr)
        self.assertIn("member-output name=visible-member", by_index.stdout)

        missing = self.run_cli("ensemble", "log", str(job_id), "does-not-exist")
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("member", missing.stderr.lower())

    def test_mps_off_forces_every_member_onto_a_nonexistent_bypass_pipe(self) -> None:
        inherited_pipe = self.base / "inherited-mps-pipe"
        inherited_pipe.mkdir()
        os.environ["CUDA_MPS_PIPE_DIRECTORY"] = str(inherited_pipe)
        manifest = self.write_manifest(
            [
                self.member(
                    "mps-off",
                    environment={"CUDA_MPS_PIPE_DIRECTORY": str(inherited_pipe)},
                )
            ],
            gpu={"device": "0", "budget": "32MiB", "mps": "off"},
        )
        _result, job_id = self.submit(manifest)
        assert job_id is not None
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        record = json.loads(
            (self.projects / "mps-off" / "member-result.json").read_text(
                encoding="utf-8"
            )
        )
        artifact_directory = Path(
            str(self.store.require_ensemble(job_id)["artifact_dir"])
        )
        bypass_pipe = artifact_directory / "mps-disabled"
        self.assertEqual(
            record["environment"]["CUDA_MPS_PIPE_DIRECTORY"], str(bypass_pipe)
        )
        self.assertNotEqual(str(bypass_pipe), str(inherited_pipe))
        self.assertFalse(bypass_pipe.exists())
        self.assertIn(
            '"event": "mps-disabled"',
            (artifact_directory / "events.jsonl").read_text(encoding="utf-8"),
        )

    def test_required_mps_is_rejected_before_queue_mutation(self) -> None:
        manifest = self.write_manifest(
            [self.member("mps-required")],
            gpu={
                "device": "0",
                "budget": "32MiB",
                "mps": "required",
            },
        )
        result, job_id = self.submit(manifest, expect_success=False)
        self.assertIsNone(job_id)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("supports only 'off'", result.stderr)
        self.assertEqual(self.store.list_jobs(include_all=True), [])

    def test_stored_required_mps_fails_closed_before_member_launch(self) -> None:
        manifest = self.write_manifest(
            [self.member("stored-mps-required")],
            gpu={"device": "0", "budget": "32MiB", "mps": "off"},
        )
        _result, job_id = self.submit(manifest)
        assert job_id is not None
        with self.store.connect() as connection:
            connection.execute(
                "UPDATE ensembles SET mps_mode = 'required' WHERE job_id = ?",
                (job_id,),
            )
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        parent = self.store.require(job_id)
        row = self.store.list_ensemble_members(job_id)[0]
        self.assertEqual(parent["state"], "failed")
        self.assertEqual(row["state"], "interrupted")
        self.assertFalse(
            (self.projects / "stored-mps-required" / "member.pid").exists()
        )
        parent_log = Path(str(parent["log_path"])).read_text(encoding="utf-8")
        self.assertIn("unsupported MPS mode", parent_log)
        self.assertIn("supervisor failed", str(row["message"]))

    def test_malformed_stored_launch_gate_fails_closed_before_member_launch(self) -> None:
        marker = "STORED_GATE_READY"
        manifest = self.write_manifest(
            [self.member("stored-malformed-gate", uvm_marker=marker)],
            gpu={"device": "0", "budget": "32MiB", "mps": "off"},
            launch_gate={
                "kind": "member-log-substring",
                "contains": marker,
                "timeout_seconds": 2,
            },
        )
        _result, job_id = self.submit(manifest)
        assert job_id is not None
        with self.store.connect() as connection:
            row = connection.execute(
                "SELECT spec_json FROM ensembles WHERE job_id = ?", (job_id,)
            ).fetchone()
            assert row is not None
            stored = json.loads(str(row["spec_json"]))
            stored["launch_gate"] = "corrupt"
            connection.execute(
                "UPDATE ensembles SET spec_json = ? WHERE job_id = ?",
                (json.dumps(stored, sort_keys=True), job_id),
            )
        worker = self.run_worker()
        self.assertEqual(worker.returncode, 0, worker.stdout + worker.stderr)

        parent = self.store.require(job_id)
        member = self.store.list_ensemble_members(job_id)[0]
        self.assertEqual(parent["state"], "failed")
        self.assertEqual(member["state"], "interrupted")
        self.assertFalse(
            (self.projects / "stored-malformed-gate" / "member.pid").exists()
        )
        parent_log = Path(str(parent["log_path"])).read_text(encoding="utf-8")
        self.assertIn("invalid stored ensemble launch gate", parent_log)


class EnsembleAtomicityTests(EnsembleFixture, unittest.TestCase):
    def test_partial_member_log_failure_rolls_back_database_and_all_files(self) -> None:
        first = self.member("atomic-first", gpu_memory="0B")
        second = self.member("atomic-second", gpu_memory="0B")
        members = [first, second]
        canonical = {
            "schema_version": "simqueue-ensemble-v1",
            "name": "atomic-ensemble",
            "max_parallel": 1,
            "gpu": {
                "device": "0",
                "budget_bytes": 0,
                "mps": "off",
                "mps_pipe_directory": None,
            },
            "diagnostics_interval": 1.0,
            "fail_fast": False,
            "allow_shared_cwd": False,
            "members": members,
        }
        stored_members = [
            {
                **member,
                "gpu_admission_bytes": 0,
                "uvm_arena_bytes": 0,
            }
            for member in members
        ]
        real_open = Path.open

        def fail_second_member_log(
            path: Path, *arguments: Any, **keywords: Any
        ) -> Any:
            if path.name.startswith("0001-") and path.suffix == ".log":
                raise PermissionError("injected second-member log failure")
            return real_open(path, *arguments, **keywords)

        with patch.object(Path, "open", new=fail_second_member_log):
            with self.assertRaisesRegex(
                PermissionError, "injected second-member log failure"
            ):
                self.store.submit_ensemble(
                    name="atomic-ensemble",
                    cwd=self.projects,
                    environment=dict(os.environ),
                    env_override_keys=[],
                    requested_cores=1,
                    runner_prefix=[
                        sys.executable,
                        str(self.launcher),
                        "_run-ensemble",
                    ],
                    schema_version="simqueue-ensemble-v1",
                    canonical_spec=canonical,
                    members=stored_members,
                    max_parallel=1,
                    gpu_device="0",
                    gpu_budget_bytes=0,
                    diagnostics_interval=1.0,
                    mps_mode="off",
                    mps_pipe_directory=None,
                    fail_fast=False,
                    allow_shared_cwd=False,
                )

        self.assertEqual(self.store.list_jobs(include_all=True), [])
        with self.store.connect() as connection:
            self.assertEqual(
                connection.execute("SELECT count(*) FROM ensembles").fetchone()[0], 0
            )
            self.assertEqual(
                connection.execute(
                    "SELECT count(*) FROM ensemble_members"
                ).fetchone()[0],
                0,
            )
        self.assertEqual(list((self.base / "state" / "logs").iterdir()), [])
        ensembles_directory = self.base / "state" / "ensembles"
        self.assertFalse(ensembles_directory.exists())

    def test_legacy_database_adds_ensemble_schema_concurrently(self) -> None:
        legacy_path = self.base / "legacy-ensemble.sqlite3"
        QueueStore(path=legacy_path)
        with closing(sqlite3.connect(legacy_path)) as connection, connection:
            connection.execute("DROP TABLE ensemble_members")
            connection.execute("DROP TABLE ensembles")
            connection.execute("PRAGMA user_version = 2")

        with ThreadPoolExecutor(max_workers=8) as pool:
            stores = list(
                pool.map(lambda _index: QueueStore(path=legacy_path), range(16))
            )
        self.assertEqual(len(stores), 16)
        with closing(sqlite3.connect(legacy_path)) as connection:
            tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                ).fetchall()
            }
            version = int(connection.execute("PRAGMA user_version").fetchone()[0])
            member_columns = {
                str(row[1])
                for row in connection.execute(
                    "PRAGMA table_info(ensemble_members)"
                ).fetchall()
            }
            ensemble_columns = {
                str(row[1])
                for row in connection.execute(
                    "PRAGMA table_info(ensembles)"
                ).fetchall()
            }
        self.assertIn("ensembles", tables)
        self.assertIn("ensemble_members", tables)
        self.assertEqual(version, SCHEMA_VERSION)
        self.assertTrue(
            {
                "launch_paused",
                "control_generation",
                "control_changed_at",
                "drained_generation",
                "drained_at",
                "supervisor_protocol",
                "supervisor_pid",
            }.issubset(ensemble_columns)
        )
        self.assertTrue(
            {
                "member_index",
                "gpu_admission_bytes",
                "uvm_arena_bytes",
                "member_token",
                "log_path",
            }.issubset(member_columns)
        )


class EnsembleCancellationTests(EnsembleFixture, unittest.TestCase):
    def test_canceling_outer_job_terminates_active_lanes_and_never_starts_queued_lane(
        self,
    ) -> None:
        if len(self.cpus) < 2:
            self.skipTest("cancellation test requires two available CPUs")
        members = [
            self.member("long-a", delay=30),
            self.member("long-b", delay=30),
            self.member("never-started", delay=0.01),
        ]
        manifest = self.write_manifest(members)
        _result, job_id = self.submit(
            manifest,
            "--max-parallel",
            "2",
            "--gpu-budget",
            "64MiB",
        )
        assert job_id is not None
        worker = subprocess.Popen(
            [sys.executable, str(self.launcher), "worker", "--drain"],
            cwd=REPOSITORY,
            env=dict(os.environ),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            deadline = time.monotonic() + 8
            active_paths = [
                self.projects / "long-a" / "member.pid",
                self.projects / "long-b" / "member.pid",
            ]
            while time.monotonic() < deadline:
                if all(path.is_file() for path in active_paths):
                    break
                time.sleep(0.05)
            else:
                self.fail("two ensemble lanes did not start")

            cancel = self.run_cli("cancel", str(job_id))
            self.assertEqual(cancel.returncode, 0, cancel.stderr)
            stdout, stderr = worker.communicate(timeout=10)
            self.assertEqual(worker.returncode, 0, stdout + stderr)
            self.assertEqual(self.store.require(job_id)["state"], "canceled")
            self.assertFalse(
                (self.projects / "never-started" / "member.pid").exists(),
                "a queued member started after cancellation",
            )
            for path in active_paths:
                pid = int(path.read_text(encoding="utf-8"))
                deadline = time.monotonic() + 3
                while time.monotonic() < deadline and process_is_live(pid):
                    time.sleep(0.05)
                self.assertFalse(process_is_live(pid), pid)
            payload = self.status(job_id)
            self.assertTrue(
                all(
                    member["state"] in {"canceled", "interrupted"}
                    for member in payload["members"]
                ),
                payload,
            )
        finally:
            if worker.poll() is None:
                worker.kill()
                worker.wait()
            for name in ("long-a", "long-b"):
                path = self.projects / name / "member.pid"
                if not path.exists():
                    continue
                pid = int(path.read_text(encoding="utf-8"))
                if process_is_live(pid):
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass


if __name__ == "__main__":
    unittest.main()
