from __future__ import annotations

import json
import io
import os
import pwd
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

from simqueue import cli, exclusive
from simqueue.config import format_cpu_list
from simqueue.store import QueueStore
from simqueue.worker import Worker


class ExclusiveFixture(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.projects = self.base / "projects"
        self.projects.mkdir()
        allowed = tuple(sorted(os.sched_getaffinity(0)))
        self.pcores = allowed[: min(2, len(allowed))]
        self.ecores = allowed[min(2, len(allowed)) : min(4, len(allowed))]
        if not self.ecores:
            self.ecores = self.pcores
        self.environment = {
            "SIMQUEUE_STATE_DIR": str(self.base / "state"),
            "SIMQUEUE_PROJECT_ROOT": str(self.projects),
            "SIMQUEUE_PCORES": format_cpu_list(self.pcores),
            "SIMQUEUE_CANCEL_GRACE_SECONDS": "0.1",
        }
        self.environment_patch = patch.dict(
            os.environ, self.environment, clear=False
        )
        self.environment_patch.start()
        os.environ.pop("SIMQUEUE_NO_AUTOSTART", None)
        self.marker = self.base / "exclusive.json"
        self.config = {
            "version": 1,
            "user": pwd.getpwuid(os.getuid()).pw_name,
            "uid": os.getuid(),
            "project_root": str(self.projects),
            "state_dir": str(self.base / "state"),
            "pcores": format_cpu_list(self.pcores),
            "ecores": format_cpu_list(self.ecores),
            "runner_unit": "simqueue-runner.service",
            "helper": "/usr/local/libexec/simqueue-shield",
        }
        self.marker.write_text(json.dumps(self.config), encoding="utf-8")

    def tearDown(self) -> None:
        self.environment_patch.stop()
        self.temporary.cleanup()


class ExclusiveConfigTests(ExclusiveFixture):
    def test_loads_valid_marker_and_rejects_malformed_marker(self) -> None:
        with patch("simqueue.exclusive.CONFIG_PATH", self.marker):
            self.assertEqual(exclusive.load_config(), self.config)

        self.marker.write_text("[]", encoding="utf-8")
        with patch("simqueue.exclusive.CONFIG_PATH", self.marker):
            with self.assertRaisesRegex(RuntimeError, "expected a JSON object"):
                exclusive.load_config()

    def test_environment_cannot_redirect_or_disable_exclusive_marker(self) -> None:
        with (
            patch("simqueue.exclusive.CONFIG_PATH", self.marker),
            patch.dict(
                os.environ,
                {
                    "SIMQUEUE_INTERNAL_TESTING": "1",
                    "SIMQUEUE_INTERNAL_TEST_EXCLUSIVE_CONFIG": str(
                        self.base / "missing"
                    ),
                },
                clear=False,
            ),
        ):
            self.assertEqual(exclusive.load_config(), self.config)

    def test_missing_marker_cannot_fall_back_with_install_footprint(self) -> None:
        missing_marker = self.base / "missing-production-marker"
        residual_unit = self.base / "simqueue-runner.service"
        residual_unit.touch()
        with (
            patch("simqueue.exclusive.CONFIG_PATH", missing_marker),
            patch("simqueue.exclusive.PRODUCTION_CONFIG_PATH", missing_marker),
            patch("simqueue.exclusive.INSTALL_FOOTPRINT", (residual_unit,)),
        ):
            with self.assertRaisesRegex(RuntimeError, "installation is incomplete"):
                exclusive.load_config()

    def test_cli_validates_marker_before_creating_redirected_state(self) -> None:
        redirected = self.base / "redirected-state"
        with (
            patch("simqueue.exclusive.CONFIG_PATH", self.marker),
            patch.dict(
                os.environ, {"SIMQUEUE_STATE_DIR": str(redirected)}, clear=False
            ),
            redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(cli.main(["list"]), 2)
        self.assertFalse(redirected.exists())

    def test_worker_refuses_unshielded_manual_execution(self) -> None:
        store = QueueStore()
        job = store.submit(
            name="must-stay-queued",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.pcores),
        )
        with (
            patch("simqueue.worker.load_config", return_value=self.config),
            patch(
                "simqueue.worker.verify_worker_context",
                side_effect=RuntimeError("outside protected cgroup"),
            ),
        ):
            self.assertEqual(Worker(drain=True).run(), 2)
        self.assertEqual(store.require(int(job["id"]))["state"], "queued")

    def test_worker_rechecks_shield_before_claiming_each_job(self) -> None:
        store = QueueStore()
        job = store.submit(
            name="must-remain-queued-after-shield-loss",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.pcores),
        )
        with (
            patch("simqueue.worker.load_config", return_value=self.config),
            patch(
                "simqueue.worker.verify_worker_context",
                side_effect=(None, RuntimeError("partition changed")),
            ),
            patch("simqueue.worker._enable_child_subreaper", return_value=True),
        ):
            self.assertEqual(Worker(drain=True).run(), 2)
        self.assertEqual(store.require(int(job["id"]))["state"], "queued")

    def test_cli_uses_installed_pcores_while_outside_partition(self) -> None:
        with (
            patch("simqueue.cli.load_config", return_value=self.config),
            patch(
                "simqueue.cli.detect_pcores",
                side_effect=RuntimeError("current shell is restricted to E-cores"),
            ),
        ):
            self.assertEqual(cli._cli_pcores(), self.pcores)

    def test_exclusive_nudge_never_falls_back_to_detached_worker(self) -> None:
        store = QueueStore()
        with (
            patch("simqueue.cli.load_config", return_value=self.config),
            patch(
                "simqueue.cli._start_exclusive_runner",
                side_effect=cli.WorkerStartError("shield failed"),
            ),
            patch("simqueue.cli.subprocess.Popen") as popen,
        ):
            with self.assertRaisesRegex(cli.WorkerStartError, "shield failed"):
                cli._nudge_worker(store)
        popen.assert_not_called()

    def test_exclusive_nudge_retries_the_drain_exit_race(self) -> None:
        store = QueueStore()
        job = store.submit(
            name="race",
            cwd=self.projects,
            argv=["/bin/true"],
            environment=dict(os.environ),
            env_override_keys=[],
            requested_cores=len(self.pcores),
        )
        starts = 0

        def start(_unit: str) -> None:
            nonlocal starts
            starts += 1
            if starts == 2:
                store.claim_next(self.pcores)

        def active_state(*_args: object, **_kwargs: object) -> str:
            return "active" if starts >= 2 else "inactive"

        with (
            patch("simqueue.cli._start_system_unit", side_effect=start),
            patch("simqueue.cli._systemctl_property", side_effect=active_state),
            patch("simqueue.cli.time.sleep"),
        ):
            mode = cli._start_exclusive_runner(store, self.config)
        self.assertEqual(mode, "exclusive-systemd")
        self.assertEqual(starts, 2)
        self.assertEqual(store.require(int(job["id"]))["state"], "starting")


if __name__ == "__main__":
    unittest.main()
