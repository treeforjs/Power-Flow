from __future__ import annotations

import importlib.machinery
import importlib.util
import os
import sys
import tempfile
import unittest
from contextlib import nullcontext
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
HELPER = REPOSITORY / "privileged" / "simqueue-shield"
LOADER = importlib.machinery.SourceFileLoader("simqueue_shield_helper", str(HELPER))
SPEC = importlib.util.spec_from_loader(LOADER.name, LOADER)
assert SPEC is not None
SHIELD = importlib.util.module_from_spec(SPEC)
sys.modules[LOADER.name] = SHIELD
LOADER.exec_module(SHIELD)


class ShieldTree:
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.cgroup = self.base / "sys" / "fs" / "cgroup"
        self.shield = self.cgroup / "simqueue.slice"
        self.system = self.cgroup / "system.slice"
        self.user = self.cgroup / "user.slice"
        self.runner = self.shield / "simqueue-runner.service"
        self.topology = self.base / "sys" / "devices"
        self.proc = self.base / "proc"
        for directory in (
            self.cgroup,
            self.shield,
            self.system,
            self.user,
            self.runner,
            self.topology / "cpu_core",
            self.topology / "cpu_atom",
            self.topology / "system" / "cpu",
            self.proc,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        self._put(self.cgroup / "cgroup.controllers", "cpuset cpu")
        self._put(self.cgroup / "cgroup.subtree_control", "cpuset cpu")
        self._put(self.cgroup / "cgroup.procs", "")
        self._put(self.shield / "cgroup.procs", "")
        self._put(self.runner / "cgroup.procs", "4242")
        self._put(self.shield / "cpuset.cpus.effective", "0-7")
        self._put(self.shield / "cpuset.cpus.exclusive.effective", "")
        self._put(self.shield / "cpuset.cpus.partition", "member")
        self._put(self.shield / "cpuset.mems.effective", "0")
        for branch in (self.system, self.user):
            self._put(branch / "cgroup.procs", "")
            self._put(branch / "cpuset.cpus.effective", "0-23")

        self._put(self.topology / "cpu_core" / "cpus", "0-7")
        self._put(self.topology / "cpu_atom" / "cpus", "8-23")
        self._put(self.topology / "system" / "cpu" / "online", "0-23")
        self._put(self.topology / "system" / "cpu" / "present", "0-23")

        self.paths = SHIELD.ShieldPaths(
            config_file=self.base / "etc" / "simqueue" / "shield.conf",
            cgroup_root=self.cgroup,
            shield_cgroup=self.shield,
            pcore_topology=self.topology / "cpu_core" / "cpus",
            ecore_topology=self.topology / "cpu_atom" / "cpus",
            cpu_online=self.topology / "system" / "cpu" / "online",
            cpu_present=self.topology / "system" / "cpu" / "present",
            proc_root=self.proc,
            lock_file=self.base / "run" / "lock" / "simqueue-shield.lock",
        )
        self.config = SHIELD.ShieldConfig(
            pcores=frozenset(range(8)), ecores=frozenset(range(8, 24))
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _put(path: Path, value: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value + ("\n" if value else ""), encoding="ascii")

    def kernel_write(self, path: Path, value: str) -> None:
        self._put(path, value)
        if path.name != "cpuset.cpus.partition":
            return
        if value == "root":
            self._put(self.shield / "cpuset.cpus.exclusive.effective", "0-7")
            outside = "8-23"
        elif value == "member":
            self._put(self.shield / "cpuset.cpus.exclusive.effective", "")
            outside = "0-23"
        else:
            return
        self._put(self.system / "cpuset.cpus.effective", outside)
        self._put(self.user / "cpuset.cpus.effective", outside)

    def controller(self, **overrides: object):
        arguments = {
            "write_control": self.kernel_write,
            "lock_factory": lambda _path: nullcontext(),
            "get_euid": lambda: 0,
            "get_pid": lambda: 4242,
        }
        arguments.update(overrides)
        return SHIELD.ShieldController(self.config, self.paths, **arguments)


class ConfigTests(unittest.TestCase):
    def test_parser_accepts_only_the_two_cpu_masks(self) -> None:
        config = SHIELD.parse_config("# fixed host topology\nPCORES=0-7\nECORES=8-23\n")
        self.assertEqual(SHIELD.format_cpu_list(config.pcores), "0-7")
        self.assertEqual(SHIELD.format_cpu_list(config.ecores), "8-23")

    def test_parser_rejects_shell_syntax_unknown_keys_and_overlap(self) -> None:
        invalid = (
            "PCORES=$(id)\nECORES=8-23\n",
            "PCORES=0-7\nECORES=8-23\nPATH=/tmp\n",
            "PCORES=0-8\nECORES=8-23\n",
            "PCORES=0,0\nECORES=8-23\n",
            "PCORES=0-999999999\nECORES=8-23\n",
        )
        for text in invalid:
            with self.subTest(text=text), self.assertRaises(SHIELD.ShieldError):
                SHIELD.parse_config(text)

    def test_production_paths_ignore_environment_overrides(self) -> None:
        with patch.dict(
            os.environ,
            {
                "SIMQUEUE_CGROUP_ROOT": "/tmp/attacker",
                "SIMQUEUE_SHIELD_CONFIG": "/tmp/attacker.conf",
            },
        ):
            paths = SHIELD.ShieldPaths.production()
        self.assertEqual(paths.config_file, Path("/etc/simqueue/shield.conf"))
        self.assertEqual(paths.shield_cgroup, Path("/sys/fs/cgroup/simqueue.slice"))
        self.assertEqual(paths.lock_file, Path("/run/simqueue/shield.lock"))


class ControllerTests(ShieldTree, unittest.TestCase):
    def test_acquire_and_repeated_acquire_are_valid_and_idempotent(self) -> None:
        controller = self.controller()
        first = controller.acquire()
        second = controller.acquire()

        self.assertTrue(first["valid"])
        self.assertTrue(second["valid"])
        self.assertEqual(second["state"], "root")
        self.assertEqual(second["configured_cpus"], "0-7")
        self.assertEqual(second["effective_cpus"], "0-7")
        self.assertEqual(second["exclusive_cpus"], "0-7")
        self.assertEqual(second["outsider_violations"], {})
        self.assertEqual(
            (self.shield / "cpuset.cpus.partition").read_text().strip(), "root"
        )

    def test_acquire_fails_closed_if_an_outsider_keeps_pcores(self) -> None:
        def broken_kernel_write(path: Path, value: str) -> None:
            self._put(path, value)
            if path.name == "cpuset.cpus.partition" and value == "root":
                self._put(
                    self.shield / "cpuset.cpus.exclusive.effective", "0-7"
                )

        controller = self.controller(write_control=broken_kernel_write)
        with patch.object(SHIELD, "VERIFY_TIMEOUT_SECONDS", 0):
            with self.assertRaisesRegex(SHIELD.ShieldError, "did not reach"):
                controller.acquire()
        self.assertEqual(
            (self.shield / "cpuset.cpus.partition").read_text().strip(), "root"
        )

    def test_release_refuses_surviving_process_and_leaves_root_partition(self) -> None:
        controller = self.controller()
        controller.acquire()
        self._put(self.runner / "cgroup.procs", "4242\n9001")

        with self.assertRaisesRegex(SHIELD.ShieldError, "9001"):
            controller.release()
        self.assertEqual(
            (self.shield / "cpuset.cpus.partition").read_text().strip(), "root"
        )

    def test_release_and_repeated_release_are_valid_and_idempotent(self) -> None:
        controller = self.controller()
        controller.acquire()
        first = controller.release()
        second = controller.release()

        self.assertTrue(first["valid"])
        self.assertTrue(second["valid"])
        self.assertEqual(second["state"], "member")
        self.assertFalse(second["populated"])
        self.assertEqual(second["survivor_pids"], [])
        self.assertEqual(second["exclusive_cpus"], "")
        self.assertEqual(second["outsider_violations"], {})

    def test_release_refuses_a_populated_member(self) -> None:
        self._put(self.runner / "cgroup.procs", "4242\n9001")
        with self.assertRaisesRegex(SHIELD.ShieldError, "9001"):
            self.controller().release()
        self.assertEqual(
            (self.shield / "cpuset.cpus.partition").read_text().strip(), "member"
        )

    def test_absent_status_and_release_are_safe(self) -> None:
        for path in sorted(self.shield.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            else:
                path.rmdir()
        self.shield.rmdir()
        controller = self.controller()

        status = controller.inspect()
        released = controller.release()
        self.assertEqual(status["state"], "absent")
        self.assertTrue(status["valid"])
        self.assertEqual(released["state"], "absent")
        self.assertTrue(released["valid"])

    def test_absent_is_safe_even_if_the_configured_topology_is_now_unavailable(self) -> None:
        for path in sorted(self.shield.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            else:
                path.rmdir()
        self.shield.rmdir()
        self._put(self.topology / "cpu_core" / "cpus", "0-6")

        controller = self.controller()
        self.assertTrue(controller.inspect()["valid"])
        self.assertTrue(controller.release()["valid"])

    def test_unprivileged_mutations_are_rejected_before_writing(self) -> None:
        controller = self.controller(get_euid=lambda: 1000)
        with self.assertRaisesRegex(SHIELD.ShieldError, "effective UID 0"):
            controller.acquire()
        with self.assertRaisesRegex(SHIELD.ShieldError, "effective UID 0"):
            controller.release()
        self.assertEqual(
            (self.shield / "cpuset.cpus.partition").read_text().strip(), "member"
        )

    def test_topology_mismatch_is_rejected_before_partition_write(self) -> None:
        self._put(self.topology / "cpu_core" / "cpus", "0-6")
        with self.assertRaisesRegex(SHIELD.ShieldError, "topology mismatch"):
            self.controller().acquire()
        self.assertEqual(
            (self.shield / "cpuset.cpus.partition").read_text().strip(), "member"
        )

    def test_status_reports_a_populated_valid_root(self) -> None:
        controller = self.controller()
        controller.acquire()
        self._put(self.runner / "cgroup.procs", "")
        nested = self.runner / "delegated.scope"
        nested.mkdir()
        self._put(nested / "cgroup.procs", "5150")
        status = controller.inspect()

        self.assertTrue(status["valid"])
        self.assertTrue(status["populated"])
        self.assertEqual(status["survivor_pids"], [5150])
        self.assertEqual(status["state"], "root")

    def test_root_rejects_direct_and_sibling_service_processes(self) -> None:
        sibling = self.shield / "unrelated.service"
        sibling.mkdir()
        self._put(sibling / "cgroup.procs", "9002")
        self._put(self.shield / "cgroup.procs", "4242\n9001")
        controller = self.controller()

        with patch.object(SHIELD, "VERIFY_TIMEOUT_SECONDS", 0):
            with self.assertRaisesRegex(
                SHIELD.ShieldError, "outside simqueue-runner.service"
            ):
                controller.acquire()

        status = controller.inspect(ignore_pids=frozenset({4242}))
        self.assertFalse(status["valid"])
        self.assertEqual(status["state"], "root")
        self.assertEqual(status["survivor_pids"], [9001, 9002])
        error = "\n".join(status["errors"])
        self.assertIn("9001 at <slice>", error)
        self.assertIn("9002 at unrelated.service", error)
        self.assertNotIn("4242", error)

    def test_status_rejects_a_populated_unshielded_member(self) -> None:
        self._put(self.runner / "cgroup.procs", "5150")
        status = self.controller().inspect()

        self.assertFalse(status["valid"])
        self.assertEqual(status["state"], "member")
        self.assertIn("unshielded simqueue.slice is populated", status["errors"])

    def test_status_rejects_isolated_partition(self) -> None:
        self._put(self.shield / "cpuset.cpus.partition", "isolated")
        self._put(self.shield / "cpuset.cpus.exclusive.effective", "0-7")
        status = self.controller().inspect(ignore_pids=frozenset({4242}))
        self.assertFalse(status["valid"])
        self.assertIn("disables scheduler load balancing", status["errors"][0])


if __name__ == "__main__":
    unittest.main()
