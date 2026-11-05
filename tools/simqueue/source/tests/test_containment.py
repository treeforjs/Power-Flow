"""SimQueue 0.7.2: MPI containment.

Open MPI 5's ``prun`` (and ``prterun --dvm search``) attaches to a PRRTE DVM
it finds through the PMIx rendezvous files under $TMPDIR (/tmp when unset);
the DVM then forks the ranks in its own cgroup, outside the job's unit.  Each
job now gets its own rendezvous directory; the worker reports processes
carrying a job's launch token outside its cgroup and MPI daemons outside
every job.  The two-rank reproduction needs an Open MPI 5 prefix with prte,
prun and mpiexec (SIMQ_TEST_OMPI_PREFIX, else the MOOSE conda environment).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from node_fixture import PYTHON_FROZEN, NodeFixture

from simqueue.containment import (
    ContainmentConfig,
    EscapeScanner,
    JobScope,
    environment_with_containment,
    is_launcher_identity,
    job_tmpdir,
    make_job_tmpdir,
    parse_containment_config,
    read_escaped,
    remove_job_tmpdir,
    write_escaped,
)

DEFAULT_OMPI = "/home/treefor/projects/runtimes/rz-spheral-20261029/envs/moose-openmpi"


def ompi_prefix() -> Path | None:
    prefix = Path(os.environ.get("SIMQ_TEST_OMPI_PREFIX") or DEFAULT_OMPI)
    if all((prefix / "bin" / name).exists() for name in ("prte", "prun", "mpiexec")):
        return prefix
    return None


def own_cgroup(pid: int | str = "self") -> str:
    for line in Path(f"/proc/{pid}/cgroup").read_text().splitlines():
        if line.startswith("0::"):
            return line[3:]
    return ""


class EnvironmentTests(unittest.TestCase):
    def test_launcher_identity_variables(self) -> None:
        for name in ("PMIX_SERVER_URI41", "PMIX_NAMESPACE", "PMIX_RANK", "PRTE_LAUNCHED", "PMI_PORT", "PMI_FD",
                     "OMPI_COMM_WORLD_RANK", "OMPI_UNIVERSE_SIZE", "OMPI_MCA_orte_hnp_uri", "PMIX_SYSTEM_TMPDIR"):
            self.assertTrue(is_launcher_identity(name), name)
        for name in ("PMIX_MCA_gds", "PRTE_MCA_plm", "OMPI_MCA_pml", "OMPI_MCA_btl", "MPIR_CVAR_CH4_NETMOD",
                     "HYDRA_LAUNCHER", "PATH", "TMPDIR", "SIMQ_JOB_ID"):
            self.assertFalse(is_launcher_identity(name), name)

    def test_inherited_identity_dropped_explicit_kept_rendezvous_set(self) -> None:
        inherited = {"PATH": "/usr/bin", "PMIX_SERVER_URI41": "outside", "PMIX_MCA_gds": "hash",
                     "TMPDIR": "/tmp", "PRTE_LAUNCHED": "1", "PMIX_NAMESPACE": "explicit-ns"}
        result = environment_with_containment(inherited, {"PMIX_NAMESPACE"}, "/tmp/simq.7.abc")
        self.assertNotIn("PMIX_SERVER_URI41", result)
        self.assertNotIn("PRTE_LAUNCHED", result)
        self.assertEqual(result["PMIX_MCA_gds"], "hash")
        self.assertEqual(result["PMIX_NAMESPACE"], "explicit-ns")
        for name in ("TMPDIR", "PMIX_SYSTEM_TMPDIR", "PMIX_SERVER_TMPDIR", "PRTE_MCA_prte_tmpdir_base",
                     "SIMQ_TMPDIR"):
            self.assertEqual(result[name], "/tmp/simq.7.abc", name)
        # -e TMPDIR=... wins; SIMQ_TMPDIR always names the job's directory
        explicit = environment_with_containment({"TMPDIR": "/scratch/x"}, {"TMPDIR"}, "/tmp/simq.7.abc")
        self.assertEqual(explicit["TMPDIR"], "/scratch/x")
        self.assertEqual(explicit["PMIX_SYSTEM_TMPDIR"], "/tmp/simq.7.abc")
        self.assertEqual(explicit["SIMQ_TMPDIR"], "/tmp/simq.7.abc")
        # containment off: only the identity filter
        self.assertNotIn("TMPDIR", environment_with_containment({"PMIX_RANK": "0"}, set(), None))

    def test_job_environment_uses_containment(self) -> None:
        from simqueue.nodeworker import job_environment

        with tempfile.TemporaryDirectory() as state:
            saved = dict(os.environ)
            os.environ["SIMQUEUE_STATE_DIR"] = state
            try:
                job = {"id": 3, "launch_token": "0123456789abcdef",
                       "env_json": json.dumps({"PATH": "/usr/bin", "PMIX_SERVER_URI41": "x", "SLURM_JOB_ID": "9"}),
                       "env_override_keys": "[]"}
                environment = job_environment(job, {"SIMQ_JOB_ID": "3"}, "/tmp/simq.3.0123456789ab")
            finally:
                os.environ.clear()
                os.environ.update(saved)
        self.assertNotIn("PMIX_SERVER_URI41", environment)
        self.assertNotIn("SLURM_JOB_ID", environment)
        self.assertEqual(environment["TMPDIR"], "/tmp/simq.3.0123456789ab")
        self.assertEqual(environment["SIMQ_JOB_ID"], "3")

    def test_config_and_paths(self) -> None:
        self.assertEqual(parse_containment_config(None), ContainmentConfig())
        config = parse_containment_config({"tmp_base": "/tmp/", "scan_seconds": 0})
        self.assertEqual((config.tmp_base, config.scan_seconds), ("/tmp", 0))
        for bad in ({"tmp_base": "relative"}, {"scan_seconds": 2}, {"enabled": "yes"}, {"other": 1}):
            with self.assertRaises(ValueError):
                parse_containment_config(bad)
        self.assertEqual(job_tmpdir(ContainmentConfig(), 12, "ab-cd0123456789ff"), "/tmp/simq.12.abcd01234567")


class DirectoryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_private_directory_made_reused_and_refused(self) -> None:
        path = self.base / "simq.1.abc"
        make_job_tmpdir(str(path))
        self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o700)
        make_job_tmpdir(str(path))  # a restarted wrapper reuses it
        path.chmod(0o755)
        with self.assertRaises(PermissionError):
            make_job_tmpdir(str(path))
        link = self.base / "simq.2.abc"
        link.symlink_to(path)
        with self.assertRaises(PermissionError):
            make_job_tmpdir(str(link))

    def test_removal_never_follows_symlinks(self) -> None:
        outside = self.base / "outside"
        outside.mkdir()
        (outside / "keep").write_text("x")
        path = self.base / "simq.3.abc"
        make_job_tmpdir(str(path))
        (path / "data").write_bytes(b"0" * 10000)
        (path / "link").symlink_to(outside)
        self.assertGreater(remove_job_tmpdir(str(path)), 0)
        self.assertFalse(path.exists())
        self.assertTrue((outside / "keep").exists())
        other = self.base / "not-a-job-dir"
        other.mkdir()
        remove_job_tmpdir(str(other))
        self.assertTrue(other.exists())


class ScannerTests(unittest.TestCase):
    """The scan against real processes of this test (no units)."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.children: list[subprocess.Popen] = []

    def tearDown(self) -> None:
        for child in self.children:
            child.kill()
            child.wait()
        self.temporary.cleanup()

    def spawn(self, argv: list[str], token: str | None = None) -> subprocess.Popen:
        env = {"PATH": "/usr/bin:/bin"}
        if token:
            env["SIMQ_LAUNCH_TOKEN"] = token
        child = subprocess.Popen(argv, env=env)
        self.children.append(child)
        return child

    def test_token_outside_its_cgroup_is_reported_inside_is_not(self) -> None:
        child = self.spawn(["/bin/sleep", "30"], token="tok-escaped")
        here = own_cgroup()
        time.sleep(0.2)
        scanner = EscapeScanner()
        elsewhere = JobScope(job_id=5, token="tok-escaped", cgroup="/user.slice/some-other.service")
        result = scanner.scan([elsewhere])
        pids = [item.pid for item in result.escaped.get(5, [])]
        self.assertEqual(pids, [child.pid])
        self.assertEqual(result.escaped[5][0].cgroup, here)
        self.assertEqual([item.pid for item in result.new], [child.pid])
        again = scanner.scan([elsewhere])
        self.assertEqual(again.new, [])  # reported once as new
        inside = JobScope(job_id=5, token="tok-escaped", cgroup=here)
        self.assertEqual(EscapeScanner().scan([inside]).escaped, {})
        ended = EscapeScanner().scan([], ended=[inside])
        self.assertTrue(ended.escaped[5][0].after_end)
        # escaped.json merges and keeps after_end
        write_escaped(self.base, result.escaped[5])
        write_escaped(self.base, ended.escaped[5])
        entries = read_escaped(self.base)
        self.assertEqual(len(entries), 1)
        self.assertTrue(entries[0]["after_end"])

    def test_outside_mpi_daemon_is_reported_with_children(self) -> None:
        daemon_bin = self.base / "prte"
        shutil.copy("/bin/sh", daemon_bin)
        daemon = self.spawn([str(daemon_bin), "-c", "sleep 30 & sleep 30"])
        time.sleep(0.3)
        result = EscapeScanner().scan([])
        found = [d for d in result.daemons if d["pid"] == daemon.pid]
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0]["comm"], "prte")
        self.assertGreaterEqual(len(found[0]["children"]), 1)
        inside = JobScope(job_id=1, token="t", cgroup=own_cgroup())
        self.assertFalse([d for d in EscapeScanner().scan([inside]).daemons if d["pid"] == daemon.pid])

    def test_scan_cost(self) -> None:
        scanner = EscapeScanner()
        scopes = [JobScope(job_id=1, token="none-such", cgroup="/x")]
        first = scanner.scan(scopes)
        later = [scanner.scan(scopes).seconds for _ in range(5)]
        receipts = os.environ.get("SIMQ_TEST_RECEIPTS")
        if receipts:
            Path(receipts, "containment-scan-cost.json").write_text(json.dumps({
                "processes": first.processes, "first_scan_cpu_s": first.seconds,
                "cached_scan_cpu_s": later}, indent=1))
        self.assertLess(first.seconds, 0.5)
        self.assertLess(min(later), 0.1)


class ContainmentUnitTests(NodeFixture, unittest.TestCase):
    """Real job units of a sandbox worker."""

    min_cpus = 1

    def config_document(self):
        document = super().config_document()
        document["containment"] = {"scan_seconds": 5}
        return document

    def test_job_gets_a_private_rendezvous_directory_removed_at_its_end(self) -> None:
        self.start_worker()
        script = self.projects / "env.sh"
        script.write_text(
            "#!/bin/sh\n"
            "echo \"TMPDIR=$TMPDIR\"\n"
            "echo \"PMIX_SYSTEM_TMPDIR=$PMIX_SYSTEM_TMPDIR\"\n"
            "echo \"PMIX_SERVER_TMPDIR=$PMIX_SERVER_TMPDIR\"\n"
            "echo \"PRTE_MCA_prte_tmpdir_base=$PRTE_MCA_prte_tmpdir_base\"\n"
            "echo \"PMIX_SERVER_URI41=${PMIX_SERVER_URI41:-unset}\"\n"
            "echo \"PMIX_MCA_gds=${PMIX_MCA_gds:-unset}\"\n"
            "stat -c 'MODE=%a OWNER=%u' \"$TMPDIR\"\n"
            "echo data > \"$TMPDIR/file\"\n")
        script.chmod(0o755)
        env = dict(self.environment, PMIX_SERVER_URI41="tcp4://127.0.0.1:1", PMIX_MCA_gds="hash")
        job_id = self.submit("-p", "pa", "-c", "1", "--export", "PMIX_SERVER_URI41,PMIX_MCA_gds",
                             "--", str(script), env=env)
        self.wait_until(lambda: self.state_of(job_id) == "succeeded", timeout=60, message="succeeded")
        text = self.log_text(job_id)
        lines = dict(line.split("=", 1) for line in text.splitlines()
                     if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", line))
        tmpdir = lines["TMPDIR"]
        job = self.job(job_id)
        self.assertEqual(tmpdir, job_tmpdir(ContainmentConfig(), job_id, str(job["launch_token"])))
        for name in ("PMIX_SYSTEM_TMPDIR", "PMIX_SERVER_TMPDIR", "PRTE_MCA_prte_tmpdir_base"):
            self.assertEqual(lines[name], tmpdir, name)
        self.assertEqual(lines["PMIX_SERVER_URI41"], "unset")
        self.assertEqual(lines["PMIX_MCA_gds"], "hash")
        self.assertIn(f"MODE=700 OWNER={os.getuid()}", text)
        self.wait_until(lambda: not Path(tmpdir).exists(), timeout=20, message="temporary directory removed")
        show = self.cli("show", str(job_id)).stdout
        self.assertIn(f"Temporary dir", show)

    def test_process_moved_out_of_the_job_is_reported(self) -> None:
        if shutil.which("systemd-run") is None:
            self.skipTest("systemd-run not found")
        self.start_worker()
        scope = f"{self.prefix}-escape"
        script = self.projects / "escape.sh"
        script.write_text(
            "#!/bin/sh\n"
            f"systemd-run --user --scope -q --unit={scope} /bin/sleep 60 &\n"
            "sleep 9\n")
        script.chmod(0o755)
        job_id = self.submit("-p", "pa", "-c", "1", "--", str(script))
        try:
            self.wait_until(lambda: self.state_of(job_id) == "running", timeout=60, message="running")
            job_dir = str(self.job(job_id)["job_dir"])
            entries = self.wait_until(lambda: read_escaped(job_dir), timeout=30, message="escape reported")
            self.assertEqual(entries[0]["comm"], "sleep")
            self.assertTrue(entries[0]["cgroup"].endswith(f"/{scope}.scope"), entries[0]["cgroup"])
            self.assertFalse(entries[0]["after_end"])
            self.assertIn("ESCAPED", self.cli("show", str(job_id)).stdout)
            self.wait_until(lambda: self.state_of(job_id) == "succeeded", timeout=60, message="succeeded")
            # still alive after the job: reported as left behind
            self.wait_until(lambda: read_escaped(job_dir) and read_escaped(job_dir)[0]["after_end"], timeout=20,
                            message="left behind")
            self.assertIn("carries the job's launch token", self.worker_log())
        finally:
            subprocess.run(["systemctl", "--user", "stop", f"{scope}.scope"], stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, check=False)


class MpiReproductionTests(NodeFixture, unittest.TestCase):
    """Two ranks: prun escapes into an outside DVM without the job's own
    rendezvous directory, and cannot reach it with it; mpiexec stays in the
    job."""

    min_cpus = 1

    def config_document(self):
        document = super().config_document()
        document["containment"] = {"scan_seconds": 5}
        return document

    def setUp(self) -> None:
        super().setUp()
        self.prefix_dir = ompi_prefix()
        if self.prefix_dir is None:
            self.skipTest("no Open MPI 5 prefix with prte, prun and mpiexec (SIMQ_TEST_OMPI_PREFIX)")
        self.dvm_pid: int | None = None

    def tearDown(self) -> None:
        if self.dvm_pid:
            try:
                os.kill(self.dvm_pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            deadline = time.monotonic() + 10
            while time.monotonic() < deadline and Path(f"/proc/{self.dvm_pid}").exists():
                time.sleep(0.1)
            if Path(f"/proc/{self.dvm_pid}").exists():
                os.kill(self.dvm_pid, signal.SIGKILL)
        super().tearDown()

    def start_outside_dvm(self) -> str:
        """A DVM where a stray one would be: rendezvous under /tmp."""
        env = {k: v for k, v in os.environ.items()
               if k not in ("TMPDIR", "PMIX_SYSTEM_TMPDIR", "PMIX_SERVER_TMPDIR", "PRTE_MCA_prte_tmpdir_base")}
        env["PATH"] = f"{self.prefix_dir}/bin:/usr/bin:/bin"
        pid_file = self.base / "dvm.pid"
        subprocess.run([str(self.prefix_dir / "bin" / "prte"), "--daemonize", "--report-pid", str(pid_file),
                        "--report-uri", str(self.base / "dvm.uri")], env=env, cwd=self.base,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30, check=False)
        self.wait_until(lambda: pid_file.exists() and pid_file.read_text().strip(), timeout=20, message="DVM pid")
        self.dvm_pid = int("".join(ch for ch in pid_file.read_text() if ch.isdigit()))
        self.wait_until(lambda: Path(f"/proc/{self.dvm_pid}").exists(), timeout=10, message="DVM running")
        return own_cgroup(self.dvm_pid)

    def ranks(self, job_id: int) -> list[str]:
        return [line.split("rank cgroup: ", 1)[1] for line in self.log_text(job_id).splitlines()
                if line.startswith("rank cgroup: ")]

    def run_job(self, command: str) -> int:
        script = self.projects / f"mpi-{len(list(self.projects.glob('mpi-*')))}.sh"
        script.write_text("#!/bin/sh\nulimit -c 0\n" + command + "\nexit 0\n")
        script.chmod(0o755)
        job_id = self.submit("-p", "pa", "-c", "1", "-e", f"PATH={self.prefix_dir}/bin:/usr/bin:/bin",
                             "--", str(script))
        self.wait_until(lambda: self.state_of(job_id) in ("succeeded", "failed"), timeout=120, message="finished")
        return job_id

    def test_two_rank_escape_reproduced_and_contained(self) -> None:
        dvm_cgroup = self.start_outside_dvm()
        self.start_worker()
        rank_script = self.projects / "rank.sh"  # prun re-splits a quoted sh -c command; a script is safe
        rank_script.write_text("#!/bin/sh\necho \"rank cgroup: $(cut -d: -f3 /proc/self/cgroup)\"\n")
        rank_script.chmod(0o755)
        rank = str(rank_script)
        unprotected = ("env -u TMPDIR -u PMIX_SYSTEM_TMPDIR -u PMIX_SERVER_TMPDIR -u PRTE_MCA_prte_tmpdir_base "
                       f"timeout 60 prun --map-by :OVERSUBSCRIBE -n 2 {rank} 2>&1")
        # 1. as 0.7.1 ran jobs (no rendezvous directory): the ranks run in the DVM's cgroup
        escaped = self.run_job(unprotected)
        unit = str(self.job(escaped)["unit_name"])
        ranks = self.ranks(escaped)
        self.assertEqual(len(ranks), 2, self.log_text(escaped))
        for cgroup in ranks:
            self.assertNotIn(unit, cgroup)
            self.assertEqual(cgroup, dvm_cgroup)
        # the worker saw the DVM outside every job
        document = self.wait_until(
            lambda: (json.loads((self.state / "containment.json").read_text())
                     if (self.state / "containment.json").exists() else None), timeout=30, message="containment.json")
        self.assertIn(self.dvm_pid, [d["pid"] for d in document["outside_mpi_daemons"]])
        self.assertIn("OUTSIDE MPI DAEMON", self.cli("partitions").stdout)
        # 2. 0.7.2: prun finds no DVM; nothing runs outside
        contained = self.run_job(f"timeout 60 prun --map-by :OVERSUBSCRIBE -n 2 {rank} 2>&1")
        self.assertEqual(self.ranks(contained), [])
        self.assertIn("no DVM", self.log_text(contained))
        # 3. mpiexec: both ranks in the job's unit, with the job's TMPDIR
        own = self.run_job(f"timeout 60 mpiexec --map-by :OVERSUBSCRIBE --bind-to none -n 2 {rank} 2>&1")
        unit = str(self.job(own)["unit_name"])
        ranks = self.ranks(own)
        self.assertEqual(len(ranks), 2, self.log_text(own))
        for cgroup in ranks:
            self.assertTrue(cgroup.endswith(f"/{unit}"), cgroup)


if __name__ == "__main__":
    unittest.main()
