"""Record a Linux/WSL FLASH build preflight without launching FLASH.

Only explicitly supplied source roots are inspected. Licensed source is never
copied. The optional MPI/HDF5 smoke check builds in a temporary directory.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile


FORTRAN_SMOKE = """program parallel_hdf5_smoke
use mpi
use hdf5
implicit none
integer :: ierr, rank
integer(hid_t) :: access, file_id
call MPI_Init(ierr)
if (ierr /= MPI_SUCCESS) stop 1
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call h5open_f(ierr)
if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
call h5pcreate_f(H5P_FILE_ACCESS_F, access, ierr)
if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
call h5pset_fapl_mpio_f(access, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
call h5fcreate_f('parallel_smoke.h5', H5F_ACC_TRUNC_F, file_id, ierr, access_prp=access)
if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
call h5fclose_f(file_id, ierr)
if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
call h5pclose_f(access, ierr)
call h5close_f(ierr)
if (rank == 0) print *, 'parallel HDF5 collective file creation passed'
call MPI_Finalize(ierr)
end program
"""


def run_command(argv, *, cwd=None, timeout=30, env_updates=None):
    try:
        environment = dict(os.environ, **env_updates) if env_updates else None
        result = subprocess.run(argv, cwd=cwd, env=environment, capture_output=True, text=True,
                                timeout=timeout, check=False)
        return {"argv": argv, "exit_code": result.returncode,
                "stdout": result.stdout, "stderr": result.stderr,
                "environment_updates": env_updates or {}}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"argv": argv, "exit_code": None, "error": str(exc)}


def probe_linux(source_root=None, hdf5_prefix=None, compile_smoke=False):
    if platform.system() != "Linux":
        raise ValueError("run in Linux or select --wsl-distro")
    names = ("gcc", "gfortran", "mpicc", "mpif90", "mpiexec", "make", "cmake", "ninja")
    commands = {name: shutil.which(name) for name in names}
    prefix = Path(hdf5_prefix).expanduser().resolve() if hdf5_prefix else None
    wrapper = next((str(path) for name in ("h5pfc", "h5pfc.openmpi")
                    if prefix and (path := prefix / "bin" / name).is_file()), None)
    wrapper = wrapper or shutil.which("h5pfc") or shutil.which("h5pfc.openmpi")
    commands["h5pfc"] = wrapper
    report = {"schema": "mhdlab-flash-environment-v1",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "platform": platform.platform(), "commands": commands,
              "checks": [], "source": None,
              "flash_run_performed": False, "licensed_source_copied": False,
              "parallel_hdf5_smoke_passed": False}
    for name, args in (("gcc", ["--version"]), ("gfortran", ["--version"]),
                       ("mpicc", ["--showme:version"]), ("mpif90", ["--showme:command"]),
                       ("cmake", ["--version"])):
        if commands[name]:
            report["checks"].append(run_command([commands[name], *args]))
    if shutil.which("dpkg-query"):
        report["checks"].append(run_command([
            "dpkg-query", "-W", "-f=${binary:Package} ${Version}\n", "*hdf5*", "*openmpi*"]))
    if wrapper:
        report["checks"].append(run_command([wrapper, "-show"]))
        report["checks"].append(run_command([wrapper, "-showconfig"]))
    if source_root:
        root = Path(source_root).expanduser().resolve()
        setup = root / "setup"
        report["source"] = {"root": str(root), "exists": root.is_dir(),
            "setup_exists": setup.is_file(),
            "setup_sha256": hashlib.sha256(setup.read_bytes()).hexdigest() if setup.is_file() else None,
            "simulation_tree_exists": (root / "source" / "Simulation").is_dir(),
            "version": None, "version_verified": False}
        # The owner must supply version/build provenance; directory names do not prove it.
    if compile_smoke and wrapper and commands["mpiexec"]:
        with tempfile.TemporaryDirectory(prefix="power-flow-hdf5-") as directory:
            work = Path(directory)
            (work / "smoke.f90").write_text(FORTRAN_SMOKE, encoding="utf-8")
            compiled = run_command([wrapper, "smoke.f90", "-o", "smoke"], cwd=work)
            report["checks"].append(compiled)
            if compiled["exit_code"] == 0:
                environment = {"LD_LIBRARY_PATH": str(prefix / "lib") +
                               (":" + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")} if prefix else None
                executed = run_command([commands["mpiexec"], "-n", "2", "./smoke"], cwd=work,
                                       env_updates=environment)
                report["checks"].append(executed)
                report["parallel_hdf5_smoke_passed"] = executed["exit_code"] == 0
    missing = [name for name in ("gcc", "gfortran", "mpicc", "mpif90", "make") if not commands[name]]
    blockers = [f"missing tool: {name}" for name in missing]
    if not wrapper:
        blockers.append("parallel HDF5 Fortran wrapper not found")
    if not report["parallel_hdf5_smoke_passed"]:
        blockers.append("two-rank MPI/HDF5 Fortran smoke not passed")
    if report["source"] is None:
        blockers.append("licensed FLASH source root not supplied")
    elif not (report["source"]["setup_exists"] and report["source"]["simulation_tree_exists"]):
        blockers.append("supplied source lacks expected FLASH setup/Simulation layout")
    blockers.append("FLASH version, site configuration and physical benchmark remain unverified")
    report["blockers"] = blockers
    report["status"] = "preflight_only"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wsl-distro", help="Windows only: e.g. Ubuntu-24.04")
    parser.add_argument("--source-root", help="licensed source path in Linux coordinates")
    parser.add_argument("--hdf5-prefix", help="optional Linux prefix with bin/h5pfc")
    parser.add_argument("--compile-smoke", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("use a fresh output file; accepted preflights are preserved")
    if args.wsl_distro:
        if os.name != "nt":
            parser.error("--wsl-distro is Windows only")
        # --exec bypasses shell expansion of source paths and dpkg format strings.
        script = Path(__file__).read_text(encoding="utf-8")
        invocation = ["wsl.exe", "--distribution", args.wsl_distro, "--exec", "python3", "-c",
                      script.rsplit("\nif __name__ == \"__main__\":", 1)[0]
                      + "\nprint(json.dumps(probe_linux(*json.loads(sys.argv[1])), allow_nan=False))",
                      json.dumps([args.source_root, args.hdf5_prefix, args.compile_smoke])]
        outcome = run_command(invocation, timeout=110)
        if outcome["exit_code"] != 0:
            raise RuntimeError(outcome.get("error") or outcome.get("stderr"))
        report = json.loads(outcome["stdout"])
        report["wsl_distro"] = args.wsl_distro
    else:
        report = probe_linux(args.source_root, args.hdf5_prefix, args.compile_smoke)
    report["probe_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(args.output), "status": report["status"],
                      "blockers": report["blockers"]}, indent=2))


if __name__ == "__main__":
    main()
