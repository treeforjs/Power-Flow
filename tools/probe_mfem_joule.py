"""Build unmodified upstream Joule and inspect a short nondimensional native run.

The separate installed MFEM library is consumed, never rebuilt or modified.
Outputs are numerical interface evidence, not a qualified physical foil case.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = (
    "miniapps/electromagnetics/joule.cpp", "miniapps/electromagnetics/joule_solver.cpp",
    "miniapps/electromagnetics/joule_solver.hpp", "miniapps/electromagnetics/cylinder-hex.mesh",
    "miniapps/common/mesh_extras.cpp", "miniapps/common/fem_extras.cpp",
    "miniapps/common/particles_extras.cpp", "miniapps/common/pfem_extras.cpp",
    "miniapps/common/dist_solver.cpp")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(argv, output, name, *, cwd=None):
    result = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, check=False)
    (output / (name + ".stdout.txt")).write_text(result.stdout, encoding="utf-8")
    (output / (name + ".stderr.txt")).write_text(result.stderr, encoding="utf-8")
    return {"argv": argv, "exit_code": result.returncode,
            "stdout_file": name + ".stdout.txt", "stderr_file": name + ".stderr.txt"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wsl-distro", help="Windows only: launch this script inside the selected WSL distro")
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--mfem-dir", type=Path, default=Path("/usr/local/lib/cmake/mfem"))
    parser.add_argument("--config", type=Path, default=ROOT / "configs/mfem_joule_smoke.json")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.wsl_distro:
        if os.name != "nt": parser.error("--wsl-distro is Windows only")
        # WSL inherits this checkout as cwd; all explicitly supplied backend paths use Linux coordinates.
        options = ["--source-root", args.source_root.as_posix(), "--build-dir", args.build_dir.as_posix(),
                   "--mfem-dir", args.mfem_dir.as_posix(), "--output-dir", args.output_dir.as_posix()]
        # Default config is selected inside Linux; explicit config must be relative to this checkout.
        if "--config" in sys.argv: options.extend(["--config", args.config.as_posix()])
        process = subprocess.run(["wsl.exe", "--distribution", args.wsl_distro, "--exec", "python3",
                                  "tools/probe_mfem_joule.py", *options])
        raise SystemExit(process.returncode)
    if platform.system() != "Linux": parser.error("run in Linux or use --wsl-distro")
    source, build, output = args.source_root.resolve(), args.build_dir.resolve(), args.output_dir.resolve()
    if output.exists(): parser.error("use a new output directory")
    config = json.loads(args.config.read_bytes())
    if config.get("schema") != "mhdlab-mfem-joule-smoke-v1" or config.get("physical_foil_case") is not False:
        parser.error("only an explicitly nonphysical smoke configuration is supported")
    source_hashes = {name: sha(source / name) for name in SOURCE_FILES}
    output.mkdir(parents=True)
    report = {"schema": "mhdlab-mfem-joule-probe-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
              "status": "upstream_dimensionless_demo", "physical_foil_case": False,
              "mfem_source_root": str(source), "source_sha256": source_hashes,
              "config": config, "config_sha256": sha(args.config), "commands": [],
              "source_modified": False, "mfem_library_rebuilt": False, "passed": False}
    (output / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    failure = None
    try:
        configure = ["cmake", "-S", str(ROOT / "native/mfem_joule"), "-B", str(build), "-G", "Ninja",
                     "-DCMAKE_BUILD_TYPE=Release", "-DMFEM_DIR=" + str(args.mfem_dir), "-DMFEM_SOURCE_ROOT=" + str(source)]
        for name, argv in (("configure", configure), ("build", ["cmake", "--build", str(build), "--parallel", "2"])):
            check = run(argv, output, name)
            report["commands"].append(check)
            if check["exit_code"]: raise RuntimeError(name + " failed; inspect saved stderr/stdout")
        executable, inspector = build / "power-flow-joule", build / "power-flow-joule-inspect"
        report["executable_sha256"] = sha(executable)
        report["inspector_sha256"] = sha(inspector)
        argv = ["mpiexec", "-n", str(config["mpi_ranks"]), str(executable),
                "--problem", config["problem"], "--mesh", str(source / config["mesh_relative_path"]),
                "--order", str(config["order"]), "--refine-serial", str(config["serial_refinements"]),
                "--refine-parallel", str(config["parallel_refinements"]), "--ode-solver", str(config["ode_solver"]),
                "--time-step", str(config["time_step_native"]), "--t-final", str(config["final_time_native"]),
                "--no-visualization", "--visit", "--outputfilename", config["collection_name"],
                "--visualization-steps", "1"]
        check = run(argv, output, "joule", cwd=output)
        report["commands"].append(check)
        if check["exit_code"]: raise RuntimeError("Joule failed; inspect saved native logs")
        cycles = sorted(int(path.name.rsplit("_", 1)[1]) for path in output.glob(config["collection_name"] + "_*") if path.is_dir())
        if not cycles: raise RuntimeError("native VisIt output cycles missing")
        report["inspections"] = []
        for cycle in (cycles[0], cycles[-1]):
            argv = ["mpiexec", "-n", str(config["mpi_ranks"]), str(inspector), "--name", config["collection_name"],
                    "--directory", str(output) + "/", "--cycle", str(cycle)]
            check = run(argv, output, "inspect_" + str(cycle), cwd=output)
            report["commands"].append(check)
            if check["exit_code"]: raise RuntimeError("native snapshot inspection failed")
            value = json.loads((output / check["stdout_file"]).read_text())
            if value.get("all_fields_finite") is not True: raise RuntimeError("nonfinite native fields")
            report["inspections"].append(value)
        final = report["inspections"][-1]
        report["passed"] = final["temperature_volume_integral_native"] > 0 and final["joule_volume_integral_native"] > 0
        if not report["passed"]: raise RuntimeError("short run did not produce positive integrated heat and temperature")
        if any(sha(source / name) != digest for name, digest in source_hashes.items()): raise RuntimeError("upstream source changed during verification")
    except (OSError, ValueError, RuntimeError) as exc:
        failure = str(exc)
        report["failure"] = failure
    report["implementation_sha256"] = {name: sha(ROOT / name) for name in (
        "tools/probe_mfem_joule.py", "native/mfem_joule/CMakeLists.txt", "native/mfem_joule/inspect_snapshot.cpp")}
    report["output_sha256"] = {str(path.relative_to(output)): sha(path) for path in sorted(output.rglob("*")) if path.is_file()}
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(output / "report.json"), "passed": report["passed"], "failure": failure}))
    if failure: raise SystemExit(1)


if __name__ == "__main__":
    main()
