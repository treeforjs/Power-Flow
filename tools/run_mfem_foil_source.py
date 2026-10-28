"""Run an explicit manufactured SI foil through MFEM's unchanged Joule kernel.

Native material-side DG trace tables are exported for the existing MFEM surface
adapter. This is not an experimentally qualified electrode calculation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REVISION = "d964264cdb9a13e94a201b6c236c7721e0c8765f"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def execute(argv, output, name, *, cwd=None):
    result = subprocess.run(argv, cwd=cwd, text=True, capture_output=True, check=False)
    (output / (name + ".stdout.txt")).write_text(result.stdout, encoding="utf-8")
    (output / (name + ".stderr.txt")).write_text(result.stderr, encoding="utf-8")
    return {"argv": argv, "exit_code": result.returncode,
            "stdout_file": name + ".stdout.txt", "stderr_file": name + ".stderr.txt"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wsl-distro")
    parser.add_argument("--source-root", type=Path, default=Path("/home/tsmit/src/mfem-4.10"))
    parser.add_argument("--build-dir", type=Path, default=Path("/home/tsmit/.local/share/power-flow/mfem-foil-gif-source"))
    parser.add_argument("--config", type=Path, default=ROOT / "configs/mfem_foil_gif_source.json")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.wsl_distro:
        if os.name != "nt":
            parser.error("--wsl-distro is for Windows")
        forwarded = ["--source-root", args.source_root.as_posix(), "--build-dir", args.build_dir.as_posix(),
                     "--output-dir", args.output_dir.as_posix()]
        if "--config" in sys.argv:
            forwarded += ["--config", args.config.as_posix()]
        raise SystemExit(subprocess.run(["wsl.exe", "-d", args.wsl_distro, "--exec", "python3",
                                        "tools/run_mfem_foil_source.py", *forwarded]).returncode)
    if os.name == "nt":
        parser.error("run in Linux or use --wsl-distro")
    output, source, build = args.output_dir.resolve(), args.source_root.resolve(), args.build_dir.resolve()
    if output.exists():
        parser.error("use a fresh output directory")
    config = json.loads(args.config.read_bytes())
    if config.get("schema") != "mhdlab-mfem-foil-gif-source-v1" or config.get("physical_experiment_prediction") is not False:
        parser.error("only the explicitly manufactured/provisional source is admitted")
    output.mkdir(parents=True)
    record = {"schema": "mhdlab-mfem-foil-gif-source-run-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
              "status": "manufactured_provisional_native_mfem", "physical_experiment_prediction": False,
              "configuration": config, "configuration_sha256": sha(args.config), "commands": [], "passed": False,
              "mfem_version": "4.10", "source_release_reference": "https://github.com/mfem/mfem/tree/v4.10",
              "mfem_source_revision": REVISION, "binary_source_build_match_independently_established": False,
              "source_files_sha256": {name: sha(source / "miniapps/electromagnetics" / name)
                                       for name in ("joule_solver.cpp", "joule_solver.hpp")}}
    (output / "configuration.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    failure = None
    try:
        commands = [("configure", ["cmake", "-S", str(ROOT / "native/mfem_foil_source"), "-B", str(build), "-G", "Ninja",
                                    "-DCMAKE_BUILD_TYPE=Release", "-DMFEM_DIR=/usr/local/lib/cmake/mfem",
                                    "-DMFEM_SOURCE_ROOT=" + str(source)]),
                    ("build", ["cmake", "--build", str(build), "--parallel", "2"])]
        for name, argv in commands:
            result = execute(argv, output, name)
            record["commands"].append(result)
            if result["exit_code"]:
                raise RuntimeError(name + " failed; inspect saved logs")
        executable = build / "power-flow-foil-joule"
        record["executable_sha256"] = sha(executable)
        options = {"width-cells": "width_cells", "thickness-cells": "thickness_cells", "current-cells": "current_cells",
                   "steps": "steps", "width-m": "width_m", "thickness-m": "thickness_m", "length-m": "current_length_m",
                   "final-time-s": "final_time_s", "sigma-s-m": "electrical_conductivity_s_m",
                   "capacity-j-m3-k": "volumetric_heat_capacity_j_m3_k", "conductivity-w-m-k": "thermal_conductivity_w_m_k",
                   "voltage-v": "voltage_difference_v", "rise-time-s": "voltage_rise_time_s"}
        argv = ["mpiexec", "-n", str(config["mpi_ranks"]), str(executable)]
        for option, key in options.items():
            argv.extend(["--" + option, str(config[key])])
        argv.extend(["--output-dir", str(output)])
        result = execute(argv, output, "native", cwd=output)
        record["commands"].append(result)
        if result["exit_code"]:
            raise RuntimeError("native Joule foil run failed; inspect saved logs")
        with (output / "diagnostics.csv").open(newline="") as stream:
            diagnostics = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
        if len(diagnostics) != config["steps"] + 1 or not all(math.isfinite(v) for row in diagnostics for v in row.values()):
            raise RuntimeError("native diagnostics have missing or nonfinite snapshots")
        expected_times = [row["time_s"] for row in diagnostics]
        if any(b <= a for a, b in zip(expected_times, expected_times[1:])):
            raise RuntimeError("native time support is not increasing")
        minimum = min(row["minimum_delta_temperature_k"] for row in diagnostics) + config["initial_temperature_k"]
        if minimum <= 0 or diagnostics[-1]["direct_electric_power_w"] <= 0:
            raise RuntimeError("native absolute temperature or heating failed its finite positive guard")
        record["diagnostics"] = {
            "initial_temperature_k": config["initial_temperature_k"], "minimum_absolute_temperature_k": minimum,
            "maximum_absolute_temperature_k": max(row["maximum_delta_temperature_k"] for row in diagnostics) + config["initial_temperature_k"],
            "final": diagnostics[-1],
            "maximum_current_end_mismatch_a": max(abs(row["current_lower_a"]-row["current_upper_a"]) for row in diagnostics),
            "maximum_projected_vs_direct_power_relative_difference": max(
                abs(row["projected_joule_power_w"]-row["direct_electric_power_w"])/row["direct_electric_power_w"]
                for row in diagnostics if row["direct_electric_power_w"] > 0),
            "maximum_thermal_source_energy_residual_j": max(
                abs(row["thermal_energy_j"]-row["cumulative_projected_joule_energy_j"]) for row in diagnostics),
            "physical_energy_convergence_qualified": False}
        width, thickness, length = config["width_m"], config["thickness_m"], config["current_length_m"]
        nx, nz = config["width_cells"], config["current_cells"]
        references = {}
        for face, sign, attribute in (("top", 1, 4), ("bottom", -1, 3)):
            table = output / (face + ".csv")
            with table.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            if len(rows) != len(expected_times)*nx*nz or any(not math.isfinite(float(row["temperature_value"])) for row in rows):
                raise RuntimeError("missing/nonfinite face samples")
            sample = json.loads((ROOT / "configs/mfem_surface_example.json").read_bytes())
            sample["configuration_status"] = "configured"
            sample["source"].update(source_kind="mfem", heating_backend="mfem_joule", mfem_version="4.10", mfem_source_revision=REVISION,
                run_id="manufactured-native-foil-" + face, source_reference="source_run.json plus native driver/source hashes and diagnostics.csv",
                unit_mapping_reference="Manufactured SI constant-coefficient Joule equations; lengths m, time s, sigma S/m, capacity J/(m3 K), conductivity W/(m K)",
                source_clock_reference="Native simulation time starting0s; no experimental trigger registration")
            sample["source"]["solver_scope"].update(stationary_material=True, conductivity_temperature_feedback=False,
                hall_induction=False, material_motion=False,
                qualification_reference="Manufactured native MFEM foil execution and one-sided boundary quadrature; physical geometry/material/drive convergence not established")
            sample["temperature"].update(value_basis="temperature_increment_kelvin", initial_temperature_k=config["initial_temperature_k"],
                zero_convention_reference="Unchanged linear Joule thermal operator initialized at deltaT0, with insulated boundaries; absoluteT is declared ambient plus deltaT once",
                location="surface_boundary", evaluation="one_sided_material_boundary_trace_patch_average", material_sample_depth_m=0,
                boundary_temperature_reference="")
            sample["face"].update(surface_region="manufactured foil " + face + " material face", normal_position_m=sign*thickness/2,
                outward_normal_sign=sign, material_boundary_attribute=attribute,
                material_selection_reference="Native material-side GetBdrFaceTransformations Elem1 trace evaluated with boundary face quadrature",
                full_width_coverage_reference="Native aligned Cartesian x/z face quadrature areas checked against every dx*dz patch",
                current_direction_representation="resolved_patches")
            sample["time_s"] = expected_times
            sample["x_edges_m"] = [-width/2 + i*width/nx for i in range(nx+1)]
            sample["z_edges_m"] = [i*length/nz for i in range(nz+1)]
            sample["table"] = {"path": table.name, "sha256": sha(table)}
            target = output / (face + "_surface.json")
            target.write_text(json.dumps(sample, indent=2) + "\n", encoding="utf-8")
            references[face] = {"configuration": target.name, "configuration_sha256": sha(target), "table_sha256": sha(table)}
        record["surface_samples"] = references
        record["passed"] = True
    except (RuntimeError, OSError, KeyError, ValueError) as exc:
        failure = str(exc); record["failure"] = failure
    record["owned_source_sha256"] = {name: sha(ROOT / name) for name in
        ("native/mfem_foil_source/CMakeLists.txt", "native/mfem_foil_source/foil_driver.cpp", "tools/run_mfem_foil_source.py")}
    record["outputs_sha256"] = {path.name: sha(path) for path in output.iterdir() if path.is_file()}
    (output / "source_run.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(output / "source_run.json"), "passed": record["passed"], "failure": failure}))
    if failure:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
