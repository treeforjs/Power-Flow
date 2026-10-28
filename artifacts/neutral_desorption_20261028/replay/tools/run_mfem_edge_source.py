"""Run a current-constrained MFEM foil/return/vacuum model on bounded C.19 data."""
from __future__ import annotations

import argparse
import bisect
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


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def command(argv, output, name, *, cwd=None):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    result = subprocess.run(argv, cwd=cwd, env=env, text=True, capture_output=True)
    (output / (name + ".stdout.txt")).write_text(result.stdout, encoding="utf-8")
    (output / (name + ".stderr.txt")).write_text(result.stderr, encoding="utf-8")
    return {"argv": argv, "exit_code": result.returncode,
            "stdout_file": name + ".stdout.txt", "stderr_file": name + ".stderr.txt"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wsl-distro")
    parser.add_argument("--config", type=Path, default=ROOT / "configs/mfem_edge_current_source.json")
    parser.add_argument("--build-dir", type=Path, default=Path("/home/tsmit/.local/share/power-flow/mfem-edge-current-source"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.wsl_distro:
        if os.name != "nt":
            parser.error("--wsl-distro is Windows only")
        forwarded = ["--build-dir", args.build_dir.as_posix(), "--output-dir", args.output_dir.as_posix()]
        if "--config" in sys.argv:
            forwarded += ["--config", args.config.as_posix()]
        raise SystemExit(subprocess.run(["wsl.exe", "-d", args.wsl_distro, "--exec", "python3",
                                        "tools/run_mfem_edge_source.py", *forwarded]).returncode)
    if os.name == "nt":
        parser.error("run in WSL/Linux")
    config = json.loads(args.config.read_bytes())
    if config.get("schema") != "mhdlab-mfem-edge-current-source-v1" or config.get("physical_experiment_prediction") is not False:
        parser.error("only the explicit provisional case is admitted")
    if config.get("initial_temperature_k") != 300.0:
        parser.error("this provisional native case requires initial_temperature_k=300; absolute snapshots use that baseline")
    output = args.output_dir.resolve()
    if output.exists():
        parser.error("use a fresh output directory")
    trace_path, provenance_path = (ROOT / config[key] for key in ("current_trace_csv", "current_trace_provenance"))
    provenance = json.loads(provenance_path.read_bytes())
    expected_hashes = [v for node in provenance.values() if isinstance(node, dict)
                       for k, v in node.items() if k == "csv_sha256"]
    canonical_trace_hash = hashlib.sha256(trace_path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    if expected_hashes not in ([sha(trace_path)], [canonical_trace_hash]):
        raise ValueError("current CSV is not bound to its recovery provenance")
    with trace_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    trace_time = [float(row["time_s"]) for row in rows]
    trace_current = [float(row["current_A"]) for row in rows]
    start, stop = float(trace_time[0]), config["stop_time_s"]
    if not start < stop <= trace_time[-1] or any(b <= a for a, b in zip(trace_time, trace_time[1:])):
        raise ValueError("run must remain inside strictly ordered finite trace support")
    for stamp in config["required_snapshot_times_s"]:
        if not start <= stamp <= stop:
            raise ValueError("required snapshot is outside supported run")
    knots = sorted(set([stamp for stamp in trace_time if start <= stamp <= stop]
                       + config["required_snapshot_times_s"] + [start, stop]))
    grid = [start]
    for a, b in zip(knots, knots[1:]):
        count = max(1, math.ceil((b-a)/config["maximum_time_step_s"]))
        grid.extend(a+(b-a)*j/count for j in range(1, count+1))
        grid[-1] = b
    def interpolate(stamp):
        i = min(bisect.bisect_right(trace_time, stamp)-1, len(trace_time)-2)
        fraction = (stamp-trace_time[i])/(trace_time[i+1]-trace_time[i])
        return (1-fraction)*trace_current[i]+fraction*trace_current[i+1]
    values = [interpolate(stamp) for stamp in grid]
    if not all(math.isfinite(value) for value in values) or max(b-a for a, b in zip(grid, grid[1:])) > config["maximum_time_step_s"] * (1 + 1e-12):
        raise ValueError("nonfinite interpolation or excessive timestep")
    output.mkdir(parents=True)
    with (output / "schedule.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["time_s", "current_A"])
        writer.writerows(zip(map(repr, grid), map(repr, values)))
    record = {"schema": "mhdlab-mfem-edge-current-run-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "provisional", "passed": False, "physical_experiment_prediction": False,
        "configuration": config, "configuration_sha256": sha(args.config),
        "drive": {"csv": str(trace_path), "csv_sha256": sha(trace_path), "provenance_sha256": sha(provenance_path),
            "canonical_lf_csv_sha256": canonical_trace_hash, "newline_normalization": "CRLF to LF for provenance identity only; numeric values unchanged",
            "shot": 17133, "figure": "C.19", "raw_oscilloscope_data": False, "time_shift_s": 0,
            "support_s": [float(trace_time[0]), float(trace_time[-1])], "run_support_s": [start, stop],
            "interpolation": "bounded piecewise linear; original knots and required40/60ns included; no padding/clamping",
            "clock_scenario": config["clock_scenario"], "initial_state_scenario": config["initial_state_scenario"]},
        "formulation": "MFEM H1 A_z MQS with sigma0 vacuum, two independent conductor gradients constrained by a2x2 current Schur complement; separate conductor-only insulated thermal FE",
        "stock_joule_kernel_reused": False, "commands": []}
    write(output / "configuration.json", config)
    failure = None
    try:
        jobs = [("configure", ["cmake", "-S", str(ROOT / "native/mfem_edge_source"), "-B", str(args.build_dir),
                  "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release", "-DMFEM_DIR=/usr/local/lib/cmake/mfem"]),
                ("build", ["cmake", "--build", str(args.build_dir), "--parallel", "2"])]
        for name, argv in jobs:
            result = command(argv, output, name); record["commands"].append(result)
            if result["exit_code"]:
                raise RuntimeError(name + " failed")
        binary = args.build_dir / "power-flow-edge-source"
        options = {"width-cells": "width_cells", "thickness-cells": "thickness_cells", "outside-cells": "outside_cells",
            "gap-cells": "gap_cells", "width-m": "foil_width_m", "thickness-m": "foil_thickness_m", "gap-m": "gap_m",
            "extent-m": "exterior_half_extent_m", "length-m": "current_length_m", "sigma-s-m": "electrical_conductivity_s_m",
            "capacity-j-m3-k": "volumetric_heat_capacity_j_m3_k", "conductivity-w-m-k": "thermal_conductivity_w_m_k"}
        argv = [str(binary), "--schedule", str(output / "schedule.csv"), "--output-dir", str(output)]
        for flag, key in options.items():
            argv += ["--" + flag, str(config[key])]
        record["binary_sha256"] = sha(binary)
        result = command(argv, output, "native", cwd=output); record["commands"].append(result)
        if result["exit_code"]:
            raise RuntimeError("native current-constrained solve failed")
        with (output / "diagnostics.csv").open(newline="") as stream:
            diagnostics = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
        if len(diagnostics) != len(grid) or not all(math.isfinite(v) for row in diagnostics for v in row.values()):
            raise RuntimeError("missing/nonfinite native diagnostics")
        record["diagnostics"] = {"final": diagnostics[-1],
            "maximum_current_constraint_error_a": max(max(abs(r["foil_current_a"]-r["target_current_a"]), abs(r["return_current_a"]+r["target_current_a"])) for r in diagnostics),
            "maximum_em_balance_residual_j": max(abs(r["em_balance_residual_j"]) for r in diagnostics),
            "maximum_thermal_balance_residual_j": max(abs(r["thermal_balance_residual_j"]) for r in diagnostics),
            "maximum_edge_to_center_joule_ratio": max(r["edge_to_center_joule_ratio"] for r in diagnostics)}
        nx, width, thickness, gap, length = (config[k] for k in ("width_cells", "foil_width_m", "foil_thickness_m", "gap_m", "current_length_m"))
        for face, sign, coordinate in (("top", 1, gap/2+thickness), ("bottom", -1, gap/2)):
            sample = json.loads((ROOT / "configs/mfem_surface_example.json").read_bytes())
            sample["configuration_status"] = "configured"
            sample["source"].update(source_kind="mfem", heating_backend="mfem_joule", mfem_version="4.10", mfem_source_revision=REVISION,
                run_id="C19-current-constrained-MFEM-"+face, source_reference="source_run.json with current, energy, mesh and surface-quadrature records",
                unit_mapping_reference="SI MFEM A_z MQS current constraints and conductor-only heat FE; no numerical vacuum conductivity or thermal material",
                source_clock_reference=config["clock_scenario"])
            sample["source"]["solver_scope"].update(stationary_material=True, conductivity_temperature_feedback=False, hall_induction=False, material_motion=False,
                qualification_reference="Provisional constant-property current-constrained foil/return/vacuum MQS, not stock Joule or qualified ALEGRA reproduction")
            sample["temperature"].update(value_basis="temperature_increment_kelvin", initial_temperature_k=config["initial_temperature_k"],
                zero_convention_reference=config["initial_state_scenario"], location="surface_boundary", evaluation="one_sided_material_boundary_trace_patch_average", material_sample_depth_m=0, boundary_temperature_reference="")
            sample["face"].update(surface_region="upper foil "+face+" material face", normal_position_m=coordinate, outward_normal_sign=sign,
                material_boundary_attribute=1 if sign<0 else 3,
                material_selection_reference="Conductor-only H1 thermal FE material-side boundary quadrature; vacuum has no thermal unknowns",
                full_width_coverage_reference="Aligned native rectangular material-face quadrature checked against every dx; invariant current direction extruded once",
                current_direction_representation="resolved_patches")
            sample["time_s"] = grid; sample["x_edges_m"] = [-width/2+width*i/nx for i in range(nx+1)]; sample["z_edges_m"] = [0., length]
            sample["table"] = {"path": face+".csv", "sha256": sha(output/(face+".csv"))}
            write(output/(face+"_surface.json"), sample)
        units = {"time_s":"s", "x_left_m":"m", "x_right_m":"m", "depth_left_m":"m", "depth_right_m":"m",
                 "temperature_k":"K", "current_density_a_m2":"A m^-2", "joule_power_w_m3":"W m^-3"}
        write(output / "foil_snapshots.json", {"schema":"mhdlab-mfem-foil-snapshots-v1", "status":"provisional", "sampling":"cell_average",
            "units":units, "csv_sha256":sha(output/"foil_snapshots.csv"), "geometry":{"foil":"upper", "depth_reference":"inward_from_upper_gap_face",
            "foil_width_m":width, "foil_thickness_m":thickness, "gap_m":gap}, "drive":record["drive"],
            "model":record["formulation"], "material_provenance":config["material_provenance"]})
        record["passed"] = True
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        failure = str(exc); record["failure"] = failure
    record["owned_source_sha256"] = {name: sha(ROOT/name) for name in ("native/mfem_edge_source/CMakeLists.txt", "native/mfem_edge_source/edge_driver.cpp", "tools/run_mfem_edge_source.py")}
    record["outputs_sha256"] = {p.name:sha(p) for p in output.iterdir() if p.is_file()}
    write(output / "source_run.json", record)
    print(json.dumps({"report":str(output/"source_run.json"), "passed":record["passed"], "failure":failure}))
    if failure:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
