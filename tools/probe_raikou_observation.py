"""Verify native population-output readback through the common observation model.

The two decaying levels and their wavelengths are manufactured interface inputs.
This does not admit a surface source into a plasma receiver or fit a shot.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Load this exact source file: a globally installed editable checkout can place
# its package finder ahead of sys.path and redirect an ordinary mhdlab import.
PROJECTOR_SOURCE = ROOT / "mhdlab" / "resolved_emission_observation.py"
_spec = importlib.util.spec_from_file_location("_power_flow_probe_resolved_emission", PROJECTOR_SOURCE)
if _spec is None or _spec.loader is None:
    raise ImportError(f"cannot load observation projector: {PROJECTOR_SOURCE}")
_projector = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _projector
_spec.loader.exec_module(_projector)
ResolvedTransition = _projector.ResolvedTransition
project_resolved_emission = _projector.project_resolved_emission


def _json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run_probe(executable: Path, output: Path, source_revision: str | None = None,
              saved_native_dir: Path | None = None, clock_origin_s: float = 0.0) -> dict:
    executable, output = executable.expanduser().resolve(), output.expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"native executable not found: {executable}")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError(f"use a fresh output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if not np.isfinite(clock_origin_s):
        raise ValueError("clock origin must be finite")
    binary_hash = hashlib.sha256(executable.read_bytes()).hexdigest()
    saved_report = None
    if saved_native_dir is not None:
        saved_native_dir = saved_native_dir.expanduser().resolve()
        saved_report = json.loads((saved_native_dir / "report.json").read_text(encoding="utf-8"))
        if saved_report.get("passed") is not True or saved_report.get("executable_sha256") != binary_hash:
            raise ValueError("saved native control must have passed with this exact executable")

    def saved_bytes(name):
        data = (saved_native_dir / name).read_bytes()
        if saved_report.get("artifact_sha256", {}).get(name) != hashlib.sha256(data).hexdigest():
            raise ValueError(f"saved native input/output hash does not match: {name}")
        return data

    def execute(operation, request):
        _json(output / f"{operation}_request.json", request)
        if saved_report is not None:
            if json.loads(saved_bytes(f"{operation}_request.json")) != request:
                raise ValueError(f"saved request differs for {operation}")
            stdout = saved_bytes(f"{operation}_stdout.json")
            stderr = saved_bytes(f"{operation}_stderr.txt")
            (output / f"{operation}_stdout.json").write_bytes(stdout)
            (output / f"{operation}_stderr.txt").write_bytes(stderr)
            value = json.loads(stdout)
            if not isinstance(value, dict):
                raise ValueError("saved response must be a JSON object")
            return value
        result = subprocess.run(
            [str(executable), operation, "--input", "-"], input=json.dumps(request, allow_nan=False),
            capture_output=True, text=True, timeout=60, check=False,
        )
        (output / f"{operation}_stdout.json").write_text(result.stdout, encoding="utf-8")
        (output / f"{operation}_stderr.txt").write_text(result.stderr, encoding="utf-8")
        if result.returncode != 0:
            raise RuntimeError(f"Raikou {operation} failed ({result.returncode}): {result.stderr.strip()}")
        value = json.loads(result.stdout)
        if not isinstance(value, dict):
            raise ValueError("Raikou response must be a JSON object")
        return value

    capabilities = execute("capabilities", {})
    if capabilities.get("cpu") is not True or "chemistry.advance" not in capabilities.get("operations", []):
        raise ValueError("selected runtime does not advertise CPU chemistry.advance")
    initial = {"X0": 0.0, "X1": 1e12, "X2": 2e12}
    request = {
        "mechanism": {
            "id": "manufactured-population-observation-control", "status": "manufactured",
            "species": [{"id": name, "elements": {"H": 2}} for name in initial],
            "reactions": [
                {"id": "decay1", "reactants": {"X1": 1}, "products": {"X0": 1},
                 "rate": {"type": "constant", "A": 2.0}, "source": "manufactured control"},
                {"id": "decay2", "reactants": {"X2": 1}, "products": {"X0": 1},
                 "rate": {"type": "constant", "A": 3.0}, "source": "manufactured control"},
            ],
        },
        "state": {"densities_m3": initial}, "dt_s": 0.1,
        "execution": {"device": "cpu", "mpi": False, "threads": 1},
    }
    native = execute("chemistry.advance", request)
    if native.get("success") is not True or not np.isclose(native.get("advanced_time_s"), 0.1, rtol=1e-12, atol=0):
        raise ValueError("native chemistry did not complete the supplied interval")
    final = native["state"]["densities_m3"]
    if set(final) != set(initial) or any(not np.isfinite(value) or value < 0 for value in final.values()):
        raise ValueError("native state readback does not match the explicit state contract")
    transitions = [ResolvedTransition("synthetic_1", "X1", 130.0, 2.0),
                   ResolvedTransition("synthetic_2", "X2", 150.0, 3.0)]
    populations = {name: np.asarray([[initial[name]], [final[name]]]) for name in ("X1", "X2")}
    projection = project_resolved_emission(
        time_s=np.asarray([clock_origin_s, clock_origin_s + 0.1]), upper_state_density_m3=populations,
        path_weights_m=np.asarray([0.004]), transitions=transitions,
        bin_edges_nm=np.linspace(120.0, 160.0, 401), gate_start_s=clock_origin_s, gate_end_s=clock_origin_s + 0.1,
        fwhm_nm=0.2, time_registration=f"manufactured t_photon = t_native + {clock_origin_s:.17g} seconds",
        population_provenance="exact X1/X2 keys from the saved native chemistry.advance response",
        geometry_provenance="manufactured single cell with a 4 mm collected ray length",
        status="manufactured", optical_model="optically_thin",
    )
    _json(output / "observation.json", projection)
    conserved = bool(np.isclose(sum(final.values()), sum(initial.values()), rtol=1e-10, atol=0))
    decayed = final["X0"] > 0 and final["X1"] < initial["X1"] and final["X2"] < initial["X2"]
    expected_lines = {
        transition.transition_id: (0.5 * (initial[transition.upper_state_id] + final[transition.upper_state_id])
                                  * 0.004 * transition.einstein_a_per_s * 6.62607015e-34 * 299792458.0
                                  / (4 * np.pi * transition.vacuum_wavelength_nm * 1e-9))
        for transition in transitions
    }
    observed_lines = {line["transition_id"]: line["integrated_radiance_w_m2_sr"] for line in projection["lines"]}
    projected = (set(observed_lines) == set(expected_lines)
                 and all(np.isclose(observed_lines[key], value, rtol=2e-14, atol=0)
                         for key, value in expected_lines.items()))
    captured = bool(np.isclose(sum(projection["bin_integrated_radiance_w_m2_sr"]),
                               sum(expected_lines.values()), rtol=2e-14, atol=0))
    report = {
        "schema": "mhdlab-raikou-observation-probe-v1", "status": "manufactured",
        "repository": "https://github.com/treeforjs/raikou",
        "source_revision_reference": source_revision,
        "binary_source_revision_match_established": False,
        "native_version": capabilities.get("version"),
        "executable": str(executable),
        "executable_sha256": binary_hash,
        "native_evolution_calls": 0 if saved_report is not None else 1,
        "saved_native_report_sha256": (hashlib.sha256((saved_native_dir / "report.json").read_bytes()).hexdigest()
                                       if saved_report is not None else None),
        "manufactured_photon_clock_origin_s": clock_origin_s,
        "projector_source": str(PROJECTOR_SOURCE),
        "projector_source_sha256": hashlib.sha256(PROJECTOR_SOURCE.read_bytes()).hexdigest(),
        "probe_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "power_flow_base_revision": "18d155b919e70b0b20165e22ce1eb4b7f48ac5a9",
        "checks": {"native_success": True, "state_ids_match": True, "particle_inventory_conserved": conserved,
                   "nontrivial_decay": bool(decayed), "analytic_interpolant_projection": bool(projected),
                   "full_window_line_energy_captured": captured},
        "passed": bool(conserved and decayed and projected and captured),
        "surface_source_receiver_accepted": False,
        "absolute_experimental_calibration_corrected": False,
        "neutral_density_or_temperature_inferred": False,
        "actual_FLASH_run": False,
        "interpretation": "Native population-output readback and observation projection only; no excitation data, physical radiative energy ledger, source admission or experimental validation.",
    }
    report["artifact_sha256"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(output.iterdir()) if path.is_file()
    }
    _json(output / "report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--source-revision", help="Separately recorded source reference; not a binary build attestation")
    parser.add_argument("--saved-native-dir", type=Path, help="Hash-checked completed probe to reproject without evolving native states again")
    parser.add_argument("--clock-origin-s", type=float, default=0.0, help="Explicit manufactured photon-clock origin")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    parser.add_argument("--output", type=Path, default=Path("runs") / f"raikou_observation_{stamp}")
    args = parser.parse_args()
    try:
        report = run_probe(args.executable, args.output, args.source_revision,
                           args.saved_native_dir, args.clock_origin_s)
    except (OSError, ValueError, RuntimeError, KeyError, subprocess.SubprocessError) as error:
        print(f"Raikou observation probe failed: {error}", file=sys.stderr)
        return 2
    print(f"{'PASS' if report['passed'] else 'FAIL'}: manufactured native population-output observation readback")
    print(f"Report: {args.output.resolve() / 'report.json'}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
