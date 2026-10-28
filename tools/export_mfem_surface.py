#!/usr/bin/env python3
"""Export explicit MFEM surface samples to the conservative planar-release format."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.meta_path[:] = [finder for finder in sys.meta_path
                    if type(finder).__module__ != "_editable_skbc_mhdlab"]
from mhdlab.mfem_surface import export_mfem_surface


def write_demo(output: Path) -> Path:
    """Create authored patch samples, never a native MFEM simulation claim."""
    output.mkdir(parents=True, exist_ok=False)
    config = json.loads((ROOT / "configs/mfem_surface_example.json").read_bytes())
    config["configuration_status"] = "manufactured"
    config["source"].update(source_kind="synthetic", mfem_version="", mfem_source_revision="",
        run_id="manufactured-mfem-surface-table", source_reference="Authored table; no native MFEM execution",
        unit_mapping_reference="Manufactured SI table fixture; no raw dimensionless values converted",
        source_clock_reference="Manufactured nonnegative simulation clock; no experimental registration")
    config["source"]["solver_scope"]["qualification_reference"] = "Manufactured interface fixture only"
    config["temperature"].update(value_basis="temperature_increment_kelvin", initial_temperature_k=300.,
        zero_convention_reference="Manufactured 300 K initial temperature plus authored delta-T; not a physical rod conversion")
    config["face"].update(surface_region="manufactured gap-facing face",
                           material_selection_reference="Manufactured material-side face selection",
                           full_width_coverage_reference="Manufactured complete x partition")
    table = output / "surface_samples.csv"
    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", "z_index", "x_index", "temperature_value"])
        for it, time in enumerate(config["time_s"]):
            for iz in range(len(config["z_edges_m"]) - 1):
                for ix in range(len(config["x_edges_m"]) - 1):
                    writer.writerow([time, iz, ix, it * (20 + 40 * ix + 80 * iz)])
    config["table"] = {"path": table.name, "sha256": hashlib.sha256(table.read_bytes()).hexdigest()}
    target = output / "configuration.json"
    target.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--config", type=Path)
    source.add_argument("--demo", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.demo:
            config = write_demo(args.output_dir / "manufactured_inputs")
            manifest = export_mfem_surface(config, args.output_dir / "export")
        else:
            manifest = export_mfem_surface(args.config, args.output_dir)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.exit(2, f"MFEM surface export: {exc}\n")
    print(json.dumps({"source_kind": manifest["source_attributes"]["source_kind"],
                      "surface_history": manifest["surface_history"],
                      "native_mfem_gridfunction_parsed": False,
                      "new_heating_simulation_performed": False}, indent=2))


if __name__ == "__main__":
    main()
