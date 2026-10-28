#!/usr/bin/env python3
"""Compare declared MFEM cell-average temperatures with recovered ALEGRA displays.

This is a literal-coordinate, conditional model comparison. It performs no
parameter fit, physical clock registration, curve smoothing or extrapolation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "mhdlab-mfem-foil-snapshots-v1"
COLUMNS = ("time_s", "x_left_m", "x_right_m", "depth_left_m", "depth_right_m",
           "temperature_k", "current_density_a_m2", "joule_power_w_m3")
UNITS = {"time_s": "s", "x_left_m": "m", "x_right_m": "m", "depth_left_m": "m",
         "depth_right_m": "m", "temperature_k": "K", "current_density_a_m2": "A m^-2",
         "joule_power_w_m3": "W m^-3"}
DATA = ROOT / "data/dissertation"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError(f"nonfinite JSON value: {value}")
    return json.loads(Path(path).read_bytes(), object_pairs_hook=pairs, parse_constant=nonfinite)


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def load_snapshots(path, metadata_path=None):
    path = Path(path)
    sidecar = Path(metadata_path) if metadata_path else path.with_suffix(".json")
    metadata = read_json(sidecar)
    if metadata.get("schema") != SCHEMA or metadata.get("sampling") != "cell_average":
        raise ValueError("MFEM snapshots require the explicit cell-average schema")
    if metadata.get("status") not in {"provisional", "manufactured", "owner_supplied"}:
        raise ValueError("MFEM snapshot status must be explicit")
    if metadata.get("units") != UNITS:
        raise ValueError("MFEM snapshot SI units must match the canonical columns")
    if metadata.get("csv_sha256") != sha(path):
        raise ValueError("MFEM snapshot CSV SHA256 mismatch")
    raw = rows(path)
    if not raw or tuple(raw[0]) != COLUMNS:
        raise ValueError("MFEM snapshot CSV requires the canonical ordered columns")
    table = np.asarray([[float(row[key]) for key in COLUMNS] for row in raw])
    if not np.isfinite(table).all() or np.any(table[:, 5] <= 0) or np.any(table[:, 7] < 0):
        raise ValueError("MFEM fields must be finite, with positive Kelvin temperature and nonnegative Joule power")
    if np.any(table[:, 2] <= table[:, 1]) or np.any(table[:, 4] <= table[:, 3]):
        raise ValueError("MFEM cells must have positive width and depth")
    groups = {}
    common = None
    for stamp in np.unique(table[:, 0]):
        cells = table[table[:, 0] == stamp]
        x = np.unique(cells[:, 1:3])
        depth = np.unique(cells[:, 3:5])
        nx, nd = len(x)-1, len(depth)-1
        if len(cells) != nx*nd or depth[0] != 0:
            raise ValueError("each MFEM snapshot must cover one complete rectangular foil from depth zero")
        fields = np.full((3, nd, nx), np.nan)
        for cell in cells:
            ix, iz = int(np.searchsorted(x, cell[1])), int(np.searchsorted(depth, cell[3]))
            if ix >= nx or iz >= nd or cell[2] != x[ix+1] or cell[4] != depth[iz+1]:
                raise ValueError("MFEM snapshots cannot contain nested, overlapping or noncontiguous cells")
            if np.isfinite(fields[0, iz, ix]):
                raise ValueError("duplicate MFEM cell")
            fields[:, iz, ix] = cell[5:8]
        if not np.isfinite(fields).all():
            raise ValueError("missing MFEM foil cell")
        if common is not None and (not np.array_equal(x, common[0]) or not np.array_equal(depth, common[1])):
            raise ValueError("comparison currently requires a common spatial grid across snapshots")
        common = (x, depth)
        groups[float(stamp)] = {"x_edges_m": x, "depth_edges_m": depth,
            "temperature_k": fields[0], "current_density_a_m2": fields[1], "joule_power_w_m3": fields[2]}
    geometry = metadata.get("geometry", {})
    if geometry.get("foil") != "upper" or geometry.get("depth_reference") != "inward_from_upper_gap_face":
        raise ValueError("comparison requires explicit upper-foil inward-depth association")
    sample = next(iter(groups.values()))
    for actual, declared in ((sample["x_edges_m"][-1]-sample["x_edges_m"][0], geometry.get("foil_width_m")),
                             (sample["depth_edges_m"][-1], geometry.get("foil_thickness_m"))):
        if not isinstance(declared, (int, float)) or not math.isclose(actual, declared, rel_tol=1e-12, abs_tol=1e-16):
            raise ValueError("MFEM grid conflicts with its declared foil geometry")
    return groups, metadata


def exact_snapshot(groups, requested):
    stamp = min(groups, key=lambda value: abs(value-requested))
    if abs(stamp-requested) > 2e-20:
        raise ValueError(f"no MFEM snapshot at {requested:g} s; temporal interpolation is unsupported")
    return stamp, groups[stamp]


def p0_indices(edges, positions):
    values = np.asarray(positions, dtype=float)
    if not np.isfinite(values).all() or np.any(values < edges[0]) or np.any(values > edges[-1]):
        raise ValueError("requested published coordinates lie outside the MFEM foil; extrapolation is unsupported")
    # Half-open cells own interior ties on their deeper/right side; final boundary belongs to last cell.
    return np.minimum(np.searchsorted(edges, values, side="right")-1, len(edges)-2)


def lateral_profile_columns(grid, fraction=.1):
    if not 0 < fraction < .5:
        raise ValueError("edge-band fraction must lie between zero and one half")
    edges = grid["x_edges_m"]
    mid = .5*(edges[:-1]+edges[1:])
    zero = .5*(edges[0]+edges[-1])
    center = np.flatnonzero((edges[:-1] <= zero) & (edges[1:] >= zero))
    nearest = np.unique([0, len(mid)-1])
    band = np.flatnonzero((mid <= edges[0]+fraction*(edges[-1]-edges[0]))
                         | (mid >= edges[-1]-fraction*(edges[-1]-edges[0])))
    if not len(band):
        band = nearest
    return center, nearest, band


def depth_comparison(grid, published_rows, edge_fraction):
    center, edge, band = lateral_profile_columns(grid, edge_fraction)
    widths = np.diff(grid["x_edges_m"])
    temperature = grid["temperature_k"]
    profiles = {"foil_middle": np.average(temperature[:, center], axis=1, weights=widths[center]),
                "foil_edge": np.average(temperature[:, edge], axis=1, weights=widths[edge])}
    output = []
    for row in published_rows:
        depth = float(row["nominal_depth_from_initial_gap_face_um"])*1e-6
        iz = int(p0_indices(grid["depth_edges_m"], [depth])[0])
        prediction = float(profiles[row["profile"]][iz])
        output.append({"profile": row["profile"], "nominal_depth_m": depth,
            "depth_readout_lower_m": float(row["nominal_depth_digitization_lower_um"])*1e-6,
            "depth_readout_upper_m": float(row["nominal_depth_digitization_upper_um"])*1e-6,
            "published_temperature_k": float(row["temperature_K"]),
            "published_readout_lower_k": float(row["temperature_digitization_lower_K"]),
            "published_readout_upper_k": float(row["temperature_digitization_upper_K"]),
            "mfem_p0_temperature_k": prediction,
            "mfem_minus_published_k": prediction-float(row["temperature_K"]),
            "mfem_edge_band_position_min_k": float(np.min(temperature[iz, band])),
            "mfem_edge_band_position_max_k": float(np.max(temperature[iz, band]))})
    return output, {"center_column_indices": center.tolist(), "nearest_edge_column_indices": edge.tolist(),
        "edge_band_column_indices": band.tolist(), "edge_band_fraction_per_side": edge_fraction,
        "center_column_intervals_m": [[float(grid["x_edges_m"][i]), float(grid["x_edges_m"][i+1])] for i in center],
        "nearest_edge_column_intervals_m": [[float(grid["x_edges_m"][i]), float(grid["x_edges_m"][i+1])] for i in edge]}


def width_comparison(grid, published_rows, depth_m=12.5e-6):
    iz = int(p0_indices(grid["depth_edges_m"], [depth_m])[0])
    output = []
    for row in published_rows:
        x = float(row["width_offset_mm"])*1e-3
        ix = int(p0_indices(grid["x_edges_m"], [x])[0])
        prediction = float(grid["temperature_k"][iz, ix])
        output.append({"x_m": x, "x_readout_lower_m": float(row["width_offset_digitization_lower_mm"])*1e-3,
            "x_readout_upper_m": float(row["width_offset_digitization_upper_mm"])*1e-3,
            "published_temperature_k": float(row["temperature_K"]),
            "published_readout_lower_k": float(row["temperature_digitization_lower_K"]),
            "published_readout_upper_k": float(row["temperature_digitization_upper_K"]),
            "mfem_p0_temperature_k": prediction, "mfem_minus_published_k": prediction-float(row["temperature_K"])})
    return output, {"nominal_depth_m": depth_m, "selected_depth_cell_index": iz,
        "selected_depth_cell_interval_m": grid["depth_edges_m"][iz:iz+2].tolist()}


def summarize(differences):
    return {"retained_points": len(differences),
        "max_absolute_temperature_difference_k": max(abs(row["mfem_minus_published_k"]) for row in differences),
        "statistical_likelihood": False, "parameter_fit": False}


def heating_contrast(grid, *, depth_m=12.5e-6):
    """State actual resolved field contrasts, without enforcing a hot edge."""
    center, edge, _ = lateral_profile_columns(grid)
    iz = int(p0_indices(grid["depth_edges_m"], [depth_m])[0])
    widths = np.diff(grid["x_edges_m"])
    result = {"nominal_depth_m": depth_m, "depth_cell_interval_m": grid["depth_edges_m"][iz:iz+2].tolist()}
    for field in ("temperature_k", "current_density_a_m2", "joule_power_w_m3"):
        values = grid[field][iz]
        if field == "current_density_a_m2":
            values = np.abs(values)
        edge_mean = float(np.average(values[edge], weights=widths[edge]))
        center_mean = float(np.average(values[center], weights=widths[center]))
        result[field] = {"nearest_edge_cell_mean": edge_mean, "central_cell_mean": center_mean,
            "edge_minus_center": edge_mean-center_mean,
            "edge_over_center": edge_mean/center_mean if center_mean != 0 else None}
    result["current_convention"] = "width-weighted mean absolute current density, not a fit to the unverified published J axis"
    return result


def make_plot(output, comparison, groups, diagnostics_path=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    colors = {"foil_edge": "tab:red", "foil_middle": "tab:blue"}
    for ax, item in zip(axes.flat[:2], comparison["depth_comparisons"]):
        _, grid = exact_snapshot(groups, item["requested_time_s"])
        c, e, band = lateral_profile_columns(grid, item["sampling"]["edge_band_fraction_per_side"])
        widths = np.diff(grid["x_edges_m"])
        for profile, columns in (("foil_edge", e), ("foil_middle", c)):
            data = sorted((row for row in item["readouts"] if row["profile"] == profile), key=lambda row: row["nominal_depth_m"])
            x = [row["nominal_depth_m"]*1e6 for row in data]
            y = [row["published_temperature_k"] for row in data]
            ax.plot(x, y, ".", ms=3, color=colors[profile], label=f"ALEGRA display: {profile.removeprefix('foil_')}")
            ax.fill_between(x, [row["published_readout_lower_k"] for row in data],
                            [row["published_readout_upper_k"] for row in data], color=colors[profile], alpha=.16)
            prediction = np.average(grid["temperature_k"][:, columns], axis=1, weights=widths[columns])
            ax.stairs(prediction, grid["depth_edges_m"]*1e6, color=colors[profile], lw=1.8, baseline=None,
                      label=f"MFEM P0: {'nearest edge cells' if profile == 'foil_edge' else 'central cells'}")
        depth_edges = grid["depth_edges_m"]*1e6
        band_low = np.min(grid["temperature_k"][:, band], axis=1)
        band_high = np.max(grid["temperature_k"][:, band], axis=1)
        ax.fill_between(depth_edges, np.r_[band_low, band_low[-1]], np.r_[band_high, band_high[-1]],
                        step="post", color="tab:red", alpha=.08,
                        label="MFEM outer-band position sensitivity")
        ax.set(xlabel="Nominal inward depth from initial gap face [µm]", ylabel="Temperature [K]",
               title=f"Figs. {'4.28' if item['requested_time_s'] == 40e-9 else '4.29'}: {item['requested_time_s']*1e9:g} ns")
        ratio = item["mfem_field_contrast"]["joule_power_w_m3"]["edge_over_center"]
        if ratio is not None:
            ax.text(.98, .97, f"MFEM q edge/middle at 12.5 µm: {ratio:.3g}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7, bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8})
        ax.legend(fontsize=7)
        ax.grid(alpha=.2)
    width = comparison["width_comparison"]
    ax = axes[1, 0]
    data = width["readouts"]
    x = [row["x_m"]*1e3 for row in data]
    ax.plot(x, [row["published_temperature_k"] for row in data], ".", ms=3, color="tab:orange", label="ALEGRA Fig.4.30, conditional eV→K")
    ax.fill_between(x, [row["published_readout_lower_k"] for row in data],
                    [row["published_readout_upper_k"] for row in data], color="tab:orange", alpha=.2)
    _, width_grid = exact_snapshot(groups, 40e-9)
    depth_index = width["sampling"]["selected_depth_cell_index"]
    ax.stairs(width_grid["temperature_k"][depth_index], width_grid["x_edges_m"]*1e3, baseline=None,
              color="black", label="MFEM P0 at nominal 12.5 µm")
    ax.set(xlabel="Distance from foil middle [mm]", ylabel="Temperature [K]", title="40 ns width lineout; conditional sampling")
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    ax = axes[1, 1]
    current = rows(DATA / "shot_17133_current.csv")
    selected = [r for r in current if float(r["time_s"]) <= 65e-9]
    ax.plot([float(r["time_s"])*1e9 for r in selected], [float(r["current_A"])*1e-3 for r in selected],
            color="black", label="Shot17133 C.19 published average current")
    for name, label in (("figure_4_28_simulation_drive.csv", "4.28 simulation inset"),
                        ("figure_4_29_simulation_drive.csv", "4.29 simulation inset")):
        data = rows(DATA / name)
        ax.plot([float(r["time_s"])*1e9 for r in data], [float(r["current_A"])*1e-3 for r in data],
                label=label, lw=1, alpha=.8)
    if diagnostics_path:
        data = rows(diagnostics_path)
        ax.plot([float(r["time_s"])*1e9 for r in data], [float(r["foil_current_a"])*1e-3 for r in data],
                "--", lw=1.2, label="MFEM imposed foil current")
    ax.set(xlabel="Literal published/native clock [ns]", ylabel="Current [kA]", xlim=(-3, 65),
           title="Drive context: no fitted retiming or amplitude repair")
    ax.legend(fontsize=7); ax.grid(alpha=.2)
    figure.suptitle("Provisional MFEM vs recovered ALEGRA displays\nReadout bounds are glyph envelopes, not measurement uncertainty", fontsize=13)
    figure.savefig(output / "mfem_alegra_comparison.png", dpi=150)
    figure.savefig(output / "mfem_alegra_comparison.pdf")
    plt.close(figure)


def compare(path, output_dir, *, metadata_path=None, diagnostics_path=None, edge_fraction=.1):
    output = Path(output_dir)
    if output.exists():
        raise ValueError("use a new output directory")
    groups, metadata = load_snapshots(path, metadata_path)
    comparisons = []
    inputs = {str(Path(path)): sha(path), str(Path(metadata_path) if metadata_path else Path(path).with_suffix('.json')):
              sha(metadata_path or Path(path).with_suffix('.json'))}
    for time, figure in ((40e-9, "4_28"), (60e-9, "4_29")):
        stamp, grid = exact_snapshot(groups, time)
        csv_path = DATA / f"figure_{figure}_temperature.csv"
        sidecar = csv_path.with_suffix(".json")
        source = read_json(sidecar)
        if source["simulation_context"]["time_s"] != time:
            raise ValueError("published temperature context has an unexpected snapshot time")
        values, sampling = depth_comparison(grid, rows(csv_path), edge_fraction)
        comparisons.append({"figure": figure.replace("_", "."), "requested_time_s": time,
            "actual_snapshot_time_s": stamp, "sampling": sampling, "readouts": values,
            "summary": summarize(values), "mfem_field_contrast": heating_contrast(grid)})
        inputs.update({str(csv_path): sha(csv_path), str(sidecar): sha(sidecar)})
    stamp, grid = exact_snapshot(groups, 40e-9)
    width_path = DATA / "figure_4_30_heating.csv"
    width_meta = read_json(width_path.with_suffix(".json"))
    conversion = width_meta["temperature_conversion"]
    data = rows(width_path)
    factor = conversion["elementary_charge_J_per_eV"] / conversion["Boltzmann_J_per_K"]
    if any(not math.isclose(float(row["temperature_K"]), float(row["temperature_eV"])*factor,
                            rel_tol=1e-10, abs_tol=1e-8) for row in data):
        raise ValueError("published conditional eV-to-K columns are inconsistent")
    width_values, width_sampling = width_comparison(grid, data)
    inputs.update({str(width_path): sha(width_path), str(width_path.with_suffix('.json')): sha(width_path.with_suffix('.json'))})
    for name in ("shot_17133_current.csv", "shot_17133_current.json", "figure_4_28_simulation_drive.csv", "figure_4_28_simulation_drive.json",
                 "figure_4_29_simulation_drive.csv", "figure_4_29_simulation_drive.json"):
        inputs[str(DATA/name)] = sha(DATA/name)
    if diagnostics_path:
        inputs[str(Path(diagnostics_path))] = sha(diagnostics_path)
    result = {"schema": "mhdlab-mfem-alegra-temperature-comparison-v1", "status": "provisional_conditional_comparison",
        "mfem_snapshot_metadata": metadata, "source_sha256": inputs,
        "depth_comparisons": comparisons,
        "width_comparison": {"figure": "4.30", "time_s": 40e-9, "sampling": width_sampling,
            "temperature_interpretation": conversion, "readouts": width_values, "summary": summarize(width_values)},
        "sampling_convention": "MFEM cell averages read as a P0 reconstruction; half-open right/deeper cell owns interior ties. No spatial/temporal interpolation or extrapolation.",
        "edge_location_policy": "Published precise transverse edge-lineout offset is unavailable. Nearest resolved edge-cell mean and a declared outer-band sensitivity are conditional choices, not matched physical probes.",
        "limits": ["ALEGRA displays are recovered simulation curves, not measured temperatures or original simulation arrays.",
            "No drive amplitude fit, clock registration, material calibration, likelihood or physical validation is performed.",
            "Current C.19 and simulation inset curves differ on literal plotted coordinates; no original ALEGRA drive file is recovered.",
            "Published coordinate/glyph bounds are readout support with unknown correlations, not stochastic errors.",
            "Fig.4.30 Kelvin values are conditional on interpreting its eV label as k_B*T.",
            "Physical phase changes, material properties, edge radius, thermal boundaries, support, exact probe locations and source clocks may differ."],
        "physical_alegra_reproduction_established": False, "tool_source_sha256": sha(__file__)}
    output.mkdir(parents=True)
    make_plot(output, result, groups, diagnostics_path)
    result["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
    (output / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshots", required=True, type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--diagnostics", type=Path)
    parser.add_argument("--edge-fraction", type=float, default=.1)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = compare(args.snapshots, args.output_dir, metadata_path=args.metadata,
                         diagnostics_path=args.diagnostics, edge_fraction=args.edge_fraction)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"MFEM/ALEGRA comparison: {exc}\n")
    print(json.dumps({"comparison": str(args.output_dir / "comparison.json"), "status": result["status"]}))


if __name__ == "__main__":
    main()
