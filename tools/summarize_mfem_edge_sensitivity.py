"""Compare immutable current-driven MFEM runs without fitting their fields."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_rows(path):
    with path.open(newline="") as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def overlap(a, b, c, d):
    return max(0.0, min(b, d) - max(a, c))


def weighted_mean(cells, field, intervals):
    total = weight = 0.0
    for row in cells:
        dx = sum(overlap(row["x_left_m"], row["x_right_m"], a, b) for a, b in intervals)
        area = dx * (row["depth_right_m"] - row["depth_left_m"])
        total += row[field] * area
        weight += area
    if weight <= 0:
        raise ValueError("band has no supported area")
    return total / weight


def project(source, edges, field):
    """Average nonoverlapping source patches onto fully supported target bins."""
    result = []
    for left, right in zip(edges, edges[1:]):
        total = weight = 0.0
        for row in source:
            dx = overlap(left, right, row["x_left_m"], row["x_right_m"])
            total += dx * row[field]
            weight += dx
        if not math.isclose(weight, right-left, rel_tol=1e-8, abs_tol=1e-14):
            raise ValueError("lineout does not cover target bin")
        result.append(total / weight)
    return result


def interior_line(cells, depth):
    """Linear reconstruction between cell-average depth centers at fixed depth."""
    rows = []
    by_x = {}
    for row in cells:
        by_x.setdefault((row["x_left_m"], row["x_right_m"]), []).append(row)
    for (left, right), column in sorted(by_x.items()):
        column.sort(key=lambda r: r["depth_left_m"])
        centers = [(r["depth_left_m"]+r["depth_right_m"])/2 for r in column]
        if not centers[0]-1e-14 <= depth <= centers[-1]+1e-14:
            raise ValueError("fixed-depth line requires unsupported depth extrapolation")
        for i in range(len(column)-1):
            if centers[i]-1e-14 <= depth <= centers[i+1]+1e-14:
                fraction = (depth-centers[i])/(centers[i+1]-centers[i])
                rows.append({"x_left_m": left, "x_right_m": right, **{
                    field: (1-fraction)*column[i][field]+fraction*column[i+1][field]
                    for field in ("temperature_k", "current_density_a_m2", "joule_power_w_m3")}})
                break
        else:
            if abs(depth-centers[-1]) < 1e-14:
                rows.append(dict(column[-1]))
            else:
                raise ValueError("missing fixed depth")
    return rows


def surface_line(directory, face, stamp, width):
    data = [row for row in numeric_rows(directory/(face+".csv")) if abs(row["time_s"]-stamp)<1e-18]
    data.sort(key=lambda r: r["x_index"])
    count = len(data)
    if not count:
        raise ValueError("missing surface time")
    return [{"x_left_m": -width/2+width*i/count, "x_right_m": -width/2+width*(i+1)/count,
             "temperature_k": 300+r["temperature_value"]} for i, r in enumerate(data)]


def difference(values, baseline, offset=0.0):
    rms = math.sqrt(sum((a-b)**2 for a, b in zip(values, baseline))/len(baseline))
    denominator = math.sqrt(sum((b-offset)**2 for b in baseline)/len(baseline))
    return {"rms_absolute_change": rms,
            "rms_change_over_baseline_rms": rms/denominator if denominator else None,
            "maximum_absolute_change": max(abs(a-b) for a, b in zip(values, baseline))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--case", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("use a fresh output directory")
    cases = [args.baseline, *args.case]
    record = {"schema": "mhdlab-mfem-edge-sensitivity-v1", "status": "provisional",
              "physical_experiment_prediction": False, "fixed_probe_depth_m": 12.5e-6,
              "band_sampling": "cell-average values weighted by exact geometric overlap; piecewise-constant spatial reconstruction",
              "lineout_sampling": "interior values reconstructed linearly between cell-average depth centers, then conservatively averaged onto baseline width bins; true surface patches projected by width overlap",
              "cases": {}}
    lineouts = []
    reference_lines = {}
    baseline_config = json.loads((args.baseline/"source_run.json").read_bytes())["configuration"]
    width = baseline_config["foil_width_m"]
    nx = baseline_config["width_cells"]
    edges = [-width/2+width*i/nx for i in range(nx+1)]
    for directory in cases:
        run = json.loads((directory/"source_run.json").read_bytes())
        if not run["passed"]:
            raise ValueError("cannot compare an unsuccessful source")
        config = run["configuration"]
        sidecar = json.loads((directory/"foil_snapshots.json").read_bytes())
        if sidecar["csv_sha256"] != sha(directory/"foil_snapshots.csv"):
            raise ValueError("snapshot hash mismatch")
        cells = numeric_rows(directory/"foil_snapshots.csv")
        item = {"run_record_sha256": sha(directory/"source_run.json"), "configuration": config,
                "drive_sha256": run["drive"]["csv_sha256"], "diagnostics": run["diagnostics"], "snapshots": {}}
        for stamp in (40e-9, 60e-9):
            selected = [r for r in cells if abs(r["time_s"]-stamp)<1e-18]
            if not selected:
                raise ValueError("missing exact40/60ns snapshot")
            edge = [(-width/2,-.4*width),(.4*width,width/2)]
            center = [(-.1*width,.1*width)]
            q_edge = weighted_mean(selected,"joule_power_w_m3",edge)
            q_center = weighted_mean(selected,"joule_power_w_m3",center)
            t_edge = weighted_mean(selected,"temperature_k",edge)
            t_center = weighted_mean(selected,"temperature_k",center)
            metrics = {"edge_mean_temperature_k":t_edge,"center_mean_temperature_k":t_center,
                "edge_mean_joule_power_w_m3":q_edge,"center_mean_joule_power_w_m3":q_center,
                "edge_to_center_joule_ratio":q_edge/q_center,
                "maximum_cell_average_temperature_k":max(r["temperature_k"] for r in selected)}
            for face in ("top","bottom"):
                face_line = surface_line(directory,face,stamp,width)
                metrics[face+"_maximum_patch_temperature_k"] = max(r["temperature_k"] for r in face_line)
                metrics[face+"_mean_patch_temperature_k"] = sum(r["temperature_k"] for r in face_line)/len(face_line)
            lines = {"interior":interior_line(selected,12.5e-6),
                     "top":surface_line(directory,"top",stamp,width),
                     "bottom":surface_line(directory,"bottom",stamp,width)}
            comparison = {}
            for location, source in lines.items():
                fields = ("temperature_k","current_density_a_m2","joule_power_w_m3") if location=="interior" else ("temperature_k",)
                for field in fields:
                    values = project(source,edges,field)
                    key = (stamp,location,field)
                    if directory==args.baseline:
                        reference_lines[key] = values
                    comparison[location+"_"+field] = difference(values,reference_lines[key],300 if field=="temperature_k" else 0)
                    for i, value in enumerate(values):
                        lineouts.append([str(directory),stamp,location,field,edges[i],edges[i+1],value,reference_lines[key][i]])
            metrics["lineout_change_from_baseline"] = comparison
            item["snapshots"][str(stamp)] = metrics
        record["cases"][str(directory)] = item
    args.output_dir.mkdir(parents=True)
    with (args.output_dir/"lineouts.csv").open("w",newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case","time_s","location","field","x_left_m","x_right_m","value","baseline_value"])
        writer.writerows(lineouts)
    record["lineouts_csv_sha256"] = sha(args.output_dir/"lineouts.csv")
    (args.output_dir/"sensitivity.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(args.output_dir/"sensitivity.json")


if __name__ == "__main__":
    main()
