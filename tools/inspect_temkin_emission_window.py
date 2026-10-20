#!/usr/bin/env python3
"""Inspect a partial source-emission window and its number/energy representation."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from mhdlab.planar_release_io import plain_metadata
from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.temkin_source_window import FrozenTemkinSourceWindow

CODE_PATHS = (
    "tools/inspect_temkin_emission_window.py", "mhdlab/temkin_emission_window.py",
    "mhdlab/temkin_source_window.py", "mhdlab/temkin_planar_io.py", "mhdlab/temkin_planar.py",
    "mhdlab/temkin_history.py", "mhdlab/temkin_release.py", "mhdlab/temkin_io.py",
    "mhdlab/thermal_release.py", "mhdlab/planar_release_io.py", "mhdlab/planar_release.py",
    "mhdlab/surface_release_history.py", "mhdlab/surface_sources.py", "mhdlab/flash_surface.py",
    "mhdlab/hdf_integrity.py", "mhdlab/constants.py",
)
CONFIG_KEYS = {"schema", "parameter_status", "parameter_provenance", "component", "window", "probabilities"}
COMPONENT_KEYS = {"reservoir", "initial_coverage_m2", "temperature_k", "parent_duration_s", "patch_area_m2"}
INDEX_KEYS = ("interval_index", "reservoir_index", "current_index", "width_index")
WINDOW_KEYS = {"start_s", "stop_s", "log_mass_tolerance"}


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    def pairs(items):
        result={}
        for key,value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key]=value
        return result
    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    return json.loads(path.read_bytes(),object_pairs_hook=pairs,parse_constant=nonfinite)


def configuration(path):
    raw=path.read_bytes()
    config=read_json(path)
    if not isinstance(config,dict) or set(config)!=CONFIG_KEYS or config["schema"]!="mhdlab-temkin-window-inspection-config-v1":
        raise ValueError("configuration requires exactly the mhdlab-temkin-window-inspection-config-v1 fields")
    if config["parameter_status"] not in {"manufactured","reconstructed","sourced","unvalidated","owner_supplied"}:
        raise ValueError("parameter_status must declare a supported input label")
    if not isinstance(config["parameter_provenance"],str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must be nonempty")
    if not isinstance(config["window"],dict) or set(config["window"])!=WINDOW_KEYS:
        raise ValueError("window requires explicit original offsets and log-mass tolerance")
    probabilities=config["probabilities"]
    if (not isinstance(probabilities,list) or not 2<=len(probabilities)<=1000
            or any(isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value)
                   or not 0<value<1 for value in probabilities)
            or any(a>=b for a,b in zip(probabilities,probabilities[1:]))):
        raise ValueError("probabilities require 2 to 1000 strictly increasing open probabilities")
    return raw,config


def source_component(args,config,source_hash):
    component=config["component"]
    if args.demo:
        if config["parameter_status"]!="manufactured":
            raise ValueError("--demo requires explicitly manufactured parameters")
        if not isinstance(component,dict) or set(component)!=COMPONENT_KEYS:
            raise ValueError("manufactured component requires a reservoir, initial density, temperature, duration and actual area")
        parameters=dict(component)
        parameters["reservoir"]=TemkinReservoir(**parameters["reservoir"])
        return parameters,{"kind":"manufactured_scalar_source_component"},None
    if not isinstance(component,dict) or set(component)!=set(INDEX_KEYS):
        raise ValueError("source component requires exactly four integer interval/reservoir/current/width indices")
    if any(isinstance(value,bool) or not isinstance(value,int) or value<0 for value in component.values()):
        raise ValueError("component indices must be nonnegative integers")
    history=load_temkin_planar_release_history(args.source)
    if args.source_run_summary is None:
        if config["parameter_status"]!="unvalidated":
            raise ValueError("without a hash-associated source summary, external parameters require status unvalidated")
    else:
        summary=read_json(args.source_run_summary)
        if (not isinstance(summary,dict) or summary.get("schema") not in {
                "mhdlab-temkin-planar-run-v1","mhdlab-joint-source-run-v1","mhdlab-coverage-source-run-v1"}
                or summary.get("parameter_status")!=config["parameter_status"]
                or summary.get("parameter_provenance")!=history.parameter_provenance
                or not isinstance(summary.get("output_sha256"),dict)
                or source_hash not in summary["output_sha256"].values()):
            raise ValueError("source summary must hash-associate the exact source, parameter status and provenance")
    t,r,z,x=(component[name] for name in INDEX_KEYS)
    if not (t<history.temperature_k.shape[0] and r<len(history.reservoirs)
            and z<history.temperature_k.shape[1] and x<history.temperature_k.shape[2]):
        raise ValueError("selected component is outside the source dimensions")
    parameters={"reservoir":history.reservoirs[r],"initial_coverage_m2":float(history.coverage_m2[t,r,z,x]),
                "temperature_k":float(history.temperature_k[t,z,x]),
                "parent_duration_s":float(history.time_edges_s[t+1]-history.time_edges_s[t]),
                "patch_area_m2":float(history.patch_areas_m2[z,x])}
    context={"kind":"selected_canonical_planar_source_component","indices":component,
             "source_parameter_provenance":history.parameter_provenance,
             "surface_attributes":plain_metadata(history.surface_attributes),
             "original_interval_start_s":float(history.time_edges_s[t]),
             "original_interval_stop_s":float(history.time_edges_s[t+1]),
             "x_bounds_m":[float(history.x_edges_m[x]),float(history.x_edges_m[x+1])],
             "z_bounds_m":[float(history.z_edges_m[z]),float(history.z_edges_m[z+1])],
             "whole_interval_released_m2":float(history.released_m2[t,r,z,x])}
    return parameters,context,history


def plot_window(rows,path,status):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure,axis=plt.subplots(figsize=(6.5,4),layout="constrained")
    axis.plot([row["log_local_delay_s"] for row in rows],[row["probability"] for row in rows],"o-")
    axis.set(xlabel="ln(local delay / 1 s)",ylabel="Conditional emitted-population probability",
             title=f"Parameter status: {status}\nSource-emission window; no photon gate")
    axis.grid(alpha=.25)
    figure.savefig(path,dpi=160)
    plt.close(figure)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--demo",action="store_true")
    group.add_argument("--source",type=Path)
    parser.add_argument("--config",type=Path)
    parser.add_argument("--source-run-summary",type=Path)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--require-represented-budget",action="store_true")
    parser.add_argument("--plot",action="store_true")
    args=parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    if args.demo and args.source_run_summary is not None:
        parser.error("--source-run-summary applies only to --source")
    if args.source is not None and args.config is None:
        parser.error("--source requires an explicit --config")
    config_path=args.config or ROOT/"configs/temkin_window_manufactured.json"
    try:
        code_hashes={name:sha(ROOT/name) for name in CODE_PATHS}
        raw,config=configuration(config_path)
        config_hash=hashlib.sha256(raw).hexdigest()
        source_hash=sha(args.source) if args.source else None
        summary_hash=sha(args.source_run_summary) if args.source_run_summary else None
        parameters,context,history=source_component(args,config,source_hash)
        window=FrozenTemkinSourceWindow(**parameters,**config["window"],status=config["parameter_status"],
                                        provenance=config["parameter_provenance"])
        manifest=window.manifest()
        if history is not None:
            partial=manifest["quantities"]["released_m2"]["value"]
            whole=context["whole_interval_released_m2"]
            if partial is not None and partial>whole and (whole==0 or math.log(partial/whole)>window.log_mass_tolerance):
                raise ValueError("partial amount exceeds the canonical parent interval budget within the declared allowance")
        budget=window.represented_budget() if window.all_budget_values_qualified else None
        if args.require_represented_budget and budget is None:
            raise ValueError("a represented number/energy budget was required but the window contains unqualified quantities")
        rows=[]
        for probability in config["probabilities"]:
            log_delay=float(window.time_law.log_quantile(probability))
            # A positive log delay remains useful even when its linear value
            # cannot satisfy the timing law's representation qualification.
            try:
                delay=float(window.time_law.quantile(probability))
                cdf=float(window.time_law.cdf(delay))
                flag="represented_local_delay"
            except ValueError:
                delay,cdf=None,None
                flag="linear_local_delay_unqualified_use_log"
            rows.append(dict(probability=probability,log_local_delay_s=log_delay,
                             local_delay_s=delay,cdf_at_local_delay=cdf,representation=flag))
        report={"schema":"mhdlab-temkin-window-inspection-v1","parameter_status":config["parameter_status"],
                "parameter_provenance":config["parameter_provenance"],"source_component":context,
                "window":manifest,"represented_budget":budget,
                "input_config_sha256":config_hash,"input_source_sha256":source_hash,
                "input_summary_sha256":summary_hash,"source_sha256":code_hashes,
                "runtime":{"python":sys.version,**{name:importlib.metadata.version(name) for name in ("numpy","scipy","h5py")}},
                "limits":["This selects emission during one frozen source interval; earlier emission may still contribute photons in a later diagnostic gate.",
                          "Original source start, interval origin and local conditional delay are separate. No absolute timestamp precision or shared diagnostic timing origin is inferred.",
                          "Supplied parameter labels and a hash-associated run summary preserve provenance; they do not independently validate physical rates or calibration.",
                          "Number/energy representation is separate from host acceptance, signed ledger compensation, feedback or particle launching.",
                          "Staged output assumes a coordinated single writer; ordinary existing paths are refused."]}
        args.output_dir.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".temkin-window-",dir=args.output_dir.parent) as stage:
            output=Path(stage)/"result"
            output.mkdir()
            (output/"input_config.json").write_bytes(raw)
            if args.source:
                shutil.copyfile(args.source,output/"source_input.h5")
                if sha(output/"source_input.h5")!=source_hash:
                    raise ValueError("source changed while copying")
            if args.source_run_summary:
                shutil.copyfile(args.source_run_summary,output/"source_run_summary.json")
            with (output/"quantiles.csv").open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            if args.plot:
                plot_window(rows,output/"source_window.png",config["parameter_status"])
            if (sha(config_path)!=config_hash or (args.source and sha(args.source)!=source_hash)
                    or (args.source_run_summary and (sha(args.source_run_summary)!=summary_hash
                         or sha(output/"source_run_summary.json")!=summary_hash))
                    or any(sha(ROOT/name)!=digest for name,digest in code_hashes.items())):
                raise ValueError("input or implementation changed during inspection")
            report["output_sha256"]={p.name:sha(p) for p in output.iterdir() if p.is_file()}
            (output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during inspection; refusing replacement")
            os.rename(output,args.output_dir)
    except (OSError,ValueError,TypeError,KeyError,OverflowError,RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir":str(args.output_dir),"parameter_status":config["parameter_status"],
                      "represented_budget_available":budget is not None,"quantiles":len(rows)}))


if __name__=="__main__":
    main()
