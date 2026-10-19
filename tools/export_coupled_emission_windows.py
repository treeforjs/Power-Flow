#!/usr/bin/env python3
"""Export fixed emission windows from explicitly bound accepted saved source bytes.

Source histories, RHS caches and width assembly are read as data. Only the new
pure polynomial window core is imported, after raw queries and input hashes
have been checked. No source evolution, RHS or assembly is repeated.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CONFIG_SCHEMA = "mhdlab-coupled-emission-windows-manufactured-config-v1"
REPORT_SCHEMA = "mhdlab-coupled-emission-windows-export-run-v1"
ALIASES = {"report", "assembly", "input_config", "reuse_map", "uniform_power", "edge_power", "middle_power", "independent_completion", "independent_report", "scientific_review"}
CASES = {"uniform_power", "edge_power", "middle_power"}
HEX = re.compile(r"[0-9a-f]{64}\Z")
INTEGER = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")
MAX_FILE_BYTES = 64*1024*1024
MAX_TOTAL_BYTES = 128*1024*1024
FORBIDDEN = ("mhdlab.thermal_width_power_reference", "mhdlab.thermal_forced_feedback_reference", "mhdlab.thermal_power_reference", "mhdlab.thermal_feedback_reference")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return digest(Path(path).read_bytes())


def _constant(token):
    raise ValueError("nonfinite JSON constant: "+token)


def _decimal(token):
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError("unsupported JSON numeric token") from exc
    projected = float(value)
    if not value.is_finite() or not math.isfinite(projected) or value != 0 and projected == 0:
        raise ValueError("JSON numeric token has unsupported finite/nonzero representation")
    return value


def _object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON key: "+key)
        value[key] = item
    return value


def _text(value, name):
    if type(value) is not str or not value.strip():
        raise ValueError(name+" must be nonempty supplied text")


def _sha(value, name):
    if type(value) is not str or not HEX.fullmatch(value):
        raise ValueError(name+" must be an explicit lower-case SHA256")


def _relative(value):
    if (type(value) is not str or not value or len(value)>512 or "\\" in value or "\x00" in value
            or PurePosixPath(value).is_absolute() or any(p in ("", ".", "..") for p in value.split("/"))
            or str(PurePosixPath(value)) != value):
        raise ValueError("source paths must be unambiguous safe relative POSIX paths")
    return value


def _fraction(record):
    if type(record) is not dict or set(record)!={"numerator", "denominator"}:
        raise ValueError("query fractions need explicit canonical numerator and denominator strings")
    for key in ("numerator", "denominator"):
        token = record[key]
        if type(token) is not str or len(token)>128 or not INTEGER.fullmatch(token) or str(int(token)) != token:
            raise ValueError("query rational integers must be bounded canonical strings")
    n,d=int(record["numerator"]),int(record["denominator"])
    if d<=0:
        raise ValueError("query denominator must be positive")
    result=Fraction(n,d)
    if (result.numerator,result.denominator)!=(n,d):
        raise ValueError("query fractions must be reduced")
    return result


def load_config(path):
    config=json.loads(Path(path).read_text(),parse_float=_decimal,parse_constant=_constant,object_pairs_hook=_object)
    keys={"schema","parameter_status","provenance","scope","source_binding","queries","quantity_contract","representation_contract","declared_counts"}
    if type(config) is not dict or set(config)!=keys or config["schema"]!=CONFIG_SCHEMA or config["parameter_status"]!="manufactured":
        raise ValueError("unsupported explicit manufactured window configuration")
    for key in ("provenance","scope"):
        _text(config[key],key)
    binding=config["source_binding"]
    if type(binding) is not dict or set(binding)!={"path_base","files","canonical_history_sha256","numerical_run_code_sha256","lineage"}:
        raise ValueError("complete frozen saved-source binding is required")
    _text(binding["path_base"],"source path base")
    files=binding["files"]
    if type(files) is not dict or set(files)!=ALIASES:
        raise ValueError("all ten fixed saved-source aliases are required")
    seen=set()
    for alias,record in files.items():
        if type(record) is not dict or set(record)!={"relative_path","sha256"}:
            raise ValueError("each source alias requires path and expected hash")
        relative=_relative(record["relative_path"])
        if relative in seen:
            raise ValueError("source aliases must identify distinct paths")
        seen.add(relative)
        _sha(record["sha256"],alias)
    canonical=binding["canonical_history_sha256"]
    if type(canonical) is not dict or set(canonical)!=CASES:
        raise ValueError("three explicit canonical history digests are required")
    for key,value in canonical.items():
        _sha(value,key)
    versions=binding["numerical_run_code_sha256"]
    expected_versions={"tools/run_thermal_width_power_reference.py", "mhdlab/thermal_width_power_reference.py", "mhdlab/thermal_forced_feedback_reference.py", "mhdlab/thermal_power_reference.py", "mhdlab/thermal_feedback_reference.py"}
    if type(versions) is not dict or set(versions)!=expected_versions:
        raise ValueError("the complete numerical source version vector is required")
    for key,value in versions.items():
        _sha(value,key)
    lineage=binding["lineage"]
    if type(lineage) is not dict or set(lineage)!={"numerical_run_cli_sha256","plot_only_public_cli_sha256","auxiliary_files"}:
        raise ValueError("numerical-versus-plot source lineage is required")
    for key in ("numerical_run_cli_sha256","plot_only_public_cli_sha256"):
        _sha(lineage[key],key)
    if lineage["numerical_run_cli_sha256"]!=versions["tools/run_thermal_width_power_reference.py"]:
        raise ValueError("numerical CLI lineage differs from source version vector")
    aux=lineage["auxiliary_files"]
    if type(aux) is not dict or set(aux)!={"executed_cli_snapshot","plot_correction"}:
        raise ValueError("the two explicit source-lineage files are required")
    for alias,record in aux.items():
        if type(record) is not dict or set(record)!={"relative_path","sha256"}:
            raise ValueError("lineage aliases require path and expected hash")
        relative=_relative(record["relative_path"])
        if relative in seen:
            raise ValueError("source and lineage aliases must be distinct")
        seen.add(relative);_sha(record["sha256"],alias)
    if aux["executed_cli_snapshot"]["sha256"]!=lineage["numerical_run_cli_sha256"]:
        raise ValueError("executed CLI bytes are not associated with numerical lineage")
    queries=config["queries"]
    expected={"complete_source_support":(Fraction(0),Fraction(1)),"early_interior":(Fraction(1,16),Fraction(5,16)),"off_knot_interior":(Fraction(23,128),Fraction(71,128))}
    if type(queries) is not list or len(queries)!=3:
        raise ValueError("exactly three declared source-emission windows are required")
    ids=set()
    for query in queries:
        if type(query) is not dict or set(query)!={"window_id","coordinate","start_fraction","stop_fraction","parameter_provenance"}:
            raise ValueError("complete exact query identity/support/provenance is required")
        key=query["window_id"]
        if type(key) is not str or key not in expected or key in ids:
            raise ValueError("window ids must be unique and explicitly declared")
        ids.add(key)
        if query["coordinate"]!="canonical_duration_fraction":
            raise ValueError("queries must use the exact canonical source-duration coordinate")
        a,b=_fraction(query["start_fraction"]),_fraction(query["stop_fraction"])
        if not 0<=a<b<=1 or (a,b)!=expected[key]:
            raise ValueError("query support differs from the bounded declared exact windows")
        _text(query["parameter_provenance"],"query provenance")
    required_quantity={"local_time_anchor":"exact_window_start","joint_local_time_temperature_moment_required":True,"source_RK_channel_comparison":"available_only_when_both_bounds_are_source_knots","time_temperature_mixture":"full_coefficients_and_joint_moments_retained","energy_reference":"exact_SI_constants_and_supplied_D_I_separate_from_activation_barrier"}
    if type(config["quantity_contract"]) is not dict or config["quantity_contract"]!=required_quantity or config["quantity_contract"].get("joint_local_time_temperature_moment_required") is not True:
        raise ValueError("the complete joint source moment/energy contract is required")
    required_rep={"arithmetic":"exact_rational_polynomial_coefficients_products_and_sums","normalized_zero_mass":"unavailable_with_explicit_reason_no_pseudocount","unrepresentable_nonzero_output":"reject_without_floor_clip_or_repair","source_residuals":"retained_not_repaired_by_interpolant"}
    if type(config["representation_contract"]) is not dict or config["representation_contract"]!=required_rep:
        raise ValueError("unsupported projection or source-repair contract")
    counts={"unique_source_window_calculations":9,"patch_window_references":18,"hypothesis_window_aggregates":6,"nonuniform_minus_uniform_contrasts":3,"source_evolution_RHS_assembly_calls":0}
    if type(config["declared_counts"]) is not dict or config["declared_counts"]!=counts or any(type(x) is not int for x in config["declared_counts"].values()):
        raise ValueError("explicit bounded calculation/reuse counts are required")
    return config


def _read_relative(base,relative):
    relative=_relative(relative);path=base
    for part in relative.split("/"):
        path=path/part
        if path.is_symlink():
            raise ValueError("source symlinks are unsupported")
    if not path.is_file() or not path.resolve().is_relative_to(base) or path.stat().st_size>MAX_FILE_BYTES:
        raise ValueError("source file is missing, external or exceeds the bounded byte limit")
    return path.read_bytes()


def _source_json(data):
    return json.loads(data.decode("utf-8"),parse_constant=_constant,object_pairs_hook=_object)


def read_bound_source(config,source_dir):
    """Read and bind immutable bytes before importing any polynomial core."""
    base=source_dir.resolve(strict=True)
    if not base.is_dir():
        raise ValueError("source-dir must be the declared accepted pass directory")
    binding=config["source_binding"];files={};aux={};total=0
    for destination,specs in ((files,binding["files"]),(aux,binding["lineage"]["auxiliary_files"])):
        for alias,spec in specs.items():
            data=_read_relative(base,spec["relative_path"]);total+=len(data)
            if total>MAX_TOTAL_BYTES or digest(data)!=spec["sha256"]:
                raise ValueError("source/lineage byte hash or bounded total differs: "+alias)
            destination[alias]=data
    report=_source_json(files["report"])
    if report.get("schema")!="mhdlab-thermal-width-power-reference-run-v1" or report.get("status")!="manufactured_reference_only":
        raise ValueError("unsupported numerical source report schema/status")
    if report.get("code_sha256")!=binding["numerical_run_code_sha256"]:
        raise ValueError("source numerical code versions differ")
    if report.get("config_sha256")!=binding["files"]["input_config"]["sha256"]:
        raise ValueError("source configuration association differs")
    if report.get("artifact_path_base")!="directory containing report.json":
        raise ValueError("unsupported source report artifact base")
    for alias in ("assembly","input_config","reuse_map",*sorted(CASES)):
        name=Path(binding["files"][alias]["relative_path"]).name
        if report.get("artifact_sha256",{}).get(name)!=digest(files[alias]):
            raise ValueError("source report artifact association differs: "+alias)
    for alias in CASES:
        history=_source_json(files[alias])
        canonical=digest(json.dumps(history,sort_keys=True,separators=(",",":"),allow_nan=False).encode())
        if canonical!=binding["canonical_history_sha256"][alias]:
            raise ValueError("canonical history association differs: "+alias)
    accepted=_source_json(files["independent_completion"])
    if accepted.get("schema")!="mitl-independent-width-power-completion-v1" or accepted.get("status")!="passed" or accepted.get("all_checks_passed") is not True:
        raise ValueError("accepted source verification is absent")
    if accepted.get("producer_report_sha256")!=digest(files["report"]) or accepted.get("report_sha256")!=digest(files["independent_report"]):
        raise ValueError("accepted source verification/report association differs")
    review=_source_json(files["scientific_review"])
    if review.get("schema")!="mitl-width-power-static-oracle-review-v1" or review.get("status")!="accepted":
        raise ValueError("source scientific acceptance is absent")
    plot=_source_json(aux["plot_correction"])
    if plot.get("numerical_run_report_sha256")!=digest(files["report"]) or plot.get("numerical_run_cli_sha256")!=binding["lineage"]["numerical_run_cli_sha256"] or plot.get("final_public_cli_sha256")!=binding["lineage"]["plot_only_public_cli_sha256"]:
        raise ValueError("plot-only version lineage differs")
    return files,aux


def _late_recheck(config,source_dir):
    """Recheck original immutable byte identities; no JSON or numerical replay."""
    base=source_dir.resolve(strict=True);binding=config["source_binding"]
    for specs in (binding["files"],binding["lineage"]["auxiliary_files"]):
        for alias,spec in specs.items():
            if digest(_read_relative(base,spec["relative_path"]))!=spec["sha256"]:
                raise ValueError("late source byte identity changed: "+alias)


def _write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def _flatten(value,prefix=""):
    if type(value) is dict and set(value)=={"exact","value","projection"}:
        return {prefix:value["value"],prefix+"__exact_numerator":value["exact"]["numerator"],prefix+"__exact_denominator":value["exact"]["denominator"],prefix+"__projection_json":json.dumps(value["projection"],sort_keys=True,separators=(",",":"),allow_nan=False)}
    if type(value) is dict:
        result={}
        for key,item in value.items():
            result.update(_flatten(item,prefix+"__"+key if prefix else key))
        return result
    if type(value) is list:
        return {prefix:json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False)}
    return {prefix:"null" if value is None else value}


def _write_csv(path,rows):
    if not rows:
        raise ValueError("required output table is empty")
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n",restval="null");writer.writeheader()
        for row in rows:
            writer.writerow({key:format(value,".17g") if type(value) is float else "true" if value is True else "false" if value is False else value for key,value in row.items()})


def _export_tables(document,output):
    tables={"source_windows":[],"patch_windows":[],"aggregate_windows":[],"width_window_contrasts":[]};pieces=[]
    for item in document["source_windows"]:
        row={k:v for k,v in item.items() if k!="pieces"};row["piece_count"]=len(item["pieces"])
        tables["source_windows"].append(_flatten(row))
        for index,piece in enumerate(item["pieces"]):
            pieces.append(_flatten({"source_history_id":item["source_history_id"],"query_id":item["query_id"],"state_id":item["state_id"],"piece_index":index,"piece":piece}))
    for key in ("patch_windows","aggregate_windows","width_window_contrasts"):
        tables[key]=[_flatten(item) for item in document[key]]
    for key,expected in (("source_windows",9),("patch_windows",18),("aggregate_windows",6),("width_window_contrasts",3)):
        if len(tables[key])!=expected:
            raise ValueError("bounded output count differs: "+key)
        _write_csv(output/(key+".csv"),tables[key])
    _write_csv(output/"polynomial_pieces.csv",pieces)
    return {**{k:len(v) for k,v in tables.items()},"polynomial_pieces":len(pieces)}


def run(config_path,source_dir,output):
    config=load_config(config_path)
    files,aux=read_bound_source(config,source_dir)
    if output.exists():
        raise ValueError("output must be a fresh directory")
    # No old source/history/RHS/assembly module is imported here or by the pure core.
    sys.path.insert(0,str(ROOT))
    from mhdlab.coupled_emission_windows import admit_saved_width_source,WindowQuery,build_coupled_emission_windows
    source=admit_saved_width_source(files,expected_sha256={key:v["sha256"] for key,v in config["source_binding"]["files"].items()},parameter_provenance=config["provenance"])
    queries=tuple(WindowQuery(window_id=q["window_id"],start_fraction=_fraction(q["start_fraction"]),stop_fraction=_fraction(q["stop_fraction"]),parameter_provenance=q["parameter_provenance"]) for q in config["queries"])
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix=output.name+".staging-",dir=output.parent))
    try:
        (stage/"input_config.json").write_bytes(config_path.read_bytes())
        copied={};binding=config["source_binding"]
        for data,specs in ((files,binding["files"]),(aux,binding["lineage"]["auxiliary_files"])):
            for alias,value in data.items():
                relative="source-inputs/"+specs[alias]["relative_path"];path=stage/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(value)
                copied[relative]=digest(value)
        # Sole saved-data calculation call; internal cache contains nine unique pairs.
        result=build_coupled_emission_windows(source,queries,parameter_status=config["parameter_status"],parameter_provenance=config["provenance"])
        document=result.to_dict()
        if any(name in sys.modules for name in FORBIDDEN):
            raise ValueError("forbidden source model import during saved-data export")
        _write_json(stage/"coupled_emission_windows.json",document)
        counts=_export_tables(document,stage)
        _late_recheck(config,source_dir)
        code={"tools/export_coupled_emission_windows.py":sha(Path(__file__)),"mhdlab/coupled_emission_windows.py":sha(ROOT/"mhdlab/coupled_emission_windows.py")}
        report={"schema":REPORT_SCHEMA,"status":"manufactured_saved_source_windows_only","config_original_path":str(config_path.resolve()),"config_sha256":sha(config_path),"code_sha256":code,"source_original_directory":str(source_dir.resolve()),"source_binding":binding,"copied_source_input_sha256":copied,"copied_source_path_base":"directory containing report.json","queries":config["queries"],"quantity_contract":config["quantity_contract"],"representation_contract":config["representation_contract"],"result_file":"coupled_emission_windows.json","result_sha256":sha(stage/"coupled_emission_windows.json"),"record_counts":counts,"cache":document["cache"],"quantity_units":document["quantity_units"],"artifact_path_base":"directory containing report.json","artifact_sha256":{str(p.relative_to(stage)):sha(p) for p in stage.rglob('*') if p.is_file()},"csv_encoding":"Nested keys use __. Exact quantity records retain numeric value plus __exact_numerator,__exact_denominator,__projection_json. Lists are complete JSON cells, booleans use true/false, and unavailable values use literal null with separate availability/reason fields. Polynomial pieces are separate rows; no value is replaced by zero.","qualifications":{"source_evolution_calls":0,"source_RHS_calls":0,"source_assembly_calls":0,"old_source_models_imported":False,"source_input_hashes_rechecked_after_window_calculations":True,"build_call_count":1,"three_source_error_layers_preserved":True,"joint_local_time_temperature_moment_required":True,"new_independent_window_reference_computed_here":False,"actual_physical_or_spectral_claim":False},"limits":["Exact algebra applies to the declared positive polynomial measure, not a continuous exact kinetic solution.","Source RK channel comparisons are unavailable at off-knot bounds; RK has no source time/temperature moment channels.","Polynomial energy defects, source RK discrepancies and independent window-reference errors remain separate; no repairs or automatic source refinement.","Source-emission windows only; no photon gate, experimental clock registration, physical material/water/monolayer, velocity launch, transport, spectrum, MHD, owner or host inference."]}
        _write_json(stage/"report.json",report)
        stage.rename(output)
    except Exception as exc:
        _write_json(stage/"failure.json",{"status":"incomplete_export_not_admitted","error_type":type(exc).__name__,"error":str(exc),"source_evolution_RHS_assembly_calls":0,"no_retry_or_fallback":True})
        raise ValueError("export failed closed; incomplete evidence retained at "+str(stage)+": "+str(exc)) from exc
    print(json.dumps({"report":str(output/"report.json"),"report_sha256":sha(output/"report.json"),"record_counts":counts,"cache":document["cache"],"status":report["status"]},indent=2))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=ROOT/"configs/coupled_emission_windows_manufactured.json")
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    try:
        run(args.config,args.source_dir,args.output_dir)
    except (ValueError,TypeError,OverflowError,OSError) as exc:
        parser.error(str(exc))


if __name__=="__main__":
    main()
