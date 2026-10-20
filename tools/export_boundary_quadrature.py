#!/usr/bin/env python3
"""Export a finite positive quadrature of accepted manufactured joint source pieces.

Read only hash-bound source/acceptance bytes. Query and input checks precede the
new pure polynomial/quadrature import. No source, RHS, assembly, old window
export or independent reference driver is imported or executed.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CONFIG_SCHEMA = "mhdlab-boundary-quadrature-manufactured-config-v1"
REPORT_SCHEMA = "mhdlab-boundary-quadrature-export-run-v1"
ALIASES = {"report", "windows", "input_config", "window_config", "independent_completion", "independent_report", "scientific_review"}
HEX = re.compile(r"[0-9a-f]{64}\Z")
INTEGER = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")
MAX_CONFIG_BYTES = 65536
MAX_FILE_BYTES = 64*1024*1024
MAX_TOTAL_BYTES = 128*1024*1024
FORBIDDEN = ("mhdlab.coupled_emission_windows", "mhdlab.thermal_width_power_reference", "mhdlab.thermal_forced_feedback_reference", "mhdlab.thermal_power_reference", "mhdlab.thermal_feedback_reference")
METHOD = {"rule":"positive_Gauss_Legendre_5_on_unit_interval", "decimal_precision":80, "nodes_per_positive_support":5, "polynomial_evaluation":"exact_rational_after_finite_decimal_node_generation", "export_projection":"binary64_fields_once_then_exact_fraction_defects", "weight_repair":"none", "cache_identity":"canonical_source_history_and_exact_restricted_piece_bounds"}
COUNTS = {"unique_piece_supports":66,"unique_node_evaluations":330,"source_query_node_references":570,"patch_node_references":1140,"source_queries":6,"aggregate_queries":4,"contrasts":2,"source_evolution_RHS_assembly_window_export_reference_driver_calls":0}
QUANTITY = {"local_time_anchor":"exact_query_start","joint_moments":["1","local_time","local_time_squared","temperature","temperature_squared","local_time_temperature"],"energy_reference":"exact_SI_constants_and_supplied_D_I_separate_from_activation_barrier","geometry":"retain_patch_range_and_apply_actual_area_once","zero_mass_normalization":"unavailable_with_explicit_reason","new_query_IVP_error":"unavailable_no_new_reference_driver"}
DEFECT = {"layers":["approximate_rule_minus_exact_polynomial","exported_binary_fields_minus_approximate_rule","exported_binary_fields_minus_exact_polynomial"],"closure":"total_equals_rule_plus_projection_exactly","inherited_errors":"source_interpolation_and_prior_query_IVP_errors_retained_separately","normalization":"never_renormalize_weights"}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return digest(Path(path).read_bytes())


def _constant(token):
    raise ValueError("nonfinite JSON constant: "+token)


def _decimal(token):
    if len(token)>256:
        raise ValueError("bounded numeric lexeme required")
    try:
        value=Decimal(token)
        projected=float(value)
    except (InvalidOperation, OverflowError, ValueError) as exc:
        raise ValueError("unsupported JSON numeric token") from exc
    if (not value.is_finite() or not math.isfinite(projected)
            or value!=0 and projected==0 or value==0 and value.is_signed()):
        raise ValueError("JSON numeric token has unsupported finite/nonzero/signed-zero representation")
    return value


def _integer(token):
    if len(token)>128 or token=="-0":
        raise ValueError("bounded canonical JSON integer required")
    return int(token)


def _object(pairs):
    value={}
    for key,item in pairs:
        if key in value:
            raise ValueError("duplicate JSON key: "+key)
        value[key]=item
    return value


def _text(value,name):
    if type(value) is not str or not value.strip() or len(value)>4096:
        raise ValueError(name+" must be bounded nonempty supplied text")


def _sha(value,name):
    if type(value) is not str or not HEX.fullmatch(value):
        raise ValueError(name+" must be an explicit lower-case SHA256")


def _relative(value):
    if (type(value) is not str or not value or len(value)>512 or "\\" in value or "\x00" in value
            or PurePosixPath(value).is_absolute() or any(p in ("", ".", "..") for p in value.split("/"))
            or str(PurePosixPath(value))!=value):
        raise ValueError("source paths must be unambiguous safe relative POSIX paths")
    return value


def _fraction(record):
    if type(record) is not dict or set(record)!={"numerator","denominator"}:
        raise ValueError("query fractions require explicit canonical integer strings")
    for key in ("numerator","denominator"):
        token=record[key]
        if type(token) is not str or len(token)>128 or not INTEGER.fullmatch(token) or str(int(token))!=token:
            raise ValueError("query rational integers must be bounded canonical strings")
    n,d=int(record["numerator"]),int(record["denominator"])
    if d<=0:
        raise ValueError("query denominator must be positive")
    value=Fraction(n,d)
    if (value.numerator,value.denominator)!=(n,d):
        raise ValueError("query fractions must be reduced")
    return value


def load_config(path,*,raw_bytes=None):
    data=Path(path).read_bytes() if raw_bytes is None else raw_bytes
    if type(data) is not bytes:
        raise ValueError("configuration must be immutable raw bytes")
    if len(data)>MAX_CONFIG_BYTES:
        raise ValueError("configuration exceeds bounded byte limit")
    config=json.loads(data.decode("utf-8"),parse_float=_decimal,parse_int=_integer,parse_constant=_constant,object_pairs_hook=_object)
    keys={"schema","parameter_status","provenance","scope","source_binding","method","queries","quantity_contract","defect_contract","declared_counts"}
    if type(config) is not dict or set(config)!=keys or config["schema"]!=CONFIG_SCHEMA or config["parameter_status"]!="manufactured":
        raise ValueError("unsupported explicit manufactured quadrature configuration")
    for key in ("provenance","scope"):
        _text(config[key],key)
    binding=config["source_binding"]
    expected={"path_base","files","accepted_source_report_sha256","accepted_window_result_sha256","original_width_config_sha256","original_window_config_sha256","numerical_run_cli_sha256","plot_only_public_cli_sha256"}
    if type(binding) is not dict or set(binding)!=expected:
        raise ValueError("complete accepted saved-source binding is required")
    _text(binding["path_base"],"source path base")
    files=binding["files"]
    if type(files) is not dict or set(files)!=ALIASES:
        raise ValueError("all seven fixed saved-source aliases are required")
    seen=set()
    for alias,record in files.items():
        if type(record) is not dict or set(record)!={"relative_path","sha256"}:
            raise ValueError("each source alias requires path and expected hash")
        relative=_relative(record["relative_path"])
        if relative in seen:
            raise ValueError("source aliases must identify distinct paths")
        seen.add(relative);_sha(record["sha256"],alias)
    for key in expected-{"path_base","files"}:
        _sha(binding[key],key)
    for key,alias in (("accepted_source_report_sha256","report"),("accepted_window_result_sha256","windows"),("original_width_config_sha256","input_config"),("original_window_config_sha256","window_config")):
        if binding[key]!=files[alias]["sha256"]:
            raise ValueError("declared source identity differs from file binding: "+key)
    method=config["method"]
    if type(method) is not dict or method!=METHOD or type(method.get("decimal_precision")) is not int or type(method.get("nodes_per_positive_support")) is not int:
        raise ValueError("only the declared finite Decimal80 positive GL5 method is supported")
    queries=config["queries"]
    expected_queries={"shifted_interior":(Fraction(3,32),Fraction(11,32)),"quarter_cell_interior":(Fraction(13,256),Fraction(99,256))}
    if type(queries) is not list or len(queries)!=2:
        raise ValueError("exactly two new declared source-emission queries are required")
    ids=set()
    for query in queries:
        if type(query) is not dict or set(query)!={"window_id","coordinate","start_fraction","stop_fraction","parameter_provenance"}:
            raise ValueError("complete query identity/support/provenance is required")
        key=query["window_id"]
        if type(key) is not str or key not in expected_queries or key in ids:
            raise ValueError("query ids must be unique and explicitly declared")
        ids.add(key)
        if query["coordinate"]!="canonical_duration_fraction":
            raise ValueError("queries require exact canonical source-duration coordinates")
        a,b=_fraction(query["start_fraction"]),_fraction(query["stop_fraction"])
        if not 0<=a<b<=1 or (a,b)!=expected_queries[key]:
            raise ValueError("query support differs from the two bounded declared windows")
        _text(query["parameter_provenance"],"query provenance")
    if type(config["quantity_contract"]) is not dict or config["quantity_contract"]!=QUANTITY:
        raise ValueError("complete joint moment, energy and geometry contract required")
    if type(config["defect_contract"]) is not dict or config["defect_contract"]!=DEFECT:
        raise ValueError("all three defect layers and no-normalization contract required")
    if type(config["declared_counts"]) is not dict or config["declared_counts"]!=COUNTS or any(type(v) is not int for v in config["declared_counts"].values()):
        raise ValueError("explicit bounded calculation/reuse counts required")
    return config


def _read_relative(base,relative):
    relative=_relative(relative);path=base
    for part in relative.split("/"):
        path=path/part
        if path.is_symlink():
            raise ValueError("source symlinks are unsupported")
    if not path.is_file() or not path.resolve().is_relative_to(base) or path.stat().st_size>MAX_FILE_BYTES:
        raise ValueError("source file is missing, external or exceeds bounded byte limit")
    data=path.read_bytes()
    if len(data)>MAX_FILE_BYTES:
        raise ValueError("source file exceeds bounded byte limit after read")
    return data


def _source_float(token):
    value=float(token)
    if not math.isfinite(value) or Decimal(token)!=0 and value==0:
        raise ValueError("saved source numeric token has unsupported representation")
    return value


def _source_json(data):
    value=json.loads(data.decode("utf-8"),parse_float=_source_float,parse_int=_integer,parse_constant=_constant,object_pairs_hook=_object)
    if type(value) is not dict:
        raise ValueError("saved source JSON must have an object root")
    return value


def read_bound_source(config,source_dir):
    """Bind raw bytes and accepted source identity before importing pure math."""
    base=source_dir.resolve(strict=True)
    if not base.is_dir():
        raise ValueError("source-dir must be the accepted source pass directory")
    files={};total=0;binding=config["source_binding"]
    for alias,spec in binding["files"].items():
        data=_read_relative(base,spec["relative_path"]);total+=len(data)
        if total>MAX_TOTAL_BYTES or digest(data)!=spec["sha256"]:
            raise ValueError("source byte hash or bounded total differs: "+alias)
        files[alias]=data
    report=_source_json(files["report"])
    if report.get("schema")!="mhdlab-coupled-emission-windows-export-run-v1" or report.get("status")!="manufactured_saved_source_windows_only":
        raise ValueError("unsupported accepted source report schema/status")
    if report.get("artifact_path_base")!="directory containing report.json":
        raise ValueError("unsupported source report artifact base")
    if report.get("result_sha256")!=digest(files["windows"]) or report.get("artifact_sha256",{}).get(report.get("result_file"))!=digest(files["windows"]):
        raise ValueError("saved polynomial result/report association differs")
    if report.get("config_sha256")!=digest(files["window_config"]):
        raise ValueError("source window configuration association differs")
    width=report.get("source_binding",{})
    if width.get("files",{}).get("input_config",{}).get("sha256")!=digest(files["input_config"]):
        raise ValueError("original supplied D/I configuration association differs")
    lineage=width.get("lineage",{})
    if lineage.get("numerical_run_cli_sha256")!=binding["numerical_run_cli_sha256"] or lineage.get("plot_only_public_cli_sha256")!=binding["plot_only_public_cli_sha256"]:
        raise ValueError("original numerical-versus-plot code lineage differs")
    accepted=_source_json(files["independent_completion"])
    if accepted.get("schema")!="mitl-independent-coupled-window-completion-v1" or accepted.get("status")!="passed" or accepted.get("all_checks_passed") is not True:
        raise ValueError("accepted source verification is absent")
    if accepted.get("producer_report_sha256")!=digest(files["report"]) or accepted.get("report_sha256")!=digest(files["independent_report"]) or accepted.get("final_static_review_sha256")!=digest(files["scientific_review"]):
        raise ValueError("accepted source verification/report/review association differs")
    independent=_source_json(files["independent_report"])
    if independent.get("schema")!="mitl-independent-coupled-window-query-report-v1" or independent.get("status")!="passed" or independent.get("all_checks_passed") is not True or independent.get("window_document_sha256")!=digest(files["windows"]):
        raise ValueError("accepted independent report/window association differs")
    review=_source_json(files["scientific_review"])
    if review.get("schema")!="mitl-coupled-window-static-oracle-review-v1" or review.get("status")!="accepted" or review.get("producer_report_sha256")!=digest(files["report"]) or review.get("window_document_sha256")!=digest(files["windows"]) or review.get("actual_report_sha256")!=digest(files["independent_report"]):
        raise ValueError("source scientific review/identity association differs")
    return files


def _late_recheck(config,source_dir):
    base=source_dir.resolve(strict=True)
    for alias,spec in config["source_binding"]["files"].items():
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
    expected={"piece_cache":66,"node_cache":330,"source_queries":6,"patch_queries":12,"aggregate_queries":4,"contrasts":2}
    tables={};references={"source_query_nodes":[],"patch_query_nodes":[]}
    for name,count in expected.items():
        records=document[name]
        if type(records) is not list or len(records)!=count:
            raise ValueError("bounded output record count differs: "+name)
        tables[name]=[_flatten(item) for item in records]
        _write_csv(output/(name+".csv"),tables[name])
    for name,destination in (("source_queries","source_query_nodes"),("patch_queries","patch_query_nodes")):
        for parent_index,item in enumerate(document[name]):
            nodes=item["summary"]["nodes"]
            if type(nodes) is not list:
                raise ValueError("actual quadrature nodes must be explicit records")
            identity={key:item[key] for key in ("source_history_id","query_id","state_id","hypothesis","patch_index") if key in item}
            for index,node in enumerate(nodes):
                references[destination].append(_flatten({"parent_query_index":parent_index,"node_reference_index":index,**identity,"node":node}))
    for name,count in (("source_query_nodes",570),("patch_query_nodes",1140)):
        if len(references[name])!=count:
            raise ValueError("bounded node-reference count differs: "+name)
        _write_csv(output/(name+".csv"),references[name])
    return {**{name:len(rows) for name,rows in tables.items()},**{name:len(rows) for name,rows in references.items()}}


def _load_core():
    """Load only the new stdlib-only core, bypassing the legacy package initializer."""
    name="_boundary_quadrature_saved_export_core"
    if name in sys.modules:
        raise ValueError("dedicated quadrature module identity already exists")
    spec=importlib.util.spec_from_file_location(name,ROOT/"mhdlab/boundary_quadrature.py")
    if spec is None or spec.loader is None:
        raise ValueError("new quadrature core is unavailable")
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        del sys.modules[name]
        raise
    return module


def run(config_path,source_dir,output):
    config_bytes=config_path.read_bytes()
    config=load_config(config_path,raw_bytes=config_bytes);files=read_bound_source(config,source_dir)
    if output.exists():
        raise ValueError("output must be a fresh directory")
    core=_load_core()
    source=core.admit_boundary_source(files,expected_sha256={k:v["sha256"] for k,v in config["source_binding"]["files"].items()},parameter_provenance=config["provenance"])
    queries=tuple(core.BoundaryQuery(window_id=q["window_id"],start_fraction=_fraction(q["start_fraction"]),stop_fraction=_fraction(q["stop_fraction"]),parameter_provenance=q["parameter_provenance"]) for q in config["queries"])
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix=output.name+".staging-",dir=output.parent))
    try:
        (stage/"input_config.json").write_bytes(config_bytes);copied={}
        for alias,data in files.items():
            relative="source-inputs/"+config["source_binding"]["files"][alias]["relative_path"]
            path=stage/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data);copied[relative]=digest(data)
        result=core.build_boundary_quadrature(source,queries,parameter_status=config["parameter_status"],parameter_provenance=config["provenance"])
        document=result.to_dict()
        if any(name in sys.modules for name in FORBIDDEN) or any(name=="mhdlab" or name.startswith("mhdlab.") for name in sys.modules):
            raise ValueError("forbidden previous source/math module import during saved-data export")
        _write_json(stage/"boundary_quadrature.json",document);counts=_export_tables(document,stage)
        _late_recheck(config,source_dir)
        if config_path.read_bytes()!=config_bytes:
            raise ValueError("late configuration byte identity changed")
        code={"tools/export_boundary_quadrature.py":sha(Path(__file__)),"mhdlab/boundary_quadrature.py":sha(ROOT/"mhdlab/boundary_quadrature.py")}
        report={"schema":REPORT_SCHEMA,"status":"manufactured_saved_source_boundary_quadrature_only","config_original_path":str(config_path.resolve()),"config_sha256":digest(config_bytes),"code_sha256":code,"source_original_directory":str(source_dir.resolve()),"source_binding":config["source_binding"],"copied_source_input_sha256":copied,"copied_source_path_base":"directory containing report.json","queries":config["queries"],"method":config["method"],"quantity_contract":config["quantity_contract"],"defect_contract":config["defect_contract"],"result_file":"boundary_quadrature.json","result_sha256":sha(stage/"boundary_quadrature.json"),"record_counts":counts,"cache":document["cache"],"quantity_units":document["quantity_units"],"artifact_path_base":"directory containing report.json","artifact_sha256":{str(p.relative_to(stage)):sha(p) for p in stage.rglob('*') if p.is_file()},"csv_encoding":"Nested keys use __. Exact quantity records retain numeric value plus __exact_numerator,__exact_denominator,__projection_json. Lists are complete JSON cells, booleans true/false and unavailable values literal null with separate availability/reason fields. No unavailable field is replaced by zero.","qualifications":{"build_call_count":1,"old_source_evolution_RHS_assembly_window_export_reference_driver_calls":0,"forbidden_previous_math_or_source_modules_imported":False,"core_import_convention":"dedicated_direct_file_module_bypasses_legacy_package_initializer","source_input_hashes_rechecked_after_calculation":True,"finite_Decimal_nodes_are_not_exact_algebraic_GL_nodes":True,"weights_renormalized":False,"three_new_defect_layers_retained":True,"inherited_source_and_IVP_errors_separate":True,"new_query_IVP_reference_computed_here":False,"position_or_velocity_sampling":False,"actual_physical_or_spectral_claim":False},"limits":["Deterministic positive quadrature of the saved polynomial source measure; finite rule generation and exported-field projection have separate exact defects.","No mass, velocity, position sampling, flight, gap, reaction, radiation or collection kernel is assigned.","Prior accepted queries do not certify new-query IVP errors; source and interpolation qualifications remain separate.","No physical water, monolayer, material, Hall/FLASH/MHD, host-budget or owner implementation inference."]}
        _write_json(stage/"report.json",report);stage.rename(output)
    except Exception as exc:
        _write_json(stage/"failure.json",{"status":"incomplete_export_not_admitted","error_type":type(exc).__name__,"error":str(exc),"old_source_evolution_RHS_assembly_window_export_reference_driver_calls":0,"no_retry_or_fallback":True})
        raise ValueError("export failed closed; incomplete evidence retained at "+str(stage)+": "+str(exc)) from exc
    print(json.dumps({"report":str(output/"report.json"),"report_sha256":sha(output/"report.json"),"record_counts":counts,"cache":document["cache"],"status":report["status"]},indent=2))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=ROOT/"configs/boundary_quadrature_manufactured.json")
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    try:
        run(args.config,args.source_dir,args.output_dir)
    except (ValueError,TypeError,OverflowError,OSError) as exc:
        parser.error(str(exc))


if __name__=="__main__":
    main()
