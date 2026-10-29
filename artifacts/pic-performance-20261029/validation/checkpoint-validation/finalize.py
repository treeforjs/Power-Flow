import pathlib,json,hashlib,datetime
out=pathlib.Path(__file__).parent
r=json.loads((out/"comparison.json").read_text());tests=json.loads((out/"fixture-tests.json").read_text())
assert r["status"]==tests["status"]=="passed" and tests["case_count"]>=146
assert r["source_reservoirs"]["all_values_exact"] and r["source_reservoirs"]["pool_metrics"]
assert r["execution"]["cpu"]==6 and r["execution"]["nice"]==15
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
header=pathlib.Path("/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/include/raikou/mfem_gap_pic.h")
guard=pathlib.Path("/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/src/mfem_pic_distributed.hpp")
report={"status":"ready","completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"tool":str(out/"compare_saved_pic.py"),"result":str(out/"comparison.json"),"tests":str(out/"fixture-tests.json"),"notes":str(out/"README.txt"),"actual_comparison_passed":True,"actual_checks":len(r["checks"]),"fixture_cases":tests["case_count"],"cli_exit_codes":[x["returncode"]for x in tests["cli"]],"particle_count":r["particle_counts"]["final"],"max_position_difference_m":r["particles"]["position_m"]["linf_absolute"],"max_proper_velocity_difference_m_s":r["particles"]["proper_velocity_m_s"]["linf_absolute"],"field_maxima":{key:r["frame"][key]["linf_absolute"]for key in ("potential_V","Ex_V_m","Ey_V_m")},"densities_identical":r["density_analysis"]["all_density_components_identical"],"runtime_guards":r["runtime_residual_guards"],"execution":r["execution"],"abi_header_sha256":sha(header),"guard_source_sha256":sha(guard),"source_candidate_or_dof_patch_edited":False,"gpu_work":False,"production_signals_or_edits":False,"limitations":r["limits"]}
(out/"handoff.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
(out/"progress.json").write_text(json.dumps({"status":"complete","handoff":str(out/"handoff.json"),"actual_comparison":"passed","fixture_cases":tests["case_count"],"cpu":6,"nice":15,"lease":r["execution"]["lease"]},indent=2)+"\n")
files=[p for p in out.iterdir() if p.is_file()and p.name!="output-manifest.json"]
(out/"output-manifest.json").write_text(json.dumps({"files":[{"path":str(p.relative_to(out)),"bytes":p.stat().st_size,"sha256":sha(p)}for p in sorted(files)]},indent=2)+"\n")
print(json.dumps({k:v for k,v in report.items()if k not in ("runtime_guards","limitations")},indent=2))
