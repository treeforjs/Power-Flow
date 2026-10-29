import json,pathlib,hashlib,datetime
out=pathlib.Path(__file__).parent
shared=pathlib.Path("/home/treefor/projects/runtimes/pic-cpu-candidate-20261029/source")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before=json.loads((out/"shared-source-before.json").read_text())
after={str(p.relative_to(shared)):sha(p) for p in shared.rglob("*") if p.is_file()}
assert before==after
results=json.loads((out/"results.json").read_text())
leased=json.loads((out/"leased-verification.json").read_text())
assert results["status"]==leased["status"]=="passed" and leased["lease_exclusion_confirmed"]
report={"status":"ready","completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"isolated_source":str(out/"source"),"patch":str(out/"dof-cache-invalidation.patch"),"patch_sha256":sha(out/"dof-cache-invalidation.patch"),"changed_files":["include/raikou/q1_scalar_gather.h","tests/test_q1_scalar_gather.cpp"],"fix":"explicit invalidate() before unversioned DOF renumbering; existing MFEM fallback; reconstruct gather to recache","automatic_unversioned_mutation_detection":False,"cpu":6,"nice":15,"one_job":True,"final_verification_lease":leased["testing_core_lease"],"lease_exclusion_confirmed":True,"final_checks_passed":len(leased["runs"]),"trace_parity_bytes":[455537,328651,31386],"shared_source_unchanged":True,"gpu_work":False,"active_run_modified":False,"deployment_performed":False,"reservation_changes":False,"results":str(out/"results.json"),"leased_verification":str(out/"leased-verification.json"),"notes":str(out/"README.txt")}
(out/"handoff.json").write_text(json.dumps(report,indent=2)+"\n")
(out/"progress.json").write_text(json.dumps({"status":"complete","handoff":str(out/"handoff.json"),"cpu":6,"nice":15,"leased_final_checks_passed":6,"shared_candidate_unchanged":True,"gpu_work":False},indent=2)+"\n")
paths=[p for p in out.rglob("*") if p.is_file() and p.name!="output-manifest.json"]
(out/"output-manifest.json").write_text(json.dumps({"files":[{"path":str(p.relative_to(out)),"bytes":p.stat().st_size,"sha256":sha(p)}for p in sorted(paths)]},indent=2)+"\n")
print(json.dumps(report,indent=2))
