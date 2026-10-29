import pathlib,json,datetime
import snapshot_occupancy as diag
import test_occupancy
out=pathlib.Path(__file__).parent
root=out.parent
with diag.state_io.cpu_lease()as execution:
 tests=test_occupancy.run_tests();summary={}
 for source,name in (("seed-current","seed"),("baseline-cpu6-run1","baseline-final"),("cpu-candidate-cpu6-run1","candidate-final")):
  result=diag.from_directory(root/source);result["execution"]=execution
  summary[name]={k:v for k,v in result.items()if k in ("global_time_s","input_checkpoint","input_checkpoint_sha256","decoder_sha256","global","crossing_flux_computed","residence_time_weighted")}
  diag.write_result(result,out/"snapshots"/name)
 report={"status":"passed","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"fixture_cases":tests["cases"],"execution":execution,"snapshots":summary,"independent_physical_times":2,"warning":"Baseline and candidate finals represent the same physical endpoint; do not pool them as extra samples or time."}
 (out/"snapshot-audit.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
 print(json.dumps({"status":"passed","fixture_cases":tests["cases"],"snapshots":[{"name":name,"time_s":r["global_time_s"],"full_weight":r["global"]["full_height"]["weighted_electrons"],"mean_energy_eV":r["global"]["full_height"]["mean_energy_eV"],"max_energy_eV":r["global"]["full_height"]["max_energy_eV"]}for name,r in summary.items()]},indent=2))
