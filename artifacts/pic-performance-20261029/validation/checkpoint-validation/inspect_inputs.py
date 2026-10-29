import pathlib,json,fcntl,os,csv
import compare_saved_pic as c
out=pathlib.Path(__file__).parent
with c.LEASE.open("a+") as lease:
 fcntl.flock(lease,fcntl.LOCK_EX)
 s=c.schema();root=out.parent.parent
 states=[c.read_state(root/name,s)for name in ["baseline-cpu6-run1","cpu-candidate-cpu6-run1","seed-current"]]
 report={"states":{p.directory.name:{"options":p.options,"ledger":p.ledger,"public":p.public(),"n":len(p.particles),"next_id":p.next_id,"step":p.step,"pool_lengths":[len(x["integral"])for x in p.pools]}for p in states}}
 frames=[c.read_csv(p.directory/("frame_"+str(p.next_frame)+".csv"),c.FRAME)for p in states[:2]]
 report["frame"]={k:c.vector_metric([r[k]for r in frames[0]],[r[k]for r in frames[1]])for k in c.FRAME}
 (out/"inspection.json").write_text(json.dumps(report,indent=2)+"\n")
 print(json.dumps({"abi":s,"states":report["states"],"frame":report["frame"]},indent=2))
