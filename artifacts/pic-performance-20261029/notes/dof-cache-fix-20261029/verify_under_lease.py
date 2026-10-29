import pathlib,subprocess,json,os,hashlib,datetime,fcntl
out=pathlib.Path(__file__).parent
logs=out/"leased-logs";logs.mkdir(exist_ok=True)
launcher="/home/treefor/projects/runtimes/arch-base-20261028/run-frozen"
lock=out.parent.parent/"testing-core.lock"
with lock.open("a+") as probe:
    try:fcntl.flock(probe,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:lease_exclusion=True
    else:raise RuntimeError("verification must execute beneath flock holding the testing-core lease")
assert sorted(os.sched_getaffinity(0))==[6] and os.getpriority(os.PRIO_PROCESS,0)==15
runs=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(name,args,reference=None):
    cmd=[launcher]+[str(a) for a in args]
    log=logs/(name+".log")
    with log.open("wb") as f:subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT,check=True)
    record={"name":name,"status":"passed","argv":cmd,"log":str(log),"sha256":sha(log)}
    if reference:
        assert log.read_bytes()==reference.read_bytes(),name+" trace changed"
        record.update({"byte_identical":True,"bytes":log.stat().st_size,"reference":str(reference)})
    runs.append(record);print(name+": passed",flush=True)
run("q1_gather_tests",[out/"build/q1_gather_tests"])
run("pic_foundation_tests",[out/"build/pic_foundation_tests"])
run("cache_edges",[out/"tests/cache_edges"])
run("pic_cpu_trace",[out/"build/pic_cpu_trace"],out.parent.parent/"cpu-audit/baseline-trace.txt")
run("serial_state_trace",[out/"tests/serial_state_trace",out/"tests/leased-checkpoint"],out.parent/"tests/baseline_serial_state_trace.txt")
run("native_field_trace",[out/"tests/field_cache_trace",out/"tests/leased-native-fixture",out/"tests/leased-native-timing.json","isolated_fix"],out.parent.parent/"cpu-audit/baseline-field-trace.txt")
assert json.loads((out/"tests/leased-native-timing.json").read_text())["checksum"]==-4019096847.1793637
report={"status":"passed","completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"affinity":sorted(os.sched_getaffinity(0)),"nice":os.getpriority(os.PRIO_PROCESS,0),"thread_environment":{k:os.environ.get(k) for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]},"testing_core_lease":str(lock),"lease_exclusion_confirmed":lease_exclusion,"one_job":True,"gpu_work":False,"builds_performed_in_final_check":False,"runs":runs}
(out/"leased-verification.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps({k:v for k,v in report.items() if k!="runs"},indent=2),flush=True)
