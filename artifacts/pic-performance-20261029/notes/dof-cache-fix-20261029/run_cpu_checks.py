import subprocess,pathlib,os,json,hashlib,datetime,shutil,time,fcntl
out=pathlib.Path(__file__).parent
review=out.parent
audit=review.parent/"cpu-audit"
rt=pathlib.Path("/home/treefor/projects/runtimes")
mf=rt/"mfem-4.10-pic-20261029"
gpu=rt/"mfem-cuda-benchmark-20261029/cuda-libs/nvidia"
source=out/"source";build=out/"build";tests=out/"tests";logs=out/"logs"
for p in (tests,logs):p.mkdir(exist_ok=True)
launcher=str(rt/"arch-base-20261028/run-frozen")
commands=[];results=[]
lease_path=review.parent/"testing-core.lock"
lease=lease_path.open("a+")
print("Waiting for testing-core lease",flush=True)
fcntl.flock(lease,fcntl.LOCK_EX)
print("Testing-core lease acquired",flush=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
execution={"affinity":sorted(os.sched_getaffinity(0)),"nice":os.getpriority(os.PRIO_PROCESS,0),"thread_environment":{k:os.environ.get(k) for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]},"compiler_jobs":1,"testing_core_lease":str(lease_path),"lease_held":True,"CUDA_VISIBLE_DEVICES":os.environ.get("CUDA_VISIBLE_DEVICES"),"gpu_work":False}
assert execution["affinity"]==[6] and execution["nice"]==15

def progress(status):
    (out/"progress.json").write_text(json.dumps({"status":status,"updated_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"execution":execution,"tests_completed":results},indent=2)+"\n")

def run(name,args,expected=0):
    cmd=[launcher]+[str(a) for a in args]
    record={"name":name,"argv":cmd,"log":str(logs/(name+".log")),"inherited_execution":execution}
    commands.append(record);progress(name)
    start=time.monotonic()
    with pathlib.Path(record["log"]).open("wb") as log:
        proc=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT)
    record.update({"returncode":proc.returncode,"wall_seconds":time.monotonic()-start})
    (out/"commands.json").write_text(json.dumps(commands,indent=2)+"\n")
    if proc.returncode!=expected:raise RuntimeError(name+" failed; see "+record["log"])
    print(name+": passed",flush=True)
    return pathlib.Path(record["log"])

def parity(name,a,b):
    aa=a.read_bytes();bb=b.read_bytes()
    result={"name":name,"byte_identical":aa==bb,"bytes":len(aa),"actual_sha256":sha(a),"reference_sha256":sha(b),"actual":str(a),"reference":str(b)}
    results.append(result)
    if aa!=bb:raise RuntimeError(name+" differs")
    print(name+": byte-identical ("+str(len(aa))+" bytes)",flush=True)

try:
    run("original_reproducer",[review/"tests/cache_edges"])
    text=(logs/"original_reproducer.log").read_text()
    assert "value_mismatches=6/6" in text and "available=1" in text
    results.append({"name":"original_bug_reproduced","mismatches":"6/6","cache_available":True})
    run("configure",["/usr/sbin/cmake","-S",source,"-B",build,"-DCMAKE_BUILD_TYPE=Release","-DMFEM_DIR="+str(mf/"install-mpi/lib/cmake/mfem"),"-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"])
    run("build_cpu_targets",["/usr/sbin/cmake","--build",build,"--parallel","1","--target","q1_gather_tests","pic_foundation_tests","pic_cpu_trace"])
    run("q1_gather_tests",[build/"q1_gather_tests"])
    results.append({"name":"q1_gather_tests","status":"passed","coverage":"renumbering invalidation, idempotence, reconstruction, values and gradients at corners/interiors; original Q1, Q2, triangles, refinement fallbacks"})
    run("pic_foundation_tests",[build/"pic_foundation_tests"])
    results.append({"name":"pic_foundation_tests","status":"passed","solver":"serial CG; distributed_fields=false"})
    trace=run("pic_cpu_trace",[build/"pic_cpu_trace"])
    parity("supplied_serial_trace",trace,audit/"baseline-trace.txt")
    common=["/usr/sbin/c++","-O2","-DNDEBUG","-std=c++20","-I"+str(source/"include"),"-I"+str(source/"src"),"-isystem",str(mf/"install-mpi/include"),"-isystem",str(mf/"tpl/hypre/include"),"-isystem",str(mf/"tpl/metis/include")]
    libs=[str(build/"libpic.a"),str(mf/"install-mpi/lib/libmfem.a"),str(mf/"tpl/hypre/lib/libHYPRE.a"),str(mf/"tpl/metis/lib/libmetis.a"),str(mf/"tpl/GKlib/lib/libGKlib.a"),"/usr/lib/libmpi.so","-Wl,-rpath,/usr/lib",str(gpu/"cu12/lib/libcudss.so.0"),str(gpu/"cuda_runtime/lib/libcudart.so.12"),"-Wl,-rpath,"+":".join(str(gpu/p) for p in ["cu12/lib","cublas/lib","cuda_runtime/lib"])]
    shutil.copy2(review/"tests/serial_state_trace.cpp",tests/"serial_state_trace.cpp")
    run("compile_serial_state_trace",common+[str(tests/"serial_state_trace.cpp"),"-o",str(tests/"serial_state_trace")]+libs)
    trace=run("serial_state_trace",[tests/"serial_state_trace",tests/"checkpoint"])
    parity("complete_serial_state_trace",trace,review/"tests/baseline_serial_state_trace.txt")
    s=(review/"tests/cache_edges.cpp").read_text()
    s=s.replace("ss=s2.GetSequence();s2.ReorderElementToDofTable();","ss=s2.GetSequence();hg.invalidate();s2.ReorderElementToDofTable();")
    s=s.replace("if(!mismatches)throw std::runtime_error(\"expected reorder mismatch was absent\");","if(mismatches || hg.available())throw std::runtime_error(\"invalidated cache did not fall back\");\n    compare(h,hg);Q1ScalarGather rebuilt(h);if(!rebuilt.available())throw std::runtime_error(\"rebuild unavailable\");compare(h,rebuilt);\n    std::cout<<\"DOF reorder fallback and rebuilt cache: value/gradient bitwise pass\\n\";")
    assert "hg.invalidate();" in s and "if(!mismatches)" not in s
    (tests/"cache_edges.cpp").write_text(s)
    run("compile_cache_edges",common+[str(tests/"cache_edges.cpp"),"-o",str(tests/"cache_edges")]+libs)
    run("cache_edges",[tests/"cache_edges"])
    results.append({"name":"cache_edges","status":"passed","coverage":"fresh coefficients; geometry mutation; signed DOFs; custom nonidentity transform; Q2 rebind; explicit invalidation before unversioned renumbering and reconstructed cache"})
    shutil.copy2(audit/"field_cache_trace.cpp",tests/"field_cache_trace.cpp")
    history="-DHISTORY_SOURCE=\""+str(source/"examples/mfem_heating_history_pic.cpp")+"\""
    run("compile_native_field_trace",common+[history,str(tests/"field_cache_trace.cpp"),"-o",str(tests/"field_cache_trace")]+libs)
    trace=run("native_field_trace",[tests/"field_cache_trace",tests/"native-fixture",tests/"native-timing.json","isolated_fix"])
    parity("native_field_time_cache_trace",trace,audit/"baseline-field-trace.txt")
    results.append({"name":"native_gather_checksum","actual":json.loads((tests/"native-timing.json").read_text())["checksum"],"expected":-4019096847.1793637,"equal":json.loads((tests/"native-timing.json").read_text())["checksum"]==-4019096847.1793637})
    assert results[-1]["equal"]
    shared=rt/"pic-cpu-candidate-20261029/source"
    before=json.loads((out/"shared-source-before.json").read_text())
    after={str(p.relative_to(shared)):sha(p) for p in shared.rglob("*") if p.is_file()}
    assert before==after
    changed=[str(p.relative_to(source)) for p in source.rglob("*") if p.is_file() and before.get(str(p.relative_to(source)))!=sha(p)]
    assert set(changed)=={"include/raikou/q1_scalar_gather.h","tests/test_q1_scalar_gather.cpp"}
    report={"status":"passed","completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"execution":execution,"results":results,"shared_source_unchanged":True,"isolated_changed_files":sorted(changed),"patch_sha256":sha(out/"dof-cache-invalidation.patch"),"fixed_library_sha256":sha(build/"libpic.a"),"limitations":["Unversioned numbering operations require explicit invalidate() before mutation; no automatic detection added.","GPU/distributed residual branch and end-to-end performance not executed.","No shared candidate edits, deployment, active-run operations, or reservation changes."]}
    (out/"results.json").write_text(json.dumps(report,indent=2)+"\n")
    progress("complete")
    print(json.dumps(report,indent=2),flush=True)
except Exception as exc:
    progress("failed: "+str(exc))
    raise
