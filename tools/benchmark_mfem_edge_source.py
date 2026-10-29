"""Preflight or bounded startup benchmark for the assembled serial MFEM driver."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def cost(configuration):
    nx,ny = configuration["width_cells"],configuration["thickness_cells"]
    ex = nx+2*configuration["outside_cells"]
    ey = 2*ny+configuration["gap_cells"]+2*configuration["outside_cells"]
    def grid(x,y):
        dofs = (x+1)*(y+1)
        nnz = (3*x+1)*(3*y+1)
        return {"cells":x*y,"h1_dofs":dofs,"q1_full_stencil_nonzeros":nnz,
                "csr_matrix_bytes":12*nnz+4*(dofs+1)}
    em,thermal = grid(ex,ey),grid(nx,ny)
    matrices = 3*em["csr_matrix_bytes"]+5*thermal["csr_matrix_bytes"]
    vectors = 8*(12*em["h1_dofs"]+8*thermal["h1_dofs"])
    return {"electromagnetic_grid":dict(width_cells=ex,normal_cells=ey,**em),
        "each_of_two_thermal_grids":thermal,
        "retained_and_working_csr_matrix_bytes":matrices,"field_and_work_vector_bytes":vectors,
        "conservative_matrix_vector_storage_estimate_bytes":matrices+vectors,
        "conservative_matrix_vector_storage_estimate_gib":(matrices+vectors)/2**30,
        "assumptions":["double values (8 bytes), int CSR indices (4 bytes), full Q1 nine-point stencil",
            "three electromagnetic CSR matrices: conductivity mass, magnetic stiffness, step operator",
            "five thermal CSR matrices: mass/stiffness for both conductors and one active step operator",
            "twelve electromagnetic and eight one-conductor-sized field/work vectors",
            "mesh objects/connectivity, MFEM assembly temporaries, Hypre copies/AMG hierarchy, solver vectors beyond this allowance, Python and file buffers are excluded"],
        "is_total_peak_memory_estimate":False,
        "field_snapshot_rows":len(configuration.get("field_snapshot_times_s",configuration["required_snapshot_times_s"]))*nx*ny}


def descendants(pid):
    pending,result = [pid],[]
    while pending:
        current = pending.pop()
        try:
            children = [int(value) for value in Path(f"/proc/{current}/task/{current}/children").read_text().split()]
        except (OSError,ValueError):
            children = []
        result.extend(children);pending.extend(children)
    return result


def linux_worker(args):
    import resource
    report_path = args.output_dir/"benchmark.json"
    report = json.loads(report_path.read_bytes())
    memory_bytes = int(args.memory_limit_gib*2**30)
    def limit_child():
        resource.setrlimit(resource.RLIMIT_AS,(memory_bytes,memory_bytes))
        resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    argv = [sys.executable,str(ROOT/"tools/run_mfem_edge_source.py"),"--config",str(args.config.resolve()),
            "--output-dir",str(args.output_dir.resolve()/"native_run")]
    started = time.perf_counter()
    native_seen = None
    timed_out = False
    with (args.output_dir/"runner.stdout.txt").open("w") as stdout,(args.output_dir/"runner.stderr.txt").open("w") as stderr:
        process = subprocess.Popen(argv,cwd=ROOT,stdout=stdout,stderr=stderr,start_new_session=True,preexec_fn=limit_child,
                                   env=dict(os.environ,OMP_NUM_THREADS="2"))
        report["linux_worker"] = {"pid":os.getpid(),"runner_pid":process.pid,"argv":argv,
            "address_space_limit_bytes":memory_bytes,"wall_time_limit_s":args.timeout_s}
        write(report_path,report)
        while process.poll() is None:
            elapsed = time.perf_counter()-started
            for pid in descendants(process.pid):
                try:
                    executable = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")[0].decode()
                except (OSError,UnicodeError):
                    continue
                if Path(executable).name=="power-flow-edge-source" and native_seen is None:
                    native_seen = elapsed
            if elapsed>=args.timeout_s:
                timed_out = True;os.killpg(process.pid,signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL);process.wait()
                break
            time.sleep(.2)
    source_path = args.output_dir/"native_run/source_run.json"
    source = json.loads(source_path.read_bytes()) if source_path.exists() else None
    commands = source.get("commands",[]) if source else []
    native = next((row for row in commands if row["stdout_file"]=="native.stdout.txt"),None)
    report["linux_worker"].update(wall_time_s=time.perf_counter()-started,exit_code=process.returncode,
        timed_out=timed_out,first_native_process_seen_after_s=native_seen,
        maximum_waited_child_resident_memory_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        native_command_wall_time_s=native.get("wall_time_s") if native else None,
        native_command_exit_code=native.get("exit_code") if native else None,
        native_time_includes_mesh_assembly_and_solves=True,time_step_solve_only_wall_time_s=None)
    native_log = args.output_dir/"native_run/native.stdout.txt"
    steps = []
    if native_log.exists():
        for line in native_log.read_text(errors="replace").splitlines():
            if line.startswith("step "):
                steps.append(int(line.split()[1]))
    report["completed_native_progress_steps"] = steps
    report["native_solve_completed"] = bool(source and source.get("passed"))
    report["outcome"] = ("native_run_completed" if report["native_solve_completed"] else
                         "bounded_wall_time_expired" if timed_out else
                         "native_initialization_or_solve_failed" if native else "build_or_runner_startup_failed")
    report["full_requested_simulation_time_measured"] = report["native_solve_completed"]
    write(report_path,report)
    print(report_path)
    return 0 if report["native_solve_completed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--mode",choices=("preflight","bounded-startup"),default="preflight")
    parser.add_argument("--memory-limit-gib",type=float,default=6)
    parser.add_argument("--timeout-s",type=float,default=60)
    parser.add_argument("--wsl-distro",default="Ubuntu-24.04")
    parser.add_argument("--worker",action="store_true",help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not (0<args.memory_limit_gib<=8 and 0<args.timeout_s<=300):
        parser.error("bounded startup permits at most8 GiB address space and300 wall seconds")
    if args.worker:
        if os.name=="nt":
            parser.error("worker is Linux only")
        raise SystemExit(linux_worker(args))
    started = time.perf_counter()
    if args.output_dir.exists():
        parser.error("use a fresh output directory")
    configuration = json.loads(args.config.read_bytes())
    estimate = cost(configuration)
    args.output_dir.mkdir(parents=True)
    report = {"schema":"mhdlab-mfem-timed-preflight-v1","created_utc":datetime.now(timezone.utc).isoformat(),
        "requested_configuration":configuration,"configuration_sha256":hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "mode":args.mode,"cost":estimate,"preflight_wall_time_s":time.perf_counter()-started,
        "native_solve_completed":False,"full_requested_simulation_time_measured":False,
        "outcome":"preflight_only","resource_note":"No unbounded allocation is permitted. The bounded mode measures limited startup, not an extrapolated full-run completion time."}
    report_path = args.output_dir/"benchmark.json";write(report_path,report)
    if args.mode=="preflight":
        print(report_path);return
    if os.name!="nt":
        raise SystemExit(linux_worker(args))
    argv = ["wsl.exe","-d",args.wsl_distro,"--exec","python3","tools/benchmark_mfem_edge_source.py",
            "--worker","--config",args.config.as_posix(),"--output-dir",args.output_dir.as_posix(),
            "--memory-limit-gib",str(args.memory_limit_gib),"--timeout-s",str(args.timeout_s)]
    launched = time.perf_counter()
    try:
        result = subprocess.run(argv,capture_output=True,timeout=args.timeout_s+45)
        exit_code = result.returncode
        for label,data in (("stdout",result.stdout),("stderr",result.stderr)):
            text = data.decode("utf-16le" if b"\0" in data[:200] else "utf-8",errors="replace")
            (args.output_dir/("wsl_launch."+label+".txt")).write_text(text,encoding="utf-8")
    except subprocess.TimeoutExpired:
        exit_code = None
    report = json.loads(report_path.read_bytes())
    report["windows_wsl_launch"] = {"argv":argv,"exit_code":exit_code,"wall_time_s":time.perf_counter()-launched}
    if "linux_worker" not in report:
        report["outcome"] = "wsl_or_worker_startup_failed_before_native_launch"
    write(report_path,report)
    print(report_path)
    if exit_code:
        raise SystemExit(exit_code)


if __name__=="__main__":
    main()
