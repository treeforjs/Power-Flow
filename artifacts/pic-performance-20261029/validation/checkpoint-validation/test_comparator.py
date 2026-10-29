"""Small isolated fixtures; no production directory is ever written."""
import pathlib,json,struct,shutil,math,fcntl,os,subprocess,sys,uuid
import compare_saved_pic as c
OUT=pathlib.Path(__file__).parent
FIX=OUT/("fixtures-"+uuid.uuid4().hex[:8])
S=c.schema()

def put_structure(data,base,fields,values):
    for f in fields:struct.pack_into("<"+{"f64":"d","i32":"i","u64":"Q","bool":"B"}[f["kind"]],data,base+f["offset"],values[f["name"]])

def checkpoint(directory,step,o,d,particles):
    data=bytearray(392+80*len(particles));struct.pack_into("<Q",data,0,c.MAGIC)
    put_structure(data,8,S["options"],o);put_structure(data,176,S["diagnostics"],d)
    struct.pack_into("<QQ",data,376,3,len(particles))
    for i,p in enumerate(particles):put_structure(data,392+80*i,S["particle"],p)
    path=directory/("checkpoint_"+str(step)+".bin");path.write_bytes(data);return path

def frame(directory,o,particles,time):
    h=o["width_m"]/o["nx"];full=[0.]*4;mid=[0.]*4
    for p in particles:
        cell=min(int(p["position_m[1]"]/h),1)*2+min(int(p["position_m[0]"]/h),1);full[cell]+=p["weight"]
        if abs(p["position_m[2]"])<.0005:mid[cell]+=p["weight"]
    lines=[",".join(c.FRAME)]
    for cell in range(4):
        x=(cell%2+.5)*h;y=(cell//2+.5)*h
        vals=[time,x,y,o["upper_potential_V"]*y/o["gap_m"],0.,-o["upper_potential_V"]/o["gap_m"],mid[cell]/(h*h*.001)/1e6,full[cell]/(h*h*o["depth_m"])/1e6,0]
        lines.append(",".join(format(v,".17g")for v in vals))
    (directory/"frame_0.csv").write_text("\n".join(lines)+"\n")

def make_base():
    FIX.mkdir(exist_ok=True)
    root=FIX/"valid"
    if root.exists():raise RuntimeError("fixtures already exist; use a new isolated output directory")
    root.mkdir()
    o={f["name"]:0 for f in S["options"]};o.update(nx=2,ny=2,width_m=.001,gap_m=.001,depth_m=.025,upper_potential_V=-1000.,short_z_m=-.0125,feed_z_m=.0125,distributed_fields=1)
    particles=[]
    for i in range(2):
        p={f["name"]:0 for f in S["particle"]};p.update(id=i+1,charge_C=-c.E_CHARGE,mass_kg=9.1093837139e-31,weight=1.)
        p["position_m[0]"]=p["position_m[1]"]=(.25+.5*i)*.001;particles.append(p)
    ds={key:0. for key in c.DIAGNOSTICS};ds.update(time_s=1e-9,particle_charge_C=-2*c.E_CHARGE,deposited_charge_C=-2*c.E_CHARGE,upper_field_charge_C=2*c.E_CHARGE,authored_charge_C=-2*c.E_CHARGE)
    dirs={name:root/name for name in ("seed","baseline","candidate")}
    for p in dirs.values():p.mkdir()
    seed_file=checkpoint(dirs["seed"],0,o,ds,particles)
    info={"production":"synthetic structural fixture","start_ns":1.,"checkpoint":seed_file.name,"checkpoint_sha256":c.sha(seed_file)}
    for name,directory in dirs.items():
        step=0 if name=="seed"else 1;d=dict(ds);d["time_s"]=1e-9 if name=="seed"else 1.001e-9
        checkpoint(directory,step,o,d,particles)
        meta=["checkpoint_"+str(step)+".bin "+str(step)+" 0"]
        for _ in range(4):meta.extend(["0","0","0"])
        meta.append(format(d["time_s"],".17g")+" "+format(o["upper_potential_V"],".17g"));(directory/"restart.meta").write_text("\n".join(meta)+"\n")
        (directory/"seed.json").write_text(json.dumps(info))
        rows=[",".join(c.HISTORY)]
        if name!="seed":
            for st,t in ((0,1e-9),(1,1.001e-9)):rows.append(",".join(format(v,".17g")for v in [t,st,2,0,0,0,0,0,0,0,0,0]))
            (directory/"progress.json").write_text(json.dumps({"phase":"completed","steps":1,"macro_particles":2,"mpi_ranks":1,"upper_potential_V":-1000.,"start_time_s":0.,"accepted_time_s":1.001e-9,"requested_end_time_s":1.001e-9,"source_reservoir_electrons":0.}))
            (directory/"run-receipt.json").write_text(json.dumps({"returncode":0,"seed":info,"profile_cpu":6,"window_ps":1.}))
            (directory/"scope.json").write_text(json.dumps({"fixture":"structural validation only, no physical solution claim"}))
            (directory/"parallel_execution.json").write_text(json.dumps({"completed":True,"all_ranks_participated":True,"mpi_ranks":1,"linear_solves_per_rank":4}))
            frame(directory,o,particles,d["time_s"])
        (directory/"history.csv").write_text("\n".join(rows)+"\n")
    return dirs

def mutate_binary(directory,offset,fmt,value):
    p=directory/"checkpoint_1.bin";data=bytearray(p.read_bytes());struct.pack_into("<"+fmt,data,offset,value);p.write_bytes(data)

def mutate_csv(directory,name,column,value,row=1):
    p=directory/name;lines=p.read_text().splitlines();keys=lines[0].split(",");cells=lines[row].split(",");cells[keys.index(column)]=str(value);lines[row]=",".join(cells);p.write_text("\n".join(lines)+"\n")

def main():
    cases=[];cli=[]
    with c.cpu_lease():
        dirs=make_base();a,b,seed=(dirs[k]for k in ("baseline","candidate","seed"))
        valid=c.compare(a,b,seed);assert valid["status"]=="passed",valid["failures"]
        cases.append({"name":"finite_identical_fixture","accepted":True})
        def case(name,mutation,accept=False):
            target=FIX/name;shutil.copytree(b,target);mutation(target)
            try:r=c.compare(a,target,seed);accepted=r["status"]=="passed";reason="numerical checks: "+str([x["name"]for x in r["failures"]])
            except (c.Invalid,ValueError,OSError,KeyError,TypeError,IndexError,ZeroDivisionError,OverflowError,struct.error)as exc:accepted=False;reason=type(exc).__name__+": "+str(exc)
            assert accepted==accept,(name,accepted,reason)
            cases.append({"name":name,"expected_accept":accept,"accepted":accepted,"reason":reason})
            return target
        for group,base in (("options",8),("diagnostics",176),("particle",392)):
            for f in S[group]:
                if f["kind"]!="f64":continue
                for label,value in (("nan",math.nan),("inf",math.inf)):
                    case(group+"_"+f["name"].replace(".","_").replace("[","_").replace("]","")+"_"+label,lambda p,offset=base+f["offset"],value=value:mutate_binary(p,offset,"d",value))
        case("invalid_boolean",lambda p:mutate_binary(p,8+72,"B",2))
        case("duplicate_id",lambda p:mutate_binary(p,392+80,"Q",1))
        case("invalid_next_id",lambda p:mutate_binary(p,376,"Q",1))
        case("count_guard",lambda p:mutate_binary(p,384,"Q",2000001))
        case("zero_mass",lambda p:mutate_binary(p,392+64,"d",0.))
        case("negative_weight",lambda p:mutate_binary(p,392+72,"d",-1.))
        case("finite_velocity_mismatch",lambda p:mutate_binary(p,392+32,"d",1e6))
        case("finite_position_mismatch",lambda p:mutate_binary(p,392+8,"d",.00026))
        case("saved_residual_guard_failure",lambda p:mutate_binary(p,176+192,"d",1.))
        case("bad_magic",lambda p:mutate_binary(p,0,"Q",0))
        case("truncated",lambda p:(p/"checkpoint_1.bin").write_bytes((p/"checkpoint_1.bin").read_bytes()[:-1]))
        case("trailing",lambda p:(p/"checkpoint_1.bin").write_bytes((p/"checkpoint_1.bin").read_bytes()+b"x"))
        for col in c.FRAME:
            case("csv_nonfinite_"+col,lambda p,col=col:mutate_csv(p,"frame_0.csv",col,"nan"))
        for col in c.HISTORY:
            case("history_nonfinite_"+col,lambda p,col=col:mutate_csv(p,"history.csv",col,"inf"))
        case("metadata_nonfinite",lambda p:(p/"restart.meta").write_text((p/"restart.meta").read_text().replace("\n0\n","\ninf\n",1)))
        case("json_nonfinite",lambda p:(p/"progress.json").write_text((p/"progress.json").read_text().replace("1.001e-09","NaN",1)))
        case("metadata_boolean_returncode",lambda p:(p/"run-receipt.json").write_text((p/"run-receipt.json").read_text().replace("\"returncode\": 0","\"returncode\": false",1)))
        case("metadata_wrong_object",lambda p:(p/"progress.json").write_text("[]"))
        bad_field=case("finite_field_mismatch",lambda p:mutate_csv(p,"frame_0.csv","Ex_V_m",1.))
        case("finite_density_mismatch",lambda p:mutate_csv(p,"frame_0.csv","electron_density_cm3",4100.))
        case("small_field_roundoff",lambda p:mutate_csv(p,"frame_0.csv","Ex_V_m",1e-8),True)
        # A conservative occupancy discontinuity caused by a sub-tolerance position shift.
        moved_a=FIX/"boundary_baseline";moved_b=FIX/"boundary_candidate";shutil.copytree(a,moved_a);shutil.copytree(b,moved_b)
        for directory,x in ((moved_a,.0005-2e-16),(moved_b,.0005+2e-16)):
            mutate_binary(directory,400,"d",x);state=c.read_state(directory,S);frame(directory,state.options,state.particles,state.ledger["time_s"])
        r=c.compare(moved_a,moved_b,seed);assert r["status"]=="passed",r["failures"]
        assert r["density_analysis"]["transverse_bin_change_count"]==1 and r["frame"]["electron_density_cm3"]["linf_absolute"]>0
        (OUT/"boundary-density-fixture-result.json").write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
        cases.append({"name":"explained_density_boundary_rebin","accepted":True,"position_linf_m":r["particles"]["position_m"]["linf_absolute"],"density_linf_cm3":r["frame"]["electron_density_cm3"]["linf_absolute"]})
    # CLI takes its own lease; do not hold another lock while launching it.
    for name,target,expected in (("pass",b,0),("numerical_failure",bad_field,1),("invalid_nonfinite",FIX/"particle_position_m_0_nan",2)):
        path=OUT/("cli_"+name+".json")
        proc=subprocess.run([sys.executable,str(OUT/"compare_saved_pic.py"),str(a),str(target),str(seed),"--output",str(path)],capture_output=True,text=True)
        assert proc.returncode==expected,(name,proc.returncode,proc.stdout,proc.stderr)
        cli.append({"name":name,"returncode":proc.returncode,"result":str(path)})
    report={"status":"passed","case_count":len(cases),"fixtures":str(FIX),"cases":cases,"cli":cli,"execution":{"cpu":6,"nice":15,"threads":1,"lease":str(c.LEASE),"gpu_work":False},"fixture_scope":"Synthetic data exercise parsing, consistency, numerical gates and density rebin handling; they do not assert a physical Poisson solution."}
    (OUT/"fixture-tests.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":"passed","cases":len(cases),"cli_exit_codes":[x["returncode"]for x in cli]}))
if __name__=="__main__":main()
