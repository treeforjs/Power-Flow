#!/usr/bin/env python3
"""Read-only comparator for the frozen PIC2 ABI. CLI work takes the CPU-test lease.

Saved checkpoints contain particle state and raw ledgers, not the nodal potential.
Public kinetic energy and particle-charge defect are reconstructed as diagnostics()
does. CSV comparisons cover the saved cell-center samples only.
"""
import argparse, contextlib, csv, dataclasses, fcntl, hashlib, json, math, os
import pathlib, re, struct, sys, datetime

E_CHARGE=1.602176634e-19
C_LIGHT=299792458.
EPS=sys.float_info.epsilon
MAGIC=0x5241494b50494332
MAX_PARTICLES=2000000
HERE=pathlib.Path(__file__).resolve().parent
LEASE=HERE.parent.parent/"testing-core.lock"
DIAGNOSTICS="time_s particle_charge_C deposited_charge_C lower_field_charge_C upper_field_charge_C gauss_defect_C particle_charge_defect_C field_energy_J particle_kinetic_energy_J lower_absorbed_charge_C upper_absorbed_charge_C lower_absorbed_energy_J upper_absorbed_energy_J lower_injected_charge_C upper_injected_charge_C injected_energy_J authored_charge_C authored_kinetic_energy_J maximum_nodal_continuity_defect_C escaped_charge_C escaped_kinetic_energy_J short_absorbed_charge_C short_absorbed_energy_J axial_escaped_charge_C poisson_residual_C_per_m".split()
HISTORY="time_s step macro_particles gauss_defect_C particle_charge_defect_C nodal_continuity_defect_C field_energy_J kinetic_energy_J source_reservoir_electrons short_absorbed_charge_C short_absorbed_energy_J axial_escaped_charge_C".split()
FRAME="time_s x_m y_m potential_V Ex_V_m Ey_V_m electron_density_cm3 height_average_electron_density_cm3 solid".split()
HISTORY_LEDGER={"gauss_defect_C":"gauss_defect_C","particle_charge_defect_C":"particle_charge_defect_C","nodal_continuity_defect_C":"maximum_nodal_continuity_defect_C","field_energy_J":"field_energy_J","kinetic_energy_J":"particle_kinetic_energy_J","short_absorbed_charge_C":"short_absorbed_charge_C","short_absorbed_energy_J":"short_absorbed_energy_J","axial_escaped_charge_C":"axial_escaped_charge_C"}

class Invalid(ValueError): pass

def require(condition,message):
    if not condition: raise Invalid(message)

def finite(value,label):
    require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value),label+": nonfinite or nonnumeric")
    return value

def strict_json(path):
    def pairs(items):
        result={}
        for k,v in items:
            require(k not in result,str(path)+": duplicate JSON key "+k);result[k]=v
        return result
    def bad(s): raise Invalid(str(path)+": nonfinite JSON constant "+s)
    value=json.loads(path.read_text(),parse_constant=bad,object_pairs_hook=pairs)
    def walk(v,key):
        if isinstance(v,float):finite(v,key)
        elif isinstance(v,dict):
            for k,x in v.items():walk(x,key+"."+k)
        elif isinstance(v,list):
            for i,x in enumerate(v):walk(x,key+"["+str(i)+"]")
    walk(value,str(path));return value

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def schema(path=HERE/"abi-schema.json"):
    s=strict_json(path)
    require((s["endian"],s["options_size"],s["diagnostics_size"],s["particle_size"])==("little",168,200,80),"unsupported frozen checkpoint ABI")
    require([f["name"]for f in s["diagnostics"]]==DIAGNOSTICS,"incomplete diagnostic schema")
    require(len(s["options"])==23 and len(s["particle"])==10,"incomplete options/particle schema")
    return s

def decode(data,offset,fields,label):
    result={}
    for f in fields:
        fmt={"f64":"d","i32":"i","u64":"Q","bool":"B"}[f["kind"]]
        v=struct.unpack_from("<"+fmt,data,offset+f["offset"])[0]
        if f["kind"]=="bool":require(v in (0,1),label+"."+f["name"]+": invalid boolean byte")
        else:finite(v,label+"."+f["name"])
        result[f["name"]]=v
    return result

def kinetic(p):
    u2=sum(p["proper_velocity_m_s["+str(k)+"]"]**2 for k in range(3))
    finite(u2,"particle proper velocity squared")
    value=p["mass_kg"]*u2/(math.sqrt(1+u2/(C_LIGHT*C_LIGHT))+1)
    return finite(value,"particle kinetic energy")

def read_csv(path,columns):
    with path.open(newline="") as f:
        reader=csv.reader(f)
        require(next(reader,None)==columns,str(path)+": unexpected columns")
        rows=[]
        for line,row in enumerate(reader,2):
            require(len(row)==len(columns),str(path)+": row width at "+str(line))
            values={}
            for name,text in zip(columns,row):
                try:v=float(text)
                except ValueError as exc:raise Invalid(str(path)+": invalid "+name+" at row "+str(line))from exc
                values[name]=finite(v,str(path)+":"+str(line)+":"+name)
            rows.append(values)
    return rows

@dataclasses.dataclass
class State:
    directory:pathlib.Path
    options:dict
    ledger:dict
    particles:list
    next_id:int
    step:int
    next_frame:int
    pools:list
    meta_time:float
    meta_voltage:float
    history:list
    seed_info:dict
    checkpoint:pathlib.Path
    files:list

    def public(self):
        if "_public_cache" in self.__dict__:return self._public_cache
        d=dict(self.ledger)
        # Match the sequential, double-precision diagnostics() operation order.
        for p in self.particles:d["particle_kinetic_energy_J"]+=kinetic(p)*p["weight"]
        d["particle_charge_defect_C"]=(d["particle_charge_C"]+d["lower_absorbed_charge_C"]+d["upper_absorbed_charge_C"]+d["short_absorbed_charge_C"]-d["lower_injected_charge_C"]-d["upper_injected_charge_C"]-d["authored_charge_C"]+d["escaped_charge_C"])
        self._public_cache=d
        return d

    def reservoir(self):
        return sum(2*sum(p["integral"])*(.005 if i<2 else .0002)/len(p["integral"])*self.options["depth_m"]/E_CHARGE for i,p in enumerate(self.pools))

def read_state(directory,s):
    directory=pathlib.Path(directory).resolve()
    meta=directory/"restart.meta"
    lines=meta.read_text().splitlines()
    require(len(lines)==14,str(meta)+": expected complete 4-pool metadata")
    head=lines[0].split();require(len(head)==3,str(meta)+": invalid header")
    name,step_text,frame_text=head
    require(re.fullmatch(r"checkpoint_[0-9]+\.bin",name) is not None,str(meta)+": invalid checkpoint basename")
    try:step=int(step_text);frame=int(frame_text)
    except ValueError as exc:raise Invalid(str(meta)+": nonintegral step/frame")from exc
    require(step>=0 and frame>=0 and name=="checkpoint_"+str(step)+".bin",str(meta)+": inconsistent checkpoint step")
    checkpoint=directory/name
    require(checkpoint.is_file(),str(checkpoint)+": missing checkpoint")
    require(checkpoint.stat().st_size<=392+80*MAX_PARTICLES,str(checkpoint)+": exceeds checkpoint resource guard")
    data=checkpoint.read_bytes()
    require(len(data)>=392,str(checkpoint)+": truncated header")
    magic=struct.unpack_from("<Q",data)[0];require(magic==MAGIC,str(checkpoint)+": bad magic")
    options=decode(data,8,s["options"],"options")
    ledger=decode(data,8+s["options_size"],s["diagnostics"],"diagnostics")
    off=8+s["options_size"]+s["diagnostics_size"]
    next_id,n=struct.unpack_from("<QQ",data,off)
    require(1<=next_id and n<=MAX_PARTICLES,str(checkpoint)+": identity/count guard")
    require(len(data)==off+16+n*s["particle_size"],str(checkpoint)+": truncated or trailing bytes")
    require(1<=options["nx"]<=4096 and 1<=options["ny"]<=4096,"options: invalid grid")
    for key in ("width_m","gap_m","depth_m"):require(options[key]>0,"options: nonpositive "+key)
    if options["finite_foils"]:
        for foil in ("lower_foil","upper_foil"):
            require(0<options[foil+".left_m"]<options[foil+".right_m"]<options["width_m"] and 0<options[foil+".bottom_m"]<options[foil+".top_m"]<options["gap_m"],"invalid conductor geometry")
        require(options["lower_foil.top_m"]<options["upper_foil.bottom_m"],"overlapping conductors")
    if options["finite_axial_extent"]:require(options["finite_foils"] and options["feed_z_m"]>options["short_z_m"] and abs(options["feed_z_m"]-options["short_z_m"]-options["depth_m"])<=1e-12*options["depth_m"],"invalid axial extent")
    require(ledger["time_s"]>=0,"negative checkpoint clock")
    for key,v in ledger.items():
        if "energy_J" in key or key in ("poisson_residual_C_per_m","maximum_nodal_continuity_defect_C"):require(v>=0,"negative ledger "+key)
    particles=[];ids=set()
    for i in range(n):
        p=decode(data,off+16+i*s["particle_size"],s["particle"],"particle["+str(i)+"]")
        require(0<p["id"]<next_id and p["id"]not in ids,"duplicate/invalid particle identity")
        ids.add(p["id"])
        require(p["mass_kg"]>0 and p["weight"]>0,"nonpositive particle mass/weight")
        finite(p["charge_C"]*p["weight"],"represented particle charge")
        finite(p["charge_C"]/p["mass_kg"],"particle charge-to-mass ratio")
        finite(kinetic(p)*p["weight"],"weighted particle kinetic energy")
        x,y,z=[p["position_m["+str(k)+"]"]for k in range(3)]
        require(0<=x<=options["width_m"] and 0<=y<=options["gap_m"],"particle outside transverse domain")
        if options["finite_axial_extent"]:require(options["short_z_m"]<=z<=options["feed_z_m"],"particle outside axial domain")
        if options["finite_foils"]:
            for foil in ("lower_foil","upper_foil"):require(not(options[foil+".left_m"]<x<options[foil+".right_m"] and options[foil+".bottom_m"]<y<options[foil+".top_m"]),"particle inside conductor")
        particles.append(p)
    pools=[]
    for i in range(4):
        values=[]
        for j in range(3):
            try:v=[float(x)for x in lines[1+3*i+j].split()]
            except ValueError as exc:raise Invalid("invalid reservoir token")from exc
            for x in v:finite(x,"reservoir")
            require(v and all(x>=0 for x in v),"empty/negative reservoir")
            values.append(v)
        require(len(values[0])==1 and len(values[1])==len(values[2]) and 1<=len(values[1])<=4096,"reservoir shape guard")
        for q,t in zip(values[1],values[2]):require((q==t==0)or(q>0 and 0<t/q<=1e6),"invalid reservoir weighted temperature")
        pools.append({"span":values[0][0],"integral":values[1],"weighted_temperature":values[2]})
    tail=lines[-1].split();require(len(tail)==2,"missing metadata clock/voltage")
    meta_time,meta_voltage=[finite(float(x),"metadata")for x in tail]
    require(meta_time==ledger["time_s"] and meta_voltage==options["upper_potential_V"],"checkpoint/metadata clock or voltage mismatch")
    history=read_csv(directory/"history.csv",HISTORY)
    last=None;unique=[]
    for row in history:
        require(row["step"].is_integer() and row["macro_particles"].is_integer() and 0<=row["macro_particles"]<=MAX_PARTICLES,"invalid history particle/step count")
        for key in ("nodal_continuity_defect_C","field_energy_J","kinetic_energy_J","source_reservoir_electrons","short_absorbed_energy_J"):require(row[key]>=0,"negative history "+key)
        if last:
            require(row["step"]>=last["step"] and row["time_s"]>=last["time_s"],"nonmonotonic history")
            if row["step"]==last["step"]:
                require(row==last,"conflicting repeated history step");continue
        unique.append(row);last=row
    seed_info=strict_json(directory/"seed.json")
    require(isinstance(seed_info,dict) and set(seed_info)=={"production","start_ns","checkpoint","checkpoint_sha256"},"invalid seed identity schema")
    finite(seed_info["start_ns"],"seed start_ns")
    require(isinstance(seed_info["production"],str) and isinstance(seed_info["checkpoint"],str) and isinstance(seed_info["checkpoint_sha256"],str) and re.fullmatch(r"[0-9a-f]{64}",seed_info["checkpoint_sha256"])is not None,"invalid seed identity values")
    return State(directory,options,ledger,particles,next_id,step,frame,pools,meta_time,meta_voltage,unique,seed_info,checkpoint,[meta,checkpoint,directory/"history.csv",directory/"seed.json"])

class Metrics:
    def __init__(self):self.n=0;self.sse=0.;self.a2=0.;self.b2=0.;self.maximum=0.;self.scale=0.;self.index=None
    def add(self,a,b,index=None):
        finite(a,"metric a");finite(b,"metric b")
        d=finite(a-b,"metric difference");self.n+=1
        self.sse+=d*d;self.a2+=a*a;self.b2+=b*b
        self.scale=max(self.scale,abs(a),abs(b))
        if abs(d)>self.maximum:self.maximum=abs(d);self.index=index
    def result(self):
        for v in (self.sse,self.a2,self.b2,self.scale):finite(v,"metric accumulator")
        den=math.sqrt(max(self.a2,self.b2))
        return {"count":self.n,"linf_absolute":self.maximum,"l2_absolute":math.sqrt(self.sse),"rms_absolute":math.sqrt(self.sse/max(1,self.n)),"reference_linf":self.scale,"a_l2":math.sqrt(self.a2),"b_l2":math.sqrt(self.b2),"relative_l2":math.sqrt(self.sse)/den if den else 0.,"relative_linf_global_scale":self.maximum/self.scale if self.scale else 0.,"max_difference_index":self.index}

def vector_metric(a,b):
    require(len(a)==len(b),"metric length mismatch")
    m=Metrics()
    for i,(x,y)in enumerate(zip(a,b)):m.add(x,y,i)
    return m.result()

POLICY={"position_atol_m":1e-12,"velocity_rtol_global":1e-10,"field_rtol_global":1e-10,"ledger_roundoff_factor":1024,"density_roundoff_factor":256,"time_roundoff_factor":64,"rationale":{"position":"Existing checkpoint regression budget; 4e-8 of a 25 um cell. Record measured differences and bin-boundary crossings separately.","velocity":"Existing global proper-velocity regression budget. Mass, charge, weights and identity remain exact.","fields":"Conservative integration regression budget, not a production solver tolerance. 1e-10 of each global component plus 256 epsilon V_scale/h for derivative cancellation. The corresponding e/m E_error window impulse is reported against the velocity budget.","ledger":"1024 epsilon times the appropriate full charge/energy scale, rather than relative error of a cancellation residual near zero.","density":"Validate each CSV against that checkpoint particle deposition; permit only roundoff or particle-bin/slab-boundary changes explained by matched, position-qualified particles. No blanket absolute density threshold.","time":"64 epsilon times the global clock scale, with a 1e-23 s floor reflecting the existing nonprogress guard. Pairwise clocks/options/counts remain exact."}}

class Comparison:
    def __init__(self):self.checks=[];self.sections={}
    def check(self,name,passed,**details):self.checks.append({"name":name,"passed":bool(passed),**details})
    def scalar(self,name,a,b,limit):
        m=vector_metric([a],[b]);self.check(name,m["linf_absolute"]<=limit,a=a,b=b,absolute_limit=limit,metrics=m)
    def vector(self,name,a,b,limit):
        m=vector_metric(a,b);self.check(name,m["linf_absolute"]<=limit,absolute_limit=limit,metrics=m);return m


def tolerance(key,states,guards):
    cache=states[0].__dict__.setdefault("_tolerance_cache",{})
    if key in cache:return cache[key]
    value=_tolerance(key,states,guards)
    cache[key]=value
    return value

def _tolerance(key,states,guards):
    d=[s.public()for s in states]
    particle_charge=max(sum(abs(p["charge_C"]*p["weight"])for p in s.particles)for s in states)
    particle_charge=max(particle_charge,max(sum(abs(v)for k,v in x.items()if k.endswith("charge_C")and k not in ("lower_field_charge_C","upper_field_charge_C"))for x in d),E_CHARGE)
    electrode_charge=max(abs(x["lower_field_charge_C"])+abs(x["upper_field_charge_C"])for x in d)
    energy=max(sum(abs(v)for k,v in x.items()if "energy_J"in k and k!="field_energy_J")for x in d)
    if key=="time_s":return max(1e-23,POLICY["time_roundoff_factor"]*EPS*max(abs(x["time_s"])for x in d))
    if key=="poisson_residual_C_per_m":return max(g["guard_limit"]for g in guards)
    if key=="field_energy_J":scale=max(abs(x[key])for x in d)
    elif "energy_J"in key:scale=energy
    elif key in ("lower_field_charge_C","upper_field_charge_C","gauss_defect_C"):scale=electrode_charge
    else:scale=particle_charge
    return POLICY["ledger_roundoff_factor"]*EPS*scale


def guard_probe(state):
    import subprocess
    command=["/home/treefor/projects/runtimes/arch-base-20261028/run-frozen",str(HERE/"residual_guard_probe"),str(state.checkpoint)]
    proc=subprocess.run(command,capture_output=True,text=True)
    require(proc.returncode==0,"unchanged Poisson residual guard failed: "+str(state.checkpoint)+" "+proc.stderr.strip())
    result=json.loads(proc.stdout)
    for key in ("rhs_norm","residual_norm","guard_limit"):finite(result[key],"RHS probe "+key)
    require(result["guard_passed"]and result["residual_norm"]==state.ledger["poisson_residual_C_per_m"],"invalid RHS probe result")
    result["command"]=command;return result


def validate_run(state,seed,comp,states,guards):
    root=state.directory;prefix=root.name
    progress=strict_json(root/"progress.json");receipt=strict_json(root/"run-receipt.json");scope=strict_json(root/"scope.json");parallel=strict_json(root/"parallel_execution.json")
    state.files.extend(root/name for name in ("progress.json","run-receipt.json","scope.json","parallel_execution.json"))
    require(all(isinstance(v,dict)for v in (progress,receipt,scope,parallel)),"invalid run metadata object")
    require(all(type(progress[k])is int for k in ("steps","macro_particles","mpi_ranks")) and type(receipt["returncode"])is int and type(receipt["profile_cpu"])is int and all(type(parallel[k])is int for k in ("mpi_ranks","linear_solves_per_rank")),"invalid run metadata integer type")
    for k in ("accepted_time_s","requested_end_time_s","start_time_s","upper_potential_V","source_reservoir_electrons"):finite(progress[k],"progress "+k)
    finite(receipt["window_ps"],"receipt window_ps")
    require(progress["phase"]=="completed" and receipt["returncode"]==0,"run did not complete successfully")
    require(receipt["seed"]==seed.seed_info==state.seed_info,"seed identity mismatch")
    require(receipt["profile_cpu"]==6 and receipt["window_ps"]>0,"unexpected replay CPU/window")
    require(progress["steps"]==state.step and progress["macro_particles"]==len(state.particles),"progress/checkpoint step or count mismatch")
    require(progress["upper_potential_V"]==state.options["upper_potential_V"],"progress/metadata voltage mismatch")
    require(parallel["completed"]is True and parallel["all_ranks_participated"]is True and parallel["mpi_ranks"]==progress["mpi_ranks"]==1,"incomplete parallel execution")
    require(parallel["linear_solves_per_rank"]==2*(state.step-seed.step)+2,"unexpected solve count")
    require(state.history and state.history[0]["step"]==seed.step and state.history[-1]["step"]==state.step,"history missing seed/final endpoints")
    time_limit=tolerance("time_s",states,guards)
    expected_seed_global=seed.seed_info["start_ns"]*1e-9
    comp.scalar(prefix+".seed_global_clock",state.history[0]["time_s"],expected_seed_global,time_limit)
    comp.scalar(prefix+".final_global_clock",state.history[-1]["time_s"],progress["start_time_s"]+state.ledger["time_s"],time_limit)
    comp.scalar(prefix+".progress_clock",progress["accepted_time_s"],state.history[-1]["time_s"],time_limit)
    comp.scalar(prefix+".requested_clock",progress["accepted_time_s"],progress["requested_end_time_s"],time_limit)
    comp.scalar(prefix+".window",(state.ledger["time_s"]-seed.ledger["time_s"]),receipt["window_ps"]*1e-12,time_limit)
    require(state.history[0]["macro_particles"]==len(seed.particles),"history initial particle count mismatch")
    require(state.history[-1]["macro_particles"]==len(state.particles),"history final particle count mismatch")
    public=state.public();seed_public=seed.public()
    for column,key in HISTORY_LEDGER.items():
        comp.scalar(prefix+".history_final."+column,state.history[-1][column],public[key],tolerance(key,states,guards))
        comp.scalar(prefix+".history_seed."+column,state.history[0][column],seed_public[key],tolerance(key,states,guards))
    comp.scalar(prefix+".reservoir_final",state.history[-1]["source_reservoir_electrons"],state.reservoir(),256*EPS*max(1.,state.reservoir()))
    comp.scalar(prefix+".reservoir_initial",state.history[0]["source_reservoir_electrons"],seed.reservoir(),256*EPS*max(1.,seed.reservoir()))
    comp.scalar(prefix+".progress_reservoir",progress["source_reservoir_electrons"],state.reservoir(),256*EPS*max(1.,state.reservoir()))
    return {"progress":progress,"receipt":receipt,"scope":scope,"parallel":parallel}


def grid_density(state):
    o=state.options;nx,ny=o["nx"],o["ny"]
    require(nx==ny and o["width_m"]==o["gap_m"],"unsupported non-square frozen snapshot convention")
    h=o["width_m"]/nx;volume=h*h*o["depth_m"]
    require(nx*ny<=1000000,"snapshot comparator resource limit: >1 million cells")
    full=[0.]*(nx*ny);mid=[0.]*(nx*ny);bins={}
    for p in state.particles:
        i=min(max(int(p["position_m[0]"]/h),0),nx-1);j=min(max(int(p["position_m[1]"]/h),0),ny-1);cell=j*nx+i
        slab=abs(p["position_m[2]"])<.0005
        full[cell]+=p["weight"]
        if slab:mid[cell]+=p["weight"]
        bins[p["id"]]=(cell,slab)
    return {"height_average_electron_density_cm3":[v/volume/1e6 for v in full],"electron_density_cm3":[v/(h*h*.001)/1e6 for v in mid]},bins


def validate_frame(state,comp,states,guards):
    path=state.directory/("frame_"+str(state.next_frame)+".csv")
    rows=read_csv(path,FRAME);state.files.append(path)
    o=state.options;nx,ny=o["nx"],o["ny"];h=o["width_m"]/nx
    require(len(rows)==nx*ny,"frame cell count mismatch")
    density,bins=grid_density(state)
    for i,row in enumerate(rows):
        x,y=((i%nx)+.5)*h,((i//nx)+.5)*h
        require(row["x_m"]==x and row["y_m"]==y,"frame coordinates/order mismatch at "+str(i))
        require(row["solid"]in (0.,1.),"invalid solid marker")
        solid=any(o[foil+".left_m"]<x<o[foil+".right_m"]and o[foil+".bottom_m"]<y<o[foil+".top_m"]for foil in ("lower_foil","upper_foil"))if o["finite_foils"]else False
        require(row["solid"]==int(solid),"solid marker inconsistent with options")
        require(row["time_s"]==rows[0]["time_s"],"mixed frame clocks")
        for column in density:require(row[column]>=0,"negative electron density")
    comp.scalar(state.directory.name+".frame_clock",rows[0]["time_s"],state.history[-1]["time_s"],tolerance("time_s",states,guards))
    for column,expected in density.items():
        actual=[r[column]for r in rows];m=vector_metric(actual,expected)
        slab=.001 if column=="electron_density_cm3"else o["depth_m"]
        one_electron_density=1/(h*h*slab*1e6)
        limit=POLICY["density_roundoff_factor"]*EPS*max(one_electron_density,m["reference_linf"])
        comp.check(state.directory.name+".frame_reconstruction."+column,m["linf_absolute"]<=limit,absolute_limit=limit,one_electron_density_cm3=one_electron_density,metrics=m)
        represented=sum(actual)*h*h*slab*1e6
        expected_count=sum(p["weight"]for p in state.particles if column!="electron_density_cm3"or abs(p["position_m[2]"])<.0005)
        comp.scalar(state.directory.name+".integrated_density_count."+column,represented,expected_count,1024*EPS*max(1.,expected_count))
    return rows,bins


def compare(a_path,b_path,seed_path,schema_path=HERE/"abi-schema.json"):
    input_paths={p for root in (a_path,b_path,seed_path) for p in pathlib.Path(root).resolve().iterdir() if p.is_file() and (p.suffix in (".bin",".csv",".meta",".json"))}
    before={p:sha(p)for p in input_paths}
    s=schema(schema_path);a,b,seed=[read_state(p,s)for p in (a_path,b_path,seed_path)]
    states=[a,b,seed];comp=Comparison()
    tracked={p:before[p]for state in states for p in state.files}
    require(seed.seed_info["checkpoint"]==seed.checkpoint.name and sha(seed.checkpoint)==seed.seed_info["checkpoint_sha256"],"seed checkpoint hash mismatch")
    require(a.options==b.options,"options mismatch")
    for key,v in seed.options.items():
        if key not in ("lower_potential_V","upper_potential_V"):require(a.options[key]==v,"static option differs from seed: "+key)
    require(a.step==b.step and a.next_id==b.next_id and a.next_frame==b.next_frame,"restart step/frame/next identity mismatch")
    require(a.ledger["time_s"]==b.ledger["time_s"],"final local clocks differ")
    require(a.step>seed.step and a.next_id>=seed.next_id,"nonprogressing replay")
    pool_metrics=[]
    for i,(pa,pb)in enumerate(zip(a.pools,b.pools)):
        require(len(pa["integral"])==len(seed.pools[i]["integral"]),"source patch count differs from seed")
        pool_metrics.append({"span":vector_metric([pa["span"]],[pb["span"]]),"integral":vector_metric(pa["integral"],pb["integral"]),"weighted_temperature":vector_metric(pa["weighted_temperature"],pb["weighted_temperature"])})
        require(pa==pb,"source reservoir mismatch in pool "+str(i))
    guards=[guard_probe(state)for state in states]
    for state,g in zip(states,guards):
        d=state.public();charge=sum(p["charge_C"]*p["weight"]for p in state.particles)
        comp.scalar(state.directory.name+".particle_charge_reconstruction",d["particle_charge_C"],charge,tolerance("particle_charge_C",states,guards))
        comp.scalar(state.directory.name+".deposition_charge",d["deposited_charge_C"],d["particle_charge_C"],tolerance("particle_charge_C",states,guards))
        comp.scalar(state.directory.name+".gauss_identity",d["gauss_defect_C"],d["particle_charge_C"]+d["lower_field_charge_C"]+d["upper_field_charge_C"],tolerance("gauss_defect_C",states,guards))
        comp.check(state.directory.name+".poisson_guard",g["guard_passed"],**g)
    meta=[validate_run(state,seed,comp,states,guards)for state in (a,b)]
    require(meta[0]["scope"]==meta[1]["scope"],"scope/model mismatch")
    for key in ("start_time_s","requested_end_time_s"):require(meta[0]["progress"][key]==meta[1]["progress"][key],"run time domain mismatch")
    ap={p["id"]:p for p in a.particles};bp={p["id"]:p for p in b.particles};sp={p["id"]:p for p in seed.particles}
    require(ap.keys()==bp.keys(),"particle identity/count mismatch")
    particle_metrics={}
    for key in ("charge_C","mass_kg","weight"):
        particle_metrics[key]=comp.vector("particles."+key,[p[key]for p in ap.values()],[bp[i][key]for i in ap],0.)
    for key in ("position_m","proper_velocity_m_s"):
        va=[p[key+"["+str(k)+"]"]for p in ap.values()for k in range(3)]
        vb=[bp[i][key+"["+str(k)+"]"]for i in ap for k in range(3)]
        limit=POLICY["position_atol_m"]if key=="position_m"else POLICY["velocity_rtol_global"]*max([abs(x)for x in va+vb],default=0.)
        particle_metrics[key]=comp.vector("particles."+key,va,vb,limit)
        particle_metrics[key+"_components"]=[vector_metric([p[key+"["+str(k)+"]"]for p in ap.values()],[bp[i][key+"["+str(k)+"]"]for i in ap])for k in range(3)]
    raw={};increments={}
    for key in DIAGNOSTICS:
        limit=tolerance(key,states,guards)
        comp.scalar("raw_diagnostics."+key,a.ledger[key],b.ledger[key],0. if key=="time_s"else limit)
        raw[key]={"baseline":a.ledger[key],"candidate":b.ledger[key],"absolute_difference":abs(a.ledger[key]-b.ledger[key]),"absolute_limit":0. if key=="time_s"else limit}
        da,db=a.public()[key]-seed.public()[key],b.public()[key]-seed.public()[key]
        comp.scalar("seed_to_final_increment."+key,da,db,limit)
        increments[key]={"baseline":da,"candidate":db,"absolute_difference":abs(da-db),"absolute_limit":limit}
    ah={int(r["step"]):r for r in a.history};bh={int(r["step"]):r for r in b.history};common=sorted(ah.keys()&bh.keys())
    require(common and common[0]==seed.step and common[-1]==a.step,"history lacks common endpoints")
    for step in common:
        for column in HISTORY:
            key=HISTORY_LEDGER.get(column)
            limit=tolerance(key,states,guards)if key else 0.
            comp.scalar("history_common."+str(step)+"."+column,ah[step][column],bh[step][column],limit)
    history_increments={}
    for column,key in HISTORY_LEDGER.items():
        da=a.history[-1][column]-a.history[0][column];db=b.history[-1][column]-b.history[0][column]
        comp.scalar("history_endpoint_increment."+column,da,db,2*tolerance(key,states,guards))
        history_increments[column]={"baseline":da,"candidate":db,"absolute_difference":abs(da-db)}
    fa,bina=validate_frame(a,comp,states,guards);fb,binb=validate_frame(b,comp,states,guards)
    frame_metrics={}
    vscale=max(abs(r["potential_V"])for r in fa+fb);h=a.options["width_m"]/a.options["nx"]
    field_error=0.
    for column in FRAME:
        va=[r[column]for r in fa];vb=[r[column]for r in fb];m=vector_metric(va,vb);frame_metrics[column]=m
        if column in ("potential_V","Ex_V_m","Ey_V_m"):
            floor=256*EPS*vscale/(h if column!="potential_V"else 1.)
            limit=POLICY["field_rtol_global"]*m["reference_linf"]+floor
            comp.check("frame."+column,m["linf_absolute"]<=limit,absolute_limit=limit,roundoff_floor=floor,metrics=m)
            if column!="potential_V":field_error=max(field_error,limit)
        elif column not in ("electron_density_cm3","height_average_electron_density_cm3"):
            comp.check("frame."+column,m["linf_absolute"]==0.,absolute_limit=0.,metrics=m)
    moved=[i for i in ap if bina[i][0]!=binb[i][0]];slab=[i for i in ap if bina[i][1]!=binb[i][1]]
    density_explanation={"transverse_bin_change_count":len(moved),"slab_membership_change_count":len(slab),"transverse_bin_change_ids":moved,"slab_membership_change_ids":slab,"basis":"Each saved density independently reconstructed from its own particles; paired identity/weight/position checks constrain bin changes. Different occupancy bins are not inferred crossing events.","all_density_components_identical":all(frame_metrics[k]["linf_absolute"]==0. for k in ("electron_density_cm3","height_average_electron_density_cm3"))}
    for state in states:
        for p in state.files:
            if p not in tracked:tracked[p]=before[p]
    # Detect a file being replaced while it is read; inputs must be completed outputs.
    require(all(sha(p)==digest for p,digest in tracked.items()),"input changed during comparison")
    failures=[x for x in comp.checks if not x["passed"]]
    births=a.next_id-seed.next_id
    report={"status":"passed"if not failures else"failed","comparison_complete":True,"policy":POLICY,"paths":{"baseline":str(a.directory),"candidate":str(b.directory),"seed":str(seed.directory)},"abi_schema_sha256":sha(schema_path),"options":{"all_fields_exact":True,"field_count":len(s["options"]),"values":a.options},"particle_counts":{"seed":len(seed.particles),"final":len(a.particles),"allocated_id_increment":births,"seed_survivors":len(ap.keys()&sp.keys()),"seed_ids_missing_at_final":len(sp.keys()-ap.keys()),"new_ids_present_at_final":len(ap.keys()-sp.keys()),"inferred_total_removed_count":len(seed.particles)+births-len(a.particles)},"particles":particle_metrics,"raw_diagnostics":raw,"seed_to_final_public_ledger_increments":increments,"history":{"baseline_rows":len(a.history),"candidate_rows":len(b.history),"common_steps":common,"baseline_only_steps":sorted(ah.keys()-bh.keys()),"candidate_only_steps":sorted(bh.keys()-ah.keys()),"baseline_values":a.history,"candidate_values":b.history,"cadence":"wall-clock reports; unmatched intermediate rows validated independently but not compared by interpolation","endpoint_ledger_increments":history_increments},"source_reservoirs":{"all_values_exact":True,"pool_metrics":pool_metrics,"pool_lengths":[len(p["integral"])for p in a.pools],"baseline_electrons":a.reservoir(),"candidate_electrons":b.reservoir(),"seed_electrons":seed.reservoir()},"frame":frame_metrics,"density_analysis":density_explanation,"runtime_residual_guards":{state.directory.name:g for state,g in zip(states,guards)},"field_error_budget_impulse_m_s":E_CHARGE/min(p["mass_kg"]for p in a.particles)*field_error*(a.ledger["time_s"]-seed.ledger["time_s"])if a.particles else 0.,"checks":comp.checks,"failures":failures,"input_manifest":[{"path":str(p),"sha256":digest}for p,digest in sorted(tracked.items())],"limits":["Checkpoint has no nodal potential or crossing-event log; field comparison is limited to final saved CSV cell centers.","Only matching history steps are compared; no unrecorded intermediate state is inferred.","Saved residual norms are checked against CPU-reconstructed RHS guards; residual error vectors and solver maximum-relative history are not saved.","This is a short-window integration regression comparison, not physical, grid or timestep convergence."]}
    return report


@contextlib.contextmanager
def cpu_lease():
    require(sorted(os.sched_getaffinity(0))==[6] and os.getpriority(os.PRIO_PROCESS,0)==15,"run with taskset -c 6 nice -n 15")
    for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):require(os.environ.get(key)=="1",key+" must be 1")
    with LEASE.open("a+")as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        yield {"cpu":6,"nice":15,"threads":1,"lease":str(LEASE),"gpu_work":False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline",type=pathlib.Path);parser.add_argument("candidate",type=pathlib.Path);parser.add_argument("seed",type=pathlib.Path);parser.add_argument("--output",type=pathlib.Path,required=True)
    args=parser.parse_args()
    try:
        with cpu_lease()as execution:report=compare(args.baseline,args.candidate,args.seed);report["execution"]=execution
    except (Invalid,ValueError,OSError,KeyError,TypeError,IndexError,ZeroDivisionError,OverflowError,struct.error)as exc:
        report={"status":"invalid","comparison_complete":False,"error":str(exc),"error_type":type(exc).__name__}
    report["created_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:report[k]for k in ("status","comparison_complete","error")if k in report}))
    return 0 if report["status"]=="passed"else 2 if report["status"]=="invalid"else 1
if __name__=="__main__":raise SystemExit(main())
