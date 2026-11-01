"""Narrow adapted frozen checkpoint reader for the preserved pause-snapshot schema.
Original comparator remains untouched. All original binary/history/reservoir checks
are retained. The extra required scope field is validated explicitly, not removed.
"""
import math,pathlib,re,struct
from historical_comparator import require,Invalid,decode,finite,kinetic,read_csv,strict_json,sha,State,MAGIC,MAX_PARTICLES,HISTORY
SCOPE='Snapshot reconstruction of paused accepted checkpoint; no physical evolution.'
PRODUCTION='/home/treefor/Documents/Codex/2026-10-28/we/outputs/electron-heating-pic-25um-full60ns-gpu-midplane110kV-height25mm'
def validate_seed(seed_info):
    require(type(seed_info) is dict and set(seed_info)=={'production','start_ns','checkpoint','checkpoint_sha256','scope'},'pause seed exact schema')
    require(seed_info['scope']==SCOPE and seed_info['production']==PRODUCTION,'pause seed scope/production')
    finite(seed_info['start_ns'],'pause seed start_ns')
    require(seed_info['start_ns']>=0,'negative seed time')
    require(isinstance(seed_info['checkpoint'],str) and re.fullmatch(r'checkpoint_[0-9]+\.bin',seed_info['checkpoint']) is not None,'seed checkpoint basename')
    require(isinstance(seed_info['checkpoint_sha256'],str) and re.fullmatch(r'[0-9a-f]{64}',seed_info['checkpoint_sha256']) is not None,'seed checksum')
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
    with checkpoint.open("rb") as f:data=f.read(393+80*MAX_PARTICLES)
    require(len(data)<=392+80*MAX_PARTICLES,"checkpoint actual read bound")
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
    validate_seed(seed_info)
    require(seed_info['checkpoint']==checkpoint.name and seed_info['checkpoint_sha256']==sha(checkpoint),'seed/checkpoint identity mismatch')
    return State(directory,options,ledger,particles,next_id,step,frame,pools,meta_time,meta_voltage,unique,seed_info,checkpoint,[meta,checkpoint,directory/"history.csv",directory/"seed.json"])
