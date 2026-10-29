#!/usr/bin/env python3
"""Weighted electron occupancy only. No crossing flux, reaction rate or time average."""
import argparse,bisect,csv,hashlib,json,math,pathlib,sys
sys.dont_write_bytecode=True
HERE=pathlib.Path(__file__).resolve().parent
DECODER=HERE.parent/"cpu-independent-review/checkpoint-validation"
sys.path.insert(0,str(DECODER))
import compare_saved_pic as state_io

ME=9.1093837139e-31
ENERGY_EDGES=[i*.25 for i in range(81)]+[30.,50.,80.,100.,200.,500.,1000.,2000.,5000.,10000.,20000.,50000.,100000.,200000.,500000.,1000000.]

class Bin:
    def __init__(self):self.weight=0.;self.w2=0.;self.energy=0.;self.macros=0
    def add(self,w,energy):self.weight+=w;self.w2+=w*w;self.energy+=w*energy;self.macros+=1
    def as_dict(self):return {"weighted_electrons":self.weight,"macro_particles":self.macros,"weight_square_sum":self.w2,"weighted_energy_eV":self.energy,"effective_macro_sample_count":self.weight*self.weight/self.w2 if self.w2 else None}

def valid_edges(edges):
    state_io.require(len(edges)>=2 and edges[0]==0.,"energy bins must start at 0 eV")
    for x in edges:state_io.finite(x,"energy edge")
    state_io.require(all(b>a for a,b in zip(edges,edges[1:])),"energy bins must strictly increase")

def vacuum_area(o,x0,x1,y0,y1):
    area=(x1-x0)*(y1-y0)
    if o["finite_foils"]:
        for foil in ("lower_foil","upper_foil"):
            area-=max(0.,min(x1,o[foil+".right_m"])-max(x0,o[foil+".left_m"]))*max(0.,min(y1,o[foil+".top_m"])-max(y0,o[foil+".bottom_m"]))
    state_io.require(area>=-128*state_io.EPS*(x1-x0)*(y1-y0),"negative vacuum area")
    return max(0.,area)

def analyze(state,grid=9,slab_m=.001,edges=None,mu_bins=8,phi_bins=16):
    edges=list(ENERGY_EDGES if edges is None else edges);valid_edges(edges)
    o=state.options
    state_io.require(type(grid)is int and 1<=grid<=128,"invalid diagnostic grid")
    state_io.require(type(mu_bins)is int and type(phi_bins)is int and mu_bins>0 and phi_bins>0,"invalid angular bins")
    state_io.finite(slab_m,"slab depth")
    state_io.require(o["finite_axial_extent"] and 0<slab_m<=o["depth_m"] and o["short_z_m"]<=-.5*slab_m<.5*slab_m<=o["feed_z_m"],"midpoint slab requires the finite axial model and valid bounds")
    hx=o["width_m"]/grid;hy=o["gap_m"]/grid
    energies=[];angles=[];cells=[];global_summaries={}
    for selection,depth in (("full_height",o["depth_m"]),("midpoint_slab",slab_m)):
        percell=[{}for _ in range(grid*grid)];angular=[{}for _ in range(grid*grid)];totals=[Bin()for _ in percell];stationary=[Bin()for _ in percell]
        global_total=Bin();global_stationary=Bin();global_bins={};global_overflow=Bin()
        min_energy=None;max_energy=None
        for p in state.particles:
            state_io.require(p["charge_C"]==-state_io.E_CHARGE and p["mass_kg"]==ME,"occupancy diagnostic supports electrons only; add an explicit species schema before ions")
            w=state_io.finite(p["weight"],"particle weight");state_io.require(w>0,"nonpositive weight")
            u=[state_io.finite(p["proper_velocity_m_s["+str(k)+"]"],"proper velocity")for k in range(3)]
            xyz=[state_io.finite(p["position_m["+str(k)+"]"],"position")for k in range(3)]
            state_io.require(0<=xyz[0]<=o["width_m"] and 0<=xyz[1]<=o["gap_m"] and o["short_z_m"]<=xyz[2]<=o["feed_z_m"],"position outside selected model")
            if selection=="midpoint_slab"and not abs(xyz[2])<.5*slab_m:continue
            energy=state_io.finite(state_io.kinetic(p)/state_io.E_CHARGE,"relativistic kinetic energy eV");state_io.require(energy>=0,"negative energy")
            ix=min(int(xyz[0]/hx),grid-1);iy=min(int(xyz[1]/hy),grid-1);cell=iy*grid+ix
            ebin=bisect.bisect_right(edges,energy)-1
            if energy==edges[-1]:ebin=len(edges)-2
            min_energy=energy if min_energy is None else min(min_energy,energy);max_energy=energy if max_energy is None else max(max_energy,energy)
            percell[cell].setdefault(ebin,Bin()).add(w,energy);global_bins.setdefault(ebin,Bin()).add(w,energy)
            totals[cell].add(w,energy);global_total.add(w,energy)
            if ebin==len(edges)-1:global_overflow.add(w,energy)
            norm=math.sqrt(sum(v*v for v in u))
            if norm==0:stationary[cell].add(w,energy);global_stationary.add(w,energy)
            else:
                mu=max(-1.,min(1.,u[2]/norm));phi=math.atan2(u[1],u[0])%(2*math.pi)
                im=min(int((mu+1)*.5*mu_bins),mu_bins-1);ip=min(int(phi/(2*math.pi)*phi_bins),phi_bins-1)
                angular[cell].setdefault((im,ip),Bin()).add(w,energy)
        volume_sum=0.;vacuum_sum=0.;count_from_density=0.;energy_from_hist=0.;count_from_angular=0.
        for cell,total in enumerate(totals):
            ix,iy=cell%grid,cell//grid;x0,x1=ix*hx,(ix+1)*hx;y0,y1=iy*hy,(iy+1)*hy
            volume=hx*hy*depth*1e6;vacuum=vacuum_area(o,x0,x1,y0,y1)*depth*1e6
            state_io.require(volume>0,"nonpositive sampling volume")
            volume_sum+=volume;vacuum_sum+=vacuum
            summary={"selection":selection,"ix":ix,"iy":iy,"x_lo_m":x0,"x_hi_m":x1,"y_lo_m":y0,"y_hi_m":y1,"nominal_volume_cm3":volume,"vacuum_volume_cm3":vacuum,**total.as_dict(),"stationary_weight":stationary[cell].weight,"mean_energy_eV":total.energy/total.weight if total.weight else None,"n_nominal_cm3":total.weight/volume,"n_vacuum_cm3":total.weight/vacuum if vacuum else None,"occupied_zero_vacuum_volume":vacuum==0. and total.weight>0}
            cells.append(summary)
            for j,v in sorted(percell[cell].items()):
                lo=edges[j];hi=edges[j+1]if j+1<len(edges)else None
                width=hi-lo if hi is not None else None
                n=v.weight/(volume*width)if width else None
                row={"selection":selection,"ix":ix,"iy":iy,"energy_lo_eV":lo,"energy_hi_eV":hi,"overflow":hi is None,**v.as_dict(),"nominal_volume_cm3":volume,"vacuum_volume_cm3":vacuum,"dn_dE_nominal_cm3_per_eV":n,"dn_dE_vacuum_cm3_per_eV":v.weight/(vacuum*width)if vacuum and width else None,"number_pdf_per_eV":v.weight/(total.weight*width)if total.weight and width else None}
                energies.append(row);count_from_density+=n*width*volume if width else v.weight;energy_from_hist+=v.energy
            solid_angle=(2./mu_bins)*(2*math.pi/phi_bins)
            for (im,ip),v in sorted(angular[cell].items()):
                angles.append({"selection":selection,"ix":ix,"iy":iy,"mu_lo":-1+2.*im/mu_bins,"mu_hi":-1+2.*(im+1)/mu_bins,"phi_lo_rad":2*math.pi*ip/phi_bins,"phi_hi_rad":2*math.pi*(ip+1)/phi_bins,"solid_angle_sr":solid_angle,**v.as_dict(),"dn_dOmega_nominal_cm3_per_sr":v.weight/(volume*solid_angle),"dn_dOmega_vacuum_cm3_per_sr":v.weight/(vacuum*solid_angle)if vacuum else None})
                count_from_angular+=v.weight
        count_from_angular+=global_stationary.weight
        count_limit=1e-11*max(1.,global_total.weight)
        energy_limit=1024*state_io.EPS*abs(global_total.energy)
        state_io.require(abs(count_from_density-global_total.weight)<=count_limit,"energy-density count normalization failed")
        state_io.require(abs(count_from_angular-global_total.weight)<=count_limit,"angular plus stationary normalization failed")
        state_io.require(abs(energy_from_hist-global_total.energy)<=energy_limit,"kinetic-energy histogram closure failed")
        global_summaries[selection]={**global_total.as_dict(),"min_energy_eV":min_energy,"max_energy_eV":max_energy,"mean_energy_eV":global_total.energy/global_total.weight if global_total.weight else None,"nominal_volume_cm3":volume_sum,"vacuum_volume_cm3":vacuum_sum,"mean_nominal_density_cm3":global_total.weight/volume_sum,"mean_vacuum_density_cm3":global_total.weight/vacuum_sum if vacuum_sum else None,"overflow":global_overflow.as_dict(),"finite_support_complete":global_overflow.weight==0.,"stationary_weight":global_stationary.weight,"energy_integrated_weight":count_from_density,"angular_plus_stationary_weight":count_from_angular,"histogram_energy_eV":energy_from_hist,"count_normalization_absolute_limit":count_limit,"energy_closure_absolute_limit_eV":energy_limit,"global_energy_bins":[{"energy_lo_eV":edges[j],"energy_hi_eV":edges[j+1]if j+1<len(edges)else None,**v.as_dict(),"dn_dE_nominal_cm3_per_eV":v.weight/(volume_sum*(edges[j+1]-edges[j]))if j+1<len(edges)else None,"dn_dE_vacuum_cm3_per_eV":v.weight/(vacuum_sum*(edges[j+1]-edges[j]))if vacuum_sum and j+1<len(edges)else None}for j,v in sorted(global_bins.items())]}
    result={"diagnostic":"instantaneous weighted electron occupancy only","crossing_flux_computed":False,"residence_time_weighted":False,"chemistry_rates_computed":False,"physical_prediction":False,"angular_reference":"global +z: mu=u_z/|u|, phi=atan2(u_y,u_x); not pitch angle to B or incidence at a wall","source_birth_provenance_available":False,"energy_edges_eV":edges,"last_finite_energy_bin_upper_inclusive":True,"slab_selection":"strict abs(z)<slab_m/2, centered at z=0","slab_m":slab_m,"full_height_m":o["depth_m"],"spatial_grid":[grid,grid],"mu_bins":mu_bins,"phi_bins":phi_bins,"sparse_histogram_rows":"Unlisted bins have zero weight; empty-cell probability density is undefined. Overflow is retained separately without a fictitious finite bin width.","global":global_summaries,"cells":cells,"energy_rows":energies,"angle_rows":angles}
    def check_numbers(x):
        if isinstance(x,float):state_io.finite(x,"derived occupancy output")
        elif isinstance(x,dict):
            for v in x.values():check_numbers(v)
        elif isinstance(x,list):
            for v in x:check_numbers(v)
    check_numbers(result);return result

def write_result(result,out):
    out.mkdir(parents=True,exist_ok=False)
    for key,name in (("energy_rows","energy_occupancy.csv"),("angle_rows","angular_occupancy.csv")):
        rows=result.pop(key)
        if rows:
            with (out/name).open("w",newline="")as f:writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/"occupancy.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")

def from_directory(directory,grid=9):
    state=state_io.read_state(directory,state_io.schema());result=analyze(state,grid=grid)
    progress=state.directory/"progress.json"
    if progress.is_file():
        p=state_io.strict_json(progress);global_time=state_io.finite(p["start_time_s"]+state.ledger["time_s"],"global clock")
        state_io.require(abs(global_time-p["accepted_time_s"])<=64*state_io.EPS*max(abs(global_time),1e-12),"progress/checkpoint clock mismatch")
        time_provenance="progress.start_time_s + checkpoint.time_s"
    else:
        state_io.require(state.seed_info["checkpoint"]==state.checkpoint.name and state.seed_info["checkpoint_sha256"]==state_io.sha(state.checkpoint),"global timestamp unavailable")
        global_time=state.seed_info["start_ns"]*1e-9;time_provenance="seed.start_ns, bound by checkpoint SHA-256"
    result.update({"global_time_s":global_time,"time_provenance":time_provenance,"input_checkpoint":str(state.checkpoint),"input_checkpoint_sha256":state_io.sha(state.checkpoint),"decoder_sha256":state_io.sha(DECODER/"compare_saved_pic.py")})
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("directory",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path,required=True);p.add_argument("--grid",type=int,default=9);a=p.parse_args()
    try:
        with state_io.cpu_lease()as execution:r=from_directory(a.directory,a.grid);r["execution"]=execution;write_result(r,a.out)
    except (ValueError,OSError,KeyError,TypeError,IndexError,OverflowError,ZeroDivisionError)as exc:print(str(exc),file=sys.stderr);return 2
    print(str(a.out/"occupancy.json"));return 0
if __name__=="__main__":raise SystemExit(main())
