#!/usr/bin/env python3
"""Reproduce bounded SS304 benchmarks and matched-input admission diagnostics.

No field solve, thermal evolution, production launch or coupled state mutation.
Run with tools/test_ss304_thermal_expansion.sh for frozen CPU7 resource policy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mhdlab.thermal_expansion import load_ss304_material, static_bar, transient_bar
from mhdlab.thermal_expansion_mapping import (
    affine_foil_geometry, projected_gap_y_m, conservative_remap_1d, qualify_field_handoff,
)
from mhdlab.thermal_expansion_source import load_matched_history

DEFAULT_RUN = "/home/treefor/projects/Power-Flow/runs/quick_heating_256x32_domain40mm_height25mm_20261029"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n")


def style():
    plt.rcParams.update({"figure.facecolor":"white", "axes.facecolor":"white", "savefig.facecolor":"white",
                         "font.size":10, "axes.grid":True, "grid.alpha":.22, "axes.spines.top":False,
                         "axes.spines.right":False, "figure.dpi":130})


def save(fig, out, name):
    fig.savefig(out/name, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run(out, source_run):
    out.mkdir(parents=True, exist_ok=False)
    style()
    material_path = ROOT/"data/materials/ss304-v1.json"
    record = json.loads(material_path.read_text()); mat=load_ss304_material(material_path)
    stress_free=300.; lo,hi=mat.admissible_temperature_K
    rho=mat.density_at_stress_free_temperature(stress_free)
    temps=np.linspace(stress_free,hi,150)
    strains=mat.strain(temps,stress_free_temperature_K=stress_free)
    gaps=[]
    for eps in strains:
        lower=affine_foil_geometry(np.eye(3)*(1+eps), center_reference_m=(0,-.0026,0),reference_density_kg_m3=rho)
        upper=affine_foil_geometry(np.eye(3)*(1+eps), center_reference_m=(0,.0026,0),reference_density_kg_m3=rho)
        gaps.append(projected_gap_y_m(lower,upper))
    transfer=conservative_remap_1d([0,.0001,.0002], [0,.00005,.00015,.0002], [rho,rho], kind="density")
    uniform_receipt=qualify_field_handoff(geometries=[lower,upper],temperatures_k=[hi,hi],
        admissible_temperature_range_k=[lo,hi],transfers=[transfer],
        mechanics_scope="uniform_affine_free_benchmark",
        support_description="Each disconnected foil expands about its own unchanged mass center; hypothetical free equilibrium benchmark, no feed/short restraint.")
    write_json(out/"uniform-geometry-interface.json",uniform_receipt)
    fig,ax=plt.subplots(1,3,figsize=(12.7,3.6))
    ax[0].plot(temps,strains*.025/2*1e6,label="axial end (+12.5 mm)")
    ax[0].plot(temps,strains*.005/2*1e6,label="width edge (+2.5 mm)")
    ax[0].set(ylabel="Free equilibrium displacement (µm)",xlabel="Uniform temperature (K)");ax[0].legend(fontsize=8)
    ax[1].plot(temps,(.005-np.array(gaps))*1e6,color="#b24720")
    ax[1].set(ylabel="Gap narrowing (µm)",xlabel="Uniform temperature (K)",title="Own-center free benchmark")
    ax[2].plot(temps,-mat.modulus(temps)*strains/1e6,label="uniaxial restrained trial")
    proof=record["proof_stress_screening"]
    ax[2].plot(temps,-np.interp(temps,proof["temperature_K"],proof["value_Pa"])/1e6,"--",label="annealed bar proof-stress proxy")
    ax[2].set(ylabel="Axial stress (MPa)",xlabel="Uniform temperature (K)");ax[2].legend(fontsize=8)
    fig.suptitle("SS304 solid benchmarks • stress-free reference 300 K • no fixture or plasticity model")
    fig.tight_layout();save(fig,out,"uniform-displacement-gap.png")

    history=load_matched_history(source_run,lower_volume="mirror_with_surface_validation")
    write_json(out/"input-manifest.json",history.manifest)
    admission=history.admissibility(lo,hi)
    face_max=np.max(np.stack([v.max(axis=1) for v in history.face_temperature_k.values()]),axis=0)
    face_min=np.min(np.stack([v.min(axis=1) for v in history.face_temperature_k.values()]),axis=0)
    bad_faces=np.flatnonzero((face_max>hi)|(face_min<lo))
    first_bad_face=float(history.face_time_s[bad_faces[0]]) if len(bad_faces) else np.inf
    bulk_bad=((history.upper_temperature_k<lo)|(history.upper_temperature_k>hi)).any(axis=(1,2))
    safe=np.flatnonzero((~np.maximum.accumulate(bulk_bad)) & (history.time_s<first_bad_face))
    if len(safe)<2:
        raise ValueError("fewer than two fully admitted source frames; no source mechanics benchmark")
    stop=int(safe[-1]); source_times=history.time_s[:stop+1]; bulk=history.upper_temperature_k[:stop+1]
    admission.update({"maximum_observed_surface_K":float(face_max.max()),
                      "source_benchmark_last_time_s":float(source_times[-1]),
                      "no_continuation_past_first_invalid_observation":True,
                      "physical_field_handoff_ready":False,
                      "mechanical_supports":"unspecified", "stress_free_temperature_K":stress_free})
    melting=record["melting"]["reported_range_K"][0]
    melt_faces=np.flatnonzero(face_max>=melting)
    admission["first_surface_possible_melt_time_s"]=float(history.face_time_s[melt_faces[0]]) if len(melt_faces) else None
    write_json(out/"admission.json",admission)
    fig,ax=plt.subplots(1,2,figsize=(11,3.8))
    for a in ax:
        a.plot(history.face_time_s*1e9,face_max,label="all directly recorded faces",color="#b24720")
        a.plot(history.time_s*1e9,history.upper_temperature_k.max(axis=(1,2)),label="upper bulk cell averages",color="#246a9e")
        a.axhline(hi,color="black",ls="--",label="elastic data limit 773.15 K")
        a.axvline(source_times[-1]*1e9,color="#398356",ls=":",label="last admitted common frame")
        a.set(xlabel="Source clock (ns)",ylabel="Maximum temperature (K)")
    ax[0].axhline(melting,color="#954a99",ls="--",label="304 possible melt onset")
    ax[0].set_title("Prescribed heating has no phase changes")
    ax[1].set(xlim=(0,22),ylim=(280,1100),title="No clipping or elastic extrapolation")
    ax[0].legend(fontsize=7);fig.tight_layout();save(fig,out,"source-admissibility.png")

    # Independent through-thickness bars at two width locations. No lateral
    # coupling, plate bending, Poisson response or mechanical short connection.
    hottest=int(np.unravel_index(np.argmax(bulk[-1]),bulk[-1].shape)[1])
    columns={"hottest_saved_column":hottest,"midwidth_column":bulk.shape[2]//2}
    xnodes=history.depth_edges_m
    summaries={}; curves={}; ref_errors={}
    def solve_column(j,ncells,max_dt):
        count=int(np.ceil((source_times[-1]-source_times[0])/max_dt))
        times=np.linspace(source_times[0],source_times[-1],count+1)
        # Conservative overlap in material coordinate for prescribed-T sampling.
        # Integral of T is NOT asserted to conserve heat; caloric EOS absent.
        edges=np.linspace(xnodes[0],xnodes[-1],ncells+1)
        mapped=conservative_remap_1d(xnodes,edges,bulk[:,:,j].T,kind="density").values.T
        temp=np.stack([np.interp(times,source_times,mapped[:,k]) for k in range(ncells)],axis=1)
        result=transient_bar(edges,times,temp,material=mat,stress_free_temperature_K=stress_free,
                             boundary="free",area_m2=(history.x_edges_m[j+1]-history.x_edges_m[j])*.025)
        return edges,temp,result
    fig,axes=plt.subplots(1,3,figsize=(13.4,3.9))
    for label,j in columns.items():
        edges,temp,result=solve_column(j,32,5e-11)
        static=static_bar(edges,temp[-1],material=mat,stress_free_temperature_K=stress_free,boundary="free")
        gap=.005+2*result["displacement_m"][:,0]
        ratio=np.abs(result["stress_Pa"])/np.interp(temp,proof["temperature_K"],proof["value_Pa"])
        first=np.flatnonzero(ratio.max(axis=1)>1)
        summaries[label]={"width_column":j,"x_center_m":float((history.x_edges_m[j]+history.x_edges_m[j+1])/2),
            "maximum_abs_displacement_m":float(np.max(abs(result["displacement_m"]))),
            "final_gap_m":float(gap[-1]),"maximum_abs_elastic_trial_stress_Pa":float(np.max(abs(result["stress_Pa"]))),
            "maximum_annealed_bar_proof_proxy_ratio":float(ratio.max()),
            "first_proof_proxy_exceedance_s":float(result["times_s"][first[0]]) if len(first) else None,
            "mass_relative_error":float(np.max(abs(result["current_mass_kg"]/result["reference_mass_kg"]-1))),
            "production_deformation_admissible":False}
        curves[label]={"time_s":result["times_s"].tolist(),"gap_m":gap.tolist(),
                       "inner_surface_displacement_m":result["displacement_m"][:,0].tolist()}
        tag="edge" if label.startswith("hottest") else "midwidth"
        axes[0].plot(edges*1e6,result["displacement_m"][-1]*1e9,label=f"{tag}: transient")
        axes[0].plot(edges*1e6,static["displacement_m"]*1e9,"--",label=f"{tag}: equilibrium")
        axes[1].plot(result["times_s"]*1e9,(.005-gap)*1e9,label=tag)
        axes[2].plot(result["times_s"]*1e9,np.max(abs(result["stress_Pa"]),axis=1)/1e6,label=tag)
        if len(first):
            for ax in axes[1:]:
                ax.axvline(result["times_s"][first[0]]*1e9,color="#aa2222",ls=":",lw=1,
                           label=f"{tag}: proof proxy exceeded")
        _,_,coarse_t=solve_column(j,32,1e-10)
        _,_,fine_x=solve_column(j,64,5e-11)
        ref_errors[label]={"time_dt_0p1_to_0p05_ns_final_displacement_Linf_m":float(np.max(abs(coarse_t["displacement_m"][-1]-result["displacement_m"][-1]))),
                           "mesh_32_to_64_final_displacement_Linf_m":float(np.max(abs(fine_x["displacement_m"][-1,::2]-result["displacement_m"][-1]))),
                           "interpretation":"solver refinement on fixed piecewise native source; not source-grid or physical convergence"}
    axes[0].set(xlabel="Depth from upper gap face (µm)",ylabel="Displacement (nm)")
    axes[1].set(xlabel="Source clock (ns)",ylabel="Symmetric gap narrowing (nm)")
    axes[2].set(xlabel="Source clock (ns)",ylabel="Maximum elastic trial stress (MPa)")
    for ax in axes:ax.legend(fontsize=7)
    fig.suptitle(f"Elastic thickness-bar trials through {source_times[-1]*1e9:.3f} ns • edge exceeds proof proxy • no plasticity")
    fig.tight_layout();save(fig,out,"transient-displacement-gap.png")
    write_json(out/"source-bar-curves.json",curves)
    speed=float(np.sqrt(mat.modulus(stress_free)/rho))
    report={"schema":"ss304-thermal-expansion-benchmarks-v1","stress_free_temperature_K":stress_free,
        "material_record_sha256":sha(material_path),"source_run_sha256":history.manifest["source_run_sha256"],
        "solid_admissible_temperature_K":[lo,hi],"uniform_max_temperature_strain":float(strains[-1]),
        "uniform_max_temperature_gap_m":gaps[-1],"uniform_reference_mass_per_foil_kg":lower.mass_kg,
        "extensional_wave_speed_m_s":speed,"wave_crossing_s":{name:dimension/speed for name,dimension in (("thickness",.0002),("width",.005),("axial_height",.025))},
        "admission":admission,"source_bar_benchmarks":summaries,"refinement":ref_errors,
        "temperature_interpolation":"linear time at admitted endpoints; piecewise cell overlap in reference depth; never interpreted as enthalpy conservation",
        "geometry_interface":"affine benchmark only; full electrode meshes, supports and moving-field solves unqualified",
        "thermal_mechanical_work":"elastic and kinetic energies are trial diagnostics; no thermoelastic caloric EOS, no heat debit, no claim of coupled energy closure",
        "physical_coupled_run_ready":False,"production_launched":False,
        "environment":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__,
                       "kernel":platform.release(),"affinity":sorted(os.sched_getaffinity(0)),"nice":os.getpriority(os.PRIO_PROCESS,0),
                       "threads":{k:os.environ.get(k) for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")}},
        "parent_port_requirements_sha256":sha(ROOT/"docs/interfaces/parent-coupling-port-requirements-20261029.json"),
        "heat_feedback_interface_sha256":sha(ROOT/"docs/interfaces/heat-feedback-port-20261029.json")}
    write_json(out/"summary.json",report)
    write_json(out/"output-sha256.json",{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps({"output":str(out),"physical_coupled_ready":False,"last_source_benchmark_ns":source_times[-1]*1e9,
                      "source_bar_benchmarks":summaries},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--source-run",default=DEFAULT_RUN)
    args=parser.parse_args();run(args.output,args.source_run)
