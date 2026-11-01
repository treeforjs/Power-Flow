"""Conditional conservative thermoelastic exchange for the reference heat port.

Positive mechanical debit removes base caloric energy. The debit includes
elastic INTERNAL energy (Phi-T Phi_T), not merely recoverable strain energy.
This adds no species evolution and commits no thermal or geometric state.
"""
from __future__ import annotations

import hashlib
import json
import numpy as np

from .solid_foil_mechanics import transient_foil, isotropic_constitutive_matrix
from .thermoelastic_energy import ss304_thermodynamic_increment
from .solid_foil_geometry import remap_extensive_xy


def _finite(a, name):
    a=np.asarray(a,dtype=float)
    if not a.size or not np.isfinite(a).all():
        raise ValueError(f"{name} must be nonempty and finite")
    return a


def von_mises_Pa(stress):
    s=_finite(stress,"stress")
    if s.shape[-1] != 6:
        raise ValueError("six engineering-conjugate stress entries required")
    return np.sqrt(.5*((s[...,0]-s[...,1])**2+(s[...,1]-s[...,2])**2+(s[...,2]-s[...,0])**2)
                   +3*np.sum(s[...,3:]**2,axis=-1))


def thermoelastic_exchange(mesh, history, *, material, stress_free_temperature_K,
                           caloric, state_revision, mesh_revision, transaction_id,
                           material_sha256, source_sha256, proof_stress_screen_Pa):
    """Prepare an auditable signed cell-J load for a completed mechanical trial.

    The mechanical history must start at the accepted state and end at the
    thermal trial stage. Summed subcycle stress work is used; a single endpoint
    trapezoid cannot replace the path. Coefficient-knot crossing is rejected.
    No native heat driver is mutated, and global acceptance belongs to caller.
    """
    if type(state_revision) is not int or state_revision < 0:
        raise ValueError("nonnegative integer state revision required")
    if any(not isinstance(x,str) or not x for x in (mesh_revision,transaction_id)):
        raise ValueError("mesh and transaction identities required")
    for value in (material_sha256,source_sha256):
        if not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
            raise ValueError("lowercase SHA256 material/source identity required")
    times=_finite(history["times_s"],"times")
    temperature=_finite(history["cell_temperature_K"],"temperature")
    mesh_digest=hashlib.sha256(np.asarray(mesh.nodes_m,dtype='<f8').tobytes()+
                               np.asarray(mesh.cells,dtype='<i8').tobytes()).hexdigest()
    if (history.get("mesh_geometry_sha256")!=mesh_digest or
            tuple(history.get("voigt_components",()))!=("xx","yy","zz","xy","yz","xz") or
            not np.array_equal(history["gauss_weights_m3"],mesh.gauss_weights_m3)):
        raise ValueError("mechanical history mesh/quadrature/Voigt identity mismatch")
    ncell=len(mesh.cells)
    if times.ndim != 1 or len(times)<2 or np.any(np.diff(times)<=0) or temperature.shape!=(len(times),ncell):
        raise ValueError("matched mechanical temperature/clock history required")
    material.check_temperature(temperature)
    if (history.get("material_provenance")!=material.provenance or
            history.get("stress_free_temperature_K")!=stress_free_temperature_K or
            not np.array_equal(history["cell_modulus_Pa"],material.modulus(temperature)) or
            not np.array_equal(history["cell_thermal_strain"],material.strain(temperature,stress_free_temperature_K=stress_free_temperature_K))):
        raise ValueError("mechanical material/stress-free-reference binding mismatch")
    rho=material.density_at_stress_free_temperature(stress_free_temperature_K)
    if not np.isclose(caloric.receipt.get("reference_density_kg_m3",np.nan),rho,rtol=1e-14,atol=0):
        raise ValueError("caloric and mechanical reference density conventions differ")
    screen=_finite(proof_stress_screen_Pa(temperature),"proof stress screen")
    if screen.shape!=temperature.shape or np.any(screen<=0):
        raise ValueError("one positive declared proof screen per cell/stage required")
    equivalent=von_mises_Pa(history["gauss_stress_Pa"])
    if np.any(equivalent>screen[...,None]):
        raise ValueError("nominal proof-stress screening exceeded; plastic response is not supported")
    for knots in (material.temperature_K,material.modulus_temperature_K):
        segments=np.searchsorted(knots,temperature,side="right")-1
        if (np.any(temperature[...,None]==knots) or
                np.any(segments!=segments[0][None,:])):
            raise ValueError("thermodynamic exchange cannot cross or land on coefficient knots")
    stress_work=_finite(history["gauss_stress_work_increment_J"],"stress work").sum(axis=0)
    strains=_finite(history["gauss_strain"],"strain")
    if stress_work.shape!=(ncell,8) or strains.shape!=(len(times),ncell,8,6):
        raise ValueError("mechanics Gauss ordering/shape mismatch")
    unit=isotropic_constitutive_matrix(1.,history["poisson_ratio"])
    elastic=strains-history["cell_thermal_strain"][...,None,None]*np.array([1.,1.,1.,0.,0.,0.])
    checked_stress=np.einsum('ab,teqb,te->teqa',unit,elastic,history['cell_modulus_Pa'],optimize=True)
    if not np.allclose(checked_stress,history['gauss_stress_Pa'],rtol=2e-12,atol=1e-6):
        raise ValueError("mechanical stress disagrees with bound constitutive state")
    checked_work=np.sum(.5*(checked_stress[1:]+checked_stress[:-1])*np.diff(strains,axis=0),axis=-1)*mesh.gauss_weights_m3
    if not np.allclose(checked_work,history['gauss_stress_work_increment_J'],rtol=2e-11,atol=1e-24):
        raise ValueError("local stress-work history differs from recomputed quadrature")
    dkin=float(history["kinetic_energy_J"][-1]-history["kinetic_energy_J"][0])
    mechanical_residual=dkin+float(stress_work.sum())
    mechanical_scale=abs(dkin)+float(np.abs(stress_work).sum())
    if not np.isfinite(mechanical_residual) or abs(mechanical_residual)>1e-22+5e-9*mechanical_scale:
        raise ValueError("mechanical stress-work/kinetic balance is not qualified")
    internal=np.zeros(ncell); phi=np.zeros(ncell); base=np.zeros(ncell); effective=np.zeros(ncell)
    for e in range(ncell):
        for q in range(8):
            v=mesh.gauss_weights_m3[e,q]
            item=ss304_thermodynamic_increment(strains[0,e,q],strains[-1,e,q],
                temperature[0,e],temperature[-1,e],material=material,
                stress_free_temperature_K=stress_free_temperature_K,
                poisson_ratio=history["poisson_ratio"],reference_volume_m3=v,
                base_energy_J_m3=caloric.energy_J_m3,base_capacity_J_m3_K=caloric.capacity_J_m3_K)
            internal[e]+=item["delta_elastic_internal_energy_J"]
            phi[e]+=item["delta_phi_J"]
            base[e]+=item["delta_base_energy_J"]
            effective[e]+=item["state1"]["fixed_strain_capacity_J_m3_K"]*v
    work=stress_work.sum(axis=1)
    debit=internal-work
    arrays={"cell_thermoelastic_debit_J":debit.tolist(),
            "cell_heat_to_thermal_J":(-debit).tolist(),
            "cell_elastic_internal_energy_change_J":internal.tolist(),
            "cell_recoverable_strain_energy_change_J":phi.tolist(),
            "cell_stress_work_J":work.tolist(),"cell_base_caloric_change_J":base.tolist(),
            "cell_fixed_strain_capacity_J_K":effective.tolist()}
    packet={"schema":"solid-thermoelastic-cell-exchange-v2","transaction_id":transaction_id,
        "state_revision":state_revision,"mesh_revision":mesh_revision,
        "mesh_geometry_sha256":mesh_digest,
        "interval_s":[float(times[0]),float(times[-1])],"stage_time_s":float(times[-1]),
        "material_sha256":material_sha256,"source_sha256":source_sha256,
        "temperature_stage":"trial_interval_end; emission remains on accepted interval-start state",
        "coordinate_basis":"reference physical xyz; complete3D cellJ, depth applied once",
        "cell_order":"z,y,x; x fastest","reference_cell_volumes_m3":mesh.reference_volumes_m3.tolist(),
        "caloric_model":dict(caloric.receipt),"energy_model":"U0(T)+Phi(epsilon,T)-T*dPhi_dT",
        "mechanical_boundary":history["boundary"],"external_support_work_J":0.,
        "support_work_basis":"immobile homogeneous constrained nodes or traction-free; no moving supports",
        "kinetic_energy_change_J":dkin,"mechanical_work_residual_J":mechanical_residual,
        "mechanical_plus_debit_residual_J":float(dkin+internal.sum()-debit.sum()),
        "species_mutation":"none","species_event_ids":[],"species_parent_ids":[],
        "deposition_class":"reversible_thermoelastic_internal_exchange; never authored external heat",
        "conditional_thermodynamic_exchange_ready":True,"physical_EOS_ready":False,
        "maximum_nominal_proof_screen_ratio":float(np.max(equivalent/screen[...,None])),
        "proof_screen_basis":"caller-supplied nominal scenario screen; no plasticity or batch calibration claimed",
        "native_thermal_mesh_deformation_ready":False,"accepted_state_mutated":False,**arrays}
    payload=json.dumps(packet,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
    packet["payload_sha256"]=hashlib.sha256(payload).hexdigest()
    return packet


def project_exchange_to_thermal_xy(mesh, packet, *, target_x_edges_m, target_y_edges_m,
                                    expected_revision, expected_mesh_revision):
    """Collapse signed 3-D cell J and conservatively project to a thermal XY mesh.

    Returned cell energy already includes the full axial extent. A 2-D FE
    assembler divides it by that extent once; it must retain the internal
    mechanical-exchange ledger class instead of calling it authored heat.
    No callbacks or accepted thermal states are changed here.
    """
    if packet.get('state_revision')!=expected_revision or packet.get('mesh_revision')!=expected_mesh_revision:
        raise ValueError('stale thermal exchange state/mesh revision')
    plain=dict(packet);digest=plain.pop('payload_sha256',None)
    expected=hashlib.sha256(json.dumps(plain,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    mesh_hash=hashlib.sha256(np.asarray(mesh.nodes_m,dtype='<f8').tobytes()+np.asarray(mesh.cells,dtype='<i8').tobytes()).hexdigest()
    if digest!=expected or packet.get('mesh_geometry_sha256')!=mesh_hash:
        raise ValueError('altered exchange payload or material mesh')
    cells=_finite(packet['cell_heat_to_thermal_J'],'cell heat')
    if cells.shape!=(len(mesh.cells),):raise ValueError('one signed heat amount per material cell required')
    xy=cells.reshape(mesh.cell_shape).sum(axis=0)
    mapped=remap_extensive_xy(mesh.x_edges_m,mesh.y_edges_m,target_x_edges_m,target_y_edges_m,xy)
    return {'cell_heat_J':mapped['values'],'reference_axial_extent_m':float(mesh.z_edges_m[-1]-mesh.z_edges_m[0]),
            'source_total_J':mapped['source_total'],'target_total_J':mapped['target_total'],
            'relative_conservation_error':mapped['relative_conservation_error'],
            'transaction_id':packet['transaction_id'],'state_revision':expected_revision,
            'stage_time_s':packet['stage_time_s'],'source_payload_sha256':digest,
            'ledger_class':'reversible_thermoelastic_internal_exchange',
            'depth_convention':'cell_heat_J is full3D; native2D FE load divides by axial extent exactly once',
            'native_adapter_mutated':False}


def coupled_cell_heat_step(mesh, initial_temperature_K, supplied_heat_J, *,
                           dt_s, material, stress_free_temperature_K, poisson_ratio,
                           boundary, caloric, state_revision, mesh_revision,
                           material_sha256, source_sha256, transaction_id,
                           proof_stress_screen_Pa,
                           initial_displacement_m=None, initial_velocity_m_s=None,
                           accepted_time_s=0.):
    """Small conditional coupled trial with local caloric cells and Q1 mechanics.

    Iterates end T and mechanical deformation to satisfy each cell's first law.
    V*C_eps is a diagonal quasi-Newton preconditioner, not the full coupled
    Jacobian (which also differentiates stress work and mechanical response).
    Only a converged actual first-law residual accepts the trial.
    Heat conduction, Joule generation, native MFEM and emissions are external;
    supplied_heat_J is the already-integrated trial heat per reference cell.
    No partial state is committed. Failure rejects the complete trial.
    """
    t0=_finite(initial_temperature_K,"initial temperature").reshape(-1)
    heat=_finite(supplied_heat_J,"supplied heat").reshape(-1)
    if t0.shape!=(len(mesh.cells),) or heat.shape!=t0.shape or not np.isfinite(dt_s) or dt_s<=0:
        raise ValueError("one T/heat per cell and positive finite dt required")
    volume=mesh.reference_volumes_m3
    t=np.asarray(caloric.temperature_K(caloric.energy_J_m3(t0)+heat/volume))
    last_norm=np.inf
    for iteration in range(20):
        result=transient_foil(mesh,[accepted_time_s,accepted_time_s+dt_s],np.stack([t0,t]),
            material=material,stress_free_temperature_K=stress_free_temperature_K,
            poisson_ratio=poisson_ratio,boundary=boundary,
            initial_displacement_m=initial_displacement_m,initial_velocity_m_s=initial_velocity_m_s)
        packet=thermoelastic_exchange(mesh,result,material=material,
            stress_free_temperature_K=stress_free_temperature_K,caloric=caloric,
            state_revision=state_revision,mesh_revision=mesh_revision,transaction_id=transaction_id,
            material_sha256=material_sha256,source_sha256=source_sha256,
            proof_stress_screen_Pa=proof_stress_screen_Pa)
        residual=np.array(packet["cell_base_caloric_change_J"])+np.array(packet["cell_thermoelastic_debit_J"])-heat
        scale=float(np.abs(heat).sum()+np.abs(packet["cell_base_caloric_change_J"]).sum())
        norm=float(np.max(np.abs(residual)))
        if max(norm,abs(float(residual.sum()))) <= 1e-18+5e-11*scale:
            packet.update({"cell_supplied_heat_J":heat.tolist(),"cell_energy_residual_J":residual.tolist(),
                "coupled_energy_residual_J":float(residual.sum()+packet["mechanical_plus_debit_residual_J"]),
                "corrector_iterations":iteration+1,"end_temperature_K":t.tolist(),
                "trial_scope":"local-cell caloric/Q1 mechanical corrector; not native conductive heat solver"})
            packet.pop("payload_sha256")
            payload=json.dumps(packet,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
            packet["payload_sha256"]=hashlib.sha256(payload).hexdigest()
            return {"mechanics":result,"exchange":packet}
        if iteration>2 and norm>last_norm*1.1:
            raise ValueError("coupled thermoelastic corrector did not contract; trial discarded")
        last_norm=norm
        t=t-residual/np.array(packet["cell_fixed_strain_capacity_J_K"])
    raise ValueError("coupled thermoelastic corrector failed; no accepted mutation")
