"""Transactional native conductive heat / Q1 thermoelastic corrector.

The accepted reference mesh is fixed. Native signed cell-J is an INTERNAL
exchange, never authored or Joule heat. Native nodal enthalpy is projected to
mechanical cells conservatively. This projection is spatial homogenization;
local heat balance remains the native nodal balance, not a false cellwise
identity after node-to-cell averaging.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import uuid
import numpy as np

from .solid_foil_mechanics import transient_foil, mesh_geometry_sha256, evaluate_foil_state
from .solid_foil_coupling import thermoelastic_exchange, von_mises_Pa
from .solid_foil_native_heat import NativeCommitUncertain
from .solid_foil_native_mapping import make_native_heat_stage_identity, project_native_heat_to_solid, signed_cell_exchange_to_native


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen(a):
    a=np.ascontiguousarray(a,dtype=float)
    return np.frombuffer(a.tobytes(),dtype=float).reshape(a.shape)


@dataclass(frozen=True)
class CoupledTrial:
    token: str
    base_digest: str
    end_displacement_m: np.ndarray
    end_velocity_m_s: np.ndarray
    receipt_json: str

    @property
    def receipt(self):
        return json.loads(self.receipt_json)


class NativeFoilCorrector:
    """One owner, one pending trial, joint preflight then no-throw native swap.

    Python stores only preallocated immutable mechanical arrays at commit. A
    checkpoint is published as one directory after both components are saved;
    its manifest binds the component bytes and accepted stage. No independent
    participant checkpoint may be treated as a joint accepted state.
    """
    def __init__(self, native, mesh, *, caloric, material,
                 stress_free_temperature_K, poisson_ratio, boundary,
                 mechanical_mesh_revision, material_sha256, source_sha256,
                 proof_stress_screen_Pa, proof_stress_screen_sha256, initial_displacement_m=None,
                 initial_velocity_m_s=None):
        self.native=native; self.mesh=mesh; self.caloric=caloric; self.material=material
        self.tref=float(stress_free_temperature_K); self.nu=float(poisson_ratio)
        self.boundary=boundary; self.mesh_revision=mechanical_mesh_revision
        self.material_sha256=material_sha256; self.source_sha256=source_sha256
        self.screen=proof_stress_screen_Pa; self.pending=None; self._indeterminate=False
        if (native.binding['source_sha256']!=source_sha256 or
            native.binding['caloric_sha256']!=caloric.receipt['source_seal_sha256'] or
            native.binding['reference_density_kg_m3']!=caloric.receipt['reference_density_kg_m3'] or
            not np.isclose(caloric.receipt['reference_density_kg_m3'],
                material.density_at_stress_free_temperature(self.tref),rtol=1e-14,atol=0)):
            raise ValueError('native/solid source, caloric seal or reference density differs')
        if (not isinstance(proof_stress_screen_sha256,str) or len(proof_stress_screen_sha256)!=64 or
                any(c not in '0123456789abcdef' for c in proof_stress_screen_sha256)):
            raise ValueError('explicit proof-screen source SHA256 required')
        self.u=frozen(np.zeros_like(mesh.nodes_m) if initial_displacement_m is None else initial_displacement_m)
        self.v=frozen(np.zeros_like(mesh.nodes_m) if initial_velocity_m_s is None else initial_velocity_m_s)
        if self.u.shape!=mesh.nodes_m.shape or self.v.shape!=mesh.nodes_m.shape or not np.isfinite(self.u).all() or not np.isfinite(self.v).all():
            raise ValueError('finite mechanical initial arrays must match mesh')
        self._config={'schema':'native-thermoelastic-joint-v3','solid_mesh_sha256':mesh_geometry_sha256(mesh),
            'mechanical_mesh_revision':mechanical_mesh_revision,'material_sha256':material_sha256,
            'source_sha256':source_sha256,'caloric':caloric.receipt,'native':native.binding,
            'stress_free_temperature_K':self.tref,'poisson_ratio':self.nu,'boundary':boundary,
            'initial_mechanical_state_sha256':digest({'u':self.u.tolist(),'v':self.v.tolist()}),
            'mechanics_temperature_max_K':773.15,'coefficient_knots':'reject landing/crossing',
            'proof_screen':'caller nominal scenario screen; unchanged v2 implementation',
            'proof_screen_sha256':proof_stress_screen_sha256,
            'proof_screen_knot_values_Pa':np.asarray(proof_stress_screen_Pa(material.temperature_K)).tolist()}
        self.config_sha256=digest(self.config)
        # Validate the actual initial mechanical state using the same guards as a step.
        initial=self._project(self.native.state(),float(np.nextafter(self.native.state()['time_s'],np.inf)),stage='accepted')
        self._validate_mechanical(self.u,self.v,initial.cell_temperature_K.reshape(-1))

    @property
    def config(self):return json.loads(json.dumps(self._config))

    def _validate_mechanical(self,u,v,temperature):
        if self.boundary=='free':fixed=np.zeros(len(u),dtype=bool)
        elif self.boundary=='short_plane_clamped':fixed=self.mesh.nodes_m[:,2]==self.mesh.z_edges_m[0]
        elif self.boundary=='both_end_clamped':fixed=np.isin(self.mesh.nodes_m[:,2],self.mesh.z_edges_m[[0,-1]])
        else:raise ValueError('unknown mechanical boundary')
        if np.any(u[fixed]!=0) or np.any(v[fixed]!=0):
            raise ValueError('accepted mechanics violates homogeneous support constraints')
        checked=evaluate_foil_state(self.mesh,u,temperature,material=self.material,
            stress_free_temperature_K=self.tref,poisson_ratio=self.nu)
        screen=np.asarray(self.screen(temperature))
        if screen.shape!=temperature.shape or np.any(screen<=0) or not np.isfinite(screen).all():
            raise ValueError('invalid nominal proof-stress screen')
        if np.any(von_mises_Pa(checked['gauss_stress_Pa'])>screen[:,None]):
            raise ValueError('nominal proof-stress screening exceeded')

    def _project(self,state,end,*,stage,native=None):
        native=self.native if native is None else native
        base=native.state()
        if (state['revision']!=base['revision']+(stage=='trial') or
            state['time_s']!=(end if stage=='trial' else base['time_s'])):
            raise ValueError('native state revision/time does not match requested mapping stage')
        identity=make_native_heat_stage_identity(self.mesh,native_x_edges_m=native.x_edges_m,
            native_y_edges_m=native.y_edges_m,native_height_m=native.height_m,native_mesh_revision=native.binding['mesh_revision'],
            solid_mesh_revision=self.mesh_revision,state_revision=base['revision'],
            interval_s=(base['time_s'],end),stage_time_s=(end if stage=='trial' else base['time_s']),stage=stage,
            material_sha256=self.material_sha256,source_sha256=self.source_sha256,
            caloric_sha256=native.binding['caloric_sha256'])
        return project_native_heat_to_solid(self.mesh,state['temperature_K'],
            native_x_edges_m=native.x_edges_m,native_y_edges_m=native.y_edges_m,
            native_height_m=native.height_m,caloric=self.caloric,
            identity=identity,expected_identity=identity,
            native_cell_base_energy_J=state['cell_energy_J'],native_stored_energy_J=state['stored_energy_J'],
            admissible_temperature_range_K=(298.2,773.15))

    def state_digest(self):
        if self._indeterminate:raise NativeCommitUncertain('joint owner requires checkpoint recovery')
        if self.native.binding!=self._config['native']:raise ValueError('native binding changed after joint admission')
        s=self.native.state()
        return digest({'config_sha256':self.config_sha256,'native_state':s['state_digest'],
            'revision':s['revision'],'time_s':s['time_s'],'u':self.u.tolist(),'v':self.v.tolist()})

    def prepare(self,dt_s,*,transaction_id,max_iterations=40):
        if self._indeterminate:raise NativeCommitUncertain('joint owner requires checkpoint recovery')
        if not isinstance(transaction_id,str) or not transaction_id or not np.isfinite(dt_s) or dt_s<=0:
            raise ValueError('explicit transaction and positive finite interval required')
        if max_iterations<1:raise ValueError('positive iteration count required')
        self.discard()
        base=self.native.state(); base_digest=self.state_digest(); end=base['time_s']+dt_s
        start=self._project(base,end,stage='accepted')
        t0=start.cell_temperature_K.reshape(-1)
        q=np.zeros((len(self.native.y_edges_m)-1,len(self.native.x_edges_m)-1))
        try:
            for iteration in range(max_iterations):
                thermal=self.native.prepare(dt_s,q,transaction_id=transaction_id,
                    joint_base_digest=base_digest,stage_time_s=end)
                projected=self._project(thermal,end,stage='trial')
                t1=projected.cell_temperature_K.reshape(-1)
                mechanics=transient_foil(self.mesh,[base['time_s'],end],np.stack([t0,t1]),
                    material=self.material,stress_free_temperature_K=self.tref,poisson_ratio=self.nu,
                    boundary=self.boundary,initial_displacement_m=self.u,initial_velocity_m_s=self.v)
                packet=thermoelastic_exchange(self.mesh,mechanics,material=self.material,
                    stress_free_temperature_K=self.tref,caloric=self.caloric,state_revision=base['revision'],
                    mesh_revision=self.mesh_revision,transaction_id=transaction_id,
                    material_sha256=self.material_sha256,source_sha256=self.source_sha256,
                    proof_stress_screen_Pa=self.screen)
                mapped=signed_cell_exchange_to_native(self.mesh,packet,native_x_edges_m=self.native.x_edges_m,
                    native_y_edges_m=self.native.y_edges_m,native_height_m=self.native.height_m,
                    identity=projected.identity,expected_identity=projected.identity)
                qnew=np.asarray(mapped['cell_heat_J']); error=qnew-q
                budget=thermal['budget']
                # Local representable heat-state energy floor, valid even at H=0.
                ct=np.maximum(self.caloric.capacity_J_m3_K(base['temperature_K'])*np.abs(base['temperature_K']),
                    self.caloric.capacity_J_m3_K(thermal['temperature_K'])*np.abs(thermal['temperature_K']))
                corner_ct=np.maximum.reduce([ct[:-1,:-1],ct[:-1,1:],ct[1:,:-1],ct[1:,1:]])
                volume=np.diff(self.native.y_edges_m)[:,None]*np.diff(self.native.x_edges_m)[None,:]*self.native.height_m
                representability=768*np.finfo(float).eps*volume*corner_ct + 5e-13*(np.abs(base['cell_energy_J'])+np.abs(thermal['cell_energy_J']))
                local_tolerance=representability+5e-10*(np.abs(qnew)+np.abs(q))
                tolerance=float(local_tolerance.sum())
                joint_residual=float(budget['thermal_delta_J']+sum(packet['cell_elastic_internal_energy_change_J'])+
                    packet['kinetic_energy_change_J'])
                # Compare each signed cell as well as sum: cancelling errors cannot pass.
                if (np.all(np.abs(error)<=local_tolerance) and abs(error.sum())<=tolerance and
                    abs(joint_residual)<=tolerance+budget['roundoff_bound_J']+abs(packet['mechanical_work_residual_J'])):
                    if budget['joule_J']!=0 or budget['authored_heat_J']!=0 or budget['wall_net_J']!=0:
                        raise ValueError('internal mechanical exchange contaminated external heat ledger')
                    receipt={'schema':'native-thermoelastic-corrector-v3','transaction_id':transaction_id,
                        'base_digest':base_digest,'config_sha256':self.config_sha256,'state_revision':base['revision'],
                        'interval_s':[base['time_s'],end],'stage_time_s':end,'corrector_iterations':iteration+1,
                        'mechanical_exchange':packet,'native_budget':budget,
                        'native_cell_heat_J':q.tolist(),'mechanical_cell_binding_residual_J':error.tolist(),
                        'max_local_binding_residual_J':float(np.max(np.abs(error))),
                        'global_internal_plus_kinetic_residual_J':joint_residual,
                        'energy_tolerance_J':tolerance+budget['roundoff_bound_J']+abs(packet['mechanical_work_residual_J']),
                        'cell_binding_tolerance_J':local_tolerance.tolist(),
                        'cell_heat_representability_floor_J':representability.tolist(),
                        'native_node_balance':'deltaH - signed mechanical load + conduction out = 0',
                        'mechanical_cell_balance':'deltaUelastic - stressWork + heatToThermal = 0',
                        'mapping_scope':'lumped nodal caloric energy to cell means; internal spatial homogenization',
                        'initial_mapping':start.receipt(),'trial_mapping':projected.receipt(),
                        'end_temperature_K':t1.tolist(),'native_end_temperature_K':thermal['temperature_K'].tolist(),
                        'accepted_state_mutated':False,'physical_EOS_ready':False,
                        'native_conductive_heat_solved':True,'native_geometry_deformed':False}
                    trial=CoupledTrial(uuid.uuid4().hex,base_digest,frozen(mechanics['displacement_m'][-1]),
                        frozen(mechanics['velocity_m_s'][-1]),json.dumps(receipt,sort_keys=True,allow_nan=False))
                    self.pending=trial
                    return trial
                q=qnew
            raise ValueError('native heat/mechanics corrector failed to converge; joint trial discarded')
        except BaseException:
            self.discard()
            raise

    def commit(self,trial):
        if trial is not self.pending or trial.base_digest!=self.state_digest() or not self.native.preflight():
            raise ValueError('stale or foreign joint trial')
        # Mechanical state is prevalidated and preallocated. Reconcile native
        # acknowledgement before replacing the accepted mechanical arrays.
        try:self.native.commit()
        except NativeCommitUncertain:
            self._indeterminate=True
            raise
        self.u=trial.end_displacement_m; self.v=trial.end_velocity_m_s
        self.pending=None

    def discard(self):
        if self._indeterminate:raise NativeCommitUncertain('joint owner requires checkpoint recovery')
        self.native.discard(); self.pending=None

    def save(self,path):
        if self._indeterminate:raise NativeCommitUncertain('joint owner requires checkpoint recovery')
        if self.pending is not None:raise ValueError('checkpoint requires accepted joint barrier')
        path=Path(path)
        if path.exists():raise ValueError('joint checkpoint destination exists')
        tmp=path.with_name(path.name+'.tmp-'+uuid.uuid4().hex)
        tmp.mkdir(parents=True)
        try:
            self.native.save(tmp/'native.bin')
            record={'schema':'native-thermoelastic-checkpoint-v3','config':self.config,
                'config_sha256':self.config_sha256,'state_digest':self.state_digest(),
                'native_sha256':file_digest(tmp/'native.bin'),'native_state':{
                    k:self.native.state()[k] for k in ('revision','time_s','state_digest')},
                'u':self.u.tolist(),'v':self.v.tolist()}
            record['manifest_sha256']=digest(record)
            (tmp/'joint.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
            os.replace(tmp,path)
        finally:
            if tmp.exists():shutil.rmtree(tmp)

    def load(self,path):
        if self.pending is not None or self.native.state()['revision']!=0:
            raise ValueError('joint restart requires fresh component without trial')
        path=Path(path); record=json.loads((path/'joint.json').read_text())
        seal=record.pop('manifest_sha256',None)
        if seal!=digest(record) or record['config_sha256']!=self.config_sha256 or record['config']!=self.config:
            raise ValueError('joint restart configuration/manifest mismatch')
        if file_digest(path/'native.bin')!=record['native_sha256']:
            raise ValueError('native checkpoint bytes changed')
        u=frozen(record['u']);v=frozen(record['v'])
        if u.shape!=self.u.shape or v.shape!=self.v.shape or not np.isfinite(u).all() or not np.isfinite(v).all():
            raise ValueError('invalid restart mechanical arrays')
        staged=self.native.fresh()
        try:
            staged.load(path/'native.bin')
            ns=staged.state()
            staged_digest=digest({'config_sha256':self.config_sha256,'native_state':ns['state_digest'],
                'revision':ns['revision'],'time_s':ns['time_s'],'u':u.tolist(),'v':v.tolist()})
            if staged_digest!=record['state_digest'] or any(ns[k]!=value for k,value in record['native_state'].items()):
                raise ValueError('joint restart state identity mismatch')
            projected=self._project(ns,float(np.nextafter(ns['time_s'],np.inf)),stage='accepted',native=staged)
            self._validate_mechanical(u,v,projected.cell_temperature_K.reshape(-1))
        except BaseException:
            staged.close()
            raise
        old=self.native
        self.native=staged;self.u=u;self.v=v
        old.close()
