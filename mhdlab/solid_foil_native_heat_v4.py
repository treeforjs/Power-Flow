"""Sealed heat/source-v2 plus signed mechanics, in a separate native v4 worker.

The source fields and event samples are explicit caller inputs held at accepted
start. This wrapper owns no second particle population or adsorbed reservoir.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from .solid_foil_native_heat import NativeConductiveHeat, NativeCommitUncertain, sha, quote, numbers

BUDGET_FIELDS=('thermal_delta_J','mechanical_exchange_J','mechanical_absolute_activity_J',
    'joule_J','authored_heat_J','wall_net_J','residual_J','roundoff_bound_J','accepted_reference_energy_J',
    'trial_reference_energy_J','nodal_residual_max_abs_J','mechanical_load_conservation_residual_J','newton_iterations',
    'source_wall_J','source_particle_K_change_J','source_electron_reference_change_J','source_image_energy_change_J',
    'source_field_work_J','neutral_K_J','neutral_internal_J','adsorbed_reference_change_J','gas_reference_J',
    'particle_charge_change_C','electrode_charge_change_C','physical_electron_number_change','new_gas_packet_count','new_path_count')
ACTION_FIELDS=('kind','layer_id','patch_index','patch_id','weight','normal_uniform','transverse_uniform','azimuth_uniform',
    'origin_x_m','origin_y_m','origin_z_m','duration_s','work_function_eV','energy_panels','action_nodes','tail_widths',
    'local_linear_vacuum_extent_m','image_force_fraction','image_energy_fraction','maximum_quantum_magnetic_energy_fraction',
    'maximum_gyroangle','maximum_relative_normal_step','maximum_substeps','PIC_applies_microscopic_image_force',
    'own_return_unit_accommodation','incoming_weighted_K_J','collection_accommodation_admitted')
LAYER_FIELDS=('layer_id','patch_index','status','weight','time_s','origin_x_m','origin_y_m','origin_z_m',
    'position_x_m','position_y_m','position_z_m','proper_vx_m_s','proper_vy_m_s','proper_vz_m_s',
    'current_weighted_K_J','work_function_eV','normal_x','normal_y','normal_z')
PATH_FIELDS=('layer_id','from_x_m','from_y_m','from_z_m','to_x_m','to_y_m','to_z_m','start_time_s','end_time_s','charge_C','field_work_J')
GAS_FIELDS=('patch_index','source_revision','weight','T_K','time_s','x_m','y_m','z_m','proper_vx_m_s','proper_vy_m_s',
    'proper_vz_m_s','support_impulse_x_N_s','support_impulse_y_N_s','support_impulse_z_N_s','mass_kg','charge_C',
    'kinetic_J','internal_J','reference_J','adsorbed_reference_change_J')
RESERVOIR_FIELDS=('adsorbed_N_hi','adsorbed_N_lo','gas_N_hi','gas_N_lo','adsorbed_reference_energy_J')
DRAW_FIELDS=('patch_index','normal','radius','angle','along','z','time')
INTEGER_ACTIONS={'kind','layer_id','patch_index','patch_id','energy_panels','action_nodes','maximum_substeps',
    'PIC_applies_microscopic_image_force','own_return_unit_accommodation','collection_accommodation_admitted'}


def electron_action(*,kind,layer_id,patch_index,patch_id,origin_m,weight,duration_s,
                    work_function_eV,local_linear_vacuum_extent_m,normal_uniform=.6,transverse_uniform=.35,
                    azimuth_uniform=.1,incoming_weighted_K_J=0.,collection_accommodation_admitted=False,
                    maximum_quantum_magnetic_energy_fraction=.2,image_energy_fraction=.03):
    """Explicit component event with numerical defaults from sealed source API.

    Weight is supplied by the caller, never inferred as a predicted event count.
    The clean-surface work function and valid local field extent are mandatory.
    """
    return dict(zip(ACTION_FIELDS,(kind,layer_id,patch_index,patch_id,weight,normal_uniform,transverse_uniform,azimuth_uniform,
        *origin_m,duration_s,work_function_eV,128,64,48.,local_linear_vacuum_extent_m,.03,image_energy_fraction,
        maximum_quantum_magnetic_energy_fraction,.02,.02,4096,False,True,incoming_weighted_K_J,collection_accommodation_admitted)))


def _row(fields,data,integers=()):
    if set(data)!=set(fields):raise ValueError('native source row requires exact explicit fields')
    result=[]
    for k in fields:
        v=data[k]
        if k in integers:
            if not isinstance(v,(int,np.integer,bool)) and (not np.isfinite(v) or float(v)!=int(v)):
                raise ValueError('integer native identity/configuration required')
            if not 0<=int(v)<2**64:raise ValueError('native unsigned identity/configuration overflow')
            result.append(str(int(v)))
        else:
            if not np.isfinite(v):raise ValueError('finite native source row required')
            result.append(format(float(v),'.17g'))
    return ' '.join(result)


class NativeSourceConductiveHeat(NativeConductiveHeat):
    """Whole native heat/source/gas/path participant; immutable v3 helpers reused."""
    def __init__(self,executable_path,*,expected_executable_sha256,nx,ny,x0_m,y0_m,width_m,thickness_m,height_m,
                 initial_temperature_K,reference_density_kg_m3,mesh_revision,caloric_sha256,source_sha256,
                 electrode_id,patches,electron_state_id,electron_parent_id,initial_time_s=0.):
        executable=Path(executable_path).resolve()
        if sha(executable)!=expected_executable_sha256:raise ValueError('native v4 executable seal mismatch')
        self._process=None;self._pending=None;self._candidate=None;self._stderr=None
        self.nx=int(nx);self.ny=int(ny);self.height_m=float(height_m)
        initial=np.ascontiguousarray(initial_temperature_K,dtype=float)
        if initial.shape!=(ny+1,nx+1) or not np.isfinite(initial).all():raise ValueError('full canonical nodal initial state required')
        self._initial=initial.copy();patches=json.loads(json.dumps(patches,allow_nan=False))
        self._constructor=dict(executable_path=executable,expected_executable_sha256=expected_executable_sha256,
            nx=nx,ny=ny,x0_m=x0_m,y0_m=y0_m,width_m=width_m,thickness_m=thickness_m,height_m=height_m,
            initial_temperature_K=self._initial,reference_density_kg_m3=reference_density_kg_m3,
            mesh_revision=mesh_revision,caloric_sha256=caloric_sha256,source_sha256=source_sha256,electrode_id=electrode_id,
            patches=patches,electron_state_id=electron_state_id,electron_parent_id=electron_parent_id,initial_time_s=initial_time_s)
        self.x_edges_m=np.linspace(x0_m,x0_m+width_m,nx+1);self.y_edges_m=np.linspace(y0_m,y0_m+thickness_m,ny+1)
        self._binding={'schema':'native-sealed-source-thermoelastic-v4','executable_path':str(executable),
            'executable_sha256':expected_executable_sha256,'mesh_revision':mesh_revision,'caloric_sha256':caloric_sha256,
            'source_sha256':source_sha256,'electrode_id':electrode_id,'reference_density_kg_m3':reference_density_kg_m3,
            'nx':nx,'ny':ny,'x0_m':x0_m,'y0_m':y0_m,'width_m':width_m,'thickness_m':thickness_m,'height_m':height_m,
            'initial_time_s':initial_time_s,'initial_nodal_sha256':hashlib.sha256(initial.astype('<f8').tobytes()).hexdigest(),
            'patches':patches,'electron_state_id':electron_state_id,'electron_parent_id':electron_parent_id,
            'source_stage':'accepted start T and total E/B, explicit component admissions',
            'source_ownership':'one native bound-layer state per electron; one canonical adsorbed pool per disjoint patch',
            'gas_reference_J':0.,'authored_heat_J':0.,'Joule_heat_J':0.,'field_solve_performed':False,
            'source_surface_approximation':'stationary fixed reference surface; support impulse diagnostic, no moving-wall/recoil work'}
        patch_rows=[]
        patch_numbers=('face','begin','end','water_enabled','water_range','initial_coverage','registry_mass_kg')
        for p in patches:
            patch_rows.append(' '.join(quote(p[k]) for k in ('id','gas_state_id','parent_id','scenario_sha256'))+' '+
                _row(patch_numbers,{k:p[k] for k in patch_numbers},{'face','water_enabled','water_range'}))
        self._stderr=tempfile.TemporaryFile(mode='w+t')
        self._process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self._stderr,text=True,bufsize=1)
        try:
            self._call('CREATE 4 '+str(nx)+' '+str(ny)+' '+numbers([x0_m,y0_m,width_m,thickness_m,height_m,
                reference_density_kg_m3,float(initial.flat[0]),initial_time_s])+' '+
                ' '.join(map(quote,(mesh_revision,caloric_sha256,source_sha256,electrode_id)))+' '+str(initial.size)+' '+numbers(initial)+' '+
                str(len(patches))+' '+' '.join(patch_rows)+' '+quote(electron_state_id)+' '+quote(electron_parent_id))
            coords=np.fromstring(self._call('COORD'),sep=' ').reshape(-1,2);xx,yy=np.meshgrid(self.x_edges_m,self.y_edges_m)
            if not np.allclose(coords,np.column_stack((xx.ravel(),yy.ravel())),rtol=2e-15,atol=0):raise ValueError('native grid mismatch')
        except BaseException:self.close();raise

    def state(self,trial=False):
        info=self._call('INFO '+str(int(trial))).split()
        if len(info)!=14 or int(info[2])!=(self.nx+1)*(self.ny+1) or int(info[3])!=self.nx*self.ny:
            raise ValueError('native v4 state-info shape mismatch')
        t,e=self._arrays(trial)
        return {'revision':int(info[0]),'time_s':float(info[1]),'temperature_K':t,'cell_energy_J':e,
            'stored_energy_J':float(info[5]),'state_digest':info[9],'cumulative_mechanical_exchange_J':float(info[6]),
            'patch_count':int(info[10]),'layer_count':int(info[11]),'gas_packet_count':int(info[12]),'path_count':int(info[13])}

    def _source_rows(self,command,trial,fields,count):
        values=self._call(command+' '+str(int(trial))).split()
        if len(values)!=count*len(fields):raise ValueError('native source getter dimensions changed')
        integer_fields={'layer_id','patch_index','status','source_revision'}
        return [{k:(int(v) if k in integer_fields else float(v)) for k,v in
            zip(fields,values[i*len(fields):(i+1)*len(fields)])} for i in range(count)]

    def source_state(self,trial=False):
        info=self.state(trial)
        cumulative=np.fromstring(self._call('CUMULATIVE '+str(int(trial))),sep=' ')
        if cumulative.shape!=(23,):raise ValueError('v2 nineteen-term plus mechanics/charge cumulative ledger required')
        return {'revision':info['revision'],'time_s':info['time_s'],'state_digest':info['state_digest'],
            'layers':self._source_rows('LAYERS',trial,LAYER_FIELDS,info['layer_count']),
            'paths':self._source_rows('PATHS',trial,PATH_FIELDS,info['path_count']),
            'gas_packets':self._source_rows('GAS',trial,GAS_FIELDS,info['gas_packet_count']),
            'reservoirs':self._source_rows('RESERVOIRS',trial,RESERVOIR_FIELDS,info['patch_count']),
            'cumulative':cumulative.tolist(),'patches':self.binding['patches']}

    def prepare(self,dt_s,cell_heat_J,*,transaction_id,joint_base_digest,stage_time_s,source_stage,source_actions,water_draws):
        self.discard();state=self.state();base=state['state_digest'];native_id=base+'/mechanical/'+str(state['revision']+1)
        q=np.ascontiguousarray(cell_heat_J,dtype=float)
        if q.shape!=(self.ny,self.nx) or not np.isfinite(q).all():raise ValueError('finite complete signed cell-J required')
        # Source request data is deep-copied by serialization before entering native code.
        stage=json.loads(json.dumps(source_stage,allow_nan=False));em=stage['field_sha256']
        stage_fields=('common_state','field_state','source_mesh_revision','source_material_revision','field_time_s',
                      'E_x_V_m','E_y_V_m','E_z_V_m','B_x_T','B_y_T','B_z_T',
                      'material_phase_admitted','planar_barrier_admitted','electron_supply_admitted')
        stage_text=_row(stage_fields,{k:stage[k] for k in stage_fields},set(stage_fields[:4]+stage_fields[-3:]))
        action_text=' '.join(_row(ACTION_FIELDS,a,INTEGER_ACTIONS) for a in source_actions)
        draw_text=' '.join(_row(DRAW_FIELDS,d,{'patch_index'}) for d in water_draws)
        response=self._call('PREP '+str(state['revision'])+' '+str(state['revision']+1)+' '+numbers([state['time_s'],dt_s,stage_time_s])+' '+
            ' '.join(map(quote,(self.binding['mesh_revision'],self.binding['caloric_sha256'],self.binding['source_sha256'],em,base,native_id)))+' '+
            stage_text+' '+str(q.size)+' '+numbers(q)+' '+str(len(source_actions))+' '+action_text+' '+str(len(water_draws))+' '+draw_text)
        self._pending=(native_id,base)
        try:
            values=response.split()
            if len(values)!=len(BUDGET_FIELDS):raise ValueError('native v4 budget field count mismatch')
            b={k:float(v) for k,v in zip(BUDGET_FIELDS,values)}
            for k in ('newton_iterations','new_gas_packet_count','new_path_count'):b[k]=int(b[k])
            candidate=self.state(True)
            ledger=np.fromstring(self._call('LEDGER'),sep=' ').reshape(7,self.ny+1,self.nx+1)
            copied=np.fromstring(self._call('MCELL'),sep=' ').reshape(q.shape)
            if not np.array_equal(copied,q):raise ValueError('native signed cell load changed')
            b['native_transaction_id']=native_id;b['caller_transaction_id']=transaction_id
            b['nodal_ledger']={k:a.tolist() for k,a in zip(('delta_H_J','mechanical_load_J','conduction_out_J','residual_J',
                'external_load_J','electron_load_J','water_load_J'),ledger)}
            if np.max(np.abs(ledger[3]))>b['roundoff_bound_J']:raise ValueError('native local nodal balance failed')
            self._candidate=candidate['state_digest']
            return dict(candidate,budget=b)
        except BaseException:self.discard();raise

    def commit(self):
        if self._pending is None:raise ValueError('no native v4 trial to commit')
        try:self._call('COMMIT '+' '.join(map(quote,self._pending)))
        except (ValueError,OSError) as error:
            try:proven=self._candidate is not None and self.state()['state_digest']==self._candidate
            except (ValueError,OSError):proven=False
            if not proven:raise NativeCommitUncertain('native v4 commit outcome unknown; recreate joint owner from checkpoint') from error
        self._pending=None;self._candidate=None
