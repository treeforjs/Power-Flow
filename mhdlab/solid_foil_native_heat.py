"""Sealed process binding to the separate v3 native MFEM exchange ABI.

The frozen MFEM archive is not position-independent. A persistent executable
links it unchanged; line-framed, roundtrip-precision RPC calls its native C ABI.
Children inherit the caller's CPU7 affinity, nice level and frozen runtime.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np

BUDGET_FIELDS=('thermal_delta_J','mechanical_exchange_J','mechanical_absolute_activity_J',
    'joule_J','authored_heat_J','wall_net_J','residual_J','roundoff_bound_J','accepted_reference_energy_J',
    'trial_reference_energy_J','nodal_residual_max_abs_J','mechanical_load_conservation_residual_J','newton_iterations')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def quote(value):
    value=str(value)
    if any(ord(c)<32 for c in value):raise ValueError('control characters are not allowed in native identities/paths')
    return json.dumps(value,ensure_ascii=False)
def numbers(array):return ' '.join(format(float(v),'.17g') for v in np.asarray(array).ravel())


class NativeCommitUncertain(RuntimeError):
    """Commit acknowledgement and state reconciliation both failed."""


class NativeConductiveHeat:
    """Fixed reference mesh, explicit material/source hashes, zero external heat.

    This ABI exposes only dedicated reversible mechanical cell-J. Its Joule,
    authored, wall/species and gas channels remain exactly zero. The executable
    is a private local numerical participant; no production process is touched.
    """
    def __init__(self,executable_path,*,expected_executable_sha256,nx,ny,x0_m,y0_m,width_m,thickness_m,height_m,
                 initial_temperature_K,reference_density_kg_m3,mesh_revision,caloric_sha256,source_sha256,
                 electrode_id,initial_time_s=0.):
        executable=Path(executable_path).resolve()
        if sha(executable)!=expected_executable_sha256:raise ValueError('native executable seal mismatch')
        self._process=None;self._pending=None;self._candidate=None;self._stderr=None
        self.nx=int(nx);self.ny=int(ny);self.height_m=float(height_m)
        initial=np.ascontiguousarray(initial_temperature_K,dtype=float)
        if initial.shape!=(ny+1,nx+1) or not np.isfinite(initial).all():
            raise ValueError('canonical full native nodal temperatures required')
        self._initial=initial.copy()
        self._constructor=dict(executable_path=executable,expected_executable_sha256=expected_executable_sha256,
            nx=nx,ny=ny,x0_m=x0_m,y0_m=y0_m,width_m=width_m,thickness_m=thickness_m,height_m=height_m,
            initial_temperature_K=self._initial,reference_density_kg_m3=reference_density_kg_m3,
            mesh_revision=mesh_revision,caloric_sha256=caloric_sha256,source_sha256=source_sha256,
            electrode_id=electrode_id,initial_time_s=initial_time_s)
        self.x_edges_m=np.linspace(x0_m,x0_m+width_m,nx+1)
        self.y_edges_m=np.linspace(y0_m,y0_m+thickness_m,ny+1)
        self._binding={'schema':'native-mfem-thermoelastic-abi-v3','executable_path':str(executable),
            'executable_sha256':expected_executable_sha256,'transport':'persistent executable / native C ABI; 17digit scalar roundtrip',
            'mesh_revision':mesh_revision,'caloric_sha256':caloric_sha256,
            'source_sha256':source_sha256,'electrode_id':electrode_id,'reference_density_kg_m3':reference_density_kg_m3,
            'nx':nx,'ny':ny,'x0_m':x0_m,'y0_m':y0_m,'width_m':width_m,'thickness_m':thickness_m,'height_m':height_m,
            'initial_time_s':initial_time_s,'initial_nodal_sha256':hashlib.sha256(initial.astype('<f8').tobytes()).hexdigest(),
            'caloric_discretization':'native lumped nodal reference enthalpy; signed h allowed',
            'external_heat':'none; zero Joule/authored/wall/species; no EM solve'}
        self._stderr=tempfile.TemporaryFile(mode='w+t')
        self._process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,
            stderr=self._stderr,text=True,bufsize=1)
        try:
            self._call('CREATE 3 '+str(nx)+' '+str(ny)+' '+numbers([x0_m,y0_m,width_m,thickness_m,height_m,
                reference_density_kg_m3,float(initial.flat[0]),initial_time_s])+' '+
                ' '.join(map(quote,(mesh_revision,caloric_sha256,source_sha256,electrode_id)))+' '+str(initial.size)+' '+numbers(initial))
            coords=np.fromstring(self._call('COORD'),sep=' ').reshape(-1,2)
            xx,yy=np.meshgrid(self.x_edges_m,self.y_edges_m)
            if not np.allclose(coords,np.column_stack((xx.ravel(),yy.ravel())),rtol=2e-15,atol=0):
                raise ValueError('native reference grid mismatch')
        except BaseException:
            self.close();raise

    @property
    def binding(self):return json.loads(json.dumps(self._binding))

    def _call(self,command):
        process=self._process
        if process is None or process.poll() is not None:raise ValueError('native worker is not running')
        process.stdin.write(command+'\n');process.stdin.flush()
        response=process.stdout.readline()
        if response.startswith('ERR '):raise ValueError(response[4:].rstrip())
        if response.rstrip()=='OK':return ''
        if response.startswith('OK '):return response[3:].rstrip()
        self._stderr.seek(0)
        raise ValueError('native worker protocol failure: '+response[:200]+' '+self._stderr.read()[-1000:])

    def _arrays(self,trial):
        t=np.fromstring(self._call('TEMP '+str(int(trial))),sep=' ').reshape(self.ny+1,self.nx+1)
        e=np.fromstring(self._call('CELL '+str(int(trial))),sep=' ').reshape(self.ny,self.nx)
        return t,e

    def state(self):
        info=self._call('INFO').split()
        if len(info)!=10 or int(info[2])!=(self.nx+1)*(self.ny+1) or int(info[3])!=self.nx*self.ny:
            raise ValueError('native state-info shape mismatch')
        t,e=self._arrays(False)
        return {'revision':int(info[0]),'time_s':float(info[1]),'temperature_K':t,'cell_energy_J':e,
            'stored_energy_J':float(info[5]),'state_digest':info[9],
            'cumulative_mechanical_exchange_J':float(info[6])}

    def prepare(self,dt_s,cell_heat_J,*,transaction_id,joint_base_digest,stage_time_s):
        self.discard();state=self.state();base=state['state_digest']
        native_id=base+'/mechanical/'+str(state['revision']+1)
        q=np.ascontiguousarray(cell_heat_J,dtype=float)
        if q.shape!=(self.ny,self.nx) or not np.isfinite(q).all():raise ValueError('finite signed native cell-J array required')
        zero_em=hashlib.sha256(json.dumps({'external_heat':'zero; no EM solve','joint_base':joint_base_digest},sort_keys=True).encode()).hexdigest()
        response=self._call('PREP '+str(state['revision'])+' '+str(state['revision']+1)+' '+
            numbers([state['time_s'],dt_s,stage_time_s])+' '+
            ' '.join(map(quote,(self.binding['mesh_revision'],self.binding['caloric_sha256'],self.binding['source_sha256'],
                               zero_em,base,native_id)))+' '+str(q.size)+' '+numbers(q))
        self._pending=(native_id,base)
        try:
            values=response.split()
            if len(values)!=len(BUDGET_FIELDS):raise ValueError('native budget ABI field count mismatch')
            b={k:float(v) for k,v in zip(BUDGET_FIELDS,values)};b['newton_iterations']=int(values[-1])
            t,e=self._arrays(True)
            ledger=np.fromstring(self._call('LEDGER'),sep=' ').reshape(4,self.ny+1,self.nx+1)
            copied=np.fromstring(self._call('MCELL'),sep=' ').reshape(q.shape)
            if not np.array_equal(copied,q):raise ValueError('native cell load changed during binding')
            b['native_transaction_id']=native_id;b['caller_transaction_id']=transaction_id
            b['nodal_ledger']={k:a.tolist() for k,a in zip(('delta_H_J','mechanical_load_J','conduction_out_J','residual_J'),ledger)}
            if np.max(np.abs(ledger[-1]))>b['roundoff_bound_J']:
                raise ValueError('native local nodal energy validation failed')
            self._candidate={'revision':state['revision']+1,'time_s':stage_time_s,'temperature_K':t.copy(),
                'cumulative_mechanical_exchange_J':state['cumulative_mechanical_exchange_J']+b['mechanical_exchange_J']}
            return {'revision':state['revision']+1,'time_s':stage_time_s,'temperature_K':t,
                    'cell_energy_J':e,'stored_energy_J':b['trial_reference_energy_J'],'budget':b}
        except BaseException:
            self.discard();raise

    def preflight(self):
        if self._pending is None:return False
        try:self._call('PREFLIGHT '+' '.join(map(quote,self._pending)));return True
        except ValueError:return False

    def commit(self):
        if self._pending is None:raise ValueError('no native trial to commit')
        try:
            self._call('COMMIT '+' '.join(map(quote,self._pending)))
        except (ValueError,OSError) as error:
            # A lost acknowledgement is not evidence of rollback. Reconcile the
            # deterministic candidate only if the native participant is readable.
            try:
                accepted=self.state();candidate=self._candidate
                proven=(candidate is not None and all(accepted[k]==candidate[k] for k in
                    ('revision','time_s','cumulative_mechanical_exchange_J')) and
                    np.array_equal(accepted['temperature_K'],candidate['temperature_K']))
            except (ValueError,OSError):proven=False
            if not proven:raise NativeCommitUncertain('native commit outcome unknown; recreate joint owner from a checkpoint') from error
        self._pending=None;self._candidate=None

    def discard(self):
        if self._pending is not None:
            self._call('DISCARD '+quote(self._pending[0]));self._pending=None;self._candidate=None

    def fresh(self):return type(self)(**self._constructor)
    def save(self,path):self._call('SAVE '+quote(path))
    def load(self,path):self._call('LOAD '+quote(path))
    def close(self):
        process=self._process
        if process is not None:
            if process.poll() is None:
                try:self._call('EXIT')
                except (ValueError,BrokenPipeError):pass
            try:process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.terminate();process.wait(timeout=5)
            for stream in (process.stdin,process.stdout):stream.close()
            self._process=None
        if self._stderr is not None:self._stderr.close();self._stderr=None
    def __del__(self):
        if hasattr(self,'_process'):self.close()
