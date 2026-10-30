"""Actual native MFEM + Q1 mechanics, not a mock thermal corrector."""
import importlib.util
import json
import os
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('native_v3_demo',ROOT/'tools/run_native_thermoelastic_v3.py')
demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)

@pytest.fixture
def factory():
    lib=os.environ.get('SS304_NATIVE_V3_EXECUTABLE')
    if not lib or not Path(lib).is_file():pytest.skip('requires separately built sealed MFEM v3 library')
    cases=[]
    def create(**kw):
        case=demo.make_case(Path(lib),**kw);cases.append(case);return case
    yield create
    for c in cases:c.native.close()


def test_actual_native_conduction_and_signed_work(factory):
    signs=[]
    for sign in (1.,-1.):
        c=factory(velocity_sign=sign);before=c.state_digest()
        trial=c.prepare(1e-10,transaction_id='actual');r=trial.receipt
        assert c.state_digest()==before
        assert r['native_conductive_heat_solved']
        b=r['native_budget'];signs.append(b['mechanical_exchange_J'])
        assert b['joule_J']==b['authored_heat_J']==b['wall_net_J']==0
        assert np.max(np.abs(b['nodal_ledger']['conduction_out_J']))>0
        assert abs(r['global_internal_plus_kinetic_residual_J'])<=r['energy_tolerance_J']
        assert r['mechanical_exchange']['maximum_nominal_proof_screen_ratio']<1
        assert min(r['native_end_temperature_K'][0])<300
        c.commit(trial);assert c.native.state()['revision']==1
    assert signs[0]<0<signs[1]


def test_rollback_failure_and_stale_trial(factory):
    c=factory();before=c.state_digest()
    trial=c.prepare(1e-10,transaction_id='discard');c.discard()
    assert c.state_digest()==before
    with pytest.raises(ValueError,match='stale or foreign'):c.commit(trial)
    with pytest.raises(ValueError,match='failed to converge'):c.prepare(1e-10,transaction_id='too-short',max_iterations=1)
    assert c.state_digest()==before
    assert not c.native.preflight()
    second=c.prepare(1e-10,transaction_id='success');c.commit(second)
    with pytest.raises(ValueError,match='stale or foreign'):c.commit(second)


def test_joint_restart_and_tamper_rejection(factory,tmp_path):
    c=factory();trial=c.prepare(1e-10,transaction_id='first');c.commit(trial);c.save(tmp_path/'checkpoint')
    d=factory();d.load(tmp_path/'checkpoint');assert d.state_digest()==c.state_digest()
    a=c.prepare(1e-10,transaction_id='second');b=d.prepare(1e-10,transaction_id='second')
    assert a.receipt==b.receipt
    c.commit(a);d.commit(b);assert c.state_digest()==d.state_digest()
    path=tmp_path/'checkpoint'/'joint.json';record=json.loads(path.read_text());record['u'][0][0]+=1e-12
    path.write_text(json.dumps(record))
    fresh=factory();before=fresh.state_digest()
    with pytest.raises(ValueError,match='manifest mismatch'):fresh.load(tmp_path/'checkpoint')
    assert fresh.state_digest()==before


def test_unchanged_stress_and_temperature_guards(factory):
    c=factory(velocity_sign=1000.);before=c.state_digest()
    with pytest.raises(ValueError,match='proof-stress|small-strain'):c.prepare(1e-10,transaction_id='overstress')
    assert c.state_digest()==before
    with pytest.raises(ValueError,match='outside solid|temperature|supports only'):
        factory(initial_shift_K=474.)


def test_native_rejects_bad_stage_and_cold_floor(factory):
    c=factory();before=c.state_digest();q=np.zeros((2,4))
    with pytest.raises(ValueError):
        c.native.prepare(1e-10,q,transaction_id='stage',joint_base_digest=before,stage_time_s=2e-10)
    assert c.state_digest()==before
    with pytest.raises(ValueError):
        c.native.prepare(1e-10,np.full((2,4),-1.),transaction_id='cold',joint_base_digest=before,stage_time_s=1e-10)
    assert c.state_digest()==before


def test_stage_and_returned_identity_are_bound(factory):
    c=factory();state=c.native.state()
    state['revision']+=1
    with pytest.raises(ValueError,match='revision/time'):
        c._project(state,1e-10,stage='accepted')
    binding=c.native.binding;binding['source_sha256']='f'*64
    assert c.native.binding['source_sha256']!=binding['source_sha256']
    config=c.config;config['boundary']='unrecognized'
    assert c.config['boundary']=='free'


def test_resealed_invalid_restart_rejected_atomically(factory,tmp_path):
    from mhdlab.solid_foil_native_coupling import digest
    c=factory();t=c.prepare(1e-10,transaction_id='checkpoint');c.commit(t);c.save(tmp_path/'checkpoint')
    p=tmp_path/'checkpoint'/'joint.json';record=json.loads(p.read_text())
    record['u'][0][0]+=1e-5
    record['state_digest']=digest({'config_sha256':record['config_sha256'],
        'native_state':record['native_state']['state_digest'],'revision':record['native_state']['revision'],
        'time_s':record['native_state']['time_s'],'u':record['u'],'v':record['v']})
    record.pop('manifest_sha256');record['manifest_sha256']=digest(record);p.write_text(json.dumps(record))
    fresh=factory();before=fresh.state_digest()
    with pytest.raises(ValueError,match='small-strain|Jacobian|proof-stress'):fresh.load(tmp_path/'checkpoint')
    assert fresh.state_digest()==before


def test_restart_binds_initial_mechanics(factory,tmp_path):
    c=factory();t=c.prepare(1e-10,transaction_id='checkpoint');c.commit(t);c.save(tmp_path/'checkpoint')
    different=factory(velocity_sign=-1.);before=different.state_digest()
    with pytest.raises(ValueError,match='configuration/manifest'):different.load(tmp_path/'checkpoint')
    assert different.state_digest()==before


def test_lost_commit_acknowledgement_reconciles(factory,monkeypatch):
    c=factory();trial=c.prepare(1e-10,transaction_id='lost-ack');original=c.native._call
    def lose_ack(command):
        value=original(command)
        if command.startswith('COMMIT '):raise ValueError('simulated lost commit acknowledgement')
        return value
    monkeypatch.setattr(c.native,'_call',lose_ack)
    c.commit(trial)
    assert c.native.state()['revision']==1
    np.testing.assert_array_equal(c.u,trial.end_displacement_m)
    assert c.pending is None


def test_unknown_commit_outcome_stops_joint_owner(factory,monkeypatch,tmp_path):
    from mhdlab.solid_foil_native_heat import NativeCommitUncertain
    c=factory();trial=c.prepare(1e-10,transaction_id='unknown');original=c.native._call;failed=False
    def lose_participant(command):
        nonlocal failed
        if command.startswith('COMMIT '):failed=True;raise ValueError('simulated unavailable participant')
        if failed and command=='INFO':raise ValueError('simulated unavailable participant')
        return original(command)
    monkeypatch.setattr(c.native,'_call',lose_participant)
    with pytest.raises(NativeCommitUncertain):c.commit(trial)
    with pytest.raises(NativeCommitUncertain):c.prepare(1e-10,transaction_id='unsafe-retry')
    with pytest.raises(NativeCommitUncertain):c.save(tmp_path/'ambiguous')
