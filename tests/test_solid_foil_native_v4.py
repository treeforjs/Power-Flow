"""Common transaction and ownership tests against the actual sealed-source worker."""
import importlib.util,json,os
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('v4_demo',ROOT/'tools/run_native_thermoelastic_v4.py')
demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)

@pytest.fixture
def factory():
    exe=os.environ.get('SS304_NATIVE_V4_EXECUTABLE')
    if not exe or not Path(exe).is_file():pytest.skip('requires separately built sealed v4 native executable')
    cases=[]
    def make(**kw):
        c=demo.make_case(Path(exe),**kw);cases.append(c);return c
    yield make
    for c in cases:c.native.close()


def prepare(c,dt=1e-15,begin=False,name='trial'):
    return c.prepare(dt,transaction_id=name,**demo.inputs(c,dt,begin=begin))


def test_combined_trial_rollback_has_no_source_or_reservoir_mutation(factory):
    c=factory();before=c.state_digest();source=c.native.source_state()
    t=prepare(c,begin=True);r=t.receipt;b=r['native_budget']
    assert c.state_digest()==before and c.native.source_state()==source
    assert b['source_wall_J']!=0 and b['adsorbed_reference_change_J']>0
    assert b['gas_reference_J']==b['joule_J']==b['authored_heat_J']==0
    assert r['source_candidate_held_at_accepted_start']
    assert abs(r['global_combined_energy_residual_J'])<=r['energy_tolerance_J']
    c.discard();assert c.state_digest()==before and c.native.source_state()==source
    with pytest.raises(ValueError,match='stale or foreign'):c.commit(t)
    accepted=prepare(c,begin=True,name='accepted');c.commit(accepted)
    assert c.native.state()['revision']==1 and c.native.state()['layer_count']==1
    assert c.native.state()['gas_packet_count']>0
    with pytest.raises(ValueError,match='stale or foreign'):c.commit(accepted)


def test_active_layer_and_common_source_checkpoint_restart(factory,tmp_path):
    c=factory();t=prepare(c,begin=True);c.commit(t)
    assert c.native.source_state()['layers'][0]['status']==0
    c.save(tmp_path/'joint');d=factory();d.load(tmp_path/'joint')
    assert d.state_digest()==c.state_digest() and d.native.source_state()==c.native.source_state()
    saved=json.loads((tmp_path/'joint/joint.json').read_text())
    assert json.loads(saved['source_provenance_json'])['source_stage']['field_time_s']==0
    a=prepare(c,1e-12,name='continue');b=prepare(d,1e-12,name='continue')
    assert a.receipt==b.receipt
    assert a.receipt['native_budget']['source_wall_J']==0  # flight does not redebit birth
    c.commit(a);d.commit(b);assert c.state_digest()==d.state_digest()
    assert c.native.source_state()['layers'][0]['status']==1
    a=prepare(c,1e-10,name='conduction-mechanics');r=a.receipt
    assert np.max(np.abs(r['native_budget']['nodal_ledger']['conduction_out_J']))>0
    assert r['native_budget']['mechanical_exchange_J']!=0
    c.commit(a)


def test_failure_rejects_all_participants(factory):
    c=factory();before=c.state_digest();source=c.native.source_state()
    with pytest.raises(ValueError,match='failed to converge'):
        c.prepare(1e-15,transaction_id='insufficient',max_iterations=1,**demo.inputs(c,1e-15,begin=True))
    assert c.state_digest()==before and c.native.source_state()==source
    args=demo.inputs(c,1e-15,begin=True);args['source_stage']['field_time_s']=1e-15
    with pytest.raises(ValueError):c.prepare(1e-15,transaction_id='wrong-stage',**args)
    assert c.state_digest()==before
    high=factory(velocity_sign=1000.)
    with pytest.raises(ValueError,match='proof-stress|small-strain'):
        high.prepare(1e-10,transaction_id='stress',**demo.inputs(high,1e-10,begin=True))
    assert high.native.state()['revision']==0 and high.native.state()['gas_packet_count']==0


def test_full_native_digest_recovers_lost_commit_ack(factory,monkeypatch):
    c=factory();t=prepare(c,begin=True);expected=t.receipt['native_source_state1']['state_digest'];original=c.native._call
    def lose(command):
        result=original(command)
        if command.startswith('COMMIT '):raise ValueError('simulated lost acknowledgement')
        return result
    monkeypatch.setattr(c.native,'_call',lose);c.commit(t)
    assert c.native.state()['state_digest']==expected
    assert c.native.state()['gas_packet_count']>0 and c.native.state()['layer_count']==1


def test_checkpoint_tamper_rejected_without_source_replacement(factory,tmp_path):
    c=factory();t=prepare(c,begin=True);c.commit(t);c.save(tmp_path/'joint')
    p=tmp_path/'joint/native.bin';data=p.read_bytes();p.write_bytes(data[:-1]+bytes([data[-1]^1]))
    d=factory();before=d.state_digest();source=d.native.source_state()
    with pytest.raises(ValueError,match='checkpoint bytes'):d.load(tmp_path/'joint')
    assert d.state_digest()==before and d.native.source_state()==source
