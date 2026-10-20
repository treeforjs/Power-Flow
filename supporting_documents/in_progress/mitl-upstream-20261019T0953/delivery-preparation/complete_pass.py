"""Root-only lease completion after compressed-payload verification."""
import datetime
import fcntl
import hashlib
import json
import os
from packet_common import WORK, REPO, OUT, BASE, contract, verify_archive, records, clean_head, load, sha, encoded, exists, write_new


def atomic(path, data):
    temporary = path.with_name(path.name+'.new')
    assert not exists(temporary)
    with temporary.open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def main():
    value = contract()
    status_path = WORK.parent/'mitl-upstream-status.json'
    with (WORK.parent/'mitl-upstream-status.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        before = status_path.read_bytes()
        status = load(status_path)
        lease = load(WORK/'lease.json')
        now = datetime.datetime.now(datetime.timezone.utc)
        assert status['active_lease_owner'] == lease['pass_id'] == WORK.name
        assert datetime.datetime.fromisoformat(status['active_editing_lease_until_utc']) > now
        assert status['last_complete_commit'] == lease['baseline_commit'] == BASE
        assert status['handoff'] == lease['prior_handoff'] and status['handoff_sha256'] == lease['prior_handoff_sha256']
        assert not exists(WORK/'completion.json') and not exists(WORK/'status-before-completion.json')
        # Every scientific milestone, receiver dependency and queue item is preserved verbatim.
        queue = json.loads(json.dumps(status['next_work']))
        changed = {'last_complete_commit','last_completed_pass_id','last_update_utc','handoff','handoff_sha256',
            'handoff_bytes','handoff_document','active_lease_owner','active_editing_lease_until_utc','reserved_paths',
            'active_completed_subtasks','active_phase','prior_handoffs','verification_history','completed_milestones',
            'npl_processability_1996_source_assessment'}
        preserved = {key: data for key, data in status.items() if key not in changed}
        clean_head(value)
        saved = load(WORK/'delivery_verification.json')
        assert saved['status'] == 'passed'
        actual = verify_archive(value)
        assert actual == saved
        records(value)
        now = datetime.datetime.now(datetime.timezone.utc)
        assert datetime.datetime.fromisoformat(status['active_editing_lease_until_utc']) > now
        prior = status.setdefault('prior_handoffs', [])
        if status['handoff'] not in prior:
            prior.append(status['handoff'])
        history = status.setdefault('verification_history', [])
        history.append({'pass_id': WORK.name, 'commit': value['approved_head'],
            'handoff': actual['archive'], 'handoff_sha256': actual['sha256'],
            'scope': 'New NPL report source identity/access/alloy assessment; exact-file packet verification only.'})
        status.update(last_complete_commit=value['approved_head'], last_completed_pass_id=WORK.name,
            last_update_utc=now.isoformat(), handoff=actual['archive'], handoff_sha256=actual['sha256'],
            handoff_bytes=actual['bytes'], handoff_document=str(OUT/'systems-mitl-npl-processability-source-10192026.md'),
            active_lease_owner=None, active_editing_lease_until_utc=None, reserved_paths=[], active_completed_subtasks=[],
            active_phase='Completed bounded NPL source assessment; existing receiver and physical-input dependencies remain pending')
        status['npl_processability_1996_source_assessment'] = {'commit': value['approved_head'],
            'public_source_sha256': value['public_sha256'], 'acquisition_evidence_sha256': value['evidence_sha256'],
            'source_assessment_summary': value['handoff_summary'], 'physical_parameters_adopted': 0,
            'raw_PDFs_distributed': False, 'actual_foil_qualified': False, 'source_or_model_replay': False}
        status.setdefault('completed_milestones', []).append(value['completion_milestone'])
        assert status['next_work'] == queue
        assert all(status[key] == data for key, data in preserved.items())
        after = encoded(status)
        result = {'schema': 'mitl-npl-processability-pass-completion-v1', 'status': 'completed',
            'utc': now.isoformat(), 'pass_id': WORK.name, 'source_commit': value['approved_head'],
            'archive': actual['archive'], 'archive_sha256': actual['sha256'], 'actual_archive_reverified': True,
            'status_before_sha256': hashlib.sha256(before).hexdigest(), 'status_after_sha256': hashlib.sha256(after).hexdigest(),
            'queue_and_prior_scientific_evidence_preserved': True, 'lease_released': True,
            'source_or_scientific_checks_repeated': False, 'physical_parameters_adopted': 0}
        with (WORK/'status-before-completion.json').open('xb') as stream:
            stream.write(before)
        atomic(status_path, after)
        write_new(WORK/'completion.json', result)
        with (WORK.parent/'mitl-upstream-journal.jsonl').open('a') as stream:
            stream.write(json.dumps(result, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
