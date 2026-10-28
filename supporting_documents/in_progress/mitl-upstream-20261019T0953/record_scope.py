"""Bind lease/bootstrap/source-discovery provenance without additional requests."""
import datetime
import hashlib
import json
from pathlib import Path

work = Path(__file__).resolve().parent
def pin(path):
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

status = json.loads((work.parent / 'mitl-upstream-status.json').read_text())
assert status['active_lease_owner'] == work.name
assert set(status['reserved_paths']) == set(json.loads((work / 'lease.json').read_text())['reserved_paths'])
record = {
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'event': 'new_primary_source_acquisition_scope', 'pass_id': work.name,
    'lease': pin(work / 'lease.json'), 'claim_script': pin(work / 'claim_lease.py'),
    'initial_claim_command': ['python', '/tmp/claim_mitl_npl_source.py'], 'claim_exit_code': 0,
    'workspace_available': True, 'branch_checked': 'research/flash-publication-foundation',
    'latest_handoff_read_before_selection': True, 'prior_handoff_hash_verified': True,
    'no_prior_unexpired_editing_lease': True, 'new_source_report_previously_unassessed': True,
    'initial_search_before_lease': {
        **pin(work / 'source-acquisition/initial-discovery-response.json'),
        'qualification': 'Read-only discovery; no filesystem or external-application mutation.'},
    'second_search_after_lease': pin(work / 'source-acquisition/official-discovery-response.json'),
    'source_request_limits': {'new_discovery_search_calls': 2, 'official_requests': 6},
    'network_readiness': {'observations_current': True, 'spec_revision': '62',
                          'network_policy_mode': 'unrestricted', 'network_policy_state': 'enforced',
                          'vpn_configured': False, 'credential_bindings': 0},
    'bootstrap': ['Copied unchanged 0810/record_command.py and claim script into newly leased scope.',
                  'Created newly reserved source-acquisition and external source-cache directories.'],
    'public_source_adoption': False, 'source_or_test_or_solver_replayed': False,
}
out = work / 'scope.json'
assert not out.exists()
out.write_text(json.dumps(record, indent=2) + '\n')
with (work.parent / 'mitl-upstream-journal.jsonl').open('a') as stream:
    stream.write(json.dumps(record) + '\n')
print(json.dumps({'scope': pin(out)}, indent=2))
