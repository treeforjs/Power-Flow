import json
import os
from pathlib import Path

p = Path('/home/treefor/projects/runtimes/pic-gpu-candidate-v2-20261029/evidence/record-tests')
assert os.sched_getaffinity(0) == {6}
assert os.getpriority(os.PRIO_PROCESS, 0) >= 15
cases = {}
for name in ('off', 'fallback', 'qualified', 'late_failure', 'outer_failure', 'unclosed'):
    j = json.loads((p / (name + '.json')).read_text())
    initial = json.loads((p / (name + '.json.initial')).read_text())
    assert initial['run_status'] == 'running' and initial['graph_status'] == 'pending'
    assert initial['run_id'] != 'old-run' and initial['solvers'] == []
    assert initial['run_id'] == j['run_id']
    assert j['outer_process_receipt_required'] is True
    assert [s['solver_id'] for s in j['solvers']] == list(range(len(j['solvers'])))
    if name == 'off':
        assert j['run_status'] == 'completed' and j['graph_status'] == 'off'
        assert j['solvers'][0]['direct_solve_calls'] == 1
        assert not j['all_solvers_graph_qualified']
    elif name == 'fallback':
        assert j['run_status'] == 'completed' and j['graph_status'] == 'fallback'
        assert [s['status'] for s in j['solvers']] == ['fallback', 'qualified']
        assert not j['all_solvers_graph_qualified'] and not j['all_solvers_requested_required']
    elif name == 'qualified':
        assert j['run_status'] == 'completed' and j['graph_status'] == 'qualified'
        assert j['all_solvers_graph_qualified'] and j['all_solvers_requested_required']
        assert all(s['closed'] and s['graph_replays'] == 1 for s in j['solvers'])
    else:
        assert j['run_status'] == 'failed' and j['graph_status'] == 'failed'
        assert not j['all_solvers_graph_qualified']
        if name == 'late_failure':
            assert j['solvers'][0]['status'] == 'failed'
            assert j['reason'] == 'Late launch failure: "quoted"\nline'
    cases[name] = 'passed'
result = {'cases': cases, 'invalid_relative_path_case': 'passed by CPU executable',
          'affinity': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
          'gpu_executed': False, 'scope': 'Record lifecycle and JSON semantics only; no CUDA stubs or GPU calls'}
(p.parent / 'record-test-results.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
