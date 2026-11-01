"""Policy must never certify no-data channels merely because their rate is unknown."""
import json
import subprocess
import sys
from pathlib import Path

r = Path(__file__).resolve().parents[1]
m = json.loads((r / 'data/tables/chemistry-v3.1-manifest.json').read_text())
p = json.loads((r / 'data/tables/chemistry-v3.1-coverage-policy.json').read_text())
assert p['model_mode'] == 'strict-production'
assert p['reachability_schema'] == 'pfc-region-reachable-support-v1'
assert 'region_support[r] > 0' in p['policy']['unbounded_regions']
assert p['admitted_unbounded_regions'] == []
placeholders = set(p['placeholder_processes'])
for region in m['region_info']:
    if region['above_grid'] or region['below_grid'] or placeholders.intersection(region['processes']):
        assert not region['bounded'], region
refused = {v['index'] for v in p['refusing_unbounded_regions']}
assert all(x['index'] in refused for x in m['region_info'] if not x['bounded'])
diagnostic = r / 'evidence/diagnostic-policy.json'
subprocess.run([sys.executable, '-I', str(r / 'tools/make_coverage_policy.py'),
                str(r / 'data/tables/chemistry-v3.1-manifest.json'), str(diagnostic),
                '--diagnostic-admit-unknown'], check=True, stdout=subprocess.DEVNULL)
diag = json.loads(diagnostic.read_text())
assert diag['model_mode'] == 'diagnostic-admit-unknown' and diag['admitted_unbounded_regions']
print('PASS strict policy refuses all unknown/placeholder/outside-grid regions; diagnostic opt-in is labeled')
