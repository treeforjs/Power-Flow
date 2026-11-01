#!/bin/bash
# tools/isa_hook_rebuild.sh -- after adding the PF_BUILD_ISA_FLAGS hook to the three build scripts: rebuild with the hook empty (b14, c4, r8, r9) and
# compare every binary with the builds of record (b13, c3, r6, r7). Runs on one leased Power-Flow core (any-core.sh).
set -u
unset PF_BUILD_ISA_FLAGS
R=/home/treefor/projects/runtimes/hall-corrective-v3.4-20261031
O=/home/treefor/projects/.workspace/outputs/hall-corrective-v3.4-20261031
F=/home/treefor/projects/runtimes/arch-base-20261028/run-frozen
cd $R
echo "$(date +%H:%M:%S) b14"; /usr/bin/python -I run_hall_v34.py build b14 test_hall_v3 test_resident test_v32 test_x87 test_v34
echo "$(date +%H:%M:%S) c4";  /usr/bin/python -I run_coordinator_parity.py build c4
cd $R/coordinator-v7.1-copy
echo "$(date +%H:%M:%S) r8";  $F /usr/bin/python -I run_v7_1.py build r8 --hall v3.4 --native v7.1
echo "$(date +%H:%M:%S) r9";  $F /usr/bin/python -I run_v7_1.py build r9 --hall v3.4 --native v7
echo "$(date +%H:%M:%S) compare"
/usr/bin/python -I - <<EOF
import json
O='$O'
def j(p): return json.load(open(p))
res = {'schema': 'hall-v3.4-isa-hook-rebuild-v1', 'PF_BUILD_ISA_FLAGS': 'unset'}
a, b = j(O + '/b13/build/build-receipt.json'), j(O + '/b14/build/build-receipt.json')
res['b14_vs_b13'] = {'status': b['status'], 'isa_flags': b.get('isa_flags'), 'binaries_identical': a['binaries_sha256'] == b['binaries_sha256'],
                     'binaries_sha256': b['binaries_sha256'],
                     'commands_identical_up_to_label': json.dumps([c['argv'] for c in a['commands']]).replace('/b13', '/b14')
                                                      == json.dumps([c['argv'] for c in b['commands']])}
a, b = j(O + '/coordinator-parity/build-c3/build-receipt.json'), j(O + '/coordinator-parity/build-c4/build-receipt.json')
res['c4_vs_c3'] = {'status': b['status'], 'isa_flags': b.get('isa_flags'), 'binary_identical': a['binary_sha256'] == b['binary_sha256'], 'binary_sha256': b['binary_sha256']}
for new, old in (('r8', 'r6'), ('r9', 'r7')):
    a, b = j(O + '/coordinator-v7.1-copy/build-%s/build-receipt.json' % old), j(O + '/coordinator-v7.1-copy/build-%s/build-receipt.json' % new)
    res['%s_vs_%s' % (new, old)] = {'status': b['status'], 'bindings': b['bindings'], 'isa_flags': b.get('isa_flags'),
                                    'binary_identical': a['binary_sha256'] == b['binary_sha256'], 'binary_sha256': b['binary_sha256']}
open(O + '/isa-hook-rebuild.json', 'w').write(json.dumps(res, indent=1) + '\n')
print(json.dumps(res, indent=1))
EOF
echo "$(date +%H:%M:%S) DONE"
