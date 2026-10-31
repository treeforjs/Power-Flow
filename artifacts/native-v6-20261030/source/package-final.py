"""Packages native v6 after verify-final.py passed: fills COORDINATOR-V6-LINK.txt and
RESULTS-006.txt from the verification result, writes the immutable snapshot, the source
archive, the handoff and READY. Usage: python3 package-final.py [LABEL]"""
import shutil, sys, tarfile
from tooling import *

label = sys.argv[1] if len(sys.argv) > 1 else 'final'
F = OUT / label
v = json.loads((F / 'verification-result.json').read_text())
manifest = json.loads((P / 'SOURCE-SHA256.json').read_text())
identity = manifest['source_identity']
need(v['status'] == 'passed' and v['source_identity'] == identity, 'verification identity')
for name, digest in manifest['files'].items():
    need(sha256(P / name) == digest, 'source changed after verification: ' + name)
B = P / 'build'
lib, adapter, v3 = B / 'libpic_native_v6.a', B / 'libpic_hall_motion_adapter_v6.a', B / 'libpic_hall_motion_v3_adapter_v6.a'
hashes = json.loads((F / 'build-artifact-hashes.json').read_text())
for p in [lib, adapter, v3]:
    need(hashes[str(p.relative_to(P))] == sha256(p), 'archive changed after verification: ' + p.name)
genesis = json.loads((F / 'ported/genesis.json').read_text())['payload']
b = v['admission_benchmark']
ms = lambda x: '%9.1f' % (1e3 * x)
rows = ['      carriers  v5: Hall v2 JSON seal + JSON record build + native JSON admission + release = total (ms)',
        '                v6: v3 adapter bind + native v3 admission (census cached / computed in admission) (ms)']
for n in [10000, 100000, 1000000]:
    a5, a6 = b['v5-%d' % n]['median_s'], b['v6-%d' % n]['median_s']
    t5 = a5['hall_v2_json_seal_s'] + a5['v2_json_record_build_s'] + a5['native_v2_json_admission_s'] + a5['v2_json_release_s']
    rows.append('      %8d  v5: %s + %s + %s + %s = %s' % (n, ms(a5['hall_v2_json_seal_s']), ms(a5['v2_json_record_build_s']),
                                                         ms(a5['native_v2_json_admission_s']), ms(a5['v2_json_release_s']), ms(t5)))
    rows.append('                v6: %s + %s / %s   (v6 JSON route admission %s; census digest once per trial %s)' % (
        ms(a6['v3_adapter_bind_warm_s']), ms(a6['native_v3_admission_cached_census_s']),
        ms(a6['native_v3_admission_cold_census_s']), ms(a6['native_v2_json_admission_s']),
        ms(a6['native_census_digest_once_per_trial_s'])))
text = (P / 'COORDINATOR-V6-LINK.template.txt').read_text()
for key, value in {'@V6_IDENTITY@': identity, '@V6_LIB_SHA@': sha256(lib), '@V6_ADAPTER_SHA@': sha256(adapter),
                   '@V6_V3_ADAPTER_SHA@': sha256(v3), '@V6_GENESIS_PREFIX@': genesis['boundary']['accepted_arena_sha256'],
                   '@BENCH_TABLE@': '\n'.join(rows)}.items():
    text = text.replace(key, value)
need('@' not in text.replace('`', ''), 'unfilled template field')
(P / 'COORDINATOR-V6-LINK.txt').write_text(text)
fl = v['v6_flows']['measurements']
lr = v['long_run_200000']
sa = v['v6_source_arrival']
kp = fl['kinetic_port']
results = f"""Native v6 (candidate006) results -- verification run '{label}'
Receipts: {F}
 (receipt-commands.json: every build/run command with lease wrapper, return code, log SHA-256;
  verification-result.json; build-artifact-hashes.json; per-suite outputs and logs)
Source identity {identity}; base native v5 {v['base_v5_identity']} (sealed, unchanged: hashed before/after).

v5 parity (no births/sinks)
 legacy: {v['legacy_checks']} legacy assertions, result/proof JSON identical to
   native v5 final evidence and to the ASan/UBSan run.
 ported candidate003 controls: {v['ported_candidate003_controls']['check_count']} checks, names identical to v5,
   allocation-failure sweep {v['ported_candidate003_controls']['allocation_failure_sweep']} (identical to v5);
   genesis checkpoint differs from v5 only in {', '.join(v['genesis_checkpoint_differences_vs_v5'])}.
 v5 own controls (test_v5_history, only the anchor-chain tag literal adapted to v6): {v['v5_history_controls']['normal']} checks
   normal, {v['v5_history_controls']['sanitized_reduced']} ASan/UBSan reduced, names identical to v5.
 identity parity trace (60 steps, forced events, K = 4): differing fields {v['identity_parity_trace']['differing_fields']}.

v6 controls (tests/test_v6_flows.cpp): {v['v6_flows']['normal']} checks normal, {v['v6_flows']['sanitized_reduced']} ASan/UBSan reduced
 (normal-only: {', '.join(v['v6_flows']['normal_only'])}).
 kinetic first-hit port vs sealed Hall helper: {kp['segments_compared']} segments, {kp['mismatches']} mismatches
   (none {kp['none']}, lower {kp['lower_foil']}, upper {kp['upper_foil']}, outer x/y {kp['outer_xy']},
   short bridge {kp['short_bridge']}, axial {kp['axial']}).
 long run (default suite, 10000 steps, K = 1024): births {fl['long_run']['births']}, sunk {fl['long_run']['sunk_carriers']},
   steps with births and sinks {fl['long_run']['steps_with_births_and_sinks']}, anchors {fl['long_run']['anchors']},
   window rows {fl['long_run']['window_rows']}, VmRSS {fl['long_run']['VmRSS_KiB_after_1000_steps']} -> {fl['long_run']['VmRSS_KiB_end']} KiB.
 long run 200000 steps: births {lr['births']}, sunk {lr['sunk_carriers']}, anchors {lr['anchors']} ({lr['anchor_bytes']} bytes),
   window rows {lr['window_rows']}, ever-packet IDs {lr['ever_packet_ids']}, maximum packets {lr['maximum_packets']},
   VmRSS after 1000 steps {lr['VmRSS_KiB_after_1000_steps']} KiB, end {lr['VmRSS_KiB_end']} KiB, {lr['seconds']:.1f} s.

Source v5 + arrival v1 sweep (tests/test_v6_source_arrival.cpp): {sa['check_count']} checks; {sa['cases']} cases,
 {sa['admitted_steps']} admitted steps with {sa['exact_stream_equalities']} exact per-stream equalities (native credit = source
 debit), {sa['births']} births ({sa['surface_births']} on the emitting face, {sa['interior_births']} interior,
 {sa['active_layer_births']} from active electron layers), {sa['multi_step_committed_cases']} multi-step cases with source and
 arrival commits; refused: source {sa['source_refused']} (typed source v5 G1 weighted-representation refusals of tiny leaves),
 arrival {sa['arrival_refused']}, exposure {sa['exposure_refused']}. ASan/UBSan reduced:
 {v['v6_source_arrival_sanitized_reduced']['check_count']} checks, {v['v6_source_arrival_sanitized_reduced']['exact_stream_equalities']} equalities.

Admission benchmark (median, see COORDINATOR-V6-LINK.txt section B.3 for the table):
""" + '\n'.join(rows) + '\n'
(P / 'RESULTS-006.txt').write_text(results)
S = P / 'snapshots' / ('candidate006-' + identity[:8])
S.mkdir(parents=True, exist_ok=False)
for name in ['source', 'inputs', 'tests', 'scripts']:
    shutil.copytree(P / name, S / name)
for name in ['SOURCE-SHA256.json', 'API-006.txt', 'COORDINATOR-V6-LINK.txt', 'COORDINATOR-V6-LINK.template.txt', 'RESULTS-006.txt',
             'tooling.py', 'build-components.py', 'verify-final.py', 'package-final.py']:
    shutil.copy2(P / name, S / name)
A = OUT / 'native-full-time-candidate006-source.tar.gz'
with tarfile.open(A, 'w:gz') as archive:
    for path in sorted(S.rglob('*')):
        archive.add(path, arcname=str(path.relative_to(S)), recursive=False)
files = {str(p.relative_to(S)): sha256(p) for p in sorted(S.rglob('*')) if p.is_file()}
(OUT / 'source-archive-manifest.json').write_text(json.dumps(files, indent=2) + '\n')
for name in ['SOURCE-SHA256.json', 'API-006.txt', 'COORDINATOR-V6-LINK.txt', 'RESULTS-006.txt']:
    shutil.copy2(P / name, OUT / name)
artifacts = {str(p): sha256(p) for p in [lib, adapter, v3] + [B / x for x in ['test_legacy', 'test_ported', 'test_history',
                                                                              'test_flows', 'test_source_arrival', 'test_bench']]}
artifacts[str(EXACT / 'build/native.o')] = sha256(EXACT / 'build/native.o')
handoff = {'status': 'ready_for_read_only_component_review', 'source_identity': identity, 'base_v5_identity': V5_IDENTITY,
           'candidate_prefix': str(P), 'immutable_snapshot': str(S), 'source_archive': str(A), 'archive_sha256': sha256(A),
           'libraries': {'core': str(lib), 'hall_adapter_v2_json': str(adapter), 'hall_adapter_v3_binary': str(v3)},
           'artifact_sha256': artifacts, 'verification': str(F / 'verification-result.json'),
           'remaining_unqualified': [
               'Hall-side flows motion provider (births with partial flight, first-hit truncation, truncated-path deposits, '
               'conductor/exterior banks, wall heat): births/sinks are admitted from a zero-field reference provider only',
               'Active source layers are credited at the layer state at tau and moved by ordinary motion (no source-owned '
               'exclusion until handoff)',
               'Carrier growth from the per-patch neutral floor rule (at least one leaf per active patch-step) is unbounded '
               'without an H7 merge/weighting policy or a quantized release',
               'Nonempty native disk restore remains deferred (typed refusal after validating anchors and window)',
               'No joint coupled step with a field solve was run'],
           'utc': utc()}
(OUT / 'handoff.json').write_text(json.dumps(handoff, indent=2) + '\n')
(OUT / 'READY').write_text('Native v6 ready for read-only component review. Source ' + identity + '.\n')
outputs = {str(p.relative_to(OUT)): sha256(p) for p in sorted(OUT.rglob('*'))
           if p.is_file() and p.name != 'output-manifest.json' and 'dev' not in p.relative_to(OUT).parts
           and 'superseded-final-nice-guard' not in p.relative_to(OUT).parts}
(OUT / 'output-manifest.json').write_text(json.dumps(outputs, indent=2) + '\n')
print(json.dumps({'source_identity': identity, 'archive_sha256': handoff['archive_sha256'], 'snapshot': str(S)}, indent=1))
