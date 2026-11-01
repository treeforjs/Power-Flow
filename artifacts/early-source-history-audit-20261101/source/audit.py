#!/usr/bin/env python3
"""Read existing ledgers and support metadata; never evaluate a source model or replay."""
import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import zipfile

ROOT = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def rows_from_bytes(data, schema):
    need(schema['byte_order'] == 'little' and schema['value'] == 'float64', 'unsupported ledger encoding')
    cols = schema['columns']
    need(len(cols) == len(set(cols)), 'duplicate ledger columns')
    width = 8 * len(cols)
    need(width > 0 and len(data) % width == 0, 'incomplete ledger row')
    return [dict(zip(cols, x)) for x in struct.iter_unpack('<' + 'd' * len(cols), data)]


def summarize(rows, genesis, summary, requested_end_s=None):
    need(rows, 'empty accepted ledger')
    need(genesis['census_rows'] == 0, 'this extractor requires the recorded empty pilot genesis')
    need(summary['fixture'] == 'pilot', 'not the supported pilot fixture')
    need(summary['genesis'] == genesis, 'summary genesis differs')
    need(summary['source_owner']['surface'] == 'benchmark', 'surface differs from reviewed benchmark')
    need(summary['source_owner']['interval_preset'] == 'benchmark', 'source interval preset differs')
    t = genesis['time_s']
    for i, row in enumerate(rows, 1):
        need(row['step'] == i, 'noncontiguous step sequence')
        need(math.isfinite(row['dt_s']) and row['dt_s'] > 0, 'invalid accepted dt')
        need(t + row['dt_s'] == row['time_s'], 'accepted clock discontinuity')
        t = row['time_s']
    need(summary['steps_completed'] == len(rows), 'summary step count differs')
    need(summary['final_time_s'] == t, 'summary final clock differs')
    if requested_end_s is not None:
        need(requested_end_s == t, 'unsupported requested endpoint: no extrapolation or partial-row integration')
    last = rows[-1]
    pairs = {
        'expected_electrons': ('arr_expected_electrons_step', 'arr_cum_expected_electrons'),
        'expected_h2o_molecules': ('arr_expected_h2o_molecules_step', 'arr_cum_expected_h2o_molecules'),
        'released_h2o_molecules': ('arr_released_h2o_molecules_step', 'arr_cum_released_h2o_molecules'),
        'electron_macroparticle_births': ('births_electron', 'arr_cum_electron_births'),
        'h2o_macroparticle_births': ('births_h2o', 'arr_cum_h2o_births'),
    }
    counts = {}
    for name, (step, cum) in pairs.items():
        values = [r[step] for r in rows]
        need(all(math.isfinite(x) and x >= 0 for x in values), 'invalid recorded count: ' + step)
        total = math.fsum(values)
        # Printed ledger values are rounded binary64 projections of exact owner sums.
        need(math.isclose(total, last[cum], rel_tol=1e-12, abs_tol=0), 'step/cumulative mismatch: ' + name)
        counts[name] = {'recorded_cumulative': last[cum], 'fsum_of_binary64_step_values': total}
    for ledger_name, summary_name in [('births_electron', 'electron_births'), ('births_h2o', 'h2o_births')]:
        need(math.fsum(r[ledger_name] for r in rows) == summary['totals'][summary_name], 'birth total differs')
    settings = dict(item.split('=', 1) for item in genesis['settings'].split(';'))
    inventory_columns = ['h2o_source_stock', 'h2o_native_neutral', 'h2o_native_ion',
                         'h2o_sunk_neutral_cum', 'h2o_sunk_ion_cum', 'h2o_total']
    return {
        'interval_s': [genesis['time_s'], t],
        'duration_s_from_clock': t - genesis['time_s'],
        'sum_recorded_dt_s': math.fsum(r['dt_s'] for r in rows),
        'accepted_rows': len(rows), 'counts': counts,
        'sampled_electrons': last['arr_cum_sampled_electrons'],
        'h2o_expected_in_nominal_macro_units': last['arr_cum_expected_h2o_molecules'] / float.fromhex(settings['w0n']),
        'recorded_inventory_accounting': {
            'initial_source_stock_molecules': summary['initial']['source_stock_molecules'],
            'binary64_ledger_projections': {name: {'first': rows[0][name], 'last': last[name],
                'min': min(r[name] for r in rows), 'max': max(r[name] for r in rows)} for name in inventory_columns},
            'all_exact_water_conservation_flags_true': all(r['h2o_total_exactly_conserved'] == 1 for r in rows),
            'separate_per_patch_remainder_or_credit_exported': False,
            'final_arrival_checkpoint_exported': bool(summary['checkpoints']),
            'rule': 'Arrival v1.1 samples n=SR(D/w0) quanta per step, capped by inventory; only n*w0 is debited. Unsampled material stays in source stock. Cumulative expectation is retained separately in arrival expected_sum, projected to the ledger.',
            'scope': 'No exact Count mantissa/exponent or per-patch credit values are exported by this pilot. These raw scalar values and exact-owner conservation flags are reported as recorded; no missing fractional state is invented.',
            'convergence_gate': 'Macro-weight/source-quantization convergence remains unqualified. No weight change or plasma-onset prediction follows from this extraction.',
        },
        'unresolved_emission_patch_steps': math.fsum(r['arr_unresolved_emission_patches_step'] for r in rows),
        'source_patch_temperature_range_K': [min(r['src_patch_T_min_K'] for r in rows), max(r['src_patch_T_max_K'] for r in rows)],
        'signed_extraction_field_range_V_m': [min(r['src_extraction_field_min_V_m'] for r in rows), max(r['src_extraction_field_max_V_m'] for r in rows)],
        'source_configuration': {
            'settings': genesis['settings'], 'owner': summary['source_owner'],
            'phi_eV': float.fromhex(settings['phi']), 'theta0': float.fromhex(settings['theta0']),
            'electron_nominal_macro_weight': float.fromhex(settings['w0e']),
            'h2o_nominal_macro_weight': float.fromhex(settings['w0n']),
            'source_owners': genesis['source_owners'], 'source_patches': genesis['source_patches'],
            'source_configuration_sha256': genesis['source_configuration_sha256'],
            'arrival_policy_sha256': genesis['arrival_policy_sha256'],
            'source_sha256_at_genesis': genesis['source_sha256'],
            'benchmark_semantics': 'constant clean-surface phi; one flat beta=1 class; wall and photo off; conditional analytic oxide-water desorption',
        },
        'binding_pair': summary['binding_pair'], 'physics_qualification': summary['physics_qualification'],
        'interpretation': 'Recorded model expectations and one sampled realization, not physical confidence bounds. Zero sampled macroparticles does not establish zero source expectation.',
    }


def cache_metadata(path):
    arrays = {}
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            with z.open(name) as f:
                magic = f.read(8)
                need(magic[:6] == b'\x93NUMPY' and magic[6] in (1, 2, 3), 'unsupported cache array')
                n = 2 if magic[6] == 1 else 4
                length = struct.unpack('<H' if n == 2 else '<I', f.read(n))[0]
                header = ast.literal_eval(f.read(length).decode())
                arrays[name] = header
                if name.endswith('_time_s.npy'):
                    need(header['descr'] == '<f8' and len(header['shape']) == 1, 'unsupported time array')
                    times = struct.unpack('<' + 'd' * header['shape'][0], f.read())
                    need(all(b > a for a, b in zip(times, times[1:])), 'cache time not strictly increasing')
                    arrays[name]['first_last_time_s'] = [times[0], times[-1]]
    return arrays


def audit(manifest):
    for path, expected in manifest['files_sha256'].items():
        need(sha(path) == expected, 'input pin differs: ' + path)
    p = manifest['paths']; pilot = Path(p['pilot'])
    receipt = load(pilot / 'run-receipt.json')
    need(receipt['status'] == 'passed' and receipt['command']['returncode'] == 0, 'pilot did not pass')
    for key in ['sealed_unchanged_after_run', 'own_sources_unchanged_after_run', 'binary_unchanged_after_run']:
        need(receipt[key] is True, 'pilot provenance guard failed: ' + key)
    for name, expected in receipt['files_sha256'].items():
        need(sha(pilot / name) == expected, 'recorded pilot output differs: ' + name)
    build = load(p['build_receipt'])
    need(build['status'] == 'passed' and build['inputs_unchanged_after_build'], 'build not qualified')
    need(build['binary_sha256'] == receipt['binary_sha256'], 'build/run binary mismatch')
    for path, expected in {**build['own_sources_sha256'], **build['dependency_sha256']}.items():
        need(sha(path) == expected, 'compiled source no longer matches receipt: ' + path)
    rows = rows_from_bytes((pilot / 'ledger.bin').read_bytes(), load(pilot / 'ledger-schema.json'))
    segment = summarize(rows, load(pilot / 'genesis.json'), load(pilot / 'run-summary.json'))
    with open(p['drive'], newline='') as f:
        drive = list(csv.DictReader(f))
    first = float(drive[0]['time_s'])
    need(segment['interval_s'][0] == first, 'pilot start differs from first supported drive knot')
    hall = Path(p['hall_v34']); prerun = load(hall / 'prerun.json')
    checkpoints = load(hall / 'checkpoints.json')
    series = [json.loads(line) for line in (hall / 'series.jsonl').read_text().splitlines()]
    need(checkpoints == prerun['checkpoints'], 'vacuum checkpoint manifests differ')
    deltas = [b['time_s'] - a['time_s'] for a, b in zip(series, series[1:])]
    need(all(x > 0 for x in deltas), 'vacuum series not strictly increasing')
    return {
        'schema': 'early-source-history-audit-v1',
        'method': 'Read-only extraction of existing qualified pilot ledger and inventory of other histories; no source evaluation, replay, interpolation, or new simulation.',
        'input_files_verified': len(manifest['files_sha256']),
        'pilot_output_receipt_hashes_verified': len(receipt['files_sha256']),
        'compiled_source_hashes_verified': len({**build['own_sources_sha256'], **build['dependency_sha256']}),
        'pilot_binary_sha256': receipt['binary_sha256'],
        'first_supported_drive_time_s': first,
        'recorded_coupled_segment': segment,
        'hall_v34_vacuum': {
            'interval_s': [prerun['start_s'], prerun['end_time_s']], 'steps': prerun['steps'],
            'aggregate_rows': len(series), 'aggregate_spacing_s_min_max': [min(deltas), max(deltas)],
            'checkpoint_times_s': [x['time_s'] for x in checkpoints],
            'first_aggregate': series[0], 'last_aggregate': series[-1],
            'per_patch_source_history_available': False,
            'cumulative_source_bound': None,
            'reason': 'Aggregate mean/max foil temperature and four absolute face-field maxima plus sparse 1 ns field/heat snapshots do not supply signed per-patch step histories, source stock/RNG, or plasma history. No interpolation bound is established.',
            'seal_scope': 'Current checkpoint files are SHA256-pinned; native owner semantic checkpoint seals were not revalidated by this metadata extractor.',
        },
        'original_temperature_cache': {
            'arrays': cache_metadata(p['original_cache']),
            'lineage': 'Earlier 2D high-resolution heating/supply cache. Distinct from corrected finite-target Hall histories; historical anode reconstruction and coarse magnetic-field limitations remain. Not a source-v7 coupled history.',
            'used_to_fill_gaps': False,
        },
        'ss304_phase_owner': {
            'path': p['ss304'],
            'status': 'Standalone v1.3 qualification is pinned for context; the audited pilot and Hall vacuum histories use solid-only Hall heat and contain no coupled SS304 phase history.',
        },
        'omitted_early_source_through_20ns': {
            'status': 'refused_unsupported_history', 'bound': None,
            'unrepresented_interval_s_after_recorded_pilot': [segment['interval_s'][1], 20e-9],
            'reason': 'No complete supported source-v7 accepted per-patch history after this 20 ps segment. A 20 ns vacuum restart resets source stock and is not continuous plasma evolution. The pilot benchmark differs from enhanced/coverage-dependent/wall source options.',
        },
        'before_first_drive_knot': 'Unsupported; no extrapolation of the drive or heating/source state.',
        'excluded_active_versions': ['hall-corrective-v3.6-20261101'],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, default=ROOT / 'INPUTS-SHA256.json')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    need(not args.output.exists(), 'refuse overwriting audit output')
    report = audit(load(args.manifest))
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
