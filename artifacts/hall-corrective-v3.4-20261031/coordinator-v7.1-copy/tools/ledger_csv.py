#!/usr/bin/env python3
"""Decimated CSV view of a coordinator v6 binary per-step ledger, plus the worst-case closure values over all steps.

  ledger_csv.py RUN_DIR OUT_CSV EVERY OUT_JSON
Every step with a birth, sink or parent event is kept in addition to every EVERY-th step and the last step.
"""
import json, sys
from pathlib import Path
import numpy as np


def main():
    run, out_csv, every, out_json = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4])
    cols = json.loads((run / 'ledger-schema.json').read_text())['columns']
    raw = np.frombuffer((run / 'ledger.bin').read_bytes(), dtype='<f8')
    rows = raw[: (len(raw) // len(cols)) * len(cols)].reshape(-1, len(cols))
    c = {n: i for i, n in enumerate(cols)}
    keep = (np.arange(len(rows)) % every == every - 1) | (rows[:, c['births_electron']] > 0) | (rows[:, c['births_h2o']] > 0) | \
           (rows[:, c['sinks']] > 0) | (rows[:, c['parent_events']] > 0)
    keep[-1] = True
    with out_csv.open('w') as h:
        h.write(','.join(cols) + '\n')
        for r in rows[keep]:
            h.write(','.join('%.17g' % x for x in r) + '\n')
    col = lambda n: rows[:, c[n]]
    worst = {
        'steps': int(len(rows)), 'first_time_s': float(col('time_s')[0]), 'last_time_s': float(col('time_s')[-1]),
        'max_abs_total_charge_change_since_genesis_C': float(np.max(np.abs(col('q_total_change_since_genesis_C')))),
        'max_gauss_over_bound': float(np.max(col('gauss_max_abs_C') / col('gauss_bound_C'))),
        'max_abs_gauss_defect_C': float(np.max(col('gauss_max_abs_C'))),
        'max_divB_over_bound': float(np.max(col('divB_max_T_m') / col('divB_bound_T_m'))),
        'max_divB_T_m': float(np.max(col('divB_max_T_m'))),
        'max_abs_step_energy_residual_J': float(np.max(np.abs(col('E_step_residual_J')))),
        'max_roundoff_bound_J': float(np.max(col('E_roundoff_bound_J'))),
        'max_abs_combined_balance_J': float(np.max(np.abs(col('E_combined_balance_J')))),
        'max_abs_reaction_defect_J': float(np.max(np.abs(col('E_reaction_defect_J')))),
        'h2o_inventory_exact_all_steps': bool(np.all(col('h2o_total_exactly_conserved') == 1)),
        'H_O_nuclei_exact_all_steps': bool(np.all(col('nuclei_exactly_conserved') == 1)),
        'max_solve_residual_over_target': float(np.max(col('solve_residual_over_target'))),
        'max_motion_continuity_C': float(np.max(col('motion_continuity_C'))),
        'max_gauss_increment_C': float(np.max(col('gauss_increment_C'))),
        'max_allocations_in_swap_window': float(np.max(col('alloc_swap_window'))),
        'cumulative_joule_J': float(np.sum(col('E_joule_J'))), 'cumulative_thermal_delta_J': float(np.sum(col('E_thermal_delta_J'))),
        'cumulative_external_loads_J': float(np.sum(col('E_external_loads_J'))), 'cumulative_wall_heat_J': float(np.sum(col('E_wall_heat_J'))),
        'cumulative_circuit_work_J': float(np.sum(col('E_circuit_work_J'))), 'cumulative_field_delta_J': float(np.sum(col('E_field_delta_J'))),
        'cumulative_plasma_work_J': float(np.sum(col('E_plasma_work_J'))),
        'births_electron': int(np.sum(col('births_electron'))), 'births_h2o': int(np.sum(col('births_h2o'))), 'sinks': int(np.sum(col('sinks'))),
        'parent_events': int(np.sum(col('parent_events'))),
        'cum_expected_electrons': float(col('arr_cum_expected_electrons')[-1]),
        'max_log10_expected_electrons_patch_step': float(np.max(col('arr_max_log10_expected_electrons_patch'))),
        'cum_expected_h2o_molecules': float(col('arr_cum_expected_h2o_molecules')[-1]),
        'cum_released_h2o_molecules': float(col('arr_cum_released_h2o_molecules')[-1]),
        'final_source_stock_molecules': float(col('h2o_source_stock')[-1]),
        'max_patch_T_K': float(np.max(col('src_patch_T_max_K'))), 'min_patch_T_K': float(np.min(col('src_patch_T_min_K'))),
        'max_extraction_field_V_m': float(np.max(col('src_extraction_field_max_V_m'))),
        'current_A_first_last': [float(col('current_A')[0]), float(col('current_A')[-1])],
        'terminal_V_first_last': [float(col('terminal_V')[0]), float(col('terminal_V')[-1])],
        'max_carriers': {k: float(np.max(col(k))) for k in ('carriers_electron', 'carriers_ion', 'carriers_h2o')},
        'unresolved_coverage_upper_events_cum': float(np.nanmax(col('cov_unresolved_upper_events_cum'))) if np.any(~np.isnan(col('cov_unresolved_upper_events_cum'))) else 0.0,
        'step_total_ms': {'median': float(np.median(col('t_step_total_ms'))), 'mean': float(np.mean(col('t_step_total_ms'))),
                          'p99': float(np.percentile(col('t_step_total_ms'), 99)), 'max': float(np.max(col('t_step_total_ms')))},
    }
    out_json.write_text(json.dumps(worst, indent=1) + '\n')
    print(json.dumps(worst, indent=1))


if __name__ == '__main__':
    main()
