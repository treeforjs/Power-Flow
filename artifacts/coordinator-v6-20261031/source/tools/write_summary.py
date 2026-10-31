#!/usr/bin/env python3
"""Compose SUMMARY.txt of the first coupled pilot from the run, qualification, closure, render and vacuum records.

  write_summary.py PILOT_DIR QUAL_DIR VACUUM_CPU_DIR VACUUM_GPU_DIR
"""
import json, sys
from pathlib import Path

T0 = 7.4978805549199997e-10


def j(p):
    return json.loads(Path(p).read_text())


def main():
    P, Q, VC, VG = map(Path, sys.argv[1:5])
    run = j(P / 'run' / 'run-summary.json')
    cl = j(P / 'ledger-closure.json')
    rs = j(P / 'figures' / 'render-summary.json')
    qs = j(Q / 'QUALIFICATION-SUMMARY.json')
    vc, vg = j(VC / 'vacuum-summary.json'), j(VG / 'vacuum-summary.json')
    deb = j(P / 'figures' / 'debye-vs-cell.json')
    L = []
    w = L.append
    w('First coupled pilot of the water-breakdown model (coordinator v6) -- %s' % 'pf-pilot-1-20261031')
    w('')
    w('Code: runtimes/pic-water-common-commit-v6-20261031 (README.txt), build %s, driver binary sha256 %s.' % (qs['build']['label'], qs['build']['binary_sha256']))
    w('Run: SimQueue job, partition powerflow, cores 6-7, GPU RTX 3080 Ti (GPU field backend gpu1 + GPU push), one process.')
    w('')
    w('1. Setup')
    w('  Geometry: two SS304 foils 5 mm wide, 0.2 mm thick, 5 mm gap, 25 mm long (z), Hall mesh 45 x 45 transverse cells of')
    w('  200 um, 4 axial foil cells of 6.25 mm + 1 air cell each end; C19 current drive (shot 17133 trace) from its first knot.')
    w('  Genesis at %.9g s (first C19 knot), empty plasma (no authored carriers), foils at 300 K.' % T0)
    w('  Sources: 8 foil-slice source owners, 32 oxide-water patches (all four faces of every slice; Temkin/Eyring, 0.87 eV,')
    w('  initial coverage 0.33, %.4g molecules in total); thermal-field emission with work function 4.3 eV.' % run['initial']['source_stock_molecules'])
    w('  Arrival: 32 electron + 32 neutral streams, w0 = 1e6 for both (surface temperature stayed near 300 K, no larger w0 needed).')
    w('  Chemistry: parent coverage v2 (Song 2021 <= 1 keV, outer3 BEB 1-10 keV, provisional RBEB > 10 keV), default refusal limit.')
    w('  dt = 10 fs (unchanged guard). Ions are moving carriers of the Hall push; screening could only arise from the Hall EM solve.')
    w('')
    w('2. Range reached')
    w('  steps %d, t = %.9g s -> %.9g s (simulated %.4f ns), stop: %s, wall %.2f h.' % (run['steps_completed'], T0, run['final_time_s'],
      run['simulated_s'] * 1e9, run['stop'], run['wall_s'] / 3600))
    if run['stop'] == 'stop_file':
        w('  Stopped by decision after 8099 steps (STOP file, 06:10 UTC) to free cores 6-7 and the GPU for the step-cost and')
        w('  adaptive-step work, which is expected to make this window take minutes.')
        w('  The run was configured for +2 ns. Pilot attempt 1 (build q2) refused at step 3849 (emission kernel field domain, README')
        w('  decision 5) and is kept as attempt-1-refused-step-3849/.')
    w('')
    w('3. Ledger closure (worst values over all steps; ledger.bin + ledger-schema.json, ledger-sample.csv)')
    for k in ('max_abs_total_charge_change_since_genesis_C', 'max_abs_gauss_defect_C', 'max_gauss_over_bound', 'max_divB_T_m', 'max_divB_over_bound',
              'max_gauss_increment_C', 'max_abs_step_energy_residual_J', 'max_roundoff_bound_J', 'max_abs_reaction_defect_J', 'max_solve_residual_over_target',
              'max_allocations_in_swap_window', 'h2o_inventory_exact_all_steps', 'H_O_nuclei_exact_all_steps'):
        w('  %-48s %s' % (k, cl[k]))
    w('  cumulative Joule heat %.6g J, thermal change %.6g J, circuit work %.6g J, field energy change %.6g J' % (
        cl['cumulative_joule_J'], cl['cumulative_thermal_delta_J'], cl['cumulative_circuit_work_J'], cl['cumulative_field_delta_J']))
    w('')
    w('4. Physics observed')
    w('  births: electrons %d, H2O packets %d; sinks %d; parent ionization events %d.' % (cl['births_electron'], cl['births_h2o'], cl['sinks'], cl['parent_events']))
    w('  expected emitted electrons (sum over the run) %.3g; largest per-patch-step expectation 10^%.1f.' % (cl['cum_expected_electrons'],
      cl['max_log10_expected_electrons_patch_step']))
    w('  expected desorbed H2O molecules %.6g, released (n w0) %.6g; source stock %.6g.' % (cl['cum_expected_h2o_molecules'], cl['cum_released_h2o_molecules'],
      cl['final_source_stock_molecules']))
    w('  patch surface temperature %.4f-%.4f K; extraction field at the patch centres up to %.4g V/m.' % (cl['min_patch_T_K'], cl['max_patch_T_K'],
      cl['max_extraction_field_V_m']))
    w('  drive current %.1f -> %.1f A; terminal voltage %.4g -> %.4g V.' % tuple(cl['current_A_first_last'] + cl['terminal_V_first_last']))
    w('  Debye length per output time: %d frames; %s' % (len(deb), 'no electrons at any output time (lambda_D undefined)' if all(d['debye_length_m'] == 'inf' for d in deb)
      else 'see debye-vs-cell.csv'))
    w('')
    w('5. Timing')
    w('  step total: median %.1f ms, mean %.1f ms, p99 %.1f ms (GPU field + GPU push).' % (cl['step_total_ms']['median'], cl['step_total_ms']['mean'],
      cl['step_total_ms']['p99']))
    for k, v in sorted(rs['phase_median_ms'].items(), key=lambda kv: -kv[1])[:12]:
        w('    %-30s %8.2f ms' % (k, v))
    pc, pg = qs['runs']['q3-pilot-cpu-20']['step_ms_median'], qs['runs']['q3-pilot-gpu-20']['step_ms_median']
    w('  CPU vs GPU backends, pilot fixture, 20 steps: %.1f ms vs %.1f ms per step (GPU run includes a CPU push parity check every step).' % (pc, pg))
    w('')
    w('6. Hall-only vacuum pre-run (no particles), Owner::prepare_vacuum / preflight_vacuum_final / publish')
    w('  CPU cpu-pcg: median %.2f ms/step; GPU gpu1: median %.2f ms/step (%d steps each).' % (vc['ms_per_step']['median'], vg['ms_per_step']['median'], vc['steps']))
    w('  dt = 20 fs refused by the field model guard: "%s" (Model::prepare requires dt <= 1e-14; not loosened).' % vc['dt_2e-14_refusal'])
    best = min(vc['ms_per_step']['median'], vg['ms_per_step']['median'])
    for ton in (10, 20, 25):
        n = (ton * 1e-9 - T0) / 1e-14
        w('  to T_on = %2d ns: %.4g steps of 10 fs -> %.1f h at %.2f ms/step' % (ton, n, n * best / 3.6e6, best))
    w('')
    for line in (P / 'LIMITATIONS.txt').read_text().splitlines():
        w(line)
    (P / 'SUMMARY.txt').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
