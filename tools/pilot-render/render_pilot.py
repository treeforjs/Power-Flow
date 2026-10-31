#!/usr/bin/env python3
"""Render a coordinator v6 run: midplane density/potential/field frames (PNG, white background) and a GIF,
per-step ledger time series, per-phase timing, the Debye-length table per output time and a JSON summary.

  render_pilot.py RUN_DIR OUT_DIR [--gif-max N]
"""
import json, struct, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
from PIL import Image

INK, INK2, GRID = '#0b0b0b', '#52514e', '#d9d8d4'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
plt.rcParams.update({'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'text.color': INK,
                     'axes.labelcolor': INK, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.grid': False,
                     'font.size': 9, 'lines.linewidth': 1.6})


def ledger(run):
    cols = json.loads((run / 'ledger-schema.json').read_text())['columns']
    raw = np.frombuffer((run / 'ledger.bin').read_bytes(), dtype='<f8')
    rows = raw[: (len(raw) // len(cols)) * len(cols)].reshape(-1, len(cols))
    return {c: rows[:, i] for i, c in enumerate(cols)}, cols


def frames(run):
    lay = json.loads((run / 'frames-layout.json').read_text())
    nx, ny = lay['nx'], lay['ny']
    rec = 2 + 6 * nx * ny + (nx + 1) * (ny + 1)
    raw = np.frombuffer((run / 'frames.bin').read_bytes(), dtype='<f8')
    n = len(raw) // rec
    out = []
    for k in range(n):
        r = raw[k * rec:(k + 1) * rec]
        o = 2
        f = {'step': int(r[0]), 'time_s': r[1]}
        for name in ('n_e', 'n_ion', 'n_h2o', 'mean_e_eV', 'Ex', 'Ey'):
            f[name] = r[o:o + nx * ny].reshape(ny, nx)
            o += nx * ny
        f['phi'] = r[o:o + (nx + 1) * (ny + 1)].reshape(ny + 1, nx + 1)
        out.append(f)
    meta = [json.loads(x) for x in (run / 'frames.jsonl').read_text().splitlines() if x.strip()]
    return lay, out, meta


def foils(ax):
    for y0 in (1.8, 7.0):
        ax.add_patch(Rectangle((2.0, y0), 5.0, 0.2, facecolor='#9a9893', edgecolor=INK2, linewidth=.6, zorder=3))


DENSITIES = (('n_e', 'electron density', 'Blues'), ('n_ion', 'H$_2$O$^+$ density', 'Oranges'), ('n_h2o', 'H$_2$O density', 'Greens'))
# Signed fields are shown in fixed display units so tick labels never change width between frames.
SIGNED = (('phi', 'quasi-static potential $\\phi$', 'kV', 1e-3), ('Ex', '$E_x$', 'MV/m', 1e-6), ('Ey', '$E_y$', 'MV/m', 1e-6))
SYMLOG_DECADES = 3


def global_scales(frames_):
    """One color scale per quantity for the whole run, so a color means the same value in every frame."""
    scales = {}
    for key, _, _ in DENSITIES:
        positive = [f[key][f[key] > 0] for f in frames_ if np.any(f[key] > 0)]
        if positive:
            vmax = max(float(p.max()) for p in positive)
            scales[key] = LogNorm(vmin=max(vmax * 1e-6, min(float(p.min()) for p in positive)), vmax=vmax)
    for key, _, _, factor in SIGNED:
        m = max(float(np.max(np.abs(f[key]))) for f in frames_) * factor or 1.0
        scales[key] = SymLogNorm(linthresh=m * 10.0 ** -SYMLOG_DECADES, linscale=0.5, vmin=-m, vmax=m, base=10)
    return scales


def fixed_axes(fig):
    """Fixed panel and colorbar positions (figure fractions); no automatic re-layout per frame."""
    left, width, gap, cbw = 0.05, 0.215, 0.115, 0.012
    bottoms, height = (0.535, 0.075), 0.36
    panels, bars = [], []
    for b in bottoms:
        for c in range(3):
            x = left + c * (width + gap)
            panels.append(fig.add_axes([x, b, width, height]))
            bars.append(fig.add_axes([x + width + 0.008, b, cbw, height]))
    return panels, bars


def density_panel(ax, cax, data, title, extent, cmap, norm):
    if norm is not None and np.any(data > 0):
        im = ax.imshow(np.where(data > 0, data, np.nan), origin='lower', extent=extent, cmap=cmap, norm=norm, interpolation='nearest')
        ax.figure.colorbar(im, cax=cax).set_label('cm$^{-3}$')
    else:
        ax.imshow(np.zeros_like(data), origin='lower', extent=extent, cmap='Greys', vmin=0, vmax=1, interpolation='nearest')
        ax.text(.5, .5, 'none in the midplane slab', transform=ax.transAxes, ha='center', va='center', color=INK2)
        if norm is not None:
            ax.figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax).set_label('cm$^{-3}$')
        else:
            cax.set_axis_off()
    ax.set_title(title)


def signed_panel(ax, cax, data, title, extent, unit, norm):
    im = ax.imshow(data, origin='lower', extent=extent, cmap='RdBu_r', norm=norm, interpolation='nearest')
    m = norm.vmax
    bar = ax.figure.colorbar(im, cax=cax, ticks=[-m, -m / 10, -m / 100, 0, m / 100, m / 10, m])
    bar.set_label(unit)
    bar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: '0' if v == 0 else '%.2g' % v))
    ax.set_title(title)


def render_frame(lay, f, meta, path, run_label, scales):
    nx, ny, dx = lay['nx'], lay['ny'], lay['dx_m'] * 1e3
    ext = [0, nx * dx, 0, ny * dx]
    ext_nodes = [-dx / 2, nx * dx + dx / 2, -dx / 2, ny * dx + dx / 2]
    fig = plt.figure(figsize=(13, 8))
    axs, bars = fixed_axes(fig)
    for (key, title, cmap), ax, cax in zip(DENSITIES, axs[:3], bars[:3]):
        density_panel(ax, cax, f[key], title, ext, cmap, scales.get(key))
    for (key, title, unit, factor), ax, cax in zip(SIGNED, axs[3:], bars[3:]):
        signed_panel(ax, cax, f[key] * factor, title, ext_nodes if key == 'phi' else ext, unit, scales[key])
    for ax in axs:
        foils(ax)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
    dbl = meta.get('debye_over_cell')
    dtxt = 'no electrons' if dbl in (None, 'inf') else '$\\lambda_D/\\Delta x$ = %.3g' % dbl
    fig.suptitle('%s   t = %.4f ns   step %d   midplane z = 0 (slab |z| <= %.3g mm)   %s   (fixed color scales; fields symmetric-log)' %
                 (run_label, f['time_s'] * 1e9, f['step'], lay['slab_m'] * 1e3 / 2, dtxt), color=INK, y=0.975)
    fig.savefig(path, dpi=90)
    plt.close(fig)


def line(ax, x, y, label, k):
    ax.plot(x, y, color=SERIES[k % len(SERIES)], label=label)


def timeseries(L, out):
    t = L['time_s'] * 1e9
    figs = {}
    # charge ledger
    fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    line(axs[0], t, L['q_total_change_since_genesis_C'], 'total charge change since genesis', 0)
    axs[0].set_ylabel('C')
    axs[0].legend(frameon=False)
    line(axs[1], t, L['gauss_max_abs_C'] / L['gauss_bound_C'], 'max |Gauss defect| / bound', 0)
    line(axs[1], t, L['divB_max_T_m'] / L['divB_bound_T_m'], 'max |div B| / bound', 1)
    axs[1].set_ylabel('ratio')
    axs[1].legend(frameon=False)
    line(axs[2], t, L['q_conductor_bank_C'], 'conductor bank', 0)
    line(axs[2], t, L['q_exterior_bank_C'], 'exterior bank', 1)
    line(axs[2], t, L['q_free_net_C'], 'free net', 2)
    axs[2].set_ylabel('C')
    axs[2].set_xlabel('t (ns)')
    axs[2].legend(frameon=False)
    fig.suptitle('charge ledger per step')
    fig.tight_layout()
    fig.savefig(out / 'ledger-charge.png', dpi=100)
    plt.close(fig)
    # energy closure
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    line(axs[0], t, np.abs(L['E_step_residual_J']), '|step energy residual|', 0)
    line(axs[0], t, L['E_roundoff_bound_J'], 'heat roundoff bound', 1)
    axs[0].set_yscale('log')
    axs[0].set_ylabel('J')
    axs[0].legend(frameon=False)
    line(axs[1], t, np.cumsum(L['E_joule_J']), 'cumulative Joule heat', 0)
    line(axs[1], t, np.cumsum(L['E_thermal_delta_J']), 'cumulative thermal change', 1)
    axs[1].set_ylabel('J')
    axs[1].set_xlabel('t (ns)')
    axs[1].legend(frameon=False)
    fig.suptitle('energy ledger per step')
    fig.tight_layout()
    fig.savefig(out / 'ledger-energy.png', dpi=100)
    plt.close(fig)
    # water inventory and arrivals
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    line(axs[0], t, L['arr_cum_released_h2o_molecules'], 'released H$_2$O molecules (n w0)', 0)
    line(axs[0], t, L['arr_cum_expected_h2o_molecules'], 'expected release (sum of quotes D)', 1)
    axs[0].set_ylabel('molecules')
    axs[0].legend(frameon=False)
    line(axs[1], t, L['carriers_h2o'], 'H$_2$O carriers', 0)
    line(axs[1], t, L['carriers_electron'], 'electron carriers', 1)
    line(axs[1], t, L['carriers_ion'], 'H$_2$O$^+$ carriers', 2)
    axs[1].set_ylabel('carriers')
    axs[1].set_xlabel('t (ns)')
    axs[1].legend(frameon=False)
    fig.suptitle('source release and census')
    fig.tight_layout()
    fig.savefig(out / 'ledger-inventory.png', dpi=100)
    plt.close(fig)
    # drive
    fig, axs = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    line(axs[0], t, L['current_A'] / 1e3, 'C19 drive current', 0)
    axs[0].set_ylabel('kA')
    line(axs[1], t, L['terminal_V'], 'terminal voltage', 0)
    axs[1].set_ylabel('V')
    axs[1].set_xlabel('t (ns)')
    for ax in axs:
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / 'drive.png', dpi=100)
    plt.close(fig)


PHASES = [('t_observe_ms', 'observe owners'), ('t_first_exposure_ms', 'first exposure'), ('t_source_field_ms', 'source field'),
          ('t_source_quotes_ms', 'desorption quotes'), ('t_emission_exposure_ms', 'emission exposures'), ('t_arrival_plan_ms', 'arrival plan'),
          ('t_emission_births_ms', 'emission births'), ('t_source_trials_ms', 'source trials'), ('t_birth_credit_ms', 'birth credit'),
          ('t_hall_admit_census_ms', 'Hall census admission'), ('t_hall_motion_ms', 'Hall motion (push)'), ('t_hall_deposit_support_ms', 'deposit support'),
          ('t_hall_field_prepare_ms', 'Hall field prepare'), ('t_hall_field_solve_ms', 'field solve'), ('t_hall_v2_seal_ms', 'Hall record seal'),
          ('t_hall_heat_ms', 'Hall heat'), ('t_native_motion_admission_ms', 'native motion admission'), ('t_second_exposure_ms', 'second exposure'),
          ('t_seal_ms', 'seal'), ('t_can_publish_ms', 'preflight (can_publish)'), ('t_swap_window_ms', 'swap window'), ('t_post_commit_ms', 'arrival commit')]


def timing(L, out, label):
    med = {name: float(np.median(L[c])) for c, name in PHASES if c in L}
    fig, ax = plt.subplots(figsize=(8, 6))
    names = list(med)[::-1]
    ax.barh(names, [med[n] for n in names], color=SERIES[0], height=.6)
    ax.set_xlabel('median ms per step')
    ax.set_title('%s: per-phase median (step total %.1f ms)' % (label, float(np.median(L['t_step_total_ms']))))
    for i, n in enumerate(names):
        ax.text(med[n], i, ' %.2f' % med[n], va='center', color=INK2, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / 'timing-phases.png', dpi=100)
    plt.close(fig)
    return med


def main():
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    gif_max = int(sys.argv[sys.argv.index('--gif-max') + 1]) if '--gif-max' in sys.argv else 120
    label = sys.argv[sys.argv.index('--label') + 1] if '--label' in sys.argv else run.name
    out.mkdir(parents=True, exist_ok=True)
    (out / 'frames').mkdir(exist_ok=True)
    L, cols = ledger(run)
    timeseries(L, out)
    med = timing(L, out, label)
    summary = {'run': str(run), 'steps': int(len(L['step'])), 'first_time_s': float(L['time_s'][0]), 'last_time_s': float(L['time_s'][-1]),
               'phase_median_ms': med, 'step_total_ms': {'median': float(np.median(L['t_step_total_ms'])), 'mean': float(np.mean(L['t_step_total_ms'])),
                                                          'p99': float(np.percentile(L['t_step_total_ms'], 99))}}
    debye = []
    if (run / 'frames.bin').exists():
        lay, fr, meta = frames(run)
        scales = global_scales(fr)
        paths = []
        for f, m in zip(fr, meta):
            p = out / 'frames' / ('frame_%05d.png' % m['frame'])
            render_frame(lay, f, m, p, label, scales)
            paths.append(p)
            debye.append({k: m[k] for k in ('frame', 'step', 'time_s', 'ne_max_cell_m3', 'Te_at_max_eV', 'debye_length_m', 'cell_size_m', 'debye_over_cell',
                                             'ne_max_slab_cm3', 'max_midplane_E_V_m', 'gap_voltage_midline_V')})
        pick = paths if len(paths) <= gif_max else [paths[int(round(i * (len(paths) - 1) / (gif_max - 1)))] for i in range(gif_max)]
        imgs = [Image.open(p).convert('RGB').quantize(colors=128, method=Image.Quantize.MEDIANCUT) for p in pick]
        imgs[0].save(out / 'pilot-midplane.gif', save_all=True, append_images=imgs[1:], duration=150, loop=0)
        summary['frames'] = len(paths)
        summary['gif_frames'] = len(pick)
        (out / 'debye-vs-cell.json').write_text(json.dumps(debye, indent=1) + '\n')
        with (out / 'debye-vs-cell.csv').open('w') as h:
            h.write('frame,step,time_ns,ne_max_cell_m3,Te_at_max_eV,debye_length_m,cell_size_m,debye_over_cell\n')
            for d in debye:
                h.write('%d,%d,%.9f,%.6g,%.6g,%s,%.6g,%s\n' % (d['frame'], d['step'], d['time_s'] * 1e9, d['ne_max_cell_m3'], d['Te_at_max_eV'],
                                                              d['debye_length_m'], d['cell_size_m'], d['debye_over_cell']))
    (out / 'render-summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'phase_median_ms'}, indent=1))


if __name__ == '__main__':
    main()
