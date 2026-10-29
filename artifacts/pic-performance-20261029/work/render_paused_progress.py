"""Render accepted native density frames; never interpolate simulation states."""
from pathlib import Path
import hashlib
import io
import json
import os
import fcntl
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from PIL import Image

W = Path('/home/treefor/Documents/Codex/2026-10-28/we')
PRODUCTION = W / 'outputs/electron-heating-pic-25um-full60ns-gpu-midplane110kV-height25mm'
SNAPSHOT = W / 'outputs/pic-profile/paused-checkpoint-snapshot/frame_27.csv'
OUT = W / 'outputs/electron-mapping-paused-26p48ns'
BROOT = OUT / 'native-b'

def render():
    OUT.mkdir(exist_ok=True)
    paths = sorted(PRODUCTION.glob('frame_*.csv'), key=lambda p: int(p.stem.split('_')[-1]))
    paths.append(SNAPSHOT)
    data, inputs = [], []
    for path in paths:
        with path.open() as f:
            names = next(f).strip().split(',')
        a = np.loadtxt(path, delimiter=',', skiprows=1, dtype=[(n, float) for n in names])
        if len(a) != 360 * 360 or any(not np.isfinite(a[n]).all() for n in names):
            raise ValueError('Invalid native frame: ' + str(path))
        if not np.all(a['time_s'] == a['time_s'][0]):
            raise ValueError('Mixed frame times')
        if np.any(a['electron_density_cm3'] < 0) or np.any(a['height_average_electron_density_cm3'] < 0):
            raise ValueError('Negative native density')
        if data and a['time_s'][0] <= data[-1]['time']:
            raise ValueError('Nonmonotone accepted frames')
        solid = a['solid'].reshape(360, 360).astype(bool)
        density = a['electron_density_cm3'].reshape(360, 360).copy()
        if np.any(density[solid] != 0):
            raise ValueError('Density in a conductor')
        represented = float(density.sum() * (25e-6)**2 * .001 * 1e6)
        bpath = BROOT / (path.stem + '_B.csv')
        with bpath.open() as f:
            bnames = next(f).strip().split(',')
        b = np.loadtxt(bpath, delimiter=',', skiprows=1, dtype=[(n, float) for n in bnames])
        if len(b) != len(a) or any(not np.isfinite(b[n]).all() for n in bnames):
            raise ValueError('Invalid native magnetic frame')
        if not np.allclose(b['time_s'], a['time_s'], rtol=0, atol=1e-22):
            raise ValueError('Magnetic/density frame time mismatch')
        if not np.array_equal(a['x_m'],b['x_m']) or not np.array_equal(a['y_m'],b['y_m']):
            raise ValueError('Magnetic/density coordinates differ')
        data.append(dict(time=float(a['time_s'][0]), density=density,
                         potential=a['potential_V'].reshape(360, 360).copy() / 1000.,
                         solid=solid, represented=represented,
                         bx=b['Bx_T'].reshape(360,360).copy(),
                         by=b['By_T'].reshape(360,360).copy(),
                         bmag=b['Bmag_T'].reshape(360,360).copy(),
                         reconstructed=(path == SNAPSHOT)))
        inputs.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                           magnetic_path=str(bpath),magnetic_sha256=hashlib.sha256(bpath.read_bytes()).hexdigest(),
                           time_ns=float(a['time_s'][0]) * 1e9,
                           represented_electrons_in_slab=represented))
    positives = [d['density'][d['density'] > 0] for d in data if np.any(d['density'] > 0)]
    if not positives:
        raise ValueError('No positive density in the accepted history')
    maximum = max(float(v.max()) for v in positives)
    minimum = min(float(v.min()) for v in positives)
    norm = LogNorm(vmin=min(1e3, minimum), vmax=maximum)
    cmap = plt.get_cmap('YlGnBu').copy(); cmap.set_bad('white')
    edges = np.linspace(-4.5, 4.5, 361)
    centres = .5 * (edges[1:] + edges[:-1])
    frames = []
    plt.rcParams.update({'figure.facecolor':'white', 'axes.facecolor':'white',
                         'savefig.facecolor':'white', 'font.size':10})
    for d in data:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios':[1.05, 1]})
        fig.subplots_adjust(left=.065, right=.88, bottom=.24, top=.79, wspace=.27)
        density = np.ma.masked_where(d['solid'] | (d['density'] <= 0), d['density'])
        magnitude = np.ma.masked_where(d['solid'], d['bmag'])
        # One immutable, uniform absolute Tesla contour family for every frame.
        levels = np.arange(2., 60.1, 2.)
        for ax in axes:
            m = ax.pcolormesh(edges, edges, density, cmap=cmap, norm=norm, shading='flat')
            visible = levels[(levels>magnitude.min()) & (levels<magnitude.max())]
            if len(visible):
                contours=ax.contour(centres, centres, magnitude, levels=visible,
                           colors='#bf6511', linewidths=.55, alpha=.8)
                labels=[v for v in visible if v in (2,10,20,30,40,50)]
                if labels:ax.clabel(contours,levels=labels,fmt='%g T',fontsize=7,inline=True)
            ax.add_patch(Rectangle((-2.5, 2.5), 5, .2, facecolor='#596777',
                                   edgecolor='#263342', linewidth=1.0, zorder=5))
            ax.add_patch(Rectangle((-2.5, -2.7), 5, .2, facecolor='#596777',
                                   edgecolor='#263342', linewidth=1.0, zorder=5))
            ax.set_aspect('equal')
            ax.set_xlabel('Width from centre [mm]')
            ax.set_ylabel('Gap-normal position [mm]')
        axes[0].set(xlim=(-4.5,4.5), ylim=(-4.5,4.5), title='Both electrodes and surrounding vacuum')
        axes[0].text(0,2.61,'Cathode',ha='center',va='center',color='white',fontsize=9,zorder=6)
        axes[0].text(0,-2.61,'Anode',ha='center',va='center',color='white',fontsize=9,zorder=6)
        axes[1].set(xlim=(2.05,2.95),ylim=(2.05,2.95),title='Cathode edge: gap, side and back faces')
        for ax, stride in zip(axes,(24,6)):
            sel=np.s_[stride//2::stride,stride//2::stride]
            bx,by=d['bx'][sel],d['by'][sel]
            mag=np.hypot(bx,by)
            valid=(mag>0) & ~d['solid'][sel]
            xx,yy=np.meshgrid(centres[stride//2::stride],centres[stride//2::stride])
            ax.quiver(xx[valid],yy[valid],bx[valid]/mag[valid],by[valid]/mag[valid],
                      color='#8c6b4f',angles='xy',scale_units='xy',scale=8 if stride==24 else 25,
                      width=.0023,alpha=.6,zorder=4)
        fig.colorbar(m, ax=axes, fraction=.035, pad=.035,
                     label='Electron number density [cm$^{-3}$]')
        fig.suptitle(f'Heating-driven electron mapping | {d["time"]*1e9:.4f} ns',
                     fontsize=15, fontweight='bold', y=.965)
        fig.text(.5,.875,'25 µm cells | 5 mm foil width × 5 mm AK gap × 25 mm target height',
                 ha='center',fontsize=11)
        note = 'Saved-checkpoint reconstruction; no further time evolution' if d['reconstructed'] else 'Accepted native frame'
        fig.text(.5,.825, note,ha='center',fontsize=9,color='#4b5563')
        if d['represented'] == 0:
            axes[0].text(.5,.02,'No represented electrons in this sampling slab',
                         transform=axes[0].transAxes,ha='center',fontsize=9,
                         bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
        fig.text(.5,.105,'Density: |z| < 0.5 mm slab. Orange: |B| contours every 2 T (fixed); arrows: B direction.',
                 ha='center',fontsize=10)
        fig.text(.5,.063,'Conditional pulse (110 kV midpoint peak); 3D trajectories / height-averaged 2D space charge.',
                 ha='center',fontsize=9,color='#4b5563')
        fig.text(.5,.025,'Mirrored anode temperature and coarse prescribed B; collisions/ions absent; sheath/grid convergence pending.',
                 ha='center',fontsize=9,color='#4b5563')
        buffer=io.BytesIO();fig.savefig(buffer,format='png',dpi=110)
        frames.append(Image.open(buffer).convert('RGB').copy())
        if d is data[-1]:
            fig.savefig(OUT/'latest-electron-map.png',dpi=150)
        plt.close(fig)
    durations=[100 if d['time'] < 24e-9 else 650 for d in data]
    durations[-1]=2500
    gif=OUT/'electron-history-with-B-through-26p48ns.gif'
    frames[0].save(gif,save_all=True,append_images=frames[1:],duration=durations,loop=0,
                   disposal=2,optimize=False)
    with Image.open(gif) as im:
        count=im.n_frames
        for i in range(count):im.seek(i);im.load()
    report=dict(gif=str(gif),frames=count,accepted_time_range_ns=[data[0]['time']*1e9,data[-1]['time']*1e9],
                final_frame_is_checkpoint_reconstruction=True,physical_advance_during_snapshot_s=0,
                density_units='cm^-3',density_sampling_depth_m=.001,background='white',
                magnetic_contour_levels_T=levels.tolist(),magnetic_arrows='direction only, normalized',
                density_scale_cm3=[norm.vmin,norm.vmax],inputs=inputs,
                cpu_affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),
                renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                not_full_60ns=True,convergence_qualified=False,fragmentation_or_ion_feedback=False)
    (OUT/'render-receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'},indent=2))

if __name__ == '__main__':
    with (W/'outputs/pic-profile/testing-core.lock').open('a') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        render()
