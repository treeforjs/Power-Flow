from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shutil
import zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/neutral_desorption_perimeter_20261028'
HEAT=ROOT/'runs/mfem_edge_current_source_four_faces'
MAPS=ROOT/'runs/neutral_desorption_current_17133_perimeter_v2'
VIEW=ROOT/'runs/neutral_desorption_current_17133_perimeter_gif'
EDGE=ROOT/'runs/neutral_desorption_current_17133_perimeter_edge_gif'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(p, relative):
    target=OUT/relative
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(p,target)
OUT.mkdir(parents=True,exist_ok=False)
for label, directory in [('overview',VIEW),('edge',EDGE)]:
    copy(directory/'neutral_transport.gif',f'neutral_{label}_equal_scale.gif')
    for name in ['preview_first.png','preview_middle.png','preview_last.png','render_report.json']:
        copy(directory/name,f'{label}/{name}')
copy(EDGE/'view_consistency.json','view_consistency.json')
copy(ROOT/'runs/mfem_alegra_comparison_v3/mfem_alegra_comparison.png','mfem_alegra_comparison.png')
copy(ROOT/'runs/mfem_alegra_comparison_v3/comparison.json','inputs/comparison.json')
for label in ['top','bottom','left','right']:
    for suffix in ['.csv','_surface.json']:
        copy(HEAT/(label+suffix),'inputs/heating/'+label+suffix)
    copy(MAPS/(label+'_release.h5'),'inputs/transport/'+label+'_release.h5')
for name in ['configuration.json','source_run.json','diagnostics.csv','schedule.csv','foil_snapshots.csv','foil_snapshots.json','four_face_verification.json']:
    copy(HEAT/name,'inputs/heating/'+name)
for name in ['case.json','metadata.json','neutral_maps.npz','particle_ledgers.json','verification.json']:
    copy(MAPS/name,'inputs/transport/'+name)
for p in (ROOT/'mhdlab').rglob('*.py'):
    copy(p,'replay/'+p.relative_to(ROOT).as_posix())
for name in ['run_neutral_desorption_map.py','render_neutral_transport_gif.py','run_mfem_edge_source.py','compare_mfem_alegra.py']:
    copy(ROOT/'tools'/name,'replay/tools/'+name)
for name in ['mfem_edge_current_source.json','mfem_surface_example.json','neutral_desorption_map_17133_perimeter.json']:
    copy(ROOT/'configs'/name,'replay/configs/'+name)
for p in (ROOT/'native/mfem_edge_source').iterdir():
    if p.is_file(): copy(p,'replay/native/mfem_edge_source/'+p.name)
for stem in ['shot_17133_current','figure_4_28_temperature','figure_4_29_temperature','figure_4_30_heating','figure_4_28_simulation_drive','figure_4_29_simulation_drive']:
    for ext in ['.csv','.json']:
        copy(ROOT/'data/dissertation'/(stem+ext),'replay/data/dissertation/'+stem+ext)
for name in ['neutral_perimeter_release.md','mfem_edge_current_source.md','mfem_alegra_comparison.md','neutral_transport_gif.md']:
    copy(ROOT/'docs'/name,'reference/'+name)

readme='''# Four-face desorption with equal physical length scales

This update corrects two limitations of the previous GIF: its y direction was
enlarged fourfold, and its source emitted only from the two broad foil faces.
Both side faces now emit from actual MFEM material-boundary temperature traces.
The current drive, material assumptions and 60 ns stop time are unchanged.

Open **neutral_edge_equal_scale.gif** for the right-edge close-up, and
**neutral_overview_equal_scale.gif** for the whole foil. Both show H, H2, H2O
and OH with equal x/y length scales, the same 65 times over 0.749788-60 ns,
the same data and the same global color ranges. The edge close-up crops the
axes to x=2.25-2.9 mm and y=2.2-3.0 mm; it changes no trajectories or densities.
The original 5.2-second playback loop retains every computed frame.

## What changed and what was checked

- A separate native material-side quadrature samples eight patches through
  each 200-micrometre-thick side. Left/right normals are physical -x/+x. Their
  local planar launch frames are explicitly permuted into physical xyz; speed,
  kinetic energy and particle weight are preserved.
- Top/bottom temperature tables, cell snapshots, current schedule and heating
  diagnostics are byte-identical to the earlier current-driven run. Every
  numeric dataset in the original broad-face release histories is unchanged.
  Tiny-tail packet allocation changes the Monte Carlo realization, while the
  physical broad-face source histories remain identical.
- The side pair has 4% of the broad pair's area. Equal per-area initial coverage
  and desorption-rate laws are assumed on all four faces. At 60 ns the sides
  supply about **9.11% of emitted particles** because their hotter pools deplete
  earlier. No side-area multiplier or x-velocity boost is introduced.
- For H emitted by the side faces, RMS displacement in x by 60 ns is about
  **150 micrometres**. The side H packets have a weighted mean age of about
  33.75 ns. Broad-face H RMS x displacement is approximately 41-56 micrometres
  depending on the face. Packet-temperature-conditioned velocity second moments
  match the flux-weighted thermal law with the normal axis on the correct face.
- 4,288,000 weighted packets are observed. Maximum particle-number residual is
  6.08e-17 of initial inventory; histogram-volume residual is 5.50e-15. The
  actual encoded overview and edge GIFs were decoded and their first, middle
  and last frames inspected. `view_consistency.json` verifies identical inputs,
  times, colors and whole-map particle totals between the views.
- Focused frame/provenance, conservative flight, numerical-tail and renderer
  tests pass. The source/frame suite covered 84 cases; the tail, flight and
  renderer suites respectively covered 8, 6 and 13 cases.

Fixed eight-decade per-species color ranges leave smaller positive densities
in the lowest color. Exact zero is white. No density floor, smoothing or time
interpolation is applied. Whole-map particle totals include regions outside
the cropped edge view.

## Physical and numerical assumptions

This remains a **provisional constant-property stationary-foil scenario**.
Current is the recovered dissertation shot 17133 Fig. C.19 trace, starting from
its first retained knot with a conditional uniform-DC magnetic state and 300 K.
The published current coordinate is not shifted; its registration against the
ALEGRA clock is unverified. The source includes a finite-width return conductor
and true nonconducting vacuum. No phase change, motion or temperature feedback
is included. Side boundary temperatures reach about 12074 K, outside the
constant-property solid model's physical range. The existing ALEGRA comparison
is included because the underlying heating snapshots are unchanged; hottest
edge cells still substantially overshoot those published profiles.

H, H2, H2O and OH are independent assumed adsorbate pools. Initial per-face
coverages are respectively 2e15, 2e15, 1.4e16 and 2e15 per m2, with common
attempt frequency 1e12/s and 0.8 eV activation/removal energies. This is not a
calibrated water-fragmentation law or an inferred absolute yield. All species
travel collisionlessly with three thermal velocity components and actual birth
times; z is invariant/periodic. No neutral force, chemistry or spectral fit is
added. The old broad-face thermal law naturally has a sqrt(2) RMS normal versus
tangential velocity ratio; its old 4x y display made that look like about 5.66.

Exhausted source pools retain machine-subnormal mathematical tails. Original
release HDF5 histories are preserved. A validated sampling-only view omits
components below the explicit 1e-250-particle representation threshold and
checks all retained release/energy components are bit-identical. In this run,
334 omitted components have a combined conservative bound of 3.34e-248
particles. Exact-decimal omission accounts are saved in transport metadata.

`inputs/heating` contains all four actual faces, current schedule and native
verification. `inputs/transport` includes full original release histories,
maps, ledgers and quantitative checks. The 205 MB packet file stays in the
original run directory and can be regenerated with the saved seed. Every
packaged file except `manifest.json` is SHA256-bound in that manifest.

## Replay

Use Python >=3.11 with numpy, scipy, h5py, matplotlib, pillow, PyYAML and imageio.
From this package's `replay` directory, using fresh output directories:

```powershell
python tools/render_neutral_transport_gif.py --input ../inputs/transport/neutral_maps.npz --metadata ../inputs/transport/metadata.json --output-dir runs/replayed_overview --fps 12 --dpi 110 --color-span-decades 8
python tools/render_neutral_transport_gif.py --input ../inputs/transport/neutral_maps.npz --metadata ../inputs/transport/metadata.json --output-dir runs/replayed_edge --fps 12 --dpi 110 --color-span-decades 8 --view-bounds-mm 2.25 2.9 2.2 3.0
python tools/run_neutral_desorption_map.py --top-surface ../inputs/heating/top_surface.json --bottom-surface ../inputs/heating/bottom_surface.json --left-surface ../inputs/heating/left_surface.json --right-surface ../inputs/heating/right_surface.json --case ../inputs/transport/case.json --output-dir runs/replayed_transport
```

Fresh heating additionally requires installed MFEM 4.10 with MPI/Hypre/Metis,
CMake/Ninja and the Ubuntu-24.04 WSL environment. The owned native driver is
included; MFEM itself is not vendored:

```powershell
python tools/run_mfem_edge_source.py --wsl-distro Ubuntu-24.04 --output-dir runs/recomputed_four_face_heating
```
'''
(OUT/'README.md').write_text(readme,encoding='utf-8')
verification=json.loads((MAPS/'verification.json').read_text())
views=json.loads((EDGE/'view_consistency.json').read_text())
record={'schema':'mhdlab-perimeter-desorption-delivery-v1','status':'provisional',
        'created_utc':datetime.now(timezone.utc).isoformat(),'time_cap_s':60e-9,
        'physical_experiment_prediction':False,'transport_verification':verification,
        'view_consistency':views,
        'files':{p.relative_to(OUT).as_posix():{'sha256':sha(p),'bytes':p.stat().st_size}
                 for p in sorted(OUT.rglob('*')) if p.is_file()}}
(OUT/'manifest.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
archive=OUT.with_suffix('.zip')
with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(OUT.rglob('*')):
        if p.is_file(): z.write(p,arcname=OUT.name+'/'+p.relative_to(OUT).as_posix())
with zipfile.ZipFile(archive) as z: assert z.testzip() is None
print(json.dumps({'package':str(OUT),'zip':str(archive),'zip_bytes':archive.stat().st_size,
                  'zip_sha256':sha(archive),'files':len(record['files'])}))
