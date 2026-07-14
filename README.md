# Power Flow MHD/CR Prototype

This project is a research scaffold for 2D reduced MHD electrode heating, Temkin
desorption, kinetic neutral expansion, BOLSIG electron kinetics, collisional
radiative populations, and CCD-based synthetic spectra.

The physics direction is intentionally ALEGRA-MHD inspired: transient
magnetics, conductivity-dependent current diffusion, Joule heating, thermal
conduction, material response, and circuit/drive coupling. AGATE is used only as
a framework-shape reference for modular Python/GPU ergonomics.

Run the example:

```powershell
python -m pip install -r requirements.txt
python MHD_heating.py --config configs/mykonos.yaml
```

The default `configs/mykonos.yaml` case is a user-initialized starter problem,
not a measured-shot reconstruction yet. It uses two 200 um thick, 5 mm wide
SS304 foils separated by a 5 mm face-to-face A-K gap, driven by a zero-voltage
parametric current pulse with an 850 kA peak and a 125 ns quarter-period rise
time. Field histories are sampled every 1 ns. The main config uses a refined
uniform mesh of roughly 20 um cells, giving about ten cells through each foil.
Use `configs/mykonos_preview.yaml` for fast coarse checks.

Measured voltage/current CSV traces can still be used later by setting
`drive.mode` back to `measured_current` or `fit_rl` and providing `trace_csv`.

The MHD package is still a reduced research model. It now solves a
resistive-induction update for `A_z` and `J_z`, enforcing the total current as a
constraint instead of assigning a prescribed current-density profile. That is a
better basis for current crowding and nonuniform Joule heating, but it is not a
full production compressible/resistive MHD model.

Create quick-look result plots after a run:

```powershell
python tools/visualize_run.py runs\run_YYYYMMDD_HHMMSS --all-gifs
```

This writes `overview_last.png`, `time_traces.png`, and optional GIFs under the
run directory's `visualization/` folder.

The reduced MHD solver uses NumPy by default and will use CuPy when
`mhd.backend: auto` or `cuda` finds a working CUDA/CuPy install. This machine
has been tested with `cupy-cuda13x[ctk]`.

A C++/pybind11 extension is included under `src/mhd_core`; it builds through
scikit-build-core. If CMake is not on PATH, prepend:

```powershell
$env:PATH = 'C:\Users\tsmit\Downloads\cmake-4.4.0-windows-x86_64\cmake-4.4.0-windows-x86_64\bin;' + $env:PATH
python -m pip install -e . --no-build-isolation
```

The bundled BOLSIG+ executable has been copied into `third_party/bolsigplus`.
The current wrapper prefers `bolsigminus.exe` and has a timeout/fallback path.
The wrapper accepts BOLSIG's Windows/Fortran EOF exit code when a parseable
output table was written, records stdout/stderr/log paths, and parses mean
energy plus per-collision rate-coefficient blocks. The included Siglo database
does not contain the full H2O/OH/H/O electron collision set, so the example can
optionally supplement missing provisional H2O rates until curated tables are
added.

Data provenance is saved under `data/`:

- `data/cross_sections/` stores the PRAB 26, 040401 refs. 75-80 reaction map plus NIST starter electron-impact CSVs.
- `data/cr_model/` stores the Yacora/Fantz/Janev/NIST/Fujimoto CR references and rovibrational table manifest.
- Rows with `file: null` are intentional placeholders for measured or digitized tables that still need to be added.
