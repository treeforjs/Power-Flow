# Analytic resistive-heating checks for FLASH

`mhdlab.heating_benchmarks.ConstantConductivityHalfSpace` supplies an analytic
reference for a controlled diffusion problem. It is neither a replacement for
FLASH nor an MITL shot model. Conductivity, permeability, volumetric heat capacity
and imposed tangential field are explicit supplied parameters.

A stationary conductor occupies inward depth y ≥ 0. Its initial magnetic field
is zero, and a constant tangential field Bs is applied at y = 0 for t > 0.
With D = 1/(μσ),

`Bx(y,t) = Bs erfc[y / (2 sqrt(Dt))]`,

`Jz(y,t) = Bs / [μ sqrt(πDt)] exp[-y²/(4Dt)]`,

`qJ = Jz²/σ`.

The positive inward depth direction fixes the current sign through
`Jz = -(∂Bx/∂y)/μ`. Fluid motion, thermal conduction, radiation, phase changes and
temperature-dependent coefficients are absent. A FLASH verification problem
must implement those same assumptions before its errors are compared.

## Generate and compare references

```bash
python tools/run_heating_benchmark.py --plot
```

The default [configuration](../configs/heating_diffusion_manufactured.json)
supplies manufactured constant coefficients and positive comparison times.
The output in `runs/heating_benchmark_reference/` contains typed cell-average
CSV profiles and SHA-linked metadata sidecars, half-space energy budgets,
pointwise adiabatic temperature diagnostics at cell centres, and an optional
figure. No FLASH calculation is performed.

An independently exported FLASH profile can use the same canonical SI contract:

```python
import json
from mhdlab.heating_benchmarks import ConstantConductivityHalfSpace
from mhdlab.heating_comparison import load_heating_profile_csv, compare_heating_profile

configuration = json.load(open("configs/heating_diffusion_manufactured.json"))
benchmark = ConstantConductivityHalfSpace(**configuration["benchmark"])
profile = load_heating_profile_csv("/path/to/flash_depth_profile.csv")
report = compare_heating_profile(profile, benchmark)
```

The required sidecar declares scalar time, exact SI units, inward depth and
cell-average semantics, source kind/run/version and the CSV hash. Every CSV row
must be one contiguous depth cell. Comparisons use the profile's actual grid
and time, with width-normalized L1/L2 errors and Linf error for B, J and Joule
power. Relative errors are explicitly undefined when the reference norm is zero.
The parameters must match the independently configured FLASH benchmark case.

## Independent energy budget

The analytic half-space integrals per unit surface area are

`U_B = (Bs²/μ) sqrt(Dt/π) (2 - sqrt(2))`,

`Q_J = (Bs²/μ) sqrt(2Dt/π)`,

`W_in = 2 (Bs²/μ) sqrt(Dt/π)`.

Thus incoming Poynting energy equals magnetic storage plus accumulated Joule
heat. Surface power is `(Bs²/μ) sqrt[D/(πt)]`. It is singular at the ideal field
step's initial instant; the API requires t > 0. These quantities distinguish
magnetic storage from heat instead of assigning all incoming energy to temperature.

The adiabatic local temperature diagnostic is

`ΔT(y,t) = Bs² / (μπCv) E1[y²/(2Dt)]`, for y > 0,

where Cv is heat capacity per volume. Its y = 0 value diverges for this idealized
step. The implementation explicitly rejects that surface-temperature query.
Use a smooth-drive benchmark or a finite source model for a finite boundary
temperature. Do not insert this singular diagnostic into the desorption workflow.

## Compare the same numerical quantities

`fields()` returns point fields. `cell_averages()` separately integrates Bx,
Jz and Joule power over each supplied depth cell. Finite-volume FLASH outputs
must be compared with these cell averages. In particular, mean(J²)/σ differs
from mean(J)²/σ across a skin layer. Temperature at a cell centre is also distinct
from a cell-averaged temperature and is not silently interchanged.

The finite computational depth must be large enough that its outer boundary
does not change the half-space solution over the comparison time. Refine normal
resolution and timestep independently; keep the magnetic boundary drive and
initial conditions fixed. A comparison on one grid does not demonstrate convergence.

After this constant-property check, add separate verification cases for thermal
conduction, variable material properties, energy-conserving phase response,
finite foil width/current crowding and the selected Hall-MHD closure. Do not use
the dissertation's unverified literal Jz axis scale as an absolute Joule-heating
reference. Published thermal profiles provide separate reconstruction benchmarks.
