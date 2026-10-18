# Linear-ramp magnetic-diffusion reference

The smooth-drive reference provides finite-volume magnetic field, current and
Joule-power averages for a stationary constant-property conductor. It also
provides magnetic storage, Joule heat and incoming Poynting energy over the
infinite half-space. This is an independent verification reference for later
FLASH work. The supplied example is manufactured and does not calibrate a
foil material or reproduce a measured shot.

Run from the repository with the configured environment:

```sh
.venv/bin/python tools/run_smooth_heating_benchmark.py \
  --config configs/smooth_heating_manufactured.json \
  --output-dir runs/smooth_heating_reference --plot
```

The output directory must be new. Generation takes place in a temporary sibling
directory, and the completed result is then moved into place. The command
refuses existing output paths, changed inputs and changed source files.
Use a dedicated path with one coordinated writer. The final existence check and
directory rename do not provide a lock against a competing writer creating an
empty destination directory in between them.

The configuration supplies every physical constant in SI: electrical
conductivity in S/m, permeability in H/m, heat capacity per volume in J/(m³ K),
and the signed boundary field rate in T/s. Time begins with zero field and the
boundary follows `B_s(t) = rate * t` for all positive times. There is no ramp
stop or inferred pulse history. Input status and provenance remain explicit.
The canonical field CSV contract requires positive times; the mathematical
reference also returns a zero initial state at time zero.

For magnetic diffusivity `D = 1/(mu*sigma)` and inward depth `y`, let
`u = y/(2*sqrt(D*t))`. The reference uses
`g(u) = exp(-u²)/sqrt(pi) - u*erfc(u)` and
`b(u) = (1+2u²)*erfc(u) - 2u*exp(-u²)/sqrt(pi)`:

- `Bx = rate*t*b(u)` and `Jz = 2*rate*t*g(u)/(mu*sqrt(D*t))`.
- The Joule-power mean integrates `Jz²/sigma` over each cell independently.
  Squaring the cell-average current loses subcell heating.
- The wall heating is finite: `DeltaT_wall = 2*rate²*t²/(mu*pi*Cv)`.
  Exported temperatures are local adiabatic rises, with no supplied absolute
  initial temperature or heat conduction.

The result contains canonical `profiles/time_*.csv` files with required,
hash-linked JSON sidecars; `energy_budget.csv`; and
`adiabatic_point_temperature.csv` with the wall and each cell center labeled
separately. `input_config.json` preserves the exact input bytes. `summary.json`
records SI units, source/configuration/output hashes, numerical qualifications,
finite-depth field fractions and model limits. `overview.png` is optional.
The finite-depth field fraction evaluates the dimensionless profile directly;
it does not divide dimensional fields whose common scale can overflow or
underflow. A zero ramp has no field ratio and is recorded as null.

Half-space energy obeys `W_in = U_B + Q_J`. With
`S = rate²*sqrt(D)*t^(5/2)/mu`, the three terms are
`U_B = (28-16*sqrt(2))*S/(15*sqrt(pi))`,
`Q_J = 16*(sqrt(2)-1)*S/(15*sqrt(pi))`, and
`W_in = 4*S/(5*sqrt(pi))`. The algebraic energy residual checks reference
generation. Export/readback checks interchange. Neither establishes a numerical
solver's convergence or energy conservation.

When a licensed FLASH setup is available, match the drive, signs, units and
constant coefficients, compare cell averages, and refine both the mesh and
time step. Move the outer material boundary sufficiently far from the diffusion
layer and check it separately; the exported energy budget extends to infinity.
The reference excludes flow, Hall terms, thermal conduction, displacement
current, radiation, phase changes and an equation of state. It cannot qualify
those parts of the eventual foil calculation.
