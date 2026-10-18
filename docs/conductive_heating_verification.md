# Periodic induction and thermal-conduction reference

This exact reference checks resistive induction, Joule heating and thermal
conduction together in a stationary periodic conductor with constant
coefficients. It generates cell averages and a spatial mean energy budget for
comparison with a later numerical calculation. It does not run FLASH or predict
a dissertation shot. The supplied configuration is manufactured; its deliberately
large thermal conductivity makes redistribution visible on the selected
timescale and is not an SS304 property recommendation.

```sh
.venv/bin/python tools/run_conductive_heating_benchmark.py \
  --config configs/conductive_heating_manufactured.json \
  --output-dir runs/conductive_heating_reference --plot
```

Use a new output path and one coordinated writer. The command builds a temporary
sibling directory and moves the completed output into place, refuses existing
paths and checks source and input hashes before publication. These checks are
not a lock against another process creating an empty destination concurrently.

Every physical coefficient is required in SI: electrical conductivity `sigma`
in S/m, permeability `mu` in H/m, volumetric heat capacity `Cv` in J/(m³ K),
thermal conductivity `kappa` in W/(m K), wavenumber `k` in rad/m and the signed
initial field `B0` in T. `Cv` is per volume; this reference does not infer it
from mass-specific constant-pressure heat capacity, density or an equation of
state. Nonnegative `kappa` includes an explicit zero-conduction limit. Time
begins at zero, with `Bx = B0*cos(k*y)` and zero temperature rise. The configured
domain contains an integer number of periods `2*pi/k`, starting at the supplied
position origin.

For `eta = 1/(mu*sigma)`, `lambda = 2*eta*k²`, `chi = kappa/Cv` and
`delta = 4*chi*k²`, the reference is

- `Bx = B0*exp(-lambda*t/2)*cos(k*y)`.
- `Jz = -dBx/dy/mu = B0*k*exp(-lambda*t/2)*sin(k*y)/mu`.
- `q = Jz²/sigma`; its cell mean integrates the square independently.
- `DeltaT = B0²/(4*mu*Cv) * [1-exp(-x)-H(x,z)*cos(2*k*y)]`, where
  `x = lambda*t`, `z = delta*t` and
  `H = x*(exp(-x)-exp(-z))/(z-x)`.

At resonance `z=x`, `H=x*exp(-x)`. At `kappa=0`,
`H=1-exp(-x)` and the temperature rise is the local adiabatic deposition.
For positive conduction and late times, temperature approaches a uniform rise.
The heat-capacity-weighted spatial mean rise equals the lost magnetic energy:
`mean(U_B) = B0²*exp(-lambda*t)/(4*mu)` and
`mean(U_thermal) = B0²*[1-exp(-lambda*t)]/(4*mu)`.
Their sum is constant because there is no net periodic energy flux. Conduction
changes the profile but not the mean.

`cell_averages.csv` contains each time, cell index, left and right SI coordinates,
and exact spatial means of `Bx`, `Jz`, `q` and `DeltaT`. Temperature is a rise,
not an absolute material temperature. `energy_budget.csv` contains spatial mean
magnetic and thermal energies in J/m³, mean Joule power in W/m³ and the algebraic
energy residual. `input_config.json` preserves the input bytes. `summary.json`
records units, parameter status and provenance, derived coefficients, source,
configuration and output SHA-256 hashes, assumptions and limits. `overview.png`
is optional. Time zero is supported by this periodic export contract.

The implementation avoids subtracting nearly equal exponentials at resonance
and uses positive formulas for the small temperature rise at a field maximum.
It combines rational raw factors and logarithmic scales before rounding final
values, including very small dimensionless times and cell widths. Both
dimensionless times at most one use a positive 16-point Gauss-Legendre integral
of a smooth bounded kernel; independent convolution and quadrature tests check
this evaluation. Positions must satisfy `abs(k*y) <= 2^20` radians. Larger
phases, unrepresentable final results and unrepresentable derived coefficients
required by the CLI are rejected explicitly. Final sub-float values round to
zero. These are numerical representation bounds, not physical cutoffs.

For a numerical host comparison, match the signs, coefficients, initial data,
periodic geometry and the finite-volume meaning of each field. Refine spatial
and temporal resolution independently and compare both the thermal profile and
the magnetic-to-thermal transfer. Squaring mean current misses subcell Joule
heating. The exact reference's residual and a successful file export do not
establish convergence of a numerical solver. The benchmark does not qualify
Hall induction, fluid motion, variable material properties, radiation, phase
change, heat-capacity closure, surface desorption or an equation of state.
