# Coverage-dependent first-order surface release

This extension supplies an independent, prescribed-temperature boundary source.
It retains a decreasing surface inventory and a barrier that changes with that
inventory. It does not solve wall cooling, gap transport or molecular kinetics.
The constant-barrier baseline and its output schema remain separate.

For each independent binding reservoir, the implemented kinetic hypothesis is

\[
\dot N=-\nu N\exp\left[-\frac{E_d}{k_BT}
              \left(1-\alpha\frac{N}{N_{\rm ref}}\right)\right].
\]

Here \(N\) and \(N_{\rm ref}\) have units m⁻², \(\nu\) is a constant supplied
attempt frequency in s⁻¹, \(E_d\) is a kinetic energy in eV, and \(k_B\) in this
expression is in eV/K. The remaining population supplies the first-order
prefactor. The reference density stays fixed when the initial population
changes. A reference monolayer must have a documented definition; it is not
automatically the starting coverage. Equivalent changes of convention require
\(\alpha_{\rm new}/N_{\rm ref,new}=\alpha_{\rm old}/N_{\rm ref,old}\).

The [primary-source audit](thermal_desorption_law_audit.md) explains the scoped
support for this kinetic form and its difference from the dissertation's
initial-population prefactor. The implemented constant attempt frequency is an
explicit modelling choice; an Eyring temperature-dependent prefactor is not
inferred. The oxide-water MD and H₂ accelerator fits in those sources do not
calibrate H₂O release from the dissertation's particular SS304 surface.

## Domain and ownership

Every `TemkinReservoir` requires the ordinary rate/energy fields, a nonnegative
`coverage_coefficient`, positive `reference_coverage_m2` and `max_coverage_m2`,
and separate nonempty `kinetic_provenance` and `energy_provenance`. Its declared
coverage domain is the closed interval from zero to `max_coverage_m2`. The
kinetic barrier must remain nonnegative over that entire domain. A population
outside it is rejected. There is no clipping, invented inventory, binding heat
or H₂O parameter default.

Each reservoir uses its **own** instantaneous coverage in its barrier. If two
gas states share one physical adsorbate layer, splitting that layer into
separate reservoirs changes this law unless the physical normalization and
shared coverage dependence are handled separately. These records do not
allocate rovibrational populations, infer a molecular mechanism or establish a
common site-occupancy model. Exact emitted state IDs and internal energies must
come from the molecular owner without duplicating the removable inventory.

`TemkinReleaseTrial` and `TemkinReleaseHistory` are distinct record types. Their
files declare a separate schema, kinetic-law identity and interval model. The
ordinary first-order history loader and packet sampler reject them. A truncated
exponential conditional event-time distribution is generally incorrect for this
coverage-dependent source. This version supplies integrated number/energy
moments and **no packet-time adapter**. Do not relabel its file as an ordinary
first-order release or pass its moments through a fabricated baseline history.

For a frozen temperature, the numerical kernel solves a separated, monotone
integral in log depletion \(s=\log(N_a/N_b)\):

\[
\int_0^s \exp\{c[1-\exp(-u)]\}\,du
 = \nu\Delta t\exp\{-E_d[1-\alpha N_a/N_{\rm ref}]/(k_BT)\},
\qquad c=\frac{E_d}{k_BT}\frac{\alpha N_a}{N_{\rm ref}}.
\]

Endpoint-scaled quadrature and logarithmic bracketing avoid subtracting nearly
equal exponential integrals. Exact ratios of the supplied float64 inputs
preserve a small positive barrier near the declared coverage limit. This
version qualifies nonlinear coupling **c ≤ 512**; greater coupling raises
`UnsupportedTemkinRangeError`. That bound is a numerical qualification, not a
physical cutoff, and there is no fallback to a frozen initial rate. Zero
coefficient or zero kinetic energy uses the analytic constant-barrier limit.
Bounds on the same integral handle tiny depletion and an unrepresentable
survivor. Representable small releases and explicit compensated inventory
residuals remain available even below the remaining inventory's spacing.

Each interval retains an inversion regime, iteration count and logarithmic time
residual. The latter is checked against a 2×10⁻¹¹ inversion tolerance in the
nonlinear-root regime. Zero outside that regime means no nonlinear inversion
was required; it is not an experimental uncertainty or proof of zero error.

## Energy and source controls

`activation_energy_ev` and `coverage_coefficient` affect the rate. They do not
specify thermodynamic heat. `desorption_energy_ev` remains a separately supplied
constant cost from the adsorbate to a defined gas reference; `internal_energy_ev`
is the emitted state's energy above that same reference. The prescribed
thermally accommodating flux-Maxwellian boundary contributes mean translation
\(2k_BT\) per released particle. Those three contributions are counted once.
The implementation reports a wall debit and does not apply it to FLASH.

This version does not implement a coverage-dependent surface internal-energy
function. Equating a fitted kinetic barrier with a differential binding cost
would require separate thermodynamic evidence and an integral over the finite
coverage change. The audit describes that additional hypothesis explicitly.

Matched uniform-width moments redistribute the local source after nonlinear
release. They preserve each time/reservoir width-integrated number and energy
moment. They represent the publication's spatial-source control. A single
Maxwellian at an averaged temperature does not preserve the original emitted
velocity mixture. A separate optional **constant-barrier kinetic sensitivity**
sets the coverage coefficient to zero while retaining initial inventory, drive
and all supplied wall-energy costs; it answers a different question from the
uniform-width control.

## Running the independent workflow

From the repository root, with the project's Python environment:

```sh
.venv/bin/python tools/run_temkin_release.py --demo \
  --constant-barrier-control --plot --output-dir runs/temkin_example
```

The example uses a manufactured hot-edge temperature history, one generic
emitted state `M:manufactured-ground` and manufactured parameters in
[`temkin_manufactured.json`](../configs/temkin_manufactured.json). It is not a
water fit or an execution of FLASH. Existing output paths are refused.

For a separately prepared canonical surface history:

```sh
.venv/bin/python tools/run_temkin_release.py \
  --surface path/to/surface.h5 --config path/to/explicit_parameters.json \
  --output-dir runs/temkin_external
```

Positive-depth material samples require an explicit `--boundary-proxy-reference`;
whole-current-averaged temperatures additionally require
`--current-uniformity-reference`. These references preserve a stated modelling
assumption; they do not validate it. Resolving temperature in the current
direction requires a separate planar extension. Do not average those
temperatures before applying this nonlinear law.

The output directory contains the canonical surface input, a distinct Temkin
release HDF5 file, extensive interval source CSV, and a JSON summary with input,
implementation and output hashes. The summary retains manufactured, sourced or
unvalidated parameter labels and independent kinetic/energy provenance. The
optional constant-barrier history has its own baseline schema. Source quantities
are per area; width integration produces particles/m and J/m, and the declared
current-direction length is applied exactly once in the CSV and summary.

Input temperature is interpolated linearly and frozen at each subinterval's
midpoint. Numerical qualification of the frozen-temperature inversion is
separate from convergence of that temperature approximation. Refine integration
subdivisions, original snapshot cadence and surface partition independently.
Stationary geometry, local thermal accommodation, independent pools and a
constant supplied wall-energy cost remain physical assumptions to test.
