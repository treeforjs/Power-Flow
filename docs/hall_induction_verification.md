# Periodic Hall-induction verification

This reference checks the Hall term's sign, dispersion, lack of dissipation and
interaction with scalar resistive damping. It is an exact solution of the
**induction subsystem with prescribed stationary ions**. It generates reference
data for a future configured FLASH verification case. No FLASH execution,
full Hall-MHD verification or physical material calibration is claimed.

## Defined equations and solution

All coefficients are supplied. Set constant permeability μ, scalar electrical
resistivity ρ, electron number density nₑ, positive elementary charge e, and
axial field B₀. Electron pressure and inertia are omitted. The equations are

\[
\mathbf J=\nabla\times\mathbf B/\mu,\qquad
\mathbf E=\rho\mathbf J+\frac{\mathbf J\times\mathbf B}{n_e e},\qquad
\partial_t\mathbf B=-\nabla\times\mathbf E.
\]

On a periodic z interval, choose signed k = helicity × 2π/λ and initial phase φ.
With γ = ρk²/μ and ω = B₀k²/(μnₑe), an exact circular mode is

\[
\mathbf B=(A e^{-\gamma t}\cos\theta,
           A e^{-\gamma t}\sin\theta,B_0),\qquad
\theta=kz+\omega t+\phi.
\]

The transverse current is −k times the transverse magnetic field divided by μ.
The electron drift relative to the prescribed ions is −J/(nₑe).
Reversing B₀ reverses the temporal Hall rotation. Reversing helicity reverses
the spatial winding and current sign; ω depends on k². Frequency ω is angular
frequency in rad/s. At zero resistivity the mode rotates without decay or heat;
at B₀ = 0 it decays without Hall rotation.

For free ions, J × B₀ would drive motion. Holding the ions stationary is
essential to this reference's scope. It cannot establish the gap's ion
magnetization, distribution functions or momentum response.

## Energy and finite-volume comparisons

The Hall term has J · Eₕ = 0. Joule power is spatially uniform for this circular
mode, ρk²A² exp(−2γt)/μ². The perturbation energy and deposited heat per volume are

\[
U_b=\frac{A^2}{2\mu}e^{-2\gamma t},\qquad
Q=\frac{A^2}{2\mu}[1-e^{-2\gamma t}].
\]

Their sum equals the initial perturbation energy. The constant background
energy B₀²/(2μ) is recorded separately. The periodic net Poynting flux is zero;
there is no ion fluid work in the prescribed-ion subsystem. Conversion from
deposited heat to temperature requires a separately qualified thermal model.

For a cell of width Δz centered at zc, multiply each linear transverse field's
point value at zc by sinc(kΔz/2), where sinc(x) = sin(x)/x. The implementation
provides exact cell means of B, J, both electric-field terms and electron drift.
It averages Joule power independently. Squaring the cell-mean current would
underestimate heating and gives zero across a complete period despite positive
local heating. Likewise, magnetic energy must be compared as a separately
averaged physical quantity, not reconstructed by squaring the cell-mean B.

The independent tests check Ampère's law, Faraday signs for both helicities and
background signs, Ohm's law, cell quadrature, periodic energy accounting, ideal
and resistive limits and representable decay tails. These checks qualify the
reference mathematics. A future numerical solver must reproduce its field and
energy evolution under mesh/time refinement before its Hall branch is qualified.

## Export and future FLASH setup

```sh
.venv/bin/python tools/run_hall_induction_benchmark.py \
  --plot --output-dir runs/hall_induction_reference
```

The explicit configuration is
[`hall_induction_manufactured.json`](../configs/hall_induction_manufactured.json).
Its values are manufactured numerical inputs. Each snapshot contains a CSV of
cell means and a sidecar declaring SI units, time, coefficients, representation
and CSV hash. The summary pins inputs, implementation and output hashes and
lists energy budgets. Existing output paths are refused.

When licensed FLASH source is installed, inspect its actual induction/Hall
options and normalization first. The reference assumes electron number density,
positive elementary charge and a resistivity in Ω m; a magnetic diffusivity
coefficient instead equals ρ/μ. Set periodic boundaries along the axial
coordinate, uniform density/coefficient/background fields, and the stated
cell-averaged helical perturbation. Suppress the fluid response through a
documented induction-only verification configuration, if available. If the host
cannot provide that subsystem, use its complete Hall-MHD eigenmode with the
fluid perturbations and dispersion relation instead of claiming this mode
solves the complete coupled equations. Record component staggering and compare
the quantities represented by the host to the corresponding analytic averages.

Time and spatial resolution must cover both the Hall rotation and resistive
decay. A successful zero-Hall diffusion test alone does not qualify this term.
This independent workspace reference does not request any Quark build or runner.
