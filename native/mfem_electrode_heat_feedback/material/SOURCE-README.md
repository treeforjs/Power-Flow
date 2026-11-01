# Conditional nominal SS304 thermal properties

Isolated thermal/material component for the evolving MFEM adapter. The source
and checks use only C++17/Python standard libraries. No MFEM, GPU, PIC, dependency,
shared checkout or production changes occur here. This is an implementation of
an explicitly selected literature scenario, not a foil calibration or a complete
EOS. Mechanical strain, density evolution, geometry and remapping belong to the
mechanics component. Its 773.15 K mechanical-data ceiling remains in force.

## Explicit choices required

```cpp
#include "include/ss304_thermal.hpp"
ss304::Model material(
  ss304::Qualification::nominal_synthesis_1800,
  ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,
  8000.0); // example literature-model reference density; not measured foil density
auto state = material.evaluate(1000.0);
double T = material.temperature_from_u_model_k(state.u_model_j_kg);
```

`evaluated_solid_1500` instead caps **all** queries and inverses at 1500 K.
`nominal_synthesis_1800` explicitly admits the source's synthesized/extrapolated
high-temperature recommendations and the phase assumptions below. Neither mode
is selected by default. Both reject nonfinite, below-300 K and above-mode-range
temperature or energy. Adjacent-float excursions are rejected, never clamped.
Domain endpoints are inclusive. There is no vapor/superheated continuation.

`ss304_c_api.cpp` exposes scalar evaluation, specific/reference-volume energy
inversion and signed heat-deposition transactions. Errors return nonzero, retain
the caller's output and expose `ss304_error()`; they never return a clipped state.
The header can be included directly by the parent MFEM component.
Heat transactions also return the signed numerical energy residual
`m*(u_new-u_old)-Q`. A deposit below temperature resolution may leave T unchanged;
its energy is retained in that explicit residual for the parent energy ledger.
The component never debits an inventory or commits a coupled trial itself.

## Caloric/reference/volume contract

- Temperature is kelvin; specific energies J/kg; conductivities W/(m K) and S/m.
- `h(300 K)=u_model(300 K)=0`. The solid enthalpy is the analytic integral of
  Mills Eq.10, not the rounded printed Eq.11 integration constant. The published
  reference-rounding discrepancy is retained in the copied provenance.
- Density is caller-supplied **mass per fixed reference volume**, invariant with
  T and phase. The example 8000 kg/m3 comes from Sjostrom's nominal 304L EOS
  reference parameter, explicitly used here as a scenario approximation. It
  does not imply 304L measurements establish this foil's density or an EOS.
- In this incompressible, fixed-geometry thermal approximation, the published
  Cp is used as the caloric heat-capacity proxy. `u_model=h` for energy increments
  at the common reference pressure; constant `p_ref/rho_ref` cancels on re-zeroing.
  The reference pressure is 101325 Pa. This is **not** a recovered physical Cv
  or `u(rho,T)` EOS. No pressure input or p(rho,T) function is supplied.
- The model includes no mechanical work because reference volume is fixed by
  this declared approximation. A moving/deforming invocation must be rejected
  by the parent unless a separate work/geometry closure is admitted; zero work
  is not permission to omit that term in moving mechanics.
- The heat adapter can use `rho_ref*u_model` J/m3 of reference volume and its
  derivative `rho_ref*C_app`. Equivalently, a cell with accepted mass `m` stores
  `m*u_model`. Depositing actual J uses `u_new=u_old+Q/m`; never divide a
  deposited energy by thickness/height twice. Bounds failures abort the trial.

## Thermal branches

Mills et al. (2004), DOI [10.2355/isijinternational.44.1661](https://doi.org/10.2355/isijinternational.44.1661),
pp.1663-1665, supplies evaluated austenitic heat capacity (Eq.10), its enthalpy
integral (Eq.11), liquid Cp=800 J/(kg K) (Eq.12) and fusion-entropy estimate
160 J/(kg K) (Eq.13). It is an evaluated correlation source using 304/316 data,
not new characterized-foil measurements. Retain distinct source uncertainty
statements, including liquid Cp +/-50 J/(kg K).

Ho and Chu (1977), CINDAS45/ADA129160, Table2 printed p.15/PDF22, supplies
estimated melting-range markers **1672 and 1727 K**. We explicitly adopt these
as nominal scenario endpoints `Ts,Tl`. They are not actual-foil phase boundaries.
The Mills generic Table3 transition convention 1723 K is not silently substituted
for this scenario. Sjostrom's liquid-model parameter 1124.3 K is not used as a
melting point. No Kim304L liquid-table conflict is ingested.

The scenario has three branches:

1. Solid, 300<=T<Ts: Mills Cp and its analytic re-zeroed integral.
2. Mushy, Ts<=T<Tl: `f=(T-Ts)/(Tl-Ts)`, sensible Cp interpolates linearly from
   `Cp_s(Ts)` to 800; latent heat `L=160*Tl=276320 J/kg` is added once as `L*f`.
   `C_app= Cp_sensible + L/(Tl-Ts)` is the exact derivative of this chosen h.
3. Liquid, Tl<=T<=1800: Cp=800 and a continuous h anchored at the mushy upper
   endpoint. No further latent contribution is added.

The linear phase fraction, linear sensible-Cp bridge and use of liquidus in L
are explicit **model approximations**. They are not phase-equilibrium data or
source-measured mushy Cp. They produce a continuous, strictly increasing energy
and an unambiguous inverse; the derivative is right-sided at Ts/Tl and jumps
there as the latent distribution begins/ends. Reversible equilibrium fraction
is assumed: no supercooling, superheating, hysteresis or phase kinetics.

## Electrical/thermal transport

The 18-row recovery in `data/ho_chu_304_recovered.csv` is visually verified from
the original Table2, 300-1800 K subset; no graph digitization is used. Resistivity
and k are piecewise linear in T. Sigma is the reciprocal of **interpolated
resistivity**, not an independently interpolated sigma. The table's endpoint
values and the bridge over the estimated melting interval are preserved.

The source explicitly **pools 304 and 304L** measurements for its AISI304
recommendation and selects a quench-annealed condition (printed pp.2,13). This
component uses that identified pooled recommendation as a nominal scenario;
it does not call it measured grade304-only data or a cold-worked foil curve.
Values above1500 K are explicitly extrapolated/synthesized by the source.
Source uncertainty is of order +/-5%, and +/-15% at/above melting, without a
coverage factor. These statements are not independent Gaussian error bars.
No new extrapolation is performed beyond the published rows. No solid transport
polynomial is extended into liquid. Liquid transport is source-estimated.

Original reports acquired with TLS verification and matched to the existing
manifest hashes: Mills `bf120d...74836`; Ho-Chu `673829...6209a` via the
[DTIC scan mirror](https://archive.org/details/DTIC_ADA129160). Exact hashes,
acquisition failures and verified pages are in the external report directory.
Raw PDFs/text/images are not included in the source package or Git repository.

## Verification and integration

Run `python3 run_checks.py /absolute/report/path`. The runner acquires the shared
CPU7 lease, pins CPU7/nice15 and one-thread environment, builds the C++ shared
library without fast-math, then verifies source anchors, Decimal caloric values,
bounded inverses, phase continuity/latent energy, derivatives, signed deposition
and domain/error behavior. It checks shared input hashes and production hold
before and after. See the result JSON for exact counts and numeric residuals.

Parent integration must bind `data/scenario.json`, source hashes and the accepted
reference mass/volume in its restart revision. Crossing 1500 K requires the
explicit nominal synthesis mode; crossing1800 K must stop. The broad thermal
scenario does not lift the mechanical773.15 K limit. Melt flow, mass loss,
vaporization enthalpy/pressure, plasma conductivity, two-temperature physics,
surface desorption/internal energies and moving-geometry work remain absent.
There is no deployment or coupled-breakdown validation claim in this packet.
