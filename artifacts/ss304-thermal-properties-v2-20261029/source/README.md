# SS304 source-conditioned thermal/vapor alternative, v2

This additive release preserves the sealed1800K release. Its **nominal1800**
mode supports298.2-1800K and reproduces all old property values bitwise on
300-1800K. A separate explicit **source_conditioned_liquid_3000** mode supports
298.2-3000K condensed caloric/transport and1800-3000K neutral-vapor alternatives.
This is a conditional material model and transaction interface, not an acquired
wide-range steel EOS or a completed coupled60ns run. All callbacks reject
nonfinite/out-of-domain states. Existing mechanical773.15K support is unchanged.

## Native use and low-temperature correction

```cpp
#include "include/ss304_thermal_v2.hpp"
ss304v2::Model material(ss304v2::Qualification::nominal_1800,
  ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000);
auto properties = material.evaluate(299.0);
double T = material.temperature_from_u(properties.u_model_j_kg);
```

Modes0/1/2 mean evaluated_solid_1500, nominal_1800 and the explicit3000K scenario.
The lower limit298.2K is a new reviewed source selection: Mills Table3a
explicitly evaluates298.2K and its Cp/H equations; Ho-Chu Table2 brackets it with
293K (R=71.3e-8 ohm m,k=14.76W/m/K) and300K rows. We evaluate the same Mills
analytic integral and interpolate transport between those rows. There is no
extension below298.2K and no change to the300K energy zero. At298.2K the specific
reference energy is-866.079237344J/kg; negative cold energy is valid. Initial300K
patches can cool physically without an artificial initial-temperature increase.

## Liquid branch above1800K: source and model choices

Kim(1975), ANL-75-55, DOI10.2172/4152287, is an **estimated304L** compilation
that the author also recommends for304; we retain that as a nominal surrogate,
not measured grade304/foil equivalence. PDF12/printed3 estimates liquid
Cp=0.190cal/(g K)=794.96J/(kg K). PDF27/printed18 gives Eq29
`k_Kim(T)=100*(0.1229+3.248e-5*T)` W/(m K), based on an elemental-metal
correlation and an assumed conductivity change on melting. Its range to3000K
is an evaluated/extrapolated recommendation, not a high-T steel measurement.
The scanned exponent is-5; OCR misreads it and is not ingested.
The SI conversion explicitly chooses the thermochemical calorie4.184J; the
report's unqualified calorie is not presented as a verified calorie-definition
measurement. The same convention is applied to latent energy.

Above1800K, `h=h_v1(1800)+794.96*(T-1800)`; this keeps the accepted reference and
phase energy continuous. The small Cp jump800->794.96 is an explicit source
switch at1800K. The endpoint itself retains the old left-sided Cp800.

The selected conditional conductivity uses Kim's **shape** anchored to the
admitted Ho-Chu1800K value: `k=28.99*k_Kim(T)/k_Kim(1800)`. The normalization
factor is1.5984429104. The59.84% difference from Kim's absolute liquid values
is a real model discrepancy and must be explored as a sensitivity; neither
source's uncertainty statement covers it automatically. We do not label the
anchored curve as Kim's original absolute data. At3000K the selected k is35.2201
W/(m K), versus Kim22.034. No held endpoint or solid transport tail is used.

Liquid electrical conductivity is a **Wiedemann-Franz model extension**:
`rho_e=L_anchor*T/k(T)`, with L_anchor set from both admitted1800K transport
values. Sigma=1/rho_e is continuous there. This is a stated constant-Lorenz
approximation, not a recovered3000K resistivity measurement. The source's
minimum-resistance/annealed-stock qualifications remain relevant.

## Thermodynamically closed neutral-vapor alternative

Kim PDF13/printed4 Eq13 gives `log10(P_atm)=6.1210-18836/T`. PDF12 estimates
vaporization enthalpy1770cal/g=7.40568MJ/kg. For the stated nominal atom fractions
Fe/Cr/Ni=.713/.202/.085, CIAAW atomic weights give M=.0553096362kg/mol. The
ideal-vapor Clausius-Clapeyron slope of Kim's pressure fit implies6.51985MJ/kg,
the stated latent is13.59% larger than this pressure-slope value. We do not silently combine them as
an exact EOS or repair either source value.

The alternative explicitly adopts the latent estimate **at T*=3000K**, the
source pressure at3000K, liquid Cp above, and a congruent monatomic ideal vapor
with R=Ru/M,Cp_g=5R/2. Placement of the source latent at3000K is a model anchor;
the report does not establish a temperature-resolved vaporization latent law.
Set `delta_Cp=Cp_g-Cp_l`, `L(T)=L*+delta_Cp*(T-T*)`, and

```
ln(Psat(T)/Pstar) = delta_Cp/R*ln(T/Tstar)
                    - (L(T)/T-Lstar/Tstar)/R
h_g(T) = h_l(Tstar)+Lstar+Cp_g*(T-Tstar)
u_g(T) = h_g(T)-R*T
u0_atom = h_g(T)-Cp_g*T                 # constant reference energy per kg
```

This enforces Kirchhoff dL/dT=Cp_g-Cp_l and
`R*T^2*dlnPsat/dT=L(T)` exactly within the declared dilute-vapor approximation.
The selected pressure differs from Kim's raw fit away from3000K: at1800K it is
0.8033Pa versus4.5949Pa; at2400K1077.95Pa versus1898.38Pa. Both curves are
reported, without asserting the alternative is the measured one. Condensate
specific-volume corrections to chemical potential are neglected: Psat/rho_ref
is at most8.81J/kg in this range for the nominal8000kg/m3, versus megajoule/kg
phase energies. This is not permission to omit mechanical work in moving MHD.

Neutral-gas thermal conductivity has an optional **parameterized hard-sphere**
alternative, `k_g=75*kB/(64*d^2)*sqrt(kB*T/(pi*m_atom))`. NBS(1981),
DOI10.6028/jres.086.011, Eq13 supplies the generic hard-sphere expression;
none of its argon data/uncertainty is transferred to steel. A positive collision
diameter must be supplied, and ionized_fraction must be exactly0. No fitted
metal-vapor diameter, charged-gas conductivity or plasma closure is invented.
The verification example derives d from the nominal reference atomic volume,
explicitly as a sphere-model assumption. k scales as d^-2. A Fourier coefficient
does not imply the vacuum plume is in the continuum regime: assess mean free
path/Knudsen number or use gas-owner particle transport. Liquid sigma/k must
never be applied to dilute or ionized vapor.

## Mass loss, packets and reference-volume energy

`NeutralVapor::net_mass_flux` implements signed Hertz-Knudsen flux
`alpha*(Psat-Pback)/sqrt(2*pi*R*T)` at a common surface/background temperature.
Alpha and background partial pressure are caller-supplied; alpha1/vacuum is an
upper kinetic scenario, not a measured accommodation coefficient. A negative
result is condensation and is never clipped. Generic number-flux formalism is
verified in the NIST Shakeel/Wei/Pomeroy desorption report, Eq6/PDF7; no cryogenic
gas/water material property is borrowed. A nonequilibrium background needs its
actual incident-packet distribution rather than this same-T formula.

Heat owner owns one canonical **metal** mass/energy inventory per patch/cell.
It is separate from adsorbedH2O inventory. The flux is only a proposal; gas owner
samples/emits or consumes actual weighted packets with caller-owned draws.
`prepare_transfer` takes the actual signed outward mass dm and actual total
packet kinetic energy K. No second mean2RT energy is added. It returns

```
E_packet = dm*u0_atom + K_actual
m_next = m_old-dm
E_next = E_old+Q_other-E_packet
T_next = bounded_inverse(E_next/m_next)
```

No extra latent-vaporization debit is added: it is already in the constant gas
reference u0_atom. Fusion remains only in condensed h. Incoming condensation
uses negative dm and negative K for the SAME packet/reference. Base states are
immutable; stale/replayed event/revision tokens and failed trials leave outputs
untouched. Only the parent preflight/commit can accept state and gas packets
together. Binding/nontranslational excited/ionization energies require separate
gas/reaction-owner state-specific debits under the same energy reference.

After mass transfer, thermal energy/capacity per reference volume are
`(m_current/V_ref)*u(T)` and `(m_current/V_ref)*C_app`, provided by
`evaluate_inventory`. Continuing to use initialrho_ref*u after mass loss would
create energy. T is recovered from E_total/m_current. The intrinsic k/sigma
returned here is NOT an effective mixture law for a partially empty cell.
Geometry owner receives removed reference volume dm/rho_ref (divide by actual
surface area to obtain recession). It owns conservative remap/shape/work, or
must explicitly admit a separate stationary-inventory approximation. Complete
depletion stops for an empty-cell handoff; no mass clamp or ghost solid is made.

## Quantified integration needs

At3000K, constant-T vacuum/alpha1 flux is41.868kg/(m2 s). Over60ns it gives
2.5121e-6kg/m2 and about18.04J/m2 evaporative cooling in this model. This is a
diagnostic extreme within the selected range, not an evolving-run prediction.
For a200um/8000kg/m3 foil, this is1.57ppm of areal mass, but about6.28e-4 of
a0.5um active layer; local inventory and surface refinement matter. The metal
column is about2.74e19atoms/m2, potentially relevant to plasma even at small
bulk mass loss. Each additional1eV/atom of excitation/ionization changes the
specific energy by1.744MJ/kg (about4.38J/m2 for that column). Rates, charge and
electronic excitation are gas/plasma-owner physics, not installed here.

Congruent evaporation is an approximation. Kim Table3 at3000K has pure-element
pressures Fe=.546,Cr=1.30,Ni=.414atm. Raoult weighting gives vapor atom fractions
about.567/.382/.051, different from bulk. A real selective Fe/Cr/Ni source needs
species inventories/chemical potentials and rates. Do not insert a fictitious
SS304 molecule into the H2O reaction network. An explicitly open escape reservoir
can balance removed mass/energy as a bounded alternative without claiming a
modeled reacting metal plume.

At3000K condensed u is2.153727855MJ/kg. `assess_energy(m,E)` returns signed
cold/hot energy margins and a range status, so the thermal owner can report
the exact unsupported trial energy `E-m*u_max` if a new solve exceeds support.
No excess heat is discarded. First reevaluate sigma(T)/thermal feedback under
the new model before interpreting the old constant-property corner trajectory.

The acquired Sjostrom2023 manuscript still supplies only published4273 model
parameters: its tabulated Tartarus electronic free energy and full numerical
coexistence/EOS grid are missing, and Ref27 directs EOS requests to the SESAME
contact. No request was sent. Mattsson/Desjarlais2004 likewise depends on4272
and tabulated temperature-dependent fits; it explicitly requires thermodynamic
integration. Their plots/parameters are not reconstructed into a wide-range
steel EOS. Above3000K, dense expanding fluid, critical behavior, ionization and
compressible p(rho,T) remain outside this concrete alternative.

`run_checks.py` uses oneCPU7/nice15 under the existing lease, frozen C++17/Python
standard libraries, noGPU. Verification checks sealed-v1 identity, cold-energy
reference, source-conditioned branches/inverses, phase identities, signed flux,
transaction conservation/rollback/replay and explicit unsupported-energy budgets.
No MFEM solve, production action or60ns coupled validation occurs here.
