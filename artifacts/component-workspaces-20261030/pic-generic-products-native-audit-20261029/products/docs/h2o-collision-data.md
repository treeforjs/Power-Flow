# Frozen gas-phase electron-water data

This release binds numerical water-vapour data to strict incident-energy support.
It is not a complete fragmentation mechanism and does not authorize a coupled run.
The frozen dataset is `data/collisions/h2o-gas-20261029.json`; the adjacent lock
records its SHA256. Original papers, full text, page renders, receipts and failures
remain outside the repository in the assigned `outputs/h2o-collision-data-20261029`.
The repository redistributes factual numerical transcriptions and code, not PDFs.

`tools/build_h2o_collision_data.py` deterministically reproduces 30 observables,
676 table points and the four-orbital NIST BEB formula. Sources and primary lineage
are retained per channel. Tables include elastic integral and momentum transfer,
ground-rotational excitation, vibration, molecular electronic excitation, neutral
OH yield, parent and fragment-ion yields, attachment, and gas-phase aggregate
inelastic/stopping data. Molecular H2O(A1B1) is not the OH(A) radical state.

Muñoz et al., PRA76,052707(2007), DOI10.1103/PhysRevA.76.052707, provides
gas-phase calculated elastic/inelastic integrals and recommended quantities to
10keV. Its measured ion-current data end at5keV. Integral elastic is not momentum
transfer; the integral inelastic and stopping moments overlap explicit excitation
and ionization. They cannot be added as independent collision events. The printed
200eV rotational2.39a0² value disagrees with the printed total; it is preserved
and flagged. The750eV stopping row implies36.46eV mean loss, outside the separate
40±2eV statement but within the authors'10% high-energy constancy description.

Song et al., JPCRD50,023103(2021), DOI10.1063/5.0035315, supplies recommended
tables with the original measurement/theory lineage. Positive fragment partials
are inclusive ion yields, with D2O proxy corrections where indicated. Atomic
O2+ means one O nucleus with charge+2, not molecular O2(+). Blank cells are not
zeros. Parent H2O+ is distinct from total ionization. Table11 O(1S) remains
quarantined for its factor100 unit/text discrepancy.

The two stretch channels use the review's explicit symmetric/asymmetric0.56/0.44
recommendation from its theory ratio at3eV. This is a declared constant split,
not independent measured branching at every energy. Main text/caption identify
the symmetric share as56%; the figure legend has inconsistent wording. The split
channels begin at0.53eV and preserve their unresolved onset gaps.

NIST's four-valence-orbital BEB model has eight electrons and no core/double-ion
term. Its analytic support is capped at10keV; NIST's generic5–20% accuracy claim
is only to about1keV. It supplies neither stable-parent branching nor a water
secondary spectrum. The orbital constants/formula were inspected on official
NIST pages; direct archival HTTP/CGI requests returned403/timeouts. Independent
60-digit formula checks are computational checks, not claimed CGI observations.

## Strict query API

`raikou/collision_data.h` is header-only and uses the actual `CollisionProcess`
type. `StrictCollisionDataset::read(path)` validates dimensions, phase, grids,
units, source bindings and support. `observable_value_SI(id,E)` returns the
declared quantity: usually m²; the stopping moment is eV m². It throws above or
below supported data, including threshold-to-first-point gaps. Only an explicitly
bound physical threshold permits a subthreshold zero. No tail holding, zero-fill,
cutoff, or covered-population renormalization is implemented.

`coverage(id,population)` reports covered/uncovered physical weight and a known-
support sigma*v contribution divided by ALL physical electron weight.
`strict_observable_rate_m3_s` throws if any positive-weight sample is uncovered.
Rates use target-rest energy and relativistic electron speed for a stationary
target. Moving targets require host-owned relative-frame treatment. A partial-yield
rate is still a yield rate, not an exclusive event probability.

`kernel_status(id,model)` exposes status and reason. The default is unselected.
`require_MCC_event_kernel(id,model,E)` can admit the explicitly selected conditional
parent model inside13.5–1000eV; it rejects missing kernels or unsupported energy.

`h2o_conditional_parent.h` implements `parent_ground_equal_pair_CM_v1`: a ground
intact H2O+ plus equal-energy, back-to-back electrons in the CM frame, with a
host-supplied isotropic pair axis and exact relativistic boost. It consumes the
parent-only rate, not the total. This is an intentionally simple conditional
differential model; it is not measured/calibrated water energy sharing or angles.
Ion mass is Mneutral−me+12.621eV/c². Ionization reference energy is12.621eV;
do not count that mass/reference energy again as gas/electrode heat. The pure
function returns momenta, kinetic energies and matched weights and mutates no
inventories. The host must implement atomic debit/staging/commit. Full measured
secondary, internal-state and fragment-recoil distributions remain absent.

`reduced_transport_processes(0.001,100)` exports only EMT plus one-way bend to
the native Raikou type, with an explicitly isotropic transport approximation.
It excludes rotation, other excitation and reactions; it is not a water swarm or
breakdown model. Native tables still have endpoint defaults, so the host MUST
admit the full energy domain before use. Do not evolve Raikou F0 alongside the
same PIC electrons. PIC remains the sole electron owner in the intended driver.

## Validation and integration

Run `bash tools/test_h2o_collision_data.sh` under the assigned CPU7 policy. It
holds `testing-core7.lock`, uses nice15/one thread, builds standard C++ executables,
checks independent source anchors/units/exponents, thresholds/gaps, multiplicity
eligibility, weighted coverage, high-energy values, formula precision, and
conditional event on-shell energy/momentum/weight closure. It performs no GPU,
MPI, EEDF, chemistry, trajectory or field evolution. Existing failures stay in
their own output attempts. No dependencies are added or upgraded.

The data port is `data/collisions/h2o-collision-port-v1.json`. Before a transaction,
bind dataset SHA/revision, state/mesh/material revision, accepted[t,t+dt], parent/
channel/event IDs, target-rest incident energy, represented weights and actual
3D volume in SI. The observation slab is not another reacting population. The
host RNG, parent debit, all species products, signed current, reference/internal/
kinetic energy and restart/rollback must commit together. These are interface
requirements, not a claim that a coupled driver exists.

The retained8.084keV electron tail is covered by several aggregate/theory
observables, but not by the full momentum-transfer/fragmentation/attachment set.
No scalar total-to-parent substitution closes that hole. Full water breakdown
and heavy/product reaction completion remain blocked. Supplied2023 hydrogen
kinetics is a method benchmark, not water coefficients. No instantaneous
H2O→2H+ +O+ +3e closure is introduced.

## Conditional attachment extension v2

The separately frozen `h2o-gas-20261029-events.json` preserves every v1 table
and adds `Hminus_Ram65_factorized_window_v1`. Its adjacent lock binds the bytes;
v1 is unchanged. `h2o_dea_kernel.h` exposes:

```cpp
conditional_hminus_Ram65_products(data, model_name, target_rest_energy_eV,
    physical_weight, neutral_mass_kg, Hminus_mass_kg, incident_unit_axis,
    u_vibration, u_lab_cosine, u_azimuth);
```

This consumes the incident electron and one H2O parent and returns one Hminus
and one OH(X,v), with equal represented weights. The admitted conditional
window is **6.25–6.75 eV**. The Hminus rate is the frozen Rawat/Song attachment
table. Ram, Prabhudesai and Krishnakumar, arXiv:1007.5144, measured the angular
and fragment-energy distributions at6.5eV with about0.5eV beam resolution.
Their Eq4 and fitted d/p amplitude0.31, phase3.139rad give the angular model.
Sample its normalized intensity in dcos(theta), with uniform azimuth. The fit
is a laboratory Hminus direction: the implementation inverts the exact CM
boost for each selected state, rather than interpreting the fit as a CM angle.

TableIII's present-experiment OH vibrational intensities for v0..4 are
1,0.82,0.57,0.32,0.25. NIST WebBook's OH(X) row/note22, tracing to Maillard,
Chauville and Mantz1976, supplies the vibrational term polynomial through v5.
Only v0..4 is used. The angular fit is integrated over fragment energy; the
model factorizes angle and vibrational state despite the source's observed
energy-dependent asymmetry. Holding the6.5eV fit within the selected window,
closing the fitted v0..4 population, and assigning no additional rotation are
explicit conditional assumptions. No measured energy dependence, unresolved
high-v tail bound, instrumental Gaussian broadening, or rotational distribution
is claimed. Do not extend this model to the8.5/12eV resonances.

The reference cost is4.35eV. For selected internal energy E_v the returned OH
state mass is `Mneutral + me + (4.35eV+E_v)/c² - MHminus`. The host must bind
this same immutable species/state mass convention, or reject the candidate;
silently substituting inconsistent nominal masses breaks energy closure.
Returned momenta and kinetic energies are in the stationary-target laboratory
frame. The parent driver owns moving-target/lab transforms, RNG variates,
availability, all product staging, current deposition and atomic commit.
Reference/internal energy is not an additional Joule or electrode heat source.

`bash tools/test_h2o_collision_events.sh` replays existing frozen data without
regenerating it. Compiles and native checks use `arch-base-20261028/run-frozen`
under CPU7/nice15/one-thread/testing-core7.lock. The first host-userspace receipt
is preserved separately. Tests verify both v1 and v2 loading/parent regression,
source constants, source angular/state measures, exact product mass shells,
energy/momentum/weights, and out-of-range/total-as-branch rejection. They do not
advance physical time or evolve electrons, chemistry, fields or trajectories.
The integration contract is `h2o-collision-port-v2.json`.

## Conditional OHminus endpoint v3

`h2o-gas-20261029-attachment.json` separately freezes the v3 addition without
rewriting v1/v2. `h2o_two_body_attachment.h` exposes
`conditional_ohminus_ground_products(data, model, E_eV, W, Mneutral, MH,
beam_unit, CM_anion_unit)`. Select `OHminus_ground_isotropic_CM_v1` explicitly.
It admits only4.5–12.5eV and returns one ground OHminus and one H(1s), both
weight W; the parent and incident electron are consumed. The host supplies
an isotropic CM axis. The state mass is
`M_OHminus=Mneutral+me+3.27eV/c²-MH`, with the same host mass-basis requirement
as v2. Exact on-shell two-body recoil and boost close momentum and kinetic plus
reference energy. No additional internal energy is assigned.

Melton1972, DOI10.1063/1.1678051, Fig2/TableIII and the pressure-control discussion
identify primary OHminus+H production; Ram2010 TableI supplies the3.27eV limit.
The rate remains the Song/Rawat-renormalized Melton curve, with its original
absence of a reliable OHminus uncertainty bound. The assumed ground product
state and isotropic angle are **a maximum-translation conditional endpoint**,
not measured differential data. This weak, historically disputed channel must
not silently become a calibrated loss spectrum. No zero is inferred above12.5eV.

Replay the v3 release using `bash tools/test_h2o_collision_events.sh --dataset
data/collisions/h2o-gas-20261029-attachment.json`. This adds411 mass-shell,
energy/momentum/weight and eligibility checks to the pinned v2/v1 regression.
All products remain staged candidates; this code advances no electron or gas
inventory. The contract is `h2o-collision-port-v3.json`.

## Conditional prompt neutral dissociation v4

The separately frozen `h2o-gas-20261029-dissociation.json` adds
`A1B1_fixed_loss_dipole_ground_v1` on **9–977.7 eV**. It uses only the molecular
`H2O:electronic_A1B1` excitation curve. Hargreaves et al.2012,
DOI10.1088/0953-4075/45/20/201001, identify the1,3B1 states as impulsively
dissociative to OH+H and report a broad loss profile and forward-peaked singlet
electron scattering. The conditional model replaces that continuum by the
source's7.49eV line, samples a Born/dipole-shaped CM electron angular law
proportional to1/q², and assigns an isotropic, ground-product breakup. These
energy/angle/internal-state assumptions are approximate, particularly near the
low-energy end, and are not measured joint differential data. The integral
curve remains unchanged, including its BEf theory lineage above17eV.

`h2o_neutral_dissociation.h` exposes:

```cpp
conditional_A1B1_dissociation_products(data, model, target_rest_energy_eV,
    W, Mneutral, MH, incident_unit_axis, u_electron_cosine, u_electron_azimuth,
    excited_rest_H_unit_axis);
```

The last axis is sampled isotropically by the host. Exact two-body excitation
creates a virtual excited water of mass Mneutral+7.49eV/c², then exact two-body
breakup and Lorentz boosts return one scattered electron, one H(1s), and one
OH(X,v0), all weight W. The parent is consumed; net electron gain is zero.
The reference cost is5.10eV; the remaining2.39eV of excitation becomes product
translation in this ground-state endpoint. **The7.49eV intermediate loss is not
another heat deposit.** Do not also sample the same molecular excitation, the
inclusive OH yield, or an overlapping aggregate inelastic/stopping event.

`h2o_common_mass_basis.h` binds all four conditional models to one reduced-model
state-mass table. Given Mwater and MH, it uses the rounded5.10/4.35/3.27eV source
references, implying H/OH electron affinities0.75/1.83eV. Those implied affinities
are **not precision atomic/molecular mass measurements**. This prevents the
DEA and neutral-dissociation kernels from assigning incompatible OH(v0) masses.
The host must explicitly use this conditional reference basis or provide a new
coherent model revision; it must not silently mix an independent species table.
State masses agree across the kernels within floating-point mass rounding.

Pinned CPU replay adds839 three-product conservation, on-shell, source-rate,
overlap-rejection and cross-kernel mass-basis checks. Run the event test script
with `--dataset data/collisions/h2o-gas-20261029-dissociation.json`. Direct
archival downloads of the Hargreaves primary PDF timed out; web-read locations,
factual binding notes and failed receipts are preserved externally. No measured
joint loss/product spectrum or complete8keV mechanism is claimed.

## High-energy outer-three-orbital conditional subset v5

`h2o-gas-20261029-outer3.json` adds weighted versions of the original NIST BEB
orbital formula over **1–10 keV**. No original data bytes or orbital constants
are changed. Murakami et al., PRA85,052704(2012), Eqs27–30, reproduce the
high-energy single-vacancy branching of Tan et al., ChemPhys29,299(1978):
1b1 and3a1 feed parent H2O+;1b2 feeds parent/OH+/H+ with weights0.08/0.70/0.22.
Only these three outer orbitals enter the mapped subset. The proton collision
dynamics from that paper are not used as electron-impact rates. Using NIST BEB
orbital cross sections as single-vacancy rates and holding those asymptotic
branches over1–10keV are explicit conditional approximations.

`h2o_outer_valence.h` exposes `conditional_outer3_ionization_rates(data,E_eV)`.
It returns parent, OHplus+H and Hplus+OH cross sections, the **positive unmapped
2a1 cross section**, and the original four-orbital sum. The mapped rates plus
unmapped2a1 equal the original BEB sum without redistributing the missing rate.
The oxygen1a1 core and multiple-ionization/Auger processes remain outside this
subset. In particular, this is not a complete8keV mechanism.

`conditional_outer3_ionization_products(data, model, E_eV, W, Mwater, MH,
beam_unit, CM_e_pair_unit, virtual_ion_rest_positive_fragment_unit)` selects
one of `outer3_BEB_Tan_parent_v1`, `outer3_BEB_Tan_OHp_H_v1`, or
`outer3_BEB_Tan_Hp_OH_v1`. Every event replaces the incident electron by two
electrons and debits one water parent, all with matched weight W. The parent
branch returns ground H2O+ with the same equal-energy, back-to-back CM electron
pair approximation as v1. Each fragment branch uses a19.83eV virtual ion,
an equal-energy CM electron pair, and an isotropic ground-product breakup.
The final chemical references are5.10+13.017eV for OHplus+H and
5.10+13.59844eV for Hplus+OH; the ionization energies are official NIST WebBook
values. All momenta are boosted exactly to the stationary-target laboratory.

These are conservative kinematic comparison models, not calibrated secondary
electron spectra, hole-band profiles or heavy-product recoil distributions.
The19.83eV value is the original BEB binding parameter; it is not a measured
joint fragment-energy spectrum. The parent branch assigns ground ions directly
and does not preserve an orbital-tagged excited-ion population.

Choose this high-energy subset instead of overlapping full-BEB/parent/yield
events. It cannot be added to them. A rate for the unmapped2a1 hole is not an
exclusive O+/H+ branch and has no product kernel. The existing1keV parent model
and this alternative overlap at the boundary; a host dispatcher must select
one prescription and report any boundary mismatch. All returned candidates
remain marked incomplete for the full water mechanism.

Kernel kinetic energies are per physical product and reference/internal gains
are per physical parent event. The host scales them by W once; reference energy
is not repeated for each daughter and is not another electrode/Joule heat term.
