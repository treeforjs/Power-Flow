# Water parent-ionization channel: coverage scope v2

Successor of coverage scope v1 (`pic-water-parent-coverage-v1-20261030`, sealed
copy in `Power-Flow/artifacts/parent-coverage-v1-20261030`). Audit links: A13
(versioned supported channels, no invented rates, no renormalized branches) and
A14 (quantify model seams, do not smooth them).

Decision recorded 2026-10-30: electron energies above 1 keV are covered by a
dataset-backed band (1-10 keV) and by a **provisional, model-derived** band
(above 10 keV) until measured parent data replace them. The model-derived band
is tagged as such in the exact ledger and in every event identity, so results
can be reprocessed when measured data become available.

All energies are the electron kinetic energy in the H2O rest frame (sealed
`pic_collisions::relative`).

## 1. Coverage map

| Target-rest energy E | Class | Rate | Kernel | Ledger category |
|---|---|---|---|---|
| E <= 12.621 eV | known zero | 0 | - | - |
| 12.621 < E < 13.5 eV | **unresolved onset gap** (v1, unchanged) | none; bounded by Munoz 2007 inelastic | - | unresolved |
| 13.5 <= E <= 1000 eV | dataset (v1, unchanged) | Song 2021 Table 12 parent H2O+ | `parent_ground_equal_pair_CM_v1` (sealed) | `dataset_song2021_parent_table` |
| 1000 < E <= 10000 eV | dataset BEB | `H2O:conditional_outer3:parent` = `outer3_BEB_Tan_parent_v1` rate | `parent_ground_BEB_SDCS_binary_encounter_v2` (BEB SDCS) | `dataset_outer3_BEB_Tan_parent` |
| 10000 < E <= 1e7 eV | **provisional, model-derived** | RBEB x Tan parent weights | `parent_ground_BEB_SDCS_binary_encounter_v2` (RBEB SDCS) | `provisional_RBEB_Tan_parent_model` |
| E > 1e7 eV | unresolved, no bound | none | - | unresolved (refused unless `admit_unbounded`) |

Edges are closed on the lower band: 1000 eV exactly is Song (v1 behaviour),
10000 eV exactly is outer3 (its closed support end). The provisional upper
support (1e7 eV) is a scope parameter bound into the identity; it lies above
the end of the Rieke-Prepejchal validation range (2.7 MeV) and above the
electron energies of the sub-MV drive.

## 2. Inputs and provenance

**Song band.** Unchanged: dataset `h2o-gas-20261029-v1` (sha256 `840663cd...e4bf`),
channel `H2O:ion_yield:H2O+`, Song et al., J. Phys. Chem. Ref. Data 50, 023103
(2021), Table 12 (Straub et al. 1998 lineage, D2O-based recommendation, 6 %).
The outer3 freeze carries this table byte-identically (checked at construction).

**BEB orbital constants** (present in the dataset; not refetched). Channel
`H2O:nist_BEB_valence` of the v1 dataset and channel
`H2O:conditional_outer3:parent` of `h2o-gas-20261029-outer3-v5` (sha256
`c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551`) carry the
same four valence orbitals (checked equal at construction); dataset source
`nist_beb`: NIST electron-impact ionization database,
https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=H2O (BEB
orbital constants of H2O; primary reference Kim & Rudd, Phys. Rev. A 50, 3954
(1994)):

| orbital | B (eV) | U (eV) | N | Q | parent weight |
|---|---|---|---|---|---|
| 2a1 | 36.88 | 70.71 | 2 | 1 | 0 |
| 1b2 | 19.83 | 48.36 | 2 | 1 | 0.08 |
| 3a1 | 15.57 | 59.52 | 2 | 1 | 1 |
| 1b1 | 12.61 | 61.91 | 2 | 1 | 1 |

a0 = 5.2918e-11 m, R = 13.6057 eV (dataset values). The oxygen 1s core is not
in the NIST orbital set (no automatic zero is implied for it; its parent
branch is not represented).

**Parent branching.** Tan, Brion, van der Leeuw & van der Wiel, Chem. Phys. 29,
299 (1978), high-energy (dipole) single-vacancy fragmentation branching, as
written by Murakami et al., Phys. Rev. A 85, 052704 (2012), Eq. (27)
(retained PDF sha256 `464c824d...edf8`):
sigma(H2O+) = 1.00 sigma(1b1) + 1.00 sigma(3a1) + 0.08 sigma(1b2); 2a1 holes
give H+ and O+ only (Eqs. 29-30), parent weight 0. These are the
`orbital_weights` of `outer3_BEB_Tan_parent_v1`. Definition used here: the
per-orbital parent fractions are energy-independent (high-energy limit) and
multiply each orbital's own cross section (BEB below 10 keV, RBEB above). The
effective parent fraction of the valence total, f(E) = sum w_i sigma_i /
sum sigma_i, is therefore weakly energy dependent: 0.6939 at 1 keV, 0.6834 at
10 keV, 0.6798 at 100 keV, 0.6775 at 1 MeV, tending to sum w_i N_i/B_i /
sum N_i/B_i = 0.6675 as E -> infinity. The Tan partial oscillator strengths
themselves are not in the dataset (only the derived branching); the full Tan
paper was not retrieved by the dataset work. No total-to-parent scaling, no
renormalization of the unmapped 2a1 rate.

## 3. Formulas

**BEB** (Kim & Rudd 1994), per orbital, t = T/B, u = U/B:
sigma = S/(t+u+1) [Q ln t/2 (1 - 1/t^2) + (2-Q)(1 - 1/t - ln t/(t+1))],
S = 4 pi a0^2 N (R/B)^2. The outer3 rate is the orbital-weighted sum evaluated
with the dataset arithmetic operation by operation (`beb_rbeb.hpp`,
`beb_weighted_sum_m2`); it equals the dataset's own implementation
(`h2o_outer_valence.h` with the orbital-weight `collision_data.h` of the
dataset worktree) bit for bit on 213 energies, and the unweighted sum equals
the bound v1 `H2O:nist_BEB_valence` bit for bit.

**RBEB** (Kim, Santos & Parente, Phys. Rev. A 62, 052710 (2000)), Q = 1:

    sigma = 4 pi a0^2 alpha^4 N / ((bt^2 + bu^2 + bb^2) 2b') *
            { 1/2 [ln(bt^2/(1-bt^2)) - bt^2 - ln(2b')] (1 - 1/t^2)
              + 1 - 1/t - ln t/(t+1) (1+2t')/(1+t'/2)^2
              + b'^2/(1+t'/2)^2 (t-1)/2 }

t' = T/mc^2, b' = B/mc^2, u' = U/mc^2, bx^2 = 1 - 1/(1+x')^2 (evaluated as
x'(x'+2)/(1+x')^2; bt^2/(1-bt^2) = t'(t'+2)). mc^2 = 510998.95069 eV
(CODATA 2022; the value of the dataset's `speed()` and equal to the bridge's
m_e c^2/e to 1e-11). alpha^4 is written (2R/mc^2)^2 with the dataset R, so the
nonrelativistic limit is exactly the dataset BEB; the implied alpha =
7.2973544e-3 differs from CODATA 2022 by 2.5e-7 (R = 13.6057 eV vs
13.605693 eV), a 1e-6 relative effect on sigma. Evaluated in x86-64 extended
precision; agrees with an independent 50-digit Decimal evaluation to <= 2e-16.

**SDCS.** dsigma/dW = (P/B) g(W/B), 0 <= w <= (t-1)/2, with
g(w) = -phi/(t+1)[1/(w+1) + 1/(t-w)] + [1/(w+1)^2 + 1/(t-w)^2] + c_b
+ L [1/(w+1)^3 + 1/(t-w)^3].
RBEB: P = the RBEB prefactor, phi = (1+2t')/(1+t'/2)^2, c_b = b'^2/(1+t'/2)^2,
L = ln(bt^2/(1-bt^2)) - bt^2 - ln(2b') (KSP 2000 RBED SDCS with the BEB dipole
substitution and the RBEB interference factor). BEB: P = S/(t+u+1), phi = 1,
c_b = 0, L = ln t (Kim-Rudd BEB SDCS, Q = 1). The dipole term is symmetric under
exchange of the two outgoing electrons, so the integral over the ejected
half-range reproduces the closed-form totals; the analytic integral is
G(w) = -phi/(t+1)[log1p(w) - log1p(-w/t)] + w/(w+1) + w/(t(t-w)) + c_b w
+ L/2 [w(w+2)/(w+1)^2 + w(2t-w)/(t^2 (t-w)^2)].

## 4. Comparisons and seams

Outer3 (nonrelativistic BEB x Tan) versus RBEB x Tan, 1-10 keV (this is the
relativistic correction alone; same orbitals, same weights):

| E (eV) | outer3_BEB_Tan_parent_v1 (m^2) | RBEB x Tan (m^2) | difference |
|---|---|---|---|
| 1000 | 5.4761e-21 | 5.4882e-21 | +0.223 % |
| 1500 | 4.0099e-21 | 4.0240e-21 | +0.349 % |
| 2000 | 3.1872e-21 | 3.2024e-21 | +0.478 % |
| 3000 | 2.2848e-21 | 2.3016e-21 | +0.736 % |
| 5000 | 1.4851e-21 | 1.5037e-21 | +1.258 % |
| 8000 | 9.9078e-22 | 1.0111e-21 | +2.048 % |
| 10000 | 8.1582e-22 | 8.3685e-22 | +2.578 % |

Seams (recorded in the scope identity, not smoothed):
- 1 keV, Song -> outer3 BEB: 4.3200e-21 -> 5.476053e-21 m^2, **+26.760 %**
  (relative jump 0.267605) between 1000 eV and the next binary64 above.
- 10 keV, outer3 BEB -> provisional RBEB: 8.158186e-22 -> 8.368496e-22 m^2,
  **+2.578 %** (0.0257791); entirely the relativistic correction.

Relativistic correction RBEB/BEB - 1 (parent-weighted; valence within 0.05 %):
+0.22 % (1 keV), +2.58 % (10 keV), +5.26 % (20 keV), +13.5 % (50 keV),
+27.7 % (100 keV), +57.6 % (200 keV), +155 % (500 keV), +335 % (1 MeV). It
vanishes linearly in T/mc^2 at low energy (1.1-1.3 T/mc^2 at 1-10 keV); at MeV
energies the nonrelativistic formula is invalid (wrong velocity, no
relativistic rise). RBEB has its minimum between 1 and 5 MeV.

Total ionization against measured data in the dataset (Munoz et al., Phys.
Rev. A 76, 052707 (2007), Table I, ion-current gross ionization, 7 %):

| E (eV) | measured (m^2) | BEB valence | RBEB valence |
|---|---|---|---|
| 1000 | 7.113e-21 | +10.95 % | +11.19 % |
| 1500 | 5.125e-21 | +13.26 % | +13.65 % |
| 2000 | 4.060e-21 | +13.89 % | +14.43 % |
| 3000 | 2.940e-21 | +13.07 % | +13.90 % |
| 4000 | 2.366e-21 | +10.56 % | +11.66 % |
| 5000 | 1.966e-21 | +10.23 % | +11.61 % |

Song 2021 total ionization at 1000 eV: 7.05e-21 m^2 (RBEB valence +12.2 %).
The measured observable is charge-weighted and includes O 1s and multiple
ionization; the model total is a valence counting cross section. No measured
total above 5 keV is in the dataset. Validation targets to obtain (not in the
dataset, no numbers taken from secondary sources): Schutten et al., J. Chem.
Phys. 44, 3924 (1966), 0.1-20 keV; Rieke & Prepejchal, Phys. Rev. A 6, 1507
(1972), 0.1-2.7 MeV.

## 5. Kernel above 1 keV: `parent_ground_BEB_SDCS_binary_encounter_v2`

Replaces the equal-pair split above 1 keV (the outer3 model's own equal-pair
product law is not used). Frame conventions as v1: input T and unit direction
in the target rest frame, output momenta in that frame, exact boost to the
laboratory by the bridge (sealed `to_lab`), ground H2O+ of mass
M - m_e + I/c^2, chemical reference I = 12.621 eV counted once.
1. orbital i in {1b2, 3a1, 1b1} with probability w_i sigma_i(T)/sum (the band's
   own rate terms; 2a1 never selected);
2. ejected energy W = B_i w from the orbital SDCS (BEB below 10 keV, RBEB
   above) by bisection inversion of the analytic integral G (one variate);
3. ejected direction from the relativistic binary-encounter relation
   cos^2 theta = W(T+2mc^2)/(T(W+2mc^2)) about the incident direction, azimuth
   2 pi u (one variate);
4. scattered electron and ion as an exact two-body system carrying the
   remaining four-momentum, back-to-back along it in their rest frame; the ion
   recoil carries only the binding-energy momentum mismatch.
The scattered electron receives T - W - I - K_ion (the B_i - I excitation of
3a1/1b2 holes is not retained as ion excitation, as in v1/outer3). Three
product-RNG variates per event (Song-band events keep the single isotropic
axis of v1). The angular law is a declared approximation; the SDCS supplies
only the energy distribution.

Consequences (from the SDCS, per event above 1 keV): mean ejected energy
19.9 eV (1 keV) to 31.5 eV (1 MeV); 65.7-68.5 % of ejected electrons start at or
below 12.621 eV (known zero for this channel); **1.72-1.74 % start in the
unresolved 12.621-13.5 eV onset gap**, so the gap ledger now receives a steady
fraction of the secondaries of every keV-MeV event.

Conservation (4000 samples each at 1 keV, 10 keV, 10 keV + 1 ulp, 100 keV,
1 MeV): energy (with I) and momentum close to <= 5.3e-19 relative in the target
frame; ion recoil <= 2.4e-6 eV; the laboratory receipts (256 epsilon) close
for targets at rest, 600 m/s and 2e5 m/s; the native arena's own closure check
accepts every event.

## 6. Native registration, identities and ledger

The sealed native arena (candidate003) validates every SSA child of a
stochastic owner against one registered channel binding, and its full-time
exposure carries one stochastic owner. The outer3 rate therefore has its own
binding entry inside a composite channel binding rather than a second native
channel (a second channel would need a second owner, splitting the competing
hazards, or a native change). Without high-energy bands the scope registers
exactly the v1 binding.

- composite dataset identity `b7561244f7b779ba8c4c9db99e9c1e439e4d78dcdea455b4287dde3e450dd76e`
  (Song dataset sha, outer3 dataset sha/channel/model/weights, RBEB model
  specification, provisional admission and support);
- kernel manifest `28fe8904a073fa80badd64af3b78dc0ca597418c5526359ed80c1f6256f0475f`
  (`kernel-binding-v2.json`: v2 kernel sources, outer3 dataset, sealed bridge
  and scheduler headers, the v1 kernel manifest `02285c3a...` and binding file
  `82672015...`), binding file sha `556567ee...8013`;
- model `parent_banded_song2021_outer3BEB_provisionalRBEB_v2`.
Both enter the subflow configuration identity, the genesis RNG identity and
every RNG/clock state hash (with the scope identity and policy).

Event ids: `conditional/native-parent/N` (Song, as v1),
`conditional/native-parent-outer3-beb/N`, `model-derived/native-parent-rbeb/N`.
High-energy event receipts add category, basis, model_derived, kernel,
orbital, ejected/scattered/ion-recoil energies and the ejection angle.

Ledger (`parent-channel-coverage-ledger-v2`): v1 fields unchanged (unresolved
measures, bins, applied total, distinct carriers, decisions) plus
`categories`, three records per exposure, per step and cumulatively:
- `dataset_song2021_parent_table` (basis dataset),
- `dataset_outer3_BEB_Tan_parent` (basis dataset),
- `provisional_RBEB_Tan_parent_model` (basis model_derived, `model_derived: true`),
each with exact applied candidate evaluations, expected physical reactions and
events, sigma-free physical and event exposures (sum w_e w_n F dt/V and
sum max(w) F dt/V; with these a replacement sigma'(E) reprocesses the
expectation per bin), executed event count and exact executed event-weight sum,
energy range, and the same quantities per energy bin (v1 edges plus 2e5, 5e5,
2e6, 5e6, 1e7 eV). Restore requires the categories to sum to the applied total,
the bins to the category totals and the identities to match; the unresolved
category is the v1 measure. Descendants of model-derived events are traceable
through the product ids in the native journal.

## 7. Limit

Policy and rule unchanged from v1 (part of the scope identity):
refuse iff unbounded pairs present and not admitted, or (U_ev > event limit AND
U(1-f) > f A). A is the applied expected physical reactions of all three
categories. Consequences: 1000.5 eV, 8 keV, 20 keV and 1 MeV electrons that v1
refused are now applied; 20 MeV (above the provisional support) refuses as
unbounded; 13 eV still refuses; with the provisional band not admitted, 20 keV
refuses as in v1; the keV-MeV applied measure enters the denominator, so a
small gap weight beside it is admitted and a comparable one is refused.

## 8. Qualification

`run_coverage_v2_checks.py` (both-cores.sh or any-core.sh, frozen launcher,
one thread per process): receipts in
`outputs/pic-water-parent-coverage-v2-20261030/attempt-NNN/`. Qualified run:
attempt-002 (cores 6-7, both-cores.sh, 108 s). attempt-001 is retained as
failed: its oracle used the decimal spelling of B instead of the binary64
value the code receives, which near threshold differs by the binary64 rounding
of B amplified by T/(T-B); the oracle now converts the binary64 inputs exactly.
- `test_parent_coverage_v2`: 118 checks, normal and UBSan, results identical;
- v3 parent regression port without bands: 442 checks equal to the sealed
  baseline, recorded event physics equal, 73/128 first events;
- v3 parent regression port with all v2 bands and the composite registration:
  440 checks equal to the baseline with the two 8 keV entries removed (8 keV is
  covered), event physics equal, 73/128;
- coverage-v1 component test port without bands: 171 checks equal to the
  sealed v1 evidence, all exact ledger receipts equal apart from the added
  categories;
- outer3 rate bit-for-bit equal to the dataset implementation (213 energies);
- BEB/RBEB within 2e-16 of a 50-digit Decimal oracle (dataset double
  arithmetic within 5e-15 near threshold);
- benchmark (64 electrons, one water reservoir, exposure 2e-17 s, median of
  7): sub-keV population v2/v1 = 1.006 with identical events (77); population
  with 32 electrons at 1 keV-1 MeV v2/v1 = 2.87 for 54 instead of 23 executed
  events (v1 excludes those pairs); cost per pair visit 3.31e-5 s (v2) vs
  3.23e-5 s (v1); v2 without bands / v1 = 1.002.

## 9. Coordinator interface changes

- `coverage_v2::HighEnergyBands::read(outer3_path, outer3_sha, song_data,
  song_sha, kernel_binding_path, song_kernel_sha, song_binding_file_sha,
  ProvisionalOptions{admit, maximum_eV})`;
  `CoverageScope(data, dataset_sha, bridge, bands, policy)`.
- Register the parent channel with `scope.arena_binding(binding_file)`
  (composite identity) at arena creation and restore.
- `ParentSubflow`, `initial_full_time_rng_sha256`,
  `full_time_configuration_sha256`: same signatures as v1 (scope v2).
- Owner id `native_water_parent_subflow/coverage-v2`; `CoverageRefusal` as v1.
- Per-step ledgers now carry `categories`; reprocessing uses the
  `provisional_RBEB_Tan_parent_model` records and the `model-derived/` events.

## 10. Open items

- Measured total ionization to validate and later replace the bands:
  Schutten et al. 1966 (0.1-20 keV), Rieke & Prepejchal 1972 (0.1-2.7 MeV);
  measured parent H2O+ partial cross sections above 1 keV (none in any
  dataset version).
- Measured ejected-electron spectra of water to qualify the SDCS kernel
  (e.g. Bolorizadeh & Rudd, Phys. Rev. A 33, 882 (1986), 50-2000 eV; Opal,
  Beaty & Peterson, Atomic Data 4, 209 (1972)).
- The 12.621-13.5 eV onset gap remains unresolved and now receives 1.72-1.74 %
  of the secondaries of every event above 1 keV.
- Oxygen 1s and multiple ionization are outside the valence BEB/RBEB set.

## 11. Files

- `source/beb_rbeb.hpp` (formulas, SDCS, sampling), `source/high_energy_bands.hpp`
  (bound inputs, rates, composite identity), `source/beb_sdcs_parent_kernel.hpp`
  (kernel and bridge), `source/parent_coverage_v2_scope.hpp`,
  `source/parent_coverage_v2_subflow.hpp`; v1 `exact_sum.hpp` is used unchanged.
- `tests/test_parent_coverage_v2.cpp`, `tests/dataset_reference_outer3.cpp`,
  `tests/bench_exposure_v2.cpp`; generated ports
  `tests/test_v3_parent_regression_port_v2a.cpp`,
  `tests/test_v3_parent_regression_port_v2b.cpp`,
  `tests/test_v1_coverage_port_v2.cpp` (`tools/make_ports.py`, pins
  `V3-TEST-SHA256`, `V1-TEST-SHA256`).
- `tools/make_kernel_binding.py` -> `kernel-binding-v2.json`;
  `tools/oracle_beb_rbeb.py`; `tools/dev_build.sh` (development builds only);
  `run_coverage_v2_checks.py` (qualification).
