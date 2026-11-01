# Water chemistry v3.2: scope, data, models and qualification

v3.2 (this directory) extends v3.1 (`pic-water-chemistry-v3.1-20261031`, unchanged) with ion and anion + H2O collisions
above 1 keV relative energy (the 23 regions v3.1 still refused), sampled bridges over the v3.1 ion gaps below 1 keV, and
model continuations of electron records whose bounds made breakdown cells refuse; section 12 lists the changes, sources, models, documented bounds, the
decisions taken and the qualification. Sections 1-11 describe v3 and v3.1 and remain valid except where section 12
says otherwise.

v3.1 extends v3 (`pic-water-chemistry-v3-20261031`, unchanged) in two coverage gaps that made every
native v8 breakdown step refuse; section 11 lists the changes, sources, bounds and the decisions taken. Sections
1-10 describe v3 and remain valid except where section 11 says otherwise.

Component for the Power-Flow coupled PIC water-breakdown simulation: the electron, ion and neutral chemistry needed
for H2O-desorption breakdown, as a sealed, tested library that native v8 links later (INTEGRATION.txt). It closes
the pilot limitations L1 (parent ionization only) and L2 (ions do not react or collide) of
`artifacts/pf-pilot-1-20261031/LIMITATIONS.txt` as far as data allow, and the parent-coverage items of audit
A13/A14 (onset gap, 1 keV and 10 keV seams). Nothing here advances physical time in a coupled run.

Tables: `data/tables/chemistry-v3.1.blob` (22 species, 20 projectile-target groups, 80 processes, 135 coverage regions),
built by `tools/build_tables.py` from `data/tables/chemistry-v3.1-spec.json` (written by
`tools/make_production_spec.py`). Every number comes from a retained source file through a documented source type;
`data/tables/PROCESSES.txt` lists every process with its sources, tags, coverage intervals, interpolation errors,
junctions and notes; `data/tables/chemistry-v3.1-manifest.json` carries the same in machine form with the sha256 of
every input file. Source indexes: `data/sources/ELECTRON-SOURCES.json`, `ION-NEUTRAL-SOURCES.json`,
`METAL-SOURCES.json` (citations, Crossref checks, retrieval times, sha256 of raw files, missing items).

## 1. Tags and coverage classes

Tags (per process and per segment, carried into every event): `measured`, `recommended` (evaluated compilations:
Song 2021, Itikawa & Mason 2005, Anicich 1993, Lennon 1988, Simon Wedlund 2019), `calculated` (NIST BEB),
`model_derived` (swarm-fitted sets, model extensions, model energy partitions).
Coverage classes per energy interval: ZERO (below the thermochemical/kinematic threshold or a source statement),
DATA, MODEL, UNRESOLVED (no data: sampled as zero; exposure and, where available, a bound recorded per region).
No gap is filled by interpolation: every free support edge carries an exact zero one binary64 outside it.

## 2. Species and thermochemistry

Rest energies and reference energies from ATcT v1.202 0 K formation enthalpies (Ruscic & Bross 2024,
DOI 10.17038/CSE/2440256; ions in the stationary-electron convention), CODATA 2022 (m_e c^2, m_u c^2) and NIST atomic
masses (1H, 16O; natural mean masses for Fe, Cr, Ni), NIST ASD first ionization energies for Fe, Cr, Ni.
Every reaction Q is the reference-energy difference, e.g. IE(H2O) = 12.6175 eV, D0(H-OH) = 5.1014 eV,
H2O+ + H2O -> H3O+ + OH: -0.999 eV. Masses satisfy m(products) - m(reactants) = Q by construction; kinematics use Q,
never mass differences.

## 3. Electrons on H2O (0 to 10 MeV)

Backbone: the MuroranIT complete set (Kawaguchi et al., Jpn. J. Appl. Phys. 55, 07LD03 (2016); LXCat), validated by
its authors against swarm data, with the recommended Song et al. 2021 (J. Phys. Chem. Ref. Data 50, 023103) values
where they are measurement-based. Itikawa & Mason 2005 is available only through its LXCat transcription (its
ionization blocks equal Song Table 12). Other LXCat sets (TRINITI, Hayashi, Phelps, Morgan) are used only as swarm
references.

| Process | 0-1 keV | 1-10 keV | above 10 keV |
|---|---|---|---|
| Elastic momentum transfer, screened angular law | MuroranIT to 100 eV (model_derived); Song Table 4 100 eV-1 keV (recommended) | IAM (NIST SRD 64) energy dependence normalized to Song at 1 keV | UNRESOLVED, bound IAM integral at 10 keV |
| Rotational J=1,2,3 (0.0046, 0.0087, 0.0176 eV) | MuroranIT | UNRESOLVED, bound Munoz 2007 rotational integral | same |
| Vibrational (010) 0.198 eV, (100)+(001) 0.466 eV | MuroranIT | UNRESOLVED, bound Munoz total inelastic | same |
| 25 electronic states 7.14-12.5 eV -> OH + H | MuroranIT | Munoz 2007 recommended-inelastic energy dependence normalized at 1 keV | RBEB energy dependence normalized at 10 keV |
| Ionization H2O+, OH+ + H, H+ + OH, O+ + H2, H2+ + O, O++ + H2 (+3e) | MuroranIT onset (model) below the Song start; Song Table 12 (Straub et al. 1998 lineage) | Munoz 2007 measured gross-ionization energy dependence (to 5 keV) normalized to each partial at 1 keV | RBEB energy dependence normalized at 5 keV |
| Ion pairs H+ + OH-, OH+ + H- | H+ and OH+ yields below their dissociative-ionization thresholds | UNRESOLVED (inseparable from the H+/OH+ yields), bound total inelastic | same |
| Dissociative attachment H- + OH, OH- + H, O- + H2 | MuroranIT | UNRESOLVED above the resonances, bound total inelastic | same |

Angular treatment: the elastic channel is sampled as sigma_m / (1 - <mu>(xi)) with the screened law
f(mu) = (1 - xi^2) / (2 (1 - xi mu)^2); xi(E) is solved from sigma_m/sigma_el (Song Tables 3/4 0.1-100 eV; Song
Table 4 / Munoz elastic integral 100 eV-1 keV; IAM 1-10 keV), so the momentum-transfer cross section is reproduced
(midpoint error <= 1.2e-3; isotropic below 0.1 eV). Inelastic channels are isotropic in the CM frame, as in the swarm
validation of the backbone set.

Kinematics (exact relativistic two-body chains in eV units): excitation leaves the target with internal energy
(lab-frame sink gamma*I); electronic states dissociate to ground OH(X) + H(1s) with all excess (state energy - D0) as
fragment kinetic energy; ionization chooses an orbital (NIST BEB 1b1/3a1/1b2/2a1 constants) by w_i sigma_i(T) with
Tan/Montenegro channel weights, samples the ejected energy from the BEB SDCS (RBEB above 10 keV) by 64-step bisection,
ejects along the binary-encounter angle, and shares the remaining four-momentum between the scattered electron and the
ion system; the orbital excess B_i - IE goes to the ion-internal (parent) or fragment-excitation sink (no fragment
kinetic energy release: kinetic fraction 0, open). Attachment, recombination and ion-molecule reactions decay an
isotropic compound in the CM with all exothermicity as product kinetic energy (ground products).

### 3.1 Onset gap and seams (audit A13/A14)

- 12.621-13.5 eV: no measured H2O+ value exists strictly inside (searched: Song, Itikawa/LXCat, Rao 1995, Straub 1998,
  Djuric 1988 -- the primaries are not retrievable here). 12.65-13.5 eV is covered by the MuroranIT swarm-fitted onset
  (model_derived); 12.6178-12.65 eV stays UNRESOLVED with bound 2.83e-24 m2 (monotone onset). Morgan has one interior
  point (13.235 eV, 1.0e-23 m2), also model-derived.
- 1 keV seam: removed for all six channels (continuous by normalization); residual normalization difference of the two
  measured data sets: Munoz gross ionization 7.113e-21 vs Song gross 7.05e-21 m2 at 1 keV (+0.9 %).
- 10 keV seam: removed (one RBEB energy dependence from 5 keV). No measured ionization above 5 keV was obtainable.
- Relative to the v2 bands, the v3 parent H2O+ cross section is 18-23 % lower from 1 keV to 10 MeV (measured gross
  ionization below BEB; measured 1 keV parent fraction 0.613 vs dipole 0.694; fractions held at 1 keV above).
- Thermochemistry cuts: source values below the ATcT thresholds are set to zero and listed (OH+ 17.5-18.1 eV and
  H+ 16.95-18.70 eV become the ion-pair channels; DEA onsets below 4.35/3.27/3.57 eV).

### 3.2 Consistency findings

- Sum of the six partial cross sections at 1 keV: 7.043e-21 m2 vs Munoz measured gross 7.113e-21 m2 (-1.0 %).
- Neutral dissociation: sum of the electronic states vs the measured OH yield (Harb et al. 2001, Song Table 10) is
  x2.8-3.3 at 15-20 eV, x1.2 at 50 eV, x0.5 at 200-300 eV (OH from H+ + OH ionization is extra there). The energy
  loss is swarm-validated; the product partition of electronic excitation (OH + H only) is an open model item.
- A1B1 of MuroranIT equals Song Tables 8/9 within 1 % above 30 eV.
- Elastic momentum transfer: MuroranIT/Song = 3.9-6.2 below 1 eV (swarm-adjusted low-energy set), 0.98 at 100 eV,
  5.6 at 1 keV (hence Song above 100 eV).

## 4. Electrons on other targets

- e + H(1s) -> 2e + H+: NIST BEB (B = U = 13.6057 eV, N = 1, Q = 0.5668; calculated), 14 eV-10 keV; 13.6058-14 eV and
  above 10 keV UNRESOLVED (bounded). Excitation/elastic of H: open.
- e + Fe/Cr/Ni -> 2e + M+: Lennon et al. 1988 recommended fits (McGuire 1977 Born basis; stated errors 40 % Fe/Cr,
  200 % Ni), threshold-10 keV; no measured ionization, excitation or elastic data obtained (Freund 1990, Shah 1993
  paywalled). Ionization only.
- e + H2O+ dissociative recombination (Rosen et al. 2000): alpha(300 K) = 4.3e-7 cm3/s, sigma ~ E^-1.24 (< 0.1 eV)
  and E^-1.55 (0.1-2 eV) as stated by the source, normalized to alpha(300 K); branching OH + H 0.20, O + H2 0.09,
  O + 2H 0.71 (printed Eq. 5 inconsistency of the source kept as the balanced Eq. 4 lineage); above 2 eV UNRESOLVED.
- e + H3O+ (Jensen et al. 2000): alpha(300 K) = 4.3e-7, sigma ~ E^-1.15 to 0.2 eV; branching H2O + H 0.25,
  OH + H2 0.14, OH + 2H 0.60, O + H2 + H 0.013; other measured alpha(300 K) 7.6e-7 (Neau 2000), 1.0e-6 (Herd 1990) kept
  as spread.
- e + OH+: Kalosi 2023 alpha(T) fit retained; no cross section; placeholder (UNRESOLVED).
- e + OH, O, H2 fragments: no data obtained (open).

## 5. Ions, negative ions, recombination

| Group | Process (data) | Model |
|---|---|---|
| H2O+ + H2O | proton transfer -> H3O+ + OH: Anicich 1993 k(300 K) 1.85e-9 cm3/s; Lishawa et al. 1990 Fig. 5 (digitized, isotopic analogue) E_cm 1.2-30 eV | constant k below the first beam point (junction -43 %); compound decay |
| | symmetric charge exchange: Lishawa Fig. 7 (digitized) E_cm 0.8-30 eV | forward CM charge exchange (fast neutral H2O to the reservoir) |
| H2O+, H3O+, OH+, H+, O+, H2+, O++, H-, O-, OH- + H2O | momentum transfer: no measured mobility or cross section obtained | Langevin polarization sigma_m = 1.105 k_L / v with alpha(H2O) = 1.501 A^3 (NIST CCCBDB), 1e-5-1e3 eV; lower estimate above ~1 eV (no ion-dipole, no repulsive core) |
| OH+ + H2O | Anicich k 2.89e-9, H2O+ + OH 0.55 / H3O+ + O 0.45 | constant k to E_cm = 1 eV |
| H2+ + H2O | Anicich k 7.3e-9, H2O+ + H2 0.53 / H3O+ + H 0.47 | constant k to E_cm = 1 eV |
| H+ + H2O | Anicich k 8.2e-9 -> H2O+ + H; Simon Wedlund 2019 electron capture 52 eV/u-130 keV/u (recommended) | constant k to E_cm = 1 eV; forward charge transfer at keV |
| O+ + H2O | Heninger et al. 1989 ICR sigma(E_cm) 0.033-11.6 eV; Anicich k 2.6e-9 | constant k below 0.033 eV (junction -24 %) |
| H- + H2O | collisional detachment 52 eV/u-130 keV/u (Simon Wedlund 2019, low confidence, model_derived) | |
| H- + H | associative detachment, Miller et al. 2011 k = 5.33e-9 cm3/s at 16 meV | constant k to E_cm = 0.1 eV |
| placeholders (never sampled; exposure recorded) | H- + H2O -> OH- + H2; O- + H2O -> OH- + OH; O- and OH- detachment in H2O; e + OH+ DR; H3O+ + OH- mutual neutralization | no rate obtained |

Three-body clustering (H3O+(H2O)n) is outside the binary-collision scope (open).

## 6. Neutral products

H, OH, O, H2 (and fast H2O from charge exchange) go to per-(cell, species) reservoirs and are released by the approved
quantized rule n = SR(D/w0) as fixed-weight packets (INTEGRATION.txt section 6); exact count, momentum and kinetic
energy stocks; released packets take the velocity and position of a weight-proportional pick among the step's
deposits. Collisions of neutral fragments with the background (elastic/diffusion, O(3P)/H reactions with H2O) have
no energy-resolved data here and are open; at 3e22 m^-3 a 1e4 m/s fragment collides about 0.1 times per ns
(estimate with a 3e-19 m2 geometric cross section), so they matter only in dense desorption plumes.

## 7. Numerics

- Tables: per group a union energy grid holding every source breakpoint, threshold and coverage edge; linear in E and
  sigma between grid points (the source law); tabulated sources reproduced bit for bit at their breakpoints; analytic
  models sampled at 64-1000 points per decade with recorded midpoint errors (<= 1e-4 typical; Lennon fits near
  threshold 6e-3).
- Majorant: sparse-table range maximum of sigma_tot v over grid intervals (log-concave interval maxima x (1 + 2^-20));
  per carrier the bound covers every target-frame energy reachable inside the truncated (7 sigma) target Maxwellian.
- RNG: Philox4x32-10 (arrival core), counter (seed, step, carrier id, candidate, word); every draw is a function of
  those only: partition- and order-independent.
- Deterministic math: own log/log1p/sincos (fdlibm-type, <= 1 ulp / 2 ulp / 3e-16 abs), only + - * / sqrt; host
  -ffp-contract=off, device --fmad=false --prec-div --prec-sqrt --ftz=false: bit-identical CPU/GPU.
- Exact ledgers: ExactSum (36 x 64-bit two's-complement limbs, unit 2^-1074) with Dekker products.

## 8. Swarm validation

`SWARM-VALIDATION.txt` (receipt qual-002; `tests/swarm_mc.cpp`, N = 3.22e22 m^-3, 300 K, 2000 electrons per point,
12 E/N points 30-3000 Td; measured data: Hasegawa et al. mobility x N and (alpha - eta)/N, Crompton et al. 1965
D_T/mu, all from LXCat LAPLACE; separate alpha/N and eta/N data were not found):
- (alpha - eta)/N, v3: within -11 % / +3 % from 198 to 3000 Td (+109 % at 152.5 Td, next to the zero crossing near
  140 Td where the relative error is not meaningful).
- Drift velocity, v3: flux drift within -10 % / +9 % from 60 to 800 Td; bulk drift within -2 % / +15 % up to 350 Td
  but +19 % to +38 % at 400-800 Td. The LXCat entry does not state which drift definition the measurement is.
- Angular treatment: the same tables with isotropic elastic scattering (v3-iso) reproduce the bulk drift within 8 %
  but overestimate (alpha - eta)/N by +12 % to +87 % above 500 Td; the screened law is kept as the production choice.
- D_T/mu at 64-188 Td: v3 +12 % to +59 % (all sets except Phelps and TRINITI are high); low-E/N transport is not
  qualified (no superelastic collisions, open item 4).
- Reference sets in the same Monte Carlo: MuroranIT as published has no data above 1 keV, so above ~1000 Td its
  electrons run away and (alpha - eta)/N is overestimated by +39 % to +127 %; TRINITI, Hayashi and Morgan deviate by
  tens of percent; Phelps agrees in (alpha - eta)/N within 6 % from 500 to 2300 Td.

## 9. Qualification

Receipt `/home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3-20261031/qual-002/receipt.json` (all stages
passed; summary in `RESULTS.txt`):
- inputs: every table input file unchanged (sha256); tables: spec regenerated and blob rebuilt byte-identical
  (blob sha256 3ca3ea495630e6b6e125c05ebf45488fd89a899544a2aa138cefb9b1dd415b0a);
- data: all 2263 tabulated source breakpoints reproduced bit for bit; majorant bounds sigma_tot v at 4e6 probes;
- ledger closure: 1,007,763 sampled events (2e6 carriers, all 20 groups, 9 target species): charge and nuclei
  conserved per event; energy closure per event within 2^-40 of its energy scale (max 6.1e-14 relative); exact ledger
  identity in = out + q + internal + defect; total defect 0.55 eV of 9.5e16 eV;
- per process: each of the 75 processes with data gives 2000 conserving events at its peak sigma*v energy;
- RNG determinism: records and exact ledgers byte-identical for 1, 2, 5 and 8 partitions in permuted order;
- CPU/GPU parity (RTX 3080 Ti, sm_86 cubin, NVRTC 12.8): 4.2e6 elementary-function values and 3 step cases x 3 launch
  geometries (including a split, out-of-order launch) x 2 repeats bit-identical, in two fresh processes;
- sanitizers: AddressSanitizer + UndefinedBehaviorSanitizer build passes every check with the same partition digest.

## 10. Open items (data not obtained or models to replace; v3.1 status in section 11)

1. Measured H2O+ inside 12.6178-13.5 eV (Rao 1995, Straub 1998 primaries); measured ionization above 5 keV
   (Schutten 1966, Rieke & Prepejchal 1972); channel fractions above 1 keV (held at 1 keV).
2. Fragment-ion kinetic energy release (kinetic fraction 0); ejected-electron spectra (BEB SDCS model).
3. Product partition of the electronic states (all OH + H; OH yield mismatch section 3.2), O + H2 / O + 2H / OH(A).
4. Superelastic rotational/vibrational collisions (thermal rotational population): swarm agreement below ~100 Td not
   claimed.
5. Elastic, rotational, vibrational and attachment above their data (bounded, unresolved); elastic above 10 keV.
6. Ion-neutral momentum transfer (no mobility vs E/N data; Langevin lower estimate); energy-dependent ion-molecule
   cross sections above E_cm ~1 eV (OH+, H2+, H+) and above 30 eV (H2O+); O++ reactions.
7. Negative-ion conversion and detachment in H2O (only H- detachment at >= 52 eV/u); ion-ion recombination rates.
8. DR energy dependence above 2 eV (H2O+) / 0.2 eV (H3O+); OH+ DR cross section; the H3O+ alpha spread (x2.3).
9. e + H excitation/elastic; e + OH, O, H2; metal excitation and elastic (only figures exist; NIST SRD 64 elastic
   needs an interactive session); measured metal ionization (Freund 1990, Shah 1993).
10. Neutral-fragment collisions with the background (no transport data obtained).

## 11. Changes in v3.1 (2026-10-31)

### 11.1 Dissociative recombination above the stated low-energy laws (gap 1)
- e + H2O+: E^-1.24 (< 0.1 eV) and E^-1.55 (0.1-2 eV) as stated by Rosen et al. 2000 and normalized to their
  alpha(300 K) = 4.3e-7 cm3/s (model_derived), joined continuously at 2.0 eV to the measured cross section of Rosen
  Fig. 4 (CRYRING), digitized from 2.16 to 19.1 eV (tag measured, figure-digitized; resonances near 4 and 11 eV
  included: 3.7e-21 m2 at 10.5 eV, about 8x the continued power law). The figure's last two bins (23 and 28 eV)
  scatter beyond a factor-2 reading uncertainty and end the table at 19.1 eV.
- e + H3O+: E^-1.15 to 0.2 eV (Jensen et al. 2000, alpha(300 K) = 4.3e-7 cm3/s), joined at 0.2 eV to Jensen Fig. 3
  (this work, x10 removed), digitized 0.21-19.9 eV (valley 0.6-5 eV and the 10 eV peak included).
- Above the tables: unresolved, bounded by the power-law continuation (E^-1.55, E^-1.15), not sampled, exposure
  recorded. Below: the group grids start at 1e-9 eV with the stated laws continued below the lowest data (Rosen data
  from 3e-4 eV), and constant sigma*v below 1e-9 eV (model_derived).
- Branching held at the 0 eV values at all energies (model).
- Digitization (`tools/digitize/`, provenance in `data/sources/rec-dr-figure-digitization-20261031/`): native-resolution
  figure images (pdfimages); axis calibration from the frame ticks, bottom/left and top/right fitted independently
  (residuals 0.3-1.5 px, A-B differences <= 0.009 decades); Jensen deskewed by 0.00401 rad. Rosen: filled-marker band
  trace (opening removes open markers and error bars), one marker-width bins. Jensen: empirical-template detection
  inside the H3O+ series band with a visual review (9 removals, 8 manual additions, all listed in the table);
  above 6 eV the overlapping previous-ASTRID points are traced with this work (the source states agreement outside
  the 0.6-5 eV valley). Reading uncertainty per point (1 sigma, quadrature of calibration, +-1 px, band scatter and a
  repeat pass): H2O+ 5-57 %, median 23 % (dense, overlapping markers); H3O+ markers 2-4 %, band 12-38 %. The sources'
  own uncertainty (~20 % at 1 sigma; up to 30 % statistical above 20 eV for H2O+) is separate.
- Weight of the low-energy extension (Maxwellian electrons on thermal ions; `dr-lowenergy-share.json` in the outputs):

| | kTe | phase space E < 1e-5 eV | DR events E < 1e-5 eV | DR events E < 1e-9 eV | DR events below the data (3e-4 eV) |
|---|---|---|---|---|---|
| H2O+ | 0.026 eV (300 K) | 5.7e-6 | 2.8e-3 | 1.3e-6 | 3.7e-2 |
| H2O+ | 1 eV | 2.4e-8 | 2.6e-4 | 1.2e-7 | 3.4e-3 |
| H2O+ | 3 eV | 4.6e-9 | 1.2e-4 | 5.6e-8 | 1.6e-3 |
| H2O+ | 10 eV | 7.5e-10 | 3.9e-5 | 1.8e-8 | 5.2e-4 |
| H3O+ | 1 eV | 2.4e-8 | 1.6e-4 | 3.5e-8 | 2.3e-3 |
| H3O+ | 10 eV | 7.5e-10 | 2.6e-5 | 5.9e-9 | 3.9e-4 |

### 11.2 Ion + H2O reactions below their tables (gap 2)
- Constant-rate continuation below the group grids (kernel flag): sigma v = sigma(E_0) v(E_0) down to zero relative
  energy for every ion + H2O reaction, the Langevin momentum transfer of all ions (unchanged model, 1e-5-1e3 eV) and
  H- + H associative detachment.
- H2O+ + H2O proton transfer: Anicich k(300 K) = 1.85e-9 cm3/s (+-15 %) up to the 300 K mean CM energy
  3/2 kT = 0.0388 eV, log-log interpolation to the first beam point (Lishawa, E_cm 1.23 eV), then the table. The
  constant-k value at 1.23 eV would be 43 % above the beam point; that difference is the stated uncertainty of the
  0.04-1.2 eV region.
- H2O+ + H2O symmetric charge exchange: constant rate anchored at the first beam point (E_cm 0.79 eV; no thermal
  measurement exists), model_derived.
- O+ + H2O: constant rate equal to the first ICR point's sigma*v (Heninger, E_cm 0.033 eV); Anicich k(300 K) is
  24 % above it (stated uncertainty).
- OH+, H2+, H+ + H2O: Anicich k(300 K) constant to E_cm = 1 eV and to zero energy; above 1 eV unresolved, bounded by
  k/v. H+ + H2O is one charge-transfer process (thermal, bounded gap, Simon Wedlund keV data).
- H- + H2O detachment: zero below 0.75 eV (endothermic); unresolved to 52 eV/u, bounded by the Langevin capture rate.

### 11.3 Coverage policy and the remaining refusing regions
`data/tables/chemistry-v3.1-coverage-policy.json`: 60 bounded regions (refuse iff B_r > max(1, 0.01 A)), 52 admitted
unbounded regions (above 10 MeV, or placeholder-only), 23 refusing unbounded regions. All 23 refusing regions are
ion + H2O (and negative ion + H2O) above 1 keV relative energy, where the v3 Langevin momentum-transfer model ends;
they refuse only when an ion moves faster than 1 keV relative to the gas.

### 11.4 Decisions taken (data handling)
- DR join: the anchor value of the stated power law at the join energy starts the digitized table, so the join is
  continuous; digitized points below the join verify it (H2O+ at 1.94 eV: 6.1e-21 vs the law 5.9e-21 m2).
- End of the H2O+ table: before the first bin whose reading uncertainty exceeds a factor 2.
- H3O+ above 6 eV: one band trace of this work and the previous ASTRID points (the source states agreement there).
- Digitization uncertainty: calibration RMS and A/B difference, +-1 px (manual readings +-2 to +-4 px), band
  scatter minus the marker's own extent, and the difference to a second, independent pass.
- Placeholders do not make a region unbounded (their exposure is still recorded).
- H+ + H2O: one process for the one reaction (forward charge-transfer kinematics at all energies).
- H- + H associative detachment above 0.1 eV: bound = measured rate / v (declared; the source reports a falling rate).

### 11.5 Qualification of v3.1
Receipt `/home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3.1-20261031/qual-001/receipt.json`, all stages
passed (summary in `RESULTS.txt`): table inputs unchanged; spec, blob, PROCESSES.txt and coverage policy regenerated
byte-identical; 1,008,311 sampled events with exact charge/nuclei balance, per-event energy closure within 2^-40 (max
6.1e-14) and the exact ledger identity; all 2498 tabulated breakpoints reproduced bit for bit; 74 processes with data
give conserving events; constant-rate continuation verified in 13 groups; records and exact ledgers identical for
1, 2, 5, 8 partitions; ASan + UBSan pass; CPU/GPU bit identity (4.2e6 math values; 3 step cases x 3 launch geometries
x 2 repeats in 2 fresh processes); swarm values identical to v3 (504 of 504).

## 12. Changes in v3.2 (2026-10-31)

Scope: every ion and anion of the registry (H+, H2+, O+, O++, OH+, H2O+, H3O+, H-, O-, OH-) colliding with H2O from
1 keV to 10 MeV (projectile kinetic energy in the gas frame), and the v3.1 gaps of the ion records below 1 keV (12.4a).
Below 1 keV the v3.1 data are kept; a process new in v3.2 is a placeholder there (sampled as zero, exposure recorded,
admitted), as the v3.1 placeholders. Fast neutral
products go to the quantized-release reservoir and are not projectiles (no fast-neutral collisions are modelled).

Decisions of 2026-10-31: data where they exist; where they do not, the dominant channels are
sampled from model estimates (tag model_derived, uncertainty stated per process) and the physics upper bound of each
estimate is documented, not used (manifest `documented_bound`, PROCESSES.txt, and the policy file's `model_estimates`
with the ratio bound / estimate); "fast heavy ions do not stop runs". The same rule was applied to four electron
records whose v3 bounds made typical breakdown cells refuse (12.4).

### 12.1 Sources (data)
- Simon Wedlund et al., A&A 630, A35 (2019) (arXiv:1902.05048): recommended polynomial fits log10 sigma = sum p_j
  (log10 v)^j and tables for H+, H, H- (Tables 2, 4, A.2, B.1) and He+, He2+ (Tables 1, 3, A.1, B.1) on H2O; fits used
  with their stated speed validity, and below it only where the source's appendix extrapolation is reproduced (3 %).
- Rudd, Kim, Madison, Gay, Rev. Mod. Phys. 64, 441 (1992): SDCS model; Table V water-vapour set (0.97, 82, 0.40,
  -0.30, 0.38, 1.04, 17.3, 0.76, 0.04, 0.64) and the universal inner-shell set (1.25, 0.50, 1.00, 1.00, 3.0, 1.10,
  1.30, 1.00, 0, 0.66) for shells with B > 2 B_min (2a1); Table I(b) binding energies 12.61, 14.73, 18.55, 32.2 eV.
- Gobet et al., Phys. Rev. A 70, 062716 (2004), Table II: partial ionization fits (H2O+, OH+, O+, H+; O2+ omitted,
  < 2 %) for 20-150 keV protons, used as product shares (held at the edge values outside).
- NIST PSTAR (ICRU 49) nuclear stopping of protons in water vapour; ZBL universal nuclear stopping (Ziegler, Biersack,
  Littmark 1985; formula from Barker & Mei 2012) for the other ions.
- Juhasz et al., Phys. Rev. A 100, 032713 (2019): O+ + H2O at 6.6 keV, total positive fragment-ion production
  (6 +- 3)e-20 m2 (validation only). UMIST Rate22: H- + H2O -> OH- + H2, k(300 K) = 4.8e-9 cm3/s (Martinez 2010).
- Searched without tabulated keV results for O+, O++, OH+, H2O+, H3O+, H2+, O-, OH- on H2O (literature helper
  findings: data/sources/ion-v32-literature-findings-20261031/FINDINGS.txt).

### 12.2 Content per projectile above 1 keV
| projectile | data | model estimates (sampled) | still bounded (not sampled) |
|---|---|---|---|
| H+ | capture sigma10 (Wedlund table to 131 keV, fit to 2.1 MeV); ionization sigma11 (fit 0.84 keV - 4.7 MeV) x Gobet shares (DATA 20-150 keV); momentum transfer PSTAR | ionization above 4.7 MeV (Rudd SDCS integral scaled at the edge); shares outside 20-150 keV; double capture sigma1-1 (low confidence) | capture above 2.1 MeV and sigma1-1 above 300 keV (edge values held; B/A <= 1.5e-4) |
| H2+, O+, OH+, H2O+ | - | capture = proton sigma10 at equal speed; ionization = He+ sigma11 at equal speed x Gobet shares; breakup = free-collision estimate; proton transfer = thermal rate / v (OH+, H2+; H2O+: last Lishawa point as constant rate) | H2+ capture above 4.2 MeV (edge held) |
| H3O+ | - | capture (H2O + H + H2O+) = proton sigma10 x exp(-a (|dE| - |dE_ref|)/(h v)), a = 7 A, dE = 6.10 eV, dE_ref = 0.98 eV (Massey 1949, Hasted 1951); ionization He+ sigma11; breakup free-collision | - |
| O++ | - | capture = He2+ sigma21 at equal speed (proton sigma10 shape above 5350 km/s); ionization = He2+ sigma22 x Gobet shares | - |
| H- | detachment sigma-10 and double detachment sigma-11 (Wedlund, low confidence; model_derived as in v3.1) | ionization = proton sigma11 at equal speed; sigma-10, sigma-11 continued as v^-2 above 20000 km/s; proton transfer = k300/v | - |
| O-, OH- | - | detachment = free-collision estimate; ionization = proton sigma11 at equal speed; O- rearrangement = Langevin k_L/v | - |
| all ions | momentum transfer ZBL (PSTAR for H+), log-log bridge from the Langevin value at 1-2 keV | | |

Estimates and their documented bounds (bound / estimate at 10, 50, 100 keV in the policy file):
- capture: classical over-barrier pi R_c^2, R_c = (2 sqrt(q) + 1)/I_T (Ryufuku, Sasaki, Watanabe, PRA 21, 745
  (1980); Niehaus, J. Phys. B 19, 2925 (1986)), 3.7e-19 m2 (q = 1), 6.0e-19 m2 (q = 2); estimates 1.9-3x below at
  10-100 keV (H3O+ further, by the adiabatic factor). Uncertainty factor 3 (H3O+: 5, an order of magnitude below 10 keV).
- target ionization: pi (R_P + R_T)^2 (Bondi radii, J. Phys. Chem. 68, 441 (1964)) + q^2 x proton sigma11 at equal
  speed. Uncertainty factor 3. Validation: at 6.6 keV O+ the He+ rule gives 1.4e-20 m2 against the measured total
  positive fragment-ion production (6 +- 3)e-20 m2 (the proton rule gives 6.1e-22 m2, which is why He+ is used).
- breakup / detachment (free-collision estimate): the most weakly bound projectile electron scattered as a free
  electron of the projectile speed by H2O (Munoz 2007 elastic + inelastic integral at E_e = (m_e/m_P) E) times the
  isotropic share 1 - I/(4 E_e) of collisions transferring at least the breakup threshold I (ATcT); zero below the
  onset E = (m_P/m_e) I/4. Bound: pi R^2 + the same free-electron term. Checks: the estimate is 1.9-3.1x above the
  H- + H2O sigma-10 (1-130 keV) and 2.6-11x above the He+ + H2O electron loss sigma12 (3000-10000 km/s, largest near
  the onset) -- the conservative side for a sampled estimate; molecular (close-collision) breakup below the onset is
  not represented.
- rearrangements (proton transfer, O- + H2O -> OH- + OH): falling 1/v forms anchored to the measured thermal rate or
  the last measured cross section (O-: Langevin capture rate, unit reaction probability); bound pi R^2 (H2O+: the v3.1
  bound, last measured value held). Uncertainty factor 10. Compound kinematics as below 1 keV.
- Ejected electrons of every ion-impact ionization: Rudd SDCS at the projectile speed (composition-rejection over the
  orbitals; power or exponential envelope, whichever has the smaller mass, exact for either), isotropic ejection.

### 12.3 Kernel and builder changes
- New process kinds ION_IONIZE (7, Rudd SDCS ejected-electron energy), BREAKUP (8, projectile fragments forward and the
  target product last; also used for H3O+ capture with the unbound H3O), CHARGE_X with n_out > 2 (capture with target
  fragmentation); deterministic exp and pow in det_math.h; ProcessRec layout unchanged (512 bytes, rudd[10] in former
  padding). Rudd trials use event words 32..511; kinematics words stay below 32.
- Builder: placeholder_below_eV (a process is a placeholder below that energy; such processes are marked
  "[placeholder here]" in region names), documented_bound, sources sum_of, rudd_analytic, branch_fraction,
  free_collision, adiabatic_scaled, close_distant_bound, rudd_total, wedlund_fit, nuclear_stopping_mt,
  zbl_momentum_transfer, equal_velocity, powerlaw_from, product_of, min_of.
- Build hook: PF_BUILD_ISA_FLAGS (whitespace-split, default empty) is appended to the C++ flags of the builds and
  recorded in the receipt.

### 12.4a Ion gaps below 1 keV (decision of 2026-10-31: drifting breakdown ions reach eV to hundreds of eV)
v3.1 left these gaps unresolved with bounds (B/A 2-18 for ion beams at 3 eV - 1 keV); v3.2 samples log-log bridges from
the last low-energy value to the upper end, model_derived, with the v3.1 bound documented (ratio in the policy file):
- OH+, H2+ charge transfer: from k(300 K) x branch / v at E_cm 1 eV to the equal-speed proton sigma10 at 1 keV (factor 5);
  proton transfer of OH+, H2+: the thermal rate continued as a constant rate from E_cm 1 eV (factor 10).
- H+ charge transfer: from k(300 K) / v at E_cm 1 eV to the first Simon Wedlund 2019 table value at 52 eV/u (factor 3).
- H2O+ proton transfer: constant rate of the last Lishawa point from that point on (factor 10); symmetric charge
  exchange: from the last Lishawa point to the equal-speed proton sigma10 at 1 keV (factor 5).
- O+ charge transfer: from the last Heninger 1989 point (E_cm 11.6 eV) to the equal-speed proton sigma10 at 1 keV
  (factor 3).
- H- detachment: from zero at the threshold (0.796 eV lab) to the first Wedlund table value (52 eV/u) as a power law in
  (E - E_th), continuous in value and log-log slope at the table start (factor 3).

### 12.4 Electron records changed (decisions of 2026-10-31; swarm re-validated, RESULTS.txt)
- e + H2O -> H2O(J=2), H2O(J=3) rotational: MuroranIT ends at 5 eV; continued to 1 keV with the J=0->1 energy
  dependence of the same set scaled at 5 eV (model_derived, factor 3). v3/v3.1 bounded 5 eV - 1 keV by the Munoz 2007
  rotational integral (211x and 1040x the continued values at 5 eV), which gave B/A of 5.5 %, 19 % and 36 % in breakdown
  cells at T_e 2, 5 and 10 eV.
- e + H2O -> e + H+ + OH- and e + OH+ + H- (ion pairs): continued above their threshold windows (18.70, 18.12 eV) with
  the Munoz 2007 total-inelastic energy dependence scaled at the window end, then the RBEB shape above 10 keV
  (model_derived, factor 10).
- e + H2O -> H- + OH, OH- + H, O- + H2 (dissociative attachment): zero above the end of the MuroranIT resonance tables
  (18.92, 14.31, 17.96 eV), model_derived: attachment proceeds through the H2O- resonances (2B1, 2A1, 2B2 near 6.5,
  8.6, 11.8 eV); above the last band there is no resonant capture channel, and anion production proceeds through the
  ion-pair records. v3/v3.1 bounded these five records by the total inelastic cross section (B/A 2 % at T_e 5 eV, 19 %
  at 10 eV).
- e + H2O -> 2e + O+ + H2, 2e + H2+ + O, 3e + O++ + H2: zero (model_derived) between the thermochemical threshold and
  the MuroranIT appearance energy (19.0, 20.7, 80 eV), where the swarm-fitted tables start at zero; v3/v3.1 bounded these
  onset slivers by the total inelastic cross section (B/A 3.7 % at T_e 20 eV). Sampling is unchanged (zero either way).

### 12.5 Coverage policy and remaining refusal regions
`data/tables/chemistry-v3.2-coverage-policy.json` (schema pfc-chem-v3.2-coverage-policy-1, bound to the blob):
162 regions; 0 refusing unbounded; 126 admitted unbounded (above 10 MeV, or placeholders only there: the v3.1 no-data
placeholders and the new processes below 1 keV); 36 bounded.

Ion beams in water vapour (tools/ion_beam_share.py; B/A of the bounded region at the beam energy; "adm" an admitted
placeholder-only region, never sampled):
| ion | 3-10 eV | 30-999 eV | 1 keV - 1 MeV | 3 MeV | 9 MeV |
|---|---|---|---|---|---|
| H+ | 0 | adm | 0 (1.4e-6 at 400 keV, 3e-6 at 1 MeV) | 6.6e-5 | 1.5e-4 |
| H2+ | adm | adm | 0 | 0 | 2e-5 |
| O+, OH+, H2O+, H3O+, O++, H-, O-, OH- | 0 or adm | adm | 0 | 0 | 0 |
(v3.1: 2-18 at 3 eV - 1 keV and unbounded above 1 keV.) The remaining bounded ion regions are H+ capture above 2.1 MeV,
H+ double capture above 300 keV and H2+ capture above 4.2 MeV (edge values held).

Representative breakdown mix (tools/breakdown_mix_share.py): Maxwellian electrons at T_e 2, 5, 10, 20 eV on
2.5e25 m^-3 H2O with n_e 1e21 m^-3; 300 K ions H3O+ 50 %, H2O+ 30 %, OH+ 10 %, H+ 5 %, O+ 4 %, H2+ and O++ 0.5 % each;
anions 10 % of n_e (OH- 60 %, O- 20 %, H- 20 %); plus either epithermal drifting ions (a fraction 0.1 or 1 of every
species at 3, 10, 30, 100, 300, 1000 eV) or keV ions (1e-3, 1e-2, 0.1 at 10, 50, 100 keV). 88 cases: none refuses;
largest B/A 4.7e-6, 7.5e-5, 5.5e-4, 5.7e-5 at T_e 2, 5, 10, 20 eV (electron onset regions: 3A2/1A2 dissociation and
the OH+ + H- ion-pair window). With the v3.1 tables all 88 cases refuse (largest B/A 0.87).
Remaining refusal regions: none refuses unbounded. Bounded regions can still refuse where exposure concentrates in them
(B_r > max(1, 0.01 A)): electron onset slivers that keep a bound (H2O+ 12.62-12.65 eV, ion-pair windows), electrons above
10 keV (rotational and vibrational records above 1 keV, attachment onset bounds), H+ above 2.1 MeV and H2+ above 4.2 MeV,
e + H2O+ / H3O+ recombination above 19 eV, H- + H associative detachment above E_cm 0.1 eV. None reaches 1e-3 of the
events in the mixes above.

### 12.6 Decisions taken (engineering and data handling)
- Rudd parameters: the water-vapour set of Rudd 1992 Table V with the universal inner-shell set for 2a1 (B > 2 B_min)
  in the kernel and the builder's quadratures (one set of numbers, checked by the builder against the header).
- Rudd sampling: power envelope (M h(w)) or exponential envelope (h_sup exp(-a (w - wc))) per orbital, whichever has the
  smaller mass; exact for either; removed the exhausted 160-trial draws of slow projectiles.
- Gobet Table II signs read from the page image (the text layer renders minus as 2); O2+ row omitted; shares held at
  the 20 and 150 keV values outside (in equal-speed proton energy for heavier projectiles).
- Wedlund sigma11 validity: Table 4's 400-30000 km/s (the text says 400-10000 km/s); recorded.
- Wedlund fits evaluated below their stated speed only where the source's appendix extrapolation is reproduced (the
  fits reproduce those rows within 3 %); above it, falling forms (scaled proton sigma10 or sigma-10 shapes, v^-2).
- Thermal records of OH+, H2+ and O+ charge transfer are zero from 1 keV by construction; the keV reaction is a separate
  forward charge-transfer record (CHARGE_X) zero below 1 keV, so kinematics follow the regime; H2O+ symmetric charge
  exchange (already CHARGE_X) is continued in place.
- Proton transfer of OH+ and H2+ above 1 keV: the thermal rate / v as for H2O+ and H- (same rule as the
  decision for H2O+ and H-; it brought every ion below 1 % B/A).
- H- double detachment below 650 km/s from 1 keV: the fit as the source's appendix extrapolation (replaces an edge hold
  that gave B/A 4.7 % at 1 keV).
- Process ids of placeholder-capable records carry no ',' or '|' (region-name separators; enforced by the builder).
- Radii: Bondi radius of the heaviest atom (O 1.52 A, H 1.20 A) for projectile and target.
- Gap bridges below 1 keV: log-log in E between the last low-energy value and the upper-end value; constant rates where
  the low-energy form already is k/v; the H- threshold rise as a power law in (E - E_th) with the table's first slope.
- Electron onset slivers of O+, H2+, O++ ionization: zero to the MuroranIT appearance energy (sampling unchanged).
- Swarm re-validation (decision of 2026-10-31): the v3 agreement bands of section 8 were single-seed and partly
  interpolation artefacts (repeated seeds scatter 0.5-3 %; interpolation between 300 and 500 Td shifted the bulk drift at
  396-497 Td by about +5 points). v3's records are corrected to six-seed means at the measured E/N, and each v3.2 mean
  must lie within the v3 multi-seed mean +- max(2 SD, 1 percentage point) (SWARM-VALIDATION.txt, both tables).
- Cost per collision against v3.1 and the swarm re-validation: RESULTS.txt.
