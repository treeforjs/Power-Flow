# Water chemistry v3: scope, data, models and qualification

Component for the Power-Flow coupled PIC water-breakdown simulation: the electron, ion and neutral chemistry needed
for H2O-desorption breakdown, as a sealed, tested library that native v8 links later (INTEGRATION.txt). It closes
the pilot limitations L1 (parent ionization only) and L2 (ions do not react or collide) of
`artifacts/pf-pilot-1-20261031/LIMITATIONS.txt` as far as data allow, and the parent-coverage items of audit
A13/A14 (onset gap, 1 keV and 10 keV seams). Nothing here advances physical time in a coupled run.

Tables: `data/tables/chemistry-v3.blob` (22 species, 20 projectile-target groups, 81 processes, 159 coverage regions),
built by `tools/build_tables.py` from `data/tables/chemistry-v3-spec.json` (written by
`tools/make_production_spec.py`). Every number comes from a retained source file through a documented source type;
`data/tables/PROCESSES.txt` lists every process with its sources, tags, coverage intervals, interpolation errors,
junctions and notes; `data/tables/chemistry-v3-manifest.json` carries the same in machine form with the sha256 of
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

## 10. Open items (data not obtained or models to replace)

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
