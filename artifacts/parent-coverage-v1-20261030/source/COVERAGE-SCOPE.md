# Water parent-ionization channel: coverage scope v1

Pilot prerequisite H4 (coupled-pilot plan, section 0 item 6 and section 2).
Audit links: A13 (versioned supported channels; no invented rates, no
renormalization of positive missing branches) and A14 (1 keV model seam).

Channel: `H2O:ion_yield:H2O+` with kernel `parent_ground_equal_pair_CM_v1`,
dataset `h2o-gas-20261029-v1` (sha256 `840663cd…e4bf`), as bound by the sealed
v3 parent subflow. All energies below are the electron kinetic energy in the
H2O rest frame, computed by the sealed `pic_collisions::relative`.

## 1. Coverage map

| Target-rest energy E | Status in scope v1 | Source | Ledger sigma |
|---|---|---|---|
| 0 <= E <= 12.621 eV | known zero (not recorded) | v1 `zero_below_eV` 12.621 = parent ionization reference of the kernel | 0 |
| 12.621 < E < 13.5 eV | **unresolved: onset gap** | Song et al. 2021 Table 12 (Straub 1998 lineage) starts at 13.5 eV (2.5e-22 m^2); no parent-only value below | unknown, 0 <= sigma <= 4.665e-21 m^2 |
| 13.5 <= E <= 1000 eV | **qualified** | Song 2021 Table 12 parent H2O+ (6 %, D2O-based recommendation) + conditional kernel `parent_ground_equal_pair_CM_v1` (data support 13.5-1000 eV) | applied |
| 1000 < E <= 10000 eV | **unresolved: above parent table** | v1 parent table ends at 1000 eV ("above1000eV" listed missing) | unknown, 0 <= sigma <= inelastic envelope (8.35e-21 m^2 at 1-1.5 keV, 1.75e-21 m^2 at 7.5-10 keV) |
| E > 10000 eV | **unresolved: no data in any version** | every water dataset (v1-v8) ends at or below 10 keV; NIST BEB nonrelativistic ceiling | unknown, no upper bound (refused by default) |

Upper bound: `H2O:munoz2007:recommended_inelastic` of the same v1 dataset
(Munoz et al. 2007, Table II total x inelastic fraction, 7.5-10000 eV). Parent
ionization is one inelastic process, so sigma_parent <= sigma_inelastic. Between
grid points the larger of the two bracketing tabulated values is used; the
recommended-data uncertainty (5-18 %) is not added. The scope refuses to
construct unless this envelope dominates the qualified parent table at every
overlapping table point. In the onset gap the envelope (4.665e-21 m^2) is about
19 times the first tabulated parent value, so gap exclusions are charged
conservatively.

Notes on specific energies:
- 13 eV lies in the onset gap.
- 38.3 eV is supported, but the equal-pair kernel gives each product electron
  (E - 12.621)/2 in the CM frame; incidents in 37.86 < E < 39.62 eV put both
  products into the onset gap (38.3 eV gives two electrons at 12.838 eV).
- 8000 eV lies in the 1-10 keV unresolved range. The plan's "> 8 keV" figure
  understates the gap: for the bound channel everything above 1000 eV is
  unresolved. 8.084 keV is the retained electron tail of the earlier
  electron-only run, and 8000 eV is a regression test energy.

### What existing data could close the gaps

- Onset gap 12.621-13.5 eV: nothing qualifies. Song total ionization also starts
  at 13.5 eV. `H2O:nist_BEB_valence` (v1) evaluates in the gap (1.9e-22 m^2 at
  13 eV) but is a four-orbital total valence model; the dataset documentation
  forbids total-to-parent substitution and the kernel header refuses BEB/total
  curves in this gap. Left unresolved.
- v8 H(1s) continuation (`H-1s-ionization-20261029.json`): atomic-H target,
  debits H(1s) not H2O, support 14-10000 eV with its own 13.59844-14 eV gap.
  Provides nothing for the water parent channel.
- v7.1 metadata (`metadata-v7.1/handoff.json`): metadata-only normalization of
  the v7 deep-valence release (commit 0f57b04). Its ten conditional models
  include the 2a1 Auger subset (1-10 keV) and carry the v5 outer3 models
  unchanged. Nothing below 1 keV; nothing for the onset gap.
- 1-10 keV: v5 `outer3_BEB_Tan_parent_v1` (dataset `h2o-gas-20261029-outer3`,
  c436020d…, identical in v7 bac4e1e0…) is an existing conditional exclusive
  parent-branch rate and kernel: NIST BEB orbitals 1b1 + 3a1 + 0.08 x 1b2 with
  Tan/Murakami branching. It is NOT admitted in scope v1 because
  (a) it is a different dataset, channel id, model and kernel identity from the
  natively registered parent channel (one binding per channel; admitting it
  needs a second registered channel or a rebinding and a header that evaluates
  `orbital_weights`, which the bound v1 `collision_data.h` BEB ignores), and
  (b) switching at 1 keV is a model seam: Song parent 4.32e-21 m^2 versus outer3
  parent 5.476e-21 m^2 at 1000 eV (+26.8 %), which the data documentation
  requires a dispatcher to select and report. Candidate for coverage scope v2.
- Above 10 keV: no rate or bound in any version.

## 2. Exclusion mechanism

The candidate query keeps the sealed bridge sigma unchanged for every energy
classified known-zero or qualified (same calls, same values, same errors). For an
unresolved energy the query returns no candidate for that pair, and a second
pass over the same population, with the same pairing order and the same sealed
kinematics, records the pair. The number of exclusions seen by the query must
equal the number recorded, and both passes must visit the same pairs, or the
exposure aborts. Other pairs, rates and branches are untouched: no
renormalization, no invented rate, no unrecorded zero.

The scope is bound to the sealed bridge by probing 80 energies (all bin
edges, neighbours one ulp either side, bin midpoints and named energies): known
zero must return 0, qualified must return the same sigma, unresolved must raise
the strict `CollisionCoverageError`.

## 3. Unresolved-coverage ledger

Per exposure (one SSA exposure; a full-time step has two):
- `pair_evaluations`: exact count of excluded (electron carrier, H2O carrier)
  pairs over all constant-population segments of the SSA loop, split by sigma
  status (tabulated, upper-bounded, unbounded);
- `distinct_excluded_pairs`, `distinct_excluded_electron_carriers`,
  `distinct_excluded_water_carriers`, with exact weight sums at first exclusion;
- exact sums over pair evaluations of electron, water and event weight;
- `physical_exposure_inv_m2` = sum w_e w_n F dt / V and `event_exposure_inv_m2`
  = sum max(w_e,w_n) F dt / V (F Moller flux, V collision-cell volume, dt the
  segment time). Multiplied by a cross section they give expected physical
  reactions and expected scheduler events;
- `tabulated_expected_*` where sigma is tabulated (none in v1 data, kept for
  kernel-domain gaps), `upper_bound_*` (tabulated or envelope sigma), lower
  bound 0, and `unbounded_*_exposure` where no upper bound exists;
- `applied`: expected physical reactions and events of the qualified candidates
  (the limit denominator);
- energy bins, `(lower, upper]` in eV with edges 0, 7.5, 12.621, 12.75, 13, 13.25,
  13.5, 1000, 1500, 2000, 3000, 5000, 8000, 10000, 20000, 50000, 1e5, 1e6, inf; each
  bin record carries the same measures, the coverage class, sigma status and the
  observed energy range;
- the exposure and cumulative limit decisions.

Exactness: counts are integers; weights are summed exactly; each exposure term
is evaluated once in extended precision and rounded once to binary64, and the
rounded terms are summed exactly. Sums use a fixed-point accumulator in units of
2^-1074 (order independent) and serialize as an odd hexadecimal mantissa and a
binary exponent plus the correctly rounded value.

State: the accepted state holds the last step's exposures and their totals, and
the cumulative totals and bins. The ledger is part of the RNG/clock payload, so
its hash enters the native stochastic completion, the full-time boundary RNG
seal, `checkpoint()`/`restore()` and `seal()`. A restored ledger must be
canonical, belong to the same scope and policy, and hash to the native accepted
completion. A full-time first exposure stages its ledger as `pending`; the
second exposure closes the step; an accepted state never carries a pending
ledger. `PreparedParentState::coverage()` returns the sealed step ledger for the
coordinator's per-step ledgers; `ParentSubflow::accepted_coverage()` returns the
accepted one.

## 4. Hard limit

Evaluated for each exposure and for the accepted history plus the candidate.
With U = upper bound of unresolved expected physical reactions, A = applied
expected physical reactions, U_ev = upper bound of unresolved expected events:

    refuse  iff  (unbounded pairs present and not admitted)
             or  (U_ev > event limit  AND  U (1 - f) > f A)

The refusal is `coverage_v1::CoverageRefusal` (derived from
`pic_collisions::Error`), carrying the criterion, the scope identity and the
refused exposure ledger. It is thrown inside the native compound or private
collision frame, so native population, both RNG streams and the ledger are
restored.

| Policy field | Default | Rationale |
|---|---|---|
| `max_exposure_unresolved_events` | 0.01 | Unresolved and qualified hazards are independent Poisson processes, so U_ev bounds the probability that the excluded pairs would have produced any event in the exposure: with at least 99 % probability the executed exposure is a realization of the full hazard. |
| `max_unresolved_reaction_fraction` | 0.01 | At most 1 % of the exposure's parent ionizations (upper bound) unresolved: below one sixth of the 6 % uncertainty of the applied Song parent table, so coverage loss stays subdominant to the applied cross-section uncertainty. |
| `max_cumulative_unresolved_events` | 1.0 | Prevents small per-exposure exclusions from accumulating unchecked: once the history's upper bound exceeds one expected event, the cumulative unresolved fraction must stay at or below 1 %. |
| `admit_unbounded` | false | Above 10 keV no data bound exists, so the fraction cannot be bounded. Setting it true records the sigma-free exposure as unbounded and lets the step proceed. |

Either criterion alone is a sufficient acceptance guarantee, so refusal needs
both. The policy is part of the scope identity, hence of the subflow
configuration identity and the genesis RNG identity; it cannot change across a
restart.

## 5. Parity with v3

With no unresolved pair, the candidate list, scheduler draws, product RNG,
event identities, product carriers and velocities are those of v3. The sealed v3
parent regression, ported by `make_v3_port.py` (16 line-local edits; check names,
order and line count unchanged) and run with the default policy, gives 442
identical check results, identical recorded event physics (event records
without the native audit) and 73/128 first events. The 13 eV, 8000 eV and 38.3 eV
cases still refuse under the default policy, now through `CoverageRefusal`.
What differs from v3 by design: the owner id
(`native_water_parent_subflow/coverage-v1`), checkpoint schema, configuration
identity (adds the scope) and the RNG/clock state hash (adds the ledger; legacy
mode now also binds the configuration). Native audits therefore differ in their
completion hashes.

## 6. Coordinator interface changes

- `coverage_v1::ParentSubflow(arena, grid_counts, binding, limits, bridge, scope,
  clock_seed, product_seed, full_time)`; `CoverageScope(data, dataset_sha256,
  bridge, policy)`.
- `initial_full_time_rng_sha256(...)` and `full_time_configuration_sha256(...)`
  take the scope; genesis must use them.
- Catch `CoverageRefusal` as a typed scoped refusal of the joint step.
- Add `PreparedParentState::coverage()` to the per-step ledgers.

## 7. Files

- `source/exact_sum.hpp`, `source/parent_coverage_scope.hpp`,
  `source/parent_coverage_subflow.hpp` (includes the sealed v3
  `native_parent_subflow.hpp` unchanged).
- `tests/test_parent_coverage.cpp`; `tests/test_v3_parent_regression_port.cpp`
  (generated by `make_v3_port.py` from the sealed test pinned in `V3-TEST-SHA256`).
- `run_coverage_checks.py`: pinned build/run of both tests, normal and UBSan,
  with the regression comparison; receipts in
  `outputs/pic-water-parent-coverage-v1-20261030/attempt-NNN/`.
