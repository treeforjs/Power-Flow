# Melt-ablation and layer proposal v1

This is a stateless, single-local-patch proposal evaluator for a future common
owner. It does not advance production or create accepted particles, source
inventories, thermal state, geometry, fields, clocks or random state. The opaque
`Proposal` exposes read-only candidate values and evidence. Its copies are
proposals, never independently evolving accepted reservoirs. Full integration
remains held.

`Evaluator` takes an explicit `Scenario`. `evaluate` reads a const common-owner
`View` and a caller-planned `Request`. It returns all candidate changes or throws;
even a late refusal leaves the supplied view, weights and earlier actions
unpublished. There is no stateful temperature cache, publisher, RNG, native
allocator or separate gas/count backend. `matches_base` checks freshness only.
Neither JSON evidence nor that boolean authorizes common publication.

## Required explicit inputs

The scenario names the mesh, fixed geometry, common reference convention,
qualification, material/reference-law hashes and bottom-to-top condensed stack.
Each material declares its entity mass, identity-preserving elemental inventory,
neutral product state and immutable product R/I reference. A condensed stock
also supplies its caloric callback/domain, fixed reference density, condensed
reference, trigger temperature and positive incremental removal cost excluding
actual K. A condensed trigger temperature is the supported trigger definition
in this version; a liquid-fraction or enthalpy rule must be converted under an
explicit admitted material law or supplied by another version.

All condensed layers in one view share one local temperature. Their caloric
energies add as `sum_i m_i*u_i(T)`. The borrowed thermal energy excludes the
separately derived bound-reference energy `sum_i m_i*R_i + sum_j U_ads_j(N_j)`.
An existing owner whose energy already includes binding/reference or mechanics
must provide a versioned, conservation-checked projection into this convention;
feeding that total again would double count it. Adsorbate sensible capacity is
not independently represented in this reduction. Different-temperature layers,
finite contact conductance and adsorption heat capacity require another explicit
closure; they cannot be silently represented by the metal caloric law.

Callbacks must be pure, immutable and admitted by the common owner under the
supplied hashes. Copying `std::function` does not prove the identity or immutability
of external captured state. Configuration hashes identify declared inputs; they
are not a trusted source factory, physics validation or native-registry proof.
The provided SS304 factory requires the existing seal and explicit qualification
and reference density. It calls only the unchanged condensed caloric model; it
does not extend the separate vapor sampler below 1800 K or any existing guard.

Water/adsorbates require their own bound-reference callback, support material,
source range and independent-desorption policy. No metal-melt gate is accepted
for them. No rate law is implemented or selected. The caller supplies actual
planned weights and velocities from its admitted source process; this component
does not reuse H2O activation energies or rates for metal or oxide.

## Local budget and references

The view is one local cell/patch's common-owner pre-removal trial. Its thermal
energy must already contain signed conduction, emission, desorption, Joule and
mechanical transfers admitted elsewhere, each exactly once. A lifetime deposited
energy or whole-foil headroom cannot replace this still-available local energy.
The view pins owner/trial/field/source hashes, revision, consumed event IDs and
the authoritative `(begin,end,dt)` clock. `end==begin+dt` is checked; `dt` is not
replaced by endpoint subtraction at a nonzero epoch.

For an exposed condensed source at its explicit trigger `T*`:

```
threshold = sum_i m_i*u_i(T*)
available = E_trial - threshold
dm = sum_p w_p*m_entity
packet = sum_p w_p*(R_product + I_product + actual_K_p)
delta_reference = reference_after - reference_before
debit = packet + delta_reference + named_other_out
incremental_cost = debit - dm*u_source(T*)
E_after = E_trial - debit
```

The evaluator requires `0 < incremental_cost <= available`, finite source
inventory and an admitted remaining thermal state. The product/reference inputs
must agree with the declared incremental removal cost. Their input consistency
comparison permits only eight machine epsilons of the sum of represented
operands; the actual energy debit and conservation ledger use exact represented
R/I/K values without that tolerance. Fusion and sensible energy already in the
sealed caloric threshold are not separately charged again. `named_other_out` is
nonnegative and has a named receiver; signed exchange must be handled explicitly
by the common pre-removal trial in this version.

The returned budget refusal quotes an exact rational ratio of available to
requested incremental cost. It does not rescale, floor, inflate or resample any
packet weight. At melt with zero excess, positive removal is refused. Actual K
is checked against the proposed relativistic proper velocity, and its actual
weighted value participates in the budget. A stationary-support momentum quote
uses `-sum(w*m*u)` with zero support work. This is not a moving-wall impulse/work
or a native mechanics receipt.

All extensive arithmetic reuses the sealed dyadic BIGNUM currency. A single
external condensed stock is encoded as mass, with zero count. A single external
adsorbate stock is encoded as the existing exact count, with zero mass. This
version requires its source count to export losslessly into the inherited
Count2-compatible source interface before and after debit. A third residual is
refused; exact gas weights do not repair a rounded source pool. It does not
install a second accepted Count2 mirror, change the existing exact backend or
qualify a general exact-source adapter.

The temperature is a derived nearest binary64 projection after exact domain
checks. Negative energy relative to the 300 K reference is valid inside the
admitted cold range. An exact sub-ULP overrun of a hot endpoint is refused. No
energy, mass, temperature or source clamp repairs a failed proposal.

## Retained inventory and exposure

Only the outermost nonzero condensed layer can be removed. Water may leave its
exposed support below metal melt; its oxide support remains. Partial oxide loss
retains oxide exposure. Complete oxide loss can expose the next retained material
only after attached adsorbates are explicitly depleted. No silent entrainment,
deletion, rebinding or permeability through a covering layer is supplied. Complete
substrate depletion is refused pending a geometry/empty-cell handoff. A returned
`dm/rho_reference` ratio is a volume quote, not an accepted moving boundary or a
recession-limit check.

There is no condensation, reaction network, changing composition, charge birth,
ionization, electron emission, rate selection, droplet distribution or complete
metal-plasma kinetics. This neutral, identity-preserving product restriction
deliberately refuses unsupported charged or noncongruent products.

## Common-owner integration obligations

Before any integration, the common owner must choose the physically admitted
layer amounts, order, composition, caloric/reference laws, trigger/removal closure
and product registry. Every unresolved foil input remains null in the successor
publication. The test's macroscopic entity masses, oxide density, heat capacity,
trigger and removal energy are manufactured arithmetic inputs, never foil values.

The owner must obtain the actual immutable native products and source release
law/RNG candidates, verify exact native weight/mass/element/charge/energy and
momentum against this plan, and bind the real field, geometry, thermal and source
states at final preflight. It must enforce reference-volume/recession limits,
transport/current and exposure changes, existing mechanics/source/material guards,
and any total energy terms outside this reduction. It must consume event IDs and
publish particles, inventories, energy, fields, geometry, clocks and RNG together
through its existing no-throw common commit. It must own restart serialization
and test rollback/restart with those real objects. No such joint publication,
native birth, restart proof or physical validation is claimed here.

## Scope of verification

Focused manufactured checks cover below-melt water release, retained oxide,
ordered exposure, positive-leaf exact debit, zero and insufficient melt budgets,
actual K, named other energy, sequential depletion of available energy, cold/hot
domain refusal, Count2 residual refusal and malformed/stale inputs. Compile-time
checks prohibit public default construction and mutable proposal access. A normal
and UBSan build run the same cases, then independent Python Fraction arithmetic
and 90-digit Decimal kinetic-energy checks reconstruct their receipts. This
checks component arithmetic and refusal behavior. It does not validate physical
layer values, source kinetics, spatial/time convergence, a native coupled solve
or production throughput.

The earlier audit and appendix remain immutable. Publication v2 uses the phrase
“The requested oxide/intermediate-layer model” and preserves its energy-reference
and exposure distinctions. The new implementation is separately qualified from
that wording-only publication successor.
