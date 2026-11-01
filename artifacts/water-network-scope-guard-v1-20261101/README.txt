Minimum water encounter admission -- isolated host reference
==========================================================
This opt-in scope gate refuses missing encounter groups before any chemistry
candidate/debit/RNG publication. It does not change cross sections, infer an EEDF,
or qualify a full reaction network. Current native runtimes and v3.2.2 are untouched.

The explicit minimum required-pair set is the actual table's 20 groups plus
e + OH, e + O and e + H2. A contract binds the table bytes, table species registry,
native registry, configuration and exact required pairs. Its canonical text and
SHA must be supplied from a separately reviewed native configuration.
Contracts written by the tests use SYNTHETIC native/configuration identities;
they must never be installed as an actual native/physical binding.

The all-species census precedes table-derived projectile/target filtering. Each
row carries stable native ID, collision-cell ID, table species ID and nonnegative
physical count. Original native row count and pre-filter scope are checked. Zero
count/absent target is an inactive encounter, never a statement that the missing
cross section is zero. A positive subnormal count is still positive. Different
cells do not produce an encounter. Every represented species must be registered.

The owner computes/retains census and stage bindings from its immutable input
snapshot, including state/registry/configuration identities, step/stage/duration
and cell volume. It must not re-sign inputs after tampering. Guard::check validates
those retained bindings. A caller-provided census alone cannot prove that its
owner supplied the complete native state; the native adapter must bind native
row count, full state digest and registration mapping from the actual owner.

An absent required table group refuses when both species have positive count in
the same collision cell. Present groups whose whole tabulated coverage domain is
unbounded also refuse; this retains e + OH+ and H3O+ + OH- refusal. A present group
with resolved/bounded intervals passes only this minimum preflight. Existing
energy/thermal-support region admission remains mandatory, including known-zero
thresholds, finite B/A limits, below-grid and above-grid refusal. This component
does not evaluate rates, thermal support, or zero cross sections at an energy.

Admission.network_complete and physical_model_validated always remain false.
Legacy parent/diagnostic modes bypass the opt-in gate and provide no full-network
admission. Unknown modes refuse. No legacy mode or policy semantics are changed.

Native adoption requirements (NOT implemented/qualified here):
1. Build this census from all registered rows before target flags discard OH/O/H2.
2. Bind the actual native registry/configuration and immutable stage snapshot.
3. Call the gate before host or device MCC candidate/debit/RNG publication.
4. Keep every existing region/energy/target-owner/commit/refusal guard.
5. Implement a GPU-resident compact all-species census and qualify it against this
   host reference; do not silently transfer all particles to host each step.
6. Bind the new explicit contract/configuration into checkpoint and restart identity.

Tests run only through one-core frozen SimQueue jobs, with bounded memory. They
exercise the actual pinned v3.2.2 blob and manufactured census records; no physical
simulation or swarm matrix is run. Source and input hashes are in INPUTS.json.
