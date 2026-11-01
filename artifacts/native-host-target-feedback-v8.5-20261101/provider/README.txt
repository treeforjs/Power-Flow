CPU native/gas target feedback reference

This provider feeds an actual open gas v1.4 candidate's neutral moment correction
into native v8.5 CellTargets before native hashes inputs, computes coverage and
majorants, and samples target velocities. It changes only eligible neutral drift
and thermal width. Native preserves density, species, order, atlas and all other
bytes. All-species water scope admission occurs first.

Construct in this order, using public factories:
1. Native chemistry configuration, optionally its minimum-water contract.
2. Native registry and matching gas configuration, with feedback explicitly true.
3. shared_ptr<gas::Owner>, then native host_target_feedback=true and the contract()
   result in host_target_contract_sha256.
4. Rebuild the final native registry/genesis with those fields.
5. For each actual native collision scope, create a fresh Provider and call
   prepare_chemistry_v3_host_targets(provider). No provider may be reused.
6. Before gas.take, require the live native scope to accept the token and call
   token.validate_candidate_host_targets(). Provider retains authentic immutable
   input identity, shared gas lifetime, private CandidateBinding and gas content
   digests. Same-digest retry and no-op mutation invalidate the capability.

The gas.take mutation intentionally invalidates the old feedback capability.
A common adapter must verify the exact admitted/debit/release records and new
post-stage ledgers, then establish a fresh binding for the next stage. These
interfaces never authorize publication. Every owner needs a separately prepared
non-failing publication capability and a common rollback policy.

This component is CPU-only. Native device attachment/dispatch refuses when
host feedback is configured. No hidden particle download or fallback is used.
Both normal and corrected targets are a Maxwellian moment model, not a measured
velocity distribution. Full-system energy closure remains false and the native
strict eventful target-energy-owner guard remains closed. Diagnostic tests
explicitly opt into the existing field-target approximation; no production
policy or unknown-reaction guard was relaxed.

The tests compare 1686 actual diagnostic reaction records and exact ledger
operands with zero correction, exercise a separately manufactured positive
moment perturbation through stage1 take and actual stage2 target injection,
then abandon that manufactured candidate. The manufactured perturbation is not
an assertion of a matching native debit or physical closure. Cold populated
continuations cover accepted history/restore before and at anchors. Public
scope/gas/feedback factory construction is checked without manual identity edits.

The external provider links the separately pinned native_gas_stage Binding
header to verify actual blob/species/process/grid operands. See DEPENDENCIES.json.
No measurements, production PIC run, measured EEDF, benchmark timing or complete
water network is qualified by this reference component.
