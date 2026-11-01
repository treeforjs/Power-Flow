Native v8.2 audit correction (candidate008-auditfix)
=================================================
This is an isolated fork of candidate008. AUDIT-SNAPSHOT.json proves before,
after and copied file hashes agreed. Originals and sealed inputs are not edited.

Changes
- Accepted v8.2 snapshots bind the full last-step published sink data, including
  absent/present state, start position and internal energy. Checkpoint restore
  reconstructs all sink tallies and validates the accepted digest even at a
  history anchor. Missing/tampered metadata cannot retain the same boundary.
- Compatible v5 mode retains the candidate008 named-shared_ptr lifetime fix for
  device-only verification/refusal records. Qualification adds a device-only
  bad-record fallback and a valid verified continuation.
- Conservative neutral release delegates to chemistry v3.1.1's actual-parcel
  helper. Positive partial-weight packets preserve phase-space spread and
  exact count/momentum/kinetic stock. Packet allocation is capped before release.
- Chemistry state codec v3-parcels-support and v8.2-auditfix2 configuration/snapshot tags
  explicitly reject old v8 restart data. No safe migration of signed legacy
  reservoirs is implied. Compatible v5 identity logic remains unchanged.
- Whole-system target energy/momentum ownership remains UNIMPLEMENTED. Default
  chemistry refuses eventful stages with target_energy_owner. Explicit
  diagnostic_field_target_approximation=true is identity-bound, serialized and
  exposed in receipts/status; full_system_energy_closure is always false.
  This guard does not implement target heating or momentum feedback.
- Long valid coverage labels restore using bounded text and exact membership
  in the pinned blob; unknown labels remain refused.
- Restore initializes the accepted identity schema before journal replay, so
  retained v8.1 motion, birth and merge rows pass only their proper v8 guards.
  The v8.1 guards are retained; schema-v5 records do not acquire v8 privileges.
- Unbounded regions are not recommended for production. Physics coverage uses
  the corrected chemistry corpus and strict production policy.

Verification
Use audit-tests.py through SimQueue, one CPU and --mem3G. GPU mode needs a
shared GPU claim; no timing claims. Source and binary identities and actual
results are recorded in the final artifact, not inferred from older reports.
Tests include sink tampering before/at anchors, diagnostic-only target-energy
scope with strict refusal/no mutation/restore, parcel codec/thermal spread and
nonrealizable stock refusal, chemistry controls, v5 flows/work controls,
nonempty restore, device-only verification and birth/sink/chemistry GPU parity.
The copied historical RESULTS-008.txt does not qualify this fork.

Reachability correction: v8.2-auditfix2 and chemistry-state-codec-v3-parcels-support
bind per-region reachable support counts, including zero physical rate. Unknown
support is a coverage refusal independent of dt; explicit diagnostics retain X=B=0.
