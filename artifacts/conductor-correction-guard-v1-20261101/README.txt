Per-exit conductor correction refusal, version 1

Problem
The first common-step draft checks whether each species' cumulative conductor
correction total changed. Signed corrections at different actual impact locations
can cancel exactly in that total. A zero total does not mean zero unrouted local
energy, momentum or internal-energy transfers.

Minimal integration
  #include "conductor_correction_guard.hpp"
  auto actual = coupled.motion_transport_receipt();
  require(bool(actual), "missing actual transport receipt");
  auto evidence = powerflow::conductor_correction_guard_v1::
      require_no_unrouted_conductor_correction(
          *actual, common.accepted_step, gas.configuration());

The planned native/common v1.4 getter returns the immutable actual gas
TransportReceipt retained during authentic native motion acceptance. Version1.3
only preserved aggregate proof and cannot reconstruct the missing individual
shares. Neither frozen bridge nor root common_step is modified here.

The root calls this predicate during preparation, before seal/publication, and
continues to preflight the live common native/gas capability and Hall/source
capabilities immediately before no-fail publication. The helper is a numerical
predicate, not publication authority. An arbitrary caller-constructed receipt
cannot establish native provenance, coverage or ownership. In particular, a
caller must not delete exits and claim that an empty vector authenticates motion.

Contract
- Require the expected transition index and nonempty canonical transport input,
  candidate-before and candidate-after SHA256 values.
- Require host census_conservative configuration and the existing gas carrier/
  exit-byte bounds. No new allocation proportional to the receipt is performed.
- Validate canonical ascending native IDs, registered neutral species, recognized
  conductor/exterior destinations and normalized base-2^32 exact digits.
- For every conductor exit, each of all five exact quantities must be zero:
  kinetic-energy correction, three momentum corrections, internal correction.
- No conversion to binary64, sum over exits, tolerance, clipping, or thermalization.
  Opposite signed shares and nonzero exact shares that round to zero still refuse.
- Nonzero exterior corrections remain admitted under the existing separate
  exterior ownership contract. This does not add them to foil heat.
- Successful Evidence counts exits only. It neither changes the gas state nor
  authorizes publication; absence of a conductor correction is not proof of all
  other conservation or physical-model requirements.

Regression
The actual frozen gas v1.7.1 owner is seeded with equal/opposite correction stocks
in two cells, and actual carriers exit both cells at conductors. Its authentic
transport receipt contains two nonzero exact shares while every cumulative
conductor total is unchanged. The old predicate therefore passes; this predicate
refuses. Rejection leaves the candidate and accepted state untouched, and the
attempt can be abandoned back to the exact checkpoint. Restore/replay reproduces
the receipt. Actual zero-conductor and nonzero-exterior controls are separate.
Operand fault controls cover all five quantities, sub-binary64 exact corrections,
unknown/colliding identities, missing hashes, noncanonical digits, wrong step and
capacity/execution/transport mode mismatches. Modified receipt controls are
explicitly synthetic operand tests, not authentic native motion claims.

All normal and ASan/UBSan compilation/execution uses one bounded SimQueue owner
core and frozen userspace. INPUTS.json pins the helper, tests, gas implementation,
forms arithmetic and toolchain. Evidence is retained per attempt. This component
adds no routing, native integration, relaxed physics guard or production claim.
