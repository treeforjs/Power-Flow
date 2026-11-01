Conservative full-census gas correction transport, CPU reference

Isolated runtime:
/home/treefor/projects/runtimes/gas-target-owner-v1.7.1-census-transport-20261101

Problem and correction

The v1.6 carrier_shares transport uses a rounded cell correction multiplied by
a rounded particle-weight fraction. Even when every particle leaves a cell,
these shares can leave exact nonzero energy, momentum or internal-energy digits
in the empty origin. SimQueue 563 reproduces the actual v1.6 case with weights
2, 6 and 4; global conservation alone did not catch the stranded correction.

New Transport::census_conservative = 2 requires the full before/after census and
all sunk native IDs, including charged particles, with an explicit conductor
or exterior destination. The census partition preserves species, weight and
native ID; every row is validated before neutral filtering. Actual post-motion
occupancy is decisive. A tiny resident whose weight disappears from a rounded
count, or a new incoming resident, prevents an evacuated-cell drain.

The method first performs the historical rounded shares in canonical ID order
using unchanged pre-motion rounded stock. It then sends all remaining exact
digits of each truly evacuated source cell/species to the actual destination
of that source's highest-ID outgoing carrier. Signed corrections are preserved;
there is no clipping, erased residual or fabricated carrier. A preexisting
nonzero correction with no actual before-census carrier refuses explicitly.

Exit correction stock is separated from the native incident kinetic/internal
energy flow. Appended exact conductor/exterior ledgers preserve E, Px, Py, Pz
and U by species, alongside the historical aggregate wall ledger. The returned
per-exit exact digits permit attribution to the authentic native surface record.
An exterior correction must not become conductor heating. Counts and native
birth/sink flows remain owned by the existing flow and native operators.

The model is still carrier_shares cell mixing. Newly born same-species carriers
in the before-motion census immediately share the old cell correction by their
weight, even though they did not experience its original collision. This is an
explicit coarse cell transport assumption, not particle-attached recoil history.
Exact conservation does not establish Maxwellian feedback realizability or the
validity of that mixing approximation. All existing feedback and chemistry
guards remain in force.

API and identity

Owner::transport_census(step, before, after, exits) -> TransportReceipt.
ExitRecord contains native_id and ExitDestination::{conductor,exterior}.
The receipt contains canonical input, candidate-before and candidate-after SHA,
number of movers/drained cells, and neutral ExactExitShare records in native ID
order. Each exact share has five arrays of 134 normalized signed digits.
These checked records are not native publication authority: the bridge must
obtain them from and retain its authentic prepared-motion capability.

Default maximum census cardinality is 1,000,000 per input. The exact exit payload
is bounded to 64 MiB before per-exit arrays are allocated; all declared exits
count toward this conservative bound, including charged ones. Limits are bound
into configuration identity, require finite positive limits, and refuse excess
without truncation. The exit byte cap may be explicitly configured up to 1 GiB.

Computation and allocation use a private State copy. Input, capacity, exactness
or allocation refusal leaves candidate, accepted state and its read binding
unchanged. Only complete success invalidates the old candidate binding and
publishes the private candidate within the still-open overall step. This late
invalidation behavior is explicit for this new API. Normal step abandon and
prepared common publication continue to own final commit/rollback.

Legacy carrier_shares remains available under its old mode and rounded rule;
the new mode refuses the old MoverRecord-only method. New mode construction
explicitly refuses device and device_verified execution; there is no fallback.
Existing global offsets are preserved and destination globals appended. The
configuration identity is versioned v1.7 and checkpoint codec is 2. Old v1.6
checkpoints/configurations do not silently acquire the new interpretation.

Qualification

CPU 563 and ASan/UBSan 564 passed:
- 135 new controls, including exact closure of all five quantities and both
  exit destinations; mixed births; unequal non-dyadic and split weights;
  actual tiny and incoming residents; permutation independence; malformed
  partitions/species/weights/phase space and capacity refusals; orphan stock;
  checkpoint integrity, restart and complete step rollback.
- 415 inherited gas-owner transaction, exactness, feedback and release guards.
- 196 real allocation-failure points each preserving candidate/accepted state
  and binding, then producing an identical retry.
- Byte-identical repeated controls. CPU 563 separately reproduces v1.6 residue.

558 refused missing explicit thread environment before computation. 560 passed
135 new controls and 415 inherited checks but its preimage fixture failed the
compiler's misleading-indentation guard. v1.7.1 only corrected that fixture's
formatting; numerical source remained unchanged. 561 was the first complete CPU
pass. Allocation sweep/tool were then added separately; 563/564 qualify the
final complete source set. Original files, logs and receipts are preserved.

Reproduce only through SimQueue/frozen userspace, 1 allocated core, 3 GiB,
CUDA_VISIBLE_DEVICES empty, OMP/OPENBLAS/MKL_NUM_THREADS=1:
python3 -I tools/check_transport_v1_1.py normal UNIQUE_LABEL
python3 -I tools/check_transport_v1_1.py sanitized UNIQUE_LABEL
The retained tools/qualify.py is inherited history; it is not the new
census-conservative mode's qualification entry point.

Source identity:
49bc2a15f2090fe34058de5cff1108462d842f57712d22e47b27b7c85fe7127a
Exact file pins and inherited baseline identities are supplied separately.

production_qualified=false
physical_model_validated=false
device_qualified=false
feedback_realizability_qualified=false
native_motion_authority_integrated=false
The coordinator's authentic native mapper and full coupled physics qualification
remain separate work.
