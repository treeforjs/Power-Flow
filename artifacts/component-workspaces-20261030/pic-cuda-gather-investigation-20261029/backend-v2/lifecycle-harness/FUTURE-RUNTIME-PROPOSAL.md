This is a proposal for later review. No runtime mode or application integration
is implemented or authorized by this document.

The current CudaBatchedForce API requires a complete independent CPU E/B/u
reference on every call and applies its original finite, componentwise and
weighted gates before a transactional velocity commit. Its measured force
pipeline must include fresh reference generation. The actual fixed-midpoint
CPU6 probe measured 20.055 ms median for that reference alone. Stored offline
references are useful qualification evidence; they cannot stand in for changing
application state or establish an integrated acceleration result.

A future design could expose two explicit modes, with no implicit fallback or
silent switch. A qualification shadow mode would keep the current per-batch
independent original provider/Boris reference and every comparison. An optimized
runtime mode would run only after an explicit qualification decision, using the
shared typed host/device math and requiring a separately reviewed API contract.
Selecting that mode would explicitly change the additional shadow-comparison
policy; it would not remove any pre-existing CPU physics or safety guard.

Both proposed modes must retain every existing scientific and resource check:
the original midpoint timing and two-drift/two-solve order; original native
frame interpolation, unit-psi/axial correction, physical charge-to-mass ratio,
proper-velocity Boris and macro-weight semantics; unchanged particle sampling,
source/reservoir accounting, boundary losses and charge deposition; original
Poisson factorization/finite/residual requirements; original gyroangle,
displacement, plasma-frequency and timestep limits; and the original charge,
Gauss/continuity, field/kinetic-energy and ledger accounting. Their tolerances
must remain at the existing values. Finite position/u/field/phi/parameter checks,
supported geometry and domain checks, positive mass/weight, integer/count and
allocation/capacity budgets, and ownership/context/thread checks remain
mandatory. An invalid batch must fail before partial particle mutation.

The optimized API would need typed snapshots tied to the exact live midpoint,
with explicit geometry/layout/unit-field/native-frame epochs, nodal ordering,
frame IDs, interpolation fraction, voltage, dt, and captured particle identity.
Stale static epochs must force a validated refresh; sizes exceeding capacity
must require bounded reconstruction. Changing inputs after evaluation must
reject commit. CPU fallback, if explicitly provided, must use the same typed
contract and remain independently qualified against the original callbacks.
No mode may invent a private context when the shared cuDSS primary is required.
Cleanup issues or retained resources disqualify the candidate process; no
borrowed primary reset/destroy is permitted.

Qualification would proceed from the real midpoint E/B/u comparisons, through
real shared-primary cuDSS/force lifetime tests, to a separately reviewed evolving
full-step replay. That replay must compare the original application against
each proposed mode over changing particles, fields, native-frame boundaries,
capacity/refresh cases and physical boundary/source events, retaining all
existing residual, timestep, finite, resource and conservation checks. Failure
cases must demonstrate transactional commit and explicit disqualification.
Only then could matched application timing include solves, drifts, packing,
transfers, synchronization, refresh/setup amortization, commit and all guards.
Shadow timings must also charge fresh independent references. Raw fixed-state
kernel/force timings must continue to be labeled as microbenchmarks.

This proposal does not authorize production resume, a shadow-policy change,
integration, or use of the unimplemented optimized mode. Those require their
own concrete review and the parent's explicit next gate.
