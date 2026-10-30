# Bounded independent exact-count review

No actionable count arithmetic, carrier-sum or lossy-export defect was found in
the reviewed seal. This is bounded count/material-interface evidence, not
qualification of the newly implemented full-time native owner.

Reviewed source identity:
`b8a9a27e6a0fcd8c54e49dd341503fed165da460180b608aa0ae024b0436fb0e`.
Archive SHA:
`883ec79a67cb8462037dc1b3cb9d5a226e3a54058e4eaa30f56636296f89d755`.
The supplied runtime handoff path was absent; the corrected handoff under
`outputs/pic-profile/reaction-transactions-20261029/native-exact-count` was used.
All probes compile the extracted sealed source. The current backend is untouched.

The implementation imports both Count2 limbs, stores an immutable exact integer
in2^-1074 units, and copies it before debit/credit. Positive decomposition takes
at most53 high bits per leaf; no floating accumulator determines its remainder.
Rate diagnostics round by exact half comparison and tie parity. Count2 export
subtracts the rounded high limb exactly and refuses an unrepresentable remainder.
The native count codec additionally requires canonical positive-carrier order.
The source is separate from the material/adsorbed reservoir and introduces no
source law, sampling rule, new dependency, count selection or physical floor.

Two independently compiled modes passed: ordinary O2 and UBSan O1 with recovery
disabled. Each passed **1,287 Python Fraction oracle assertions over47 targeted
cases**, with identical JSON results and no UBSan diagnostic. These cover:

- Exact imports, negative low limbs and normal/subnormal transitions.
- Nearest-even ties on both sides of one and the lowest normal exponents.
- Minimum positive units, a dense2098-bit value, and a40-leaf representation.
- Credit/debit at the explicit range boundary and unchanged values on refusal.
- Invalid nonfinite, zero/negative debit, negative credit,43-input-leaf and
  overdraw admissions.
- Exact Count2 re-import when representable;21 observed snapshot export refusals
  when a third residual or scalar overflow prevents a valid Count2 conversion.
- Exact material atom, mass and packet-energy currency computed from every
  positive leaf, including states that cannot be exported to Count2. The energy
  operand K here is explicitly manufactured; this is arithmetic port evidence,
  not native kinetic/source-energy qualification.

The47 cases admit32 initial values and42 operations, with14 operation refusals.
The existing368 native/core checks,5120 source-weight debit matrix and root's
442 parent-subflow controls were not rerun. The copied sealed arena code/diff
was inspected for count codecs, validation, census access and diagnostic use.
The whole native motion/source implementation, allocator failures, all fenv modes
and physical weight/source-rate regimes remain outside this review.

Three contract details should remain explicit during integration:

1. The domain is the stated2098-bit integer bound, not a promise that its scalar
   rounded diagnostic is finite. Exact values above DBL_MAX can retain positive
   leaves while `rounded_for_rate()` and Count2 export refuse overflow. Native
   physics representability checks remain separate. `view().weight` is a rounded
   diagnostic; use actual native positive leaves for selection/current/mass/E.
2. At hi=DBL_MAX, the inherited Count2 normalization predicate computes an infinite
   upward spacing. Consequently `{DBL_MAX,-DBL_MAX}` imports as exact zero and
   exports normalized zero. This does not lose count and is also admitted by the
   heat Count2 predicate. Tightening this unphysical input edge requires an
   explicit shared Count2 admission-contract decision; it is not classified as
   a new gas-backend correctness defect.
3. Exact gas arithmetic cannot repair a lossy debit in another owner's source
   Count2 reservoir. Common source preflight must independently prove the actual
   reservoir delta equals the sum of actual newborn native weights. Do not round,
   truncate or back-export an extended gas value simply to fit Count2. The
   spatial material port can consume all native positive leaves as exact atoms,
   mass and K/reference energy without a Count2 conversion. Its pinned v6.2
   fixture parser uses the old pair codec; production integration must read the
   new typed count codec and keep its new schema/source identity.

Concrete boundary reproductions are in `inputs.json` and `probe-normal.json`:
`Count2-third-residual`, `Count2-max`, `dense-max2098`,
`dense-one-minus-min`, `RNE-even-low/high`, `separated`, and
`Count2-max-negative-cancellation`. `oracle-normal.json` and `oracle-ubsan.json`
retain every exact integer, carrier sum, refusal and material-currency check.

CPU7/nice15/thread1, existing exclusive testing-core7.lock; maximum compute lease
3.138s. No backend/source edits, dependency installs, competing count backend,
GPU work, production signals, policy changes, deployment or production resume.
Production PID278902/start5288339 remained stopped T and cgroup frozen1.
