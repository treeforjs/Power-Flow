Gas target owner v1.3: CPU publication and motion-input milestone

Frozen normal347/h2 and ASan/UBSan348/s2 passed221 checks each. This includes
18 preparation-allocation failure sites, allocation-free accepted publication,
move/foreign/abandon/destroyed-owner/same-digest-retry capability checks,
nonempty restore/continuation, all125 inherited input/feedback tests, four original
v1 feedback/stage reproductions and nine v1.1 motion/exact-input reproductions.
The isolated intermediatev1.2 passed125 checks in normal340 and sanitizer341.

Numerical and API changes from CPU milestone5d0564b:
- Census IDs are unique and every row has valid species, phase space and weight.
  Motion survivors plus sinks must exactly partition the before inventory with
  unchanged species/binary weight. Births/debits/releases/merge-split are separate.
- Raw mover batches reject duplicate IDs, inconsistent cell-count denominators
  and RN(exact moved count) above declared RN(exact source count).
- Exact accumulation refuses bad/nonfinite/out-of-grid terms and rounded overflow.
- prepare_commit hashes, validates and compacts before freezing the candidate.
  Its move-only token binds owner lifetime and trial generation. can_publish and
  try_publish are noexcept; publication allocates nothing and launches no GPU work.
  The serial common coordinator must prepare/preflight every owner before any
  publication. This API alone neither coordinates owners nor provides concurrency.
- Configuration identity changes explicitly; no automatic checkpoint migration.

No share-transport equation, chemistry rate, particle count, precision or physical
closure is changed. Rounded carrier-share transport remains unqualified as a
full kinetic closure. Native admitted-event/debit export, atomic joint integration,
source/thermal/field coupling and full water breakdown are not qualified.
Native strict target-owner refusal is retained. GPU349/g1 is queued, not passed.
Modified device code is included for review, not qualified deployment.

Prior343 failed in the TEST allocation interceptor due to GCC inlining's
mismatched-new-delete diagnostic. The interceptor moved into a separate translation
unit with per-thread controls; no warning was disabled and no numerical code changed.
Prior sanitizer344 passed221; final source uses the new interceptor in347/348.
Original failed/succeeded job logs are retained. No production speed claim is made.
