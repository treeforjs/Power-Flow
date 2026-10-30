Authoritative collision exposure clock

The new prepare_exposure(ExposureClock, candidates) entry point retains the
actual accepted physical interval and frozen physical stage while using an
explicit local operator exposure. The full-time composition owner must supply
the authoritative common dt and dt/2 for each frozen chemistry stage. End must
equal begin + common_dt; subtracting rounded absolute timestamps is not used to
recover dt. The common endpoint must still advance in binary64.

Legacy prepare(Clock, candidates) retains its original end-begin behavior and
random decisions. Both methods share the unchanged candidate, propensity,
waiting-time and selection body; candidate construction is unchanged. Proposal
exposes the duration actually used. RNG checkpoint format remains unchanged.

47 existing checks and 3912 explicit-duration checks pass in the pinned runtime.
These include exact legacy decisions/RNG parity, physical-origin independence,
nonzero-epoch duration-boundary decisions, invalid input rejection and restart.
An initially incorrect TEST expectation for the direction of timestamp rounding
is preserved in attempt-001; the numerical implementation was unchanged by its
correction. Independent review is separate and recorded in the handoff.

This is an integration primitive, not water physics or physical evolution.
The native completion record and common coordinator must bind the same explicit
exposure operand; the old native v6.2 schema has not been reinterpreted.
Run run_checks.sh on the exclusive CPU 7 lease using frozen userspace. No GPU,
production process, dependency, guard or old sealed source is changed.
