Gas target guard correction: CPU milestone, 2026-11-01

The corrected CPU implementation passed 85 focused checks in both optimized
and AddressSanitizer/UndefinedBehaviorSanitizer builds, frozen SimQueue job324.
Four original unsafe behaviors were reproduced against the captured preimages:
negative variance accepted/clamped, early commit, duplicate stage and replayed
step number. The corrected version refuses all four.

The checks include atomic whole-span feedback refusal, empty-population residuals,
nonfinite/overflow inputs, ordered stages and complete commit, grid/volume IDs,
nonempty checkpoint continuation, complete partitioned inputs, duplicate/missing
partial rejection and foreign-history rejection. See source/GUARD-SCOPE.txt.

CPU/sanitizer evidence only. GPU compilation/behavior/parity is pending job325.
The modified device sources are retained for review, not qualified deployment.
Full coupled target ownership, energy/momentum closure, transport realizability,
material/circuit/field coupling and a physical breakdown run remain unqualified.
Native's strict coupled-target-owner guard is unchanged. The original active v1
worker was not modified. Historical 'production' preset names are not approval.

Source is an isolated coherent snapshot plus corrections. Baseline bytes and
hashes are included; adoption/guards.patch is exact-base only, not applied live.
The complete source and tests can be reconstructed at the original runtime path
with the pinned dependencies listed in the receipts. The new configuration tag
refuses silent old-checkpoint migration. CPU commands (via frozen SimQueue):
 python -I tools/qualify.py normal LABEL
 python -I tools/qualify.py sanitized LABEL
GPU command additionally requires a scheduled GPU allocation:
 python -I tools/qualify.py device LABEL
Use unique output labels. No timing or full simulation claim is made here.
