# Audit closure status (A01-A19)

Updated 2026-10-30. "Reproducer" = the audit's own counterexample executed against the corrected code; "new tests" = tests written
with the fix. Evidence paths are relative to this repository unless they start with `OUT/` (= /home/treefor/sim-workspace/powerflow/outputs)
or `RT/` (= /home/treefor/projects/runtimes).

| ID | Corrected component | Reproducer / tests | Independent review | Connected to common driver | Status |
|---|---|---|---|---|---|
| A01 | native candidate003 (c045c42d...) | New tests drive the audited route: 147/147 normal and sanitized; link-time traps record 0 solver and 0 advance calls; 368 legacy and 442 parent checks unchanged. Re-run after reboot identical. | Read-only review (no re-run) | v3 links the library | CLOSED (component) |
| A02 | Hall 8dddec0 (field), source v4 c65ac51 (source), v2 predicate, v3 wiring | New tests only; v3 constructs 0 owners | v2/v3 reviews say open | Hall declaration only; source v4 not wired | PARTIAL (Stage 1 T2/T3 target) |
| A03 | candidate003 | Sibling, reentry, destroyed owner, 48-case allocation sweep; no stale test; outer publish not exercised | Read-only | Declaration only | PARTIAL |
| A04 | candidate003 + Hall motion adapter | First RNG binding live; rest by pure guards; typed motion refuses JSON; no positive Hall receipt executed | Read-only | Declaration only | PARTIAL |
| A05 | candidate003 full-time grammar | Negative test only; nonempty restart and split/full parity unqualified | Read-only | No | FIX DRAFTED (deferred until after the pilot) |
| A06 | candidate003 typed restart codec | 14 corrupt-genesis mutations, duplicate keys, depth, option prefix, energy overflow; event/byte/collection caps, charge overflow, nonempty replay untested | Read-only | No | PARTIAL |
| A07 | source v4 c65ac51 | Audit counterexamples ported (artifacts/audit-reproducer-ports-20261030/source-v4): exact or refused, normal and sanitized identical, Fraction oracle 13/13; stock drain at extreme rate gives an untyped error | Port + oracle | Leaves not yet bound to native gas credits | CLOSED (component) |
| A08 | source v5 (artifacts/source-owner-v5-20261030) | G1 and G2 now refused with typed scoped refusals (minimal reproducer NOT-REPRODUCED); audit inputs still refused with v4-identical log scopes; exposure scalars admitted only within 2^-40 of the exact finite-log value; weighted packet values admitted only with zero weight-induced error; stock drain and energy-closure failures typed; 677 checks normal+sanitized byte-identical; Fraction oracles 30/30, 20/20, count oracle 49/49 | Same-author oracles only | No | FIX VERIFIED (component); births below ~1e-267 refused conservatively |
| A09 | Hall 8dddec0 `Model::admit` / Owner bootstrap | Audit counterexamples ported (artifacts/audit-reproducer-ports-20261030/hall-field): inconsistent charge and 50 T/m divB rejected at admit, prepare, bootstrap and restore; Gauss allowance equals the solver's pre-existing 2e-13 stopping rule | Port + bound analysis | Census built by the port, not a native arena | PARTIAL (Stage 1 T6 target) |
| A10 | Hall 8dddec0 authoritative compensated midpoint | Audit fixture: midpoint equals compensated average on all components (9,758 E edges, 1,583 B faces differ from high-word); prepare_closed_proposal executed; particle side closed by hash binding only | Port | Source consumer absent | PARTIAL |
| A11 | anchors-v1 (draft API) | Nothing implemented; cap 818 steps (~8 ps at 10 fs) | None | No | NOT STARTED (bounded in-memory history required before the pilot) |
| A12 | none (v4 layer map only) | v3 hard-codes surface source, collection and active layers off | - | No | NOT STARTED (pilot blockers H2, H3, H6) |
| A13 | water-data v7.1 metadata, v8 H(1s) continuation (raikou bf67687) | 581 component checks; 2a1 remainder, fragments, low energy open | None | No | SCIENTIFIC SCOPE (open) |
| A14 | streaming primitive v1, wrapper v1 | Audit fixture (46,726 queries) replayed byte-identically; 1 keV seam and macroparticle noise open | Same-author oracle + independent-process replay | No | PARTIAL |
| A15 | none | v4 still uses the 4.3 eV clean-surface law; no phase, return or held-T tests | - | No | NOT STARTED |
| A16 | melt proposal 1c11656 (synthetic only), SS304 expansion | Melt: 74 checks normal/UBSan, re-run identical after reboot; no single mass/energy owner; Hall thermal 298-1500 K; Count2 export conflict with source v4 | Self-review | No | NOT STARTED (core) |
| A17 | none | C19 binding, RGA calibration, slit/volume separation open | - | - | DATA SCOPE (open) |
| A18 | dispatcher v1 | 172 injected checks, 0 real children, not deployed; melt runner replaced by a versioned runner with explicit guards | None | No | FIX DRAFTED |
| A19 | - | No joint, pilot or convergence run | - | - | NOT STARTED |

Host changes: reboots at 12:09 and 14:31 MDT; host identity, frozen root (246,673 entries), launcher and the paused 26.4787 ns
checkpoint verified unchanged; every corrective suite re-run with identical scientific outputs
(artifacts/post-reboot-requalification-20261030).

Stage 1 (artifacts/stage1-joint-step-20261030): first real all-owner joint step published; T2 28/28, T4-T7 pass; T1 passes off
Hall node planes; T3 sampled sweep finds std::terminate in sealed JSON destructors; Hall node-plane guard defect found.

Next work and gate order: artifacts/coupled-pilot-plan-20261030/PLAN.md.

Re-run notes: native `verify-final.py` needs the CPU7 lease held externally and overwrites its evidence directory (back it up and
restore); do not re-run `package-final.py` (rewrites API-003.txt); Hall `run_components.py` takes the CPU6 lease itself and needs a new
label; v3 and melt runners take the CPU7 lease themselves (exit 75 when busy); the streaming primitive overwrites its fixed output
directory (back up), the wrapper honors `PIC_STREAMING_WRAPPER_OUTPUT`; inside the frozen sandbox use the physical workspace path, not
/home/treefor/sim-workspace.
