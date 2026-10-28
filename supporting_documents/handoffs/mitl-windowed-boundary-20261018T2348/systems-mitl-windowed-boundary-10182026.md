# MITL boundary launches from partial source windows

Power-Flow is pinned to `e9cd0425282de254e12b7683a197d40963ecd1eb` on `research/flash-publication-foundation`.
Local commits separate canonical source conditioning, paired packet sampling and
the export interface. `windowed-boundary.patch` contains this pass's changes.
The previous packet is
`/workspace/handoffs/mitl-emission-windows-20261018T2248.tar.gz`.

The conditioner retains each original interval/pool/current/width association,
actual patch area, initial population, exact supplied state and frozen local
temperature. Explicit selections use offsets inside original intervals; repeated
and overlapping selections are rejected. An authoritative immutable canonical
history supports revalidation of original densities, temperatures, geometry,
metadata and exact-zero skips. No new history is evolved. Positive partial
components remain present even when a density underflows; physical counts and
energies keep independent representation qualifications.

The paired sampler changes only launch width position. Component, source time,
state, temperature, current/normal position, mass, weight, velocity and momentum
remain the same in each pair. Mandatory packets retain every positive component;
exact represented-count CDF refinement avoids discarding tiny components through
probability underflow. Per-component weights may differ and retain signed number
closure. Expected translation, sampled translation, supplied internal/removal
costs and actual weighted wall debit are reported separately. No velocity or
inventory repair forces the sampled energy to its expectation.

Original interval origin, window start, local delay and log delay stay separate.
Best-effort global times have explicit precision flags. Local position offsets
remain available alongside projected positions. Projection/quantization guards
reject unsupported positive packet weights, speeds, Cartesian velocities and
momenta. Independent review found and closed coarse-subnormal velocity and
momentum distortions before freeze; the original focused evidence remains.

All 134 new checks passed in the portable subset. File counts and exact
commands/exits/log hashes appear in `verification.json`; earlier passing suites,
the previous source evolution and saved reference oracles were not rerun.
The source subset contains byte-identical public files and four named packaging
wrappers. The compressed archive and every manifest payload are verified by
the recorded delivery check. Commands after freezing remain in the local pass's
`commands.jsonl`.

After the bounded implementation, one preserved manufactured source HDF was
conditioned on two partial windows, retaining 48 components in 512 packets. An
independent NPZ/HDF readback checked original associations, explicit variates,
velocities, momenta, weighted number and expected/sampled energies without
calling a history loader, source evolution, scalar window or packet sampler.
Ten exact-rational fixed linear full-width launch controls match. A manufactured
partial aperture distinguishes the paired launch positions. These are numerical
source controls; they do not establish spectral distinguishability.

The error envelopes are estimates, not certified bounds. Source-emission windows
are separate from photon recording gates. Unequal weights and actual sampled
energy require explicit downstream support and ledger acceptance;
`host_budget_accepted` remains false. No gap, collision, radiation, molecular
state data, thermal feedback, charged trajectory or actual FLASH run was added.
Quark integration remains an owner handoff. Nothing was pushed or published.

The selected next task is deterministic window-conditioned timing/joint moments
and weighted sampling qualification. Actual response/collection/loss/field
definitions remain owner inputs. Physical water coverage, absolute VUV
normalization, original diagnostics and installed FLASH/workstation setup remain
unresolved. Read `source/docs/temkin_windowed_boundary.md` for the exact API.
