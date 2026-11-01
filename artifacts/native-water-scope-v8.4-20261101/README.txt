Native v8.4 minimum water encounter guard -- CPU reference controls passed

CPU396 and ASan/UBSan400 each passed55 actual native scope controls, plus the
inherited target-energy/restore and7 cold-support controls. See QUALIFICATION.json.
This is a targeted host-component qualification, not full-native or production
qualification. No GPU census, complete reaction network or physical validation
is claimed; no deployed runtime binding was changed.

The new opt-in mode checks authentic native all-species inventory before target
filtering, refuses missing e+OH/O/H2 and wholly unsupported present encounters,
preserves existing regional and energy-owner guards, binds both collision stages
and restores persisted census/history roots. Device chemistry and device native
lanes explicitly refuse; there is no hidden host census transfer. Legacy mode
chemistry semantics are retained under an intentionally new schema identity.

Source/forms/tests and exact-base adoption.diff are preserved here. The full
runtime inputs and frozen dependencies remain read in place locally; their hashes
are in SOURCE-SHA256, BASE-SNAPSHOT and per-run receipts. Inputs containing source
library data, generated blobs, binaries and large checkpoints are intentionally
not republished. Reproduction uses the local baseline inputs and the frozen
one-core SimQueue commands preserved in the evidence. Never apply this patch to
an active writer or assume old v8.3 checkpoints migrate to the new schema.

397 was canceled after its arena compiler hit the inherited3GiB address-space
limit.400 used an explicit7GiB compiler cap within an8GiB scheduler reservation.
This changed build resources only; source, numerical flags and physics guards
remained pinned. The earlier successful45-control source/receipt is preserved
separately and is not substituted for the final55-control evidence.
