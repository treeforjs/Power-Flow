Exact common conservation ledger file codec

Job 683 passes 46 normal and 46 UBSan controls, with deterministic
replay and identical build results. The actual v9.2 neutral common ledger is
encoded and decoded after each of two accepted steps, preserving canonical
bytes, exact inventories, endpoint energy, accumulated bounds and roots.

The bounded binary format retains all exact accumulator digits and checks
canonical digits before arithmetic, species order/catalogue, clock/root binding,
byte/length caps, trailing data, and the cumulative energy/nuclei/charge algebra.
Tests reject malformed input and rehashed inconsistent accounting. Labelled
manufactured positive/negative residuals below binary64 range survive the
roundtrip and are refused against a zero bound. No low digits are discarded.

This qualifies serialization and arithmetic admission only. It does not restore
actual owners or grant a native/gas continuation. Root must verify its coherent
checkpoint envelope, reconstruct genuine owners, recompute their actual stocks
and configuration, and bind the conservation root through the accepted common
receipt. The production run remains paused and no full breakdown is claimed.
