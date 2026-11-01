Checkpoint file-envelope component, not common owner restore

SimQueue651 passed normal and ASan/UBSan controls. The exact nested file set is hashed with a64KiB buffer; nested owner manifest.json files are included. Hashing does not read a whole checkpoint file into memory. Bytes/file,total bytes,file count and manifest caps remain explicit. A coherent caller-owned immutable temporary directory is fsynced and renamed, then its parent is fsynced. Committed corruption refuses; incomplete temporary directories do not supersede committed history.

Single-writer contract: temporary/final paths are exclusively owned by the caller and serialized bytes stay immutable through commit. This implementation does not provide concurrent-writer no-replace arbitration. The parent directory must already be durable. Tests exercise logical interruption points, not actual power loss or injected failing fsync syscalls.

This is only file integrity and durability. It cannot authorize native/gas continuation or prove cross-owner stocks, clocks, charge, field or thermal consistency. A separate concrete accepted-owner loader and persisted native/gas attestation remain required. No physical simulation or runtime speedup claim.

651 builds/runs only in frozen userspace under SimQueue. Source includes a reproducible launcher. c1/c2 preserve harness include/formatting failures; no numerical guard was relaxed.
