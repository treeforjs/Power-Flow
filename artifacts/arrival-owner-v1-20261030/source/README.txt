Arrival / count / RNG / ID owner v1 (pilot prerequisite H6), 2026-10-30

Normative contract: API-v1.txt. Plain directory, nothing committed.

include/pf_arrival/philox_core.h        Philox4x32-10, stream layout v1, exact conversions,
                                        exact Bernoulli, stochastic rounding, marks (host + NVRTC)
include/pf_arrival/canonical.hpp        canonical binary encoding, SHA-256
include/pf_arrival/arrival_owner.hpp    policy pf-arrival-policy-v1, plans, IDs, ledger, checkpoint
include/pf_arrival/source_v5_adapter.hpp  binding to heat/source owner v5 (read in place)
device/arrival_parity.cu                NVRTC parity program
tests/test_core.cpp                     known answers, layout, conversions, Bernoulli, statistics,
                                        exact debit sweeps, partition independence, IDs, ledger,
                                        checkpoint, refusals
tests/test_source_v5.cpp                against the sealed source owner v5 (neutral leaves through
                                        Owner::prepare, electron births through Owner::begin_source)
tests/test_native_census.cpp            births as native v5 genesis packets: ID = base + rank
tests/gpu_parity.cpp                    host/device bit parity through the shared GPU layer
tests/oracle.py                         independent exact oracle (Python integers/Fractions)
scripts/run_v1.py                       leased build/run driver writing receipts

Run: taskset -c 6,7 python3 -I scripts/run_v1.py LABEL [build|core|source|native|gpu|oracle ...]
Builds: core-lease-20261030/both-cores.sh + arch-base-20261028/run-frozen; CPU runs: any-core.sh
--wait (core 6 or 7, nice 15, one thread, no GPU); GPU runs additionally hold testing-gpu.lock and
testing-gpu-<UUID>.lock. Receipts: /home/treefor/projects/.workspace/outputs/pic-arrival-owner-v1-20261030/LABEL/.

Results: see RESULTS.txt.
