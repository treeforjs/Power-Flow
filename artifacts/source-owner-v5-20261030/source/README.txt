Heat/source owner v5 workspace (2026-10-30)

native/mfem_electrode_heat_source_v5/  owner v5 headers, API-v5.txt, focused tests, oracles
audit-port/                            source-v4 audit port re-targeted to v5 (A07/A08/A02 rows),
                                       A08 minimal reproducer, compile-reject cases, port oracle
scripts/run_v5.py                      configure/build (both-cores lease) and runs (any-core lease)
                                       inside the frozen launcher; writes receipts
build/                                 CMake build tree (Release; normal and ASan/UBSan targets)

Prerequisites are read in place, unmodified: sealed v3/v2 sources of the v4 worktree
(runtimes/pic-electrode-heat-feedback-20261029, commit c65ac51) and the ExactParticleCount
header (SHA256 46e01fe2...). Evidence and receipts:
/home/treefor/sim-workspace/powerflow/outputs/pic-electrode-heat-source-v5-20261030/ (workspace alias)
