Water chemistry v3.1 for the Power-Flow coupled PIC water-breakdown simulation (component, 2026-10-31)
Successor of pic-water-chemistry-v3-20261031 (unchanged); changes: SCOPE.md section 11, INTEGRATION.txt section A.

Plain directory; nothing committed. Native v8 links it later: INTEGRATION.txt. Data, models, coverage and open
items: SCOPE.md. Per-process table: data/tables/PROCESSES.txt. Qualification results: RESULTS.txt and
SWARM-VALIDATION.txt (copied from the receipt).

Layout
  include/pf_chem_v3/det_math.h         deterministic elementary functions (host + NVRTC)
  include/pf_chem_v3/chem_rng.h         chemistry stream layout v3 on the arrival owner Philox core (read in place)
  include/pf_chem_v3/chem_kinematics.h  relativistic collision kinematics (eV units)
  include/pf_chem_v3/chem_tables.h      table blob records and evaluation, majorant
  include/pf_chem_v3/chem_event.h       per-carrier null-collision stage (shared by CPU and GPU)
  include/pf_chem_v3/chem_host.hpp      blob loader, exact ledgers, CPU step, debits, quantized release
  device/chem_kernels.cu                NVRTC program
  data/sources/                         retained raw sources (one directory per source) + source indexes
  data/extracted/                       transcribed tables (one JSON per table, with locators and tags)
  data/tables/                          production spec, blob, manifest, checks, PROCESSES.txt, coverage policy
  tools/make_production_spec.py         writes the production spec from the sources
  tools/make_coverage_policy.py         manifest -> recommended coverage policy (admitted-region list)
  tools/digitize/                       figure digitization of the measured DR cross sections
  tools/build_tables.py                 spec -> blob + manifest + checks (deterministic; rebuilt in qualification)
  tools/make_lxcat_spec.py              LXCat complete set -> comparison spec (swarm references)
  tools/report_processes.py             manifest -> PROCESSES.txt
  tools/swarm_table.py                  swarm validation table
  tools/run_jobs.py                     single-thread jobs on leased cores (6, 7 and the lent cores 0, 1, 3)
  tools/dev_build.sh, dev_build_gpu.sh  builds (frozen toolchain)
  tools/research/                       acquisition and transcription helpers used for data/sources and data/extracted
  tests/test_core.cpp                   CPU qualification (math, RNG, data, 1e6-event ledgers, per-process kinematics,
                                        partitions, debits, release)
  tests/gpu_parity.cpp                  CPU/GPU bit identity
  tests/swarm_mc.cpp                    0D swarm Monte Carlo
  run_checks.py                         qualification driver writing receipts

Run (cores 6-7 and lent cores 0-1; builds and CPU stages wait for the core leases; the GPU stage takes the GPU locks for minutes):
  cd /home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031
  taskset -c 6,7 nice -n 15 python3 -I run_checks.py LABEL [inputs tables build core sanitized gpu swarm]
Receipts: /home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3.1-20261031/LABEL/receipt.json
Rebuild the tables after a spec change:
  python3 -I tools/make_production_spec.py data/tables/chemistry-v3.1-spec.json
  python3 -I tools/build_tables.py data/tables/chemistry-v3.1-spec.json data/tables/chemistry-v3.1
  python3 -I tools/report_processes.py data/tables/chemistry-v3.1-manifest.json data/tables/PROCESSES.txt
  python3 -I tools/make_coverage_policy.py data/tables/chemistry-v3.1-manifest.json data/tables/chemistry-v3.1-coverage-policy.json
