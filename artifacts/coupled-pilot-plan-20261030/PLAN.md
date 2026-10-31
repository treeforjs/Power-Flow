# Coordinator v4: first real all-owner step, short coupled pilot, any-N execution

Design review of 2026-10-30, read-only (nothing built or run for this document).
All timings are estimates; the first measurement will be the Stage 1 profile.

Paths: `RT` = /home/treefor/projects/runtimes, `OUT` = /home/treefor/sim-workspace/powerflow/outputs.
`N` = native candidate003 (`RT/pic-reaction-transactions-20261029/native-full-time-candidate003`, source identity c045c42d…),
`H` = Hall corrective owner (`RT/hall-coupled-transport-20261029/source`, commit 8dddec0, `development/hall-field-correction-20261030`),
`V3` = `RT/pic-water-common-commit-v3-20261030`, `S4` = source owner v4 (`RT/pic-electrode-heat-feedback-20261029`, commit c65ac51).

## 0. Findings that change the plan

1. **NEXT-INTEGRATION.txt describes the v2 header, not v3.** V3 `closed_parent_joint_step.hpp` already uses
   `hc::prepare_closed_proposal`, `m::bind_hall_total_field_motion(proposal.motion())` (138-149), compares the motion digest with
   `h::total_field_sha256(midpoint.E(), midpoint.B())` (189-190) and calls `owner.preflight_final(..., *boundary)` last (195).
   The real gap: nothing constructs real owners; no genesis/A09 protocol; species and frame identities are authored; no fault or
   mutation hooks inside `can_publish`; the outer reservation is not released after publication.
2. **Geometry mismatch.** The native003 fixture uses z in [0, 0.025] (`tests/full_time_fixture.hpp:35-47`); Hall uses z in
   [-depth/2, +depth/2] (`test_owner.cpp:8-11`, `axial_mesh.hpp:17-24`, `MeshSpec{45,4,1,.009}` at `field_contract.hpp:105`).
   The joint fixture must use Hall's `geometry()`.
3. **Native history cap blocks even a 1-2 ns single-process pilot.** 4096 rows / 16 MiB (`mfem_pic_arena.cpp:33-35`, journal
   append 566-578); every motion row stores full per-carrier endpoint and path JSON (`mfem_pic_full_time.inc:1077-1084`). The coupled
   dt guard is <= 10 fs (`field_model.hpp:73`, `common_commit_contract.hpp:82`), so 1 ns needs >= 1e5 steps; the cap allows 818.
   The bounded-history (anchor) part of A11 cannot be deferred; only disk restart can.
4. **Native full-time is closed to births and sinks.** Motion requires `sources == []`, `boundaries == []` (inc:888-890), coverage
   `all_native_carriers_no_sources_no_boundaries` (inc:977-978) and strictly interior endpoints (inc:218-231, 990-991). Packets are
   created only at genesis by `create_full_time` (inc:1445+). No birth, collection or escape API exists.
5. **Two electrode thermal owners.** Hall `Owner` holds a `hall_driver::Thermal` supplying sigma(T) and Joule heating (`owner.hpp:43,62`);
   source v4 privately keeps a heat-v2 Coupler driven by `mfem::Coefficient& joule` (`source_owner.hpp:172`). Coupling both would count
   Joule heat twice or give emission the wrong temperature (A16).
6. **Parent channel aborts on unsupported energies.** Refusal where sigma > 0 but no kernel exists: the 13 eV, 38.3 eV onset gap and
   8000 eV cases (`water_parent_bridge.hpp:91-92`). Gap electrons above ~8 keV in a cell with H2O would abort the step.
7. **The qualified GPU path solves a different operator.** `CudssPoisson` (`RT/mfem-gpu-pic-20261029/source/src/cudss_poisson.hpp:12-66`)
   caches an FP64 SPD factorization of a constant electrostatic Poisson matrix (IR_N_STEPS=2). Full-time mode forbids native Poisson (A01)
   and `create_full_time` requires `!distributed_fields` (inc:1449). Hall's step is a matrix-free implicit edge solve
   `M + dt/2 (S+Y) + dt^2/4 C^T H C` (`field_model.hpp:79-88`) with Jacobi PCG to residual <= max(1e-30, 2e-13 ||b||) (`foil_em.hpp:56-63`);
   S depends on sigma(T), so a cached factorization goes stale as the foils heat. Reuse the cuDSS wrapper and guard pattern, not the solve.
8. **Both owners are serial today.** Hall is built with `/usr/bin/c++` (+ libmpi only for `MPI_Init` in tests) with no ParMesh/hypre path.
   The earlier 8-rank production run kept all particles on root and scattered CSR row blocks into a HypreParMatrix
   (`mfem_pic_distributed.hpp:17-97`); only that row-scatter pattern is reusable.
9. **Libraries here:** hypre 2.33.0 and MFEM 4.10 without CUDA; cuDSS 0.6.0 with `cudss_distributed_interface.h` and
   `libcudss_commlayer_{openmpi,nccl}.so.0`.
10. **Exact-count export conflict.** The melt proposal needs two-limb export of adsorbate stocks (`ablation_proposal.hpp:127,158`);
    source v4 keeps exact stock without a two-limb export (`exact_particle_count.h:42-45,76`). Deferred: the pilot has no melt.

## 1. Stage 1: first real all-owner step (serial, N = 1, sealed sources unchanged)

New directories only: `RT/pic-water-common-commit-v4-20261030/` and `OUT/pic-water-common-commit-v4-20261030/`.

- `source/genesis_v4.hpp` (joint genesis and A09 admission), `source/identity_v4.hpp` (species and frame identities derived from real
  objects), `source/closed_parent_joint_step_v4.hpp` (v3 flow + release of the outer reservation after publication + fault points inside
  `can_publish` + `friend struct TrialTestAccess` defined only in tests).
- `tests/fixture_v4.hpp`, `tests/test_stage1_joint.cpp`, `tests/alloc_hook.cpp`, `tests/declaration_v4.cpp`, `run_stage1.py`.
- V3 `native_parent_subflow.hpp` and `common_commit_contract.hpp` are included unchanged via `-I V3/source`.

Fixture (existing files): Hall trace `H/development/hall-finite-foil-em-20261029/data/shot_17133_current.csv` (C19);
water data `RT/pic-water-exact-count-subflow-20261029/inputs/h2o-gas.json`; binding `OUT/pic-water-exact-count-subflow-20261029/run-binding.json`;
states/channels from V3 `test_native_parent_subflow.cpp:13-21`; Hall `geometry(Mesh(45,4,.009,1))`; collision grid {45,45,5}.
Population P1: one electron and one H2O+ packet co-located at (.0045, .0045, 0), counts {1,0}, different proper velocities
(e.g. 100 eV electron, thermal ion), so every node starts exactly neutral and the DC seed satisfies absolute Gauss.
Variant P2 adds a high-weight H2O packet in the same cell for >= 1 parent event.

Genesis and A09 admission:
1. `auto model = make_shared<const hc::Model>(MeshSpec{}, trace)` (`field_model.hpp:39-58`); `seed = model->dc_seed()` at
   `drive.start() = 7.4978805549199997e-10` (`Model::admit` rejects earlier times).
2. Predict the native census: IDs from 1 over `positive_carriers`, charge = charge_e x pic_e (inc:1462-1485).
3. `ledger = declared_dc(*model, seed)` (`test_owner.cpp:12-17` recipe) with `ledger.particles = hc::charge_of_census(mesh, o, predicted)`.
4. Provisional `hc::Owner(model, seed, ledger, o, predicted)` (`owner.hpp:89-90`): read `field_sha256()`, `heat_sha256()`, destroy.
5. `ParentSubflow::full_time_configuration_sha256(...)`, `initial_full_time_rng_sha256(...)` (`native_parent_subflow.hpp:332-342`);
   `config = closed_parent_configuration_sha256(chem_cfg, model->configuration_sha256(), inputs)` (V3:17-29).
6. Genesis binding: field, heat (step 4); `mesh = h::motion_mesh_sha256(model->mesh(), o)`; `material = h::material_binding_sha256()`;
   frame = derived digest of schema + shift {-w/2, -g/2, 0} (the shift native checks at inc:955-957) + motion mesh;
   `source = closed_parent_inactive_source_sha256(config)` (V3:43-47); rng (step 5); genesis clock digest.
7. `arena = m::PicArena::create_full_time(o, states, channels, packets, genesis)` (`mfem_pic_arena.h:141-143`).
8. Require `h::canonical_census_sha256(arena->rate_population(cap)) == ...(predicted)`.
9. Real `owner = make_unique<hc::Owner>(model, seed, ledger_from_actual, o, arena->rate_population(cap))`; its constructor admits absolute
   Gauss and divB against the actual census (`owner.hpp:48-66`); require field/heat hashes equal to genesis.
10. `ParentSubflow chem(*arena, counts, binding, limits, bridge, s1, s2, true)`; require `owner_values_match(observe(...))`.

One step (V3 sequence, v4 adds checks): `observe_closed_accepted` -> `arena.begin_full_time_step` -> first exposure (`chem.prepare_on`) ->
`fork_transport_trial`, `native_census` -> `h::MotionBinding` (staged arena prefix; species digest derived from `ArenaStateSpec`) ->
`hc::prepare_closed_proposal` (`root_api.hpp:40-60`, up to 16 iterations to a bitwise-stable midpoint) -> `m::FullTimeMotionBinding` ->
`motion_scope.prepare_total_field_motion(binding, m::bind_hall_total_field_motion(proposal.motion()))` (typed route inc:1211-1244;
admission inc:883-1088) -> `accept_transport_trial` -> second exposure -> `private_step.seal({field_after, heat_after, source, rng})` ->
`chem.seal(first, second, native)`.

Commit: `can_publish` (arena and chem `can_publish`; `candidate_full_time_boundary()`; exposures' transport trial == motion seal; fresh
observations with RNG recomputed; `trial_field_sha256 == h::total_field_sha256(midpoint E/B)` and != `midpoint().seal()`;
`accepts_closed_step`; last `owner.preflight_final(...)`), then `arena.try_publish` -> `chem.publish_prevalidated` ->
`owner->publish_prevalidated(move(proposal))`, no I/O or callbacks between swaps; destroy the private step and prove a second step can
be prepared.

Test matrix ("byte-identical" = native `arena.save()` bytes + `state_prefix()` + boundary JSON, Hall `owner.save()` bytes + field/heat
hashes, `chem.checkpoint().dump()`):
- T1 positive: one step at t0+dt; field, heat, RNG change, source unchanged; zero calls to wrapped `GapPic::solve_fields/advance/
  advance_empty/advance_observed`; field/motion/heat balance within the Hall bound; count of edges where compensated midpoint differs
  from the high-word average recorded; two fresh processes give byte-identical receipts (N-parity reference).
- T2 mutation after sealing (test-access hook + twin fixtures): field (other clock), heat (dt/2), source (refused by `preflight_final`),
  RNG (twin chemistry; `chem.restore()` between prepare and publish), geometry (twin with `nextafter` foil coordinate; compile-reject for
  any non-const geometry accessor), config (species/frame/caps), native token from a twin arena. Each refuses with all owners
  byte-identical; the unmutated trial still publishes afterwards.
- T3 faults: throw points at every hash/allocation in prepare and `can_publish`; global `operator new` countdown sweep; zero allocations
  between first and last swap of `try_publish`.
- T4 single use: second `try_publish`, moved-from trial, reused native token, consumed Hall proposal and stale chemistry generation refuse.
- T5 rollback: declaration test that v4 never exposes native `try_rollback`; failure atomicity: after every T2/T3 refusal a fresh trial
  reproduces T1 receipts byte for byte.
- T6 A09: lone electron with DC seed refuses (absolute Gauss); ~50 T/m divB seed refuses; ledger gross mismatch refuses; 1-ulp ion
  move after provisional hash refuses; genesis time != drive start refuses.
- T7 A01/A03: legacy entry points refuse with zero callbacks; second `begin_full_time_step`, sibling and reentry refuse.

Closure evidence: A02 (T2, T3, live re-observation), A04 (typed motion only, RNG chain, source digest), A09 (T6, actual-census bootstrap),
A10 (motion digest = actual midpoint digest with nonzero compensation), A01 (T7), A03 partial. A05, A11, A12 not closed by Stage 1.

Build: frozen launcher, `/usr/bin/mpic++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off -fno-fast-math
-ffunction-sections -fdata-sections -I v4/source -I V3/source -DPIC_REACTION_NATIVE_SOURCE_SHA256="c045…"`, `-isystem` in V3 order
(N/source/include, I, I/include, N/inputs, H/development/hall-field-correction-20261030, H/development, MFEM install-mpi/include,
hypre/include, metis/include). Link the test, `alloc_hook.cpp`, `N/source/src/mfem_pic_hall_motion_adapter.cpp`, Hall
`hall-canonical-driver-20261029/{native_arena_motion.cpp, kinetic_boundary.cpp}` (compiled separately with Hall's flags, no -Werror),
plus native-v6.2 `mfem_gap_pic.cpp` for `pic_boris` and `pic_kinetic_energy_J` (link gap found by the field reproducer port),
`N/build/libpic_native_candidate003.a` (d2568af9…), libmfem, libHYPRE, libmetis, libGKlib, the four `-Wl,--wrap=` symbols,
`-Wl,--gc-sections -lcrypto`. Pins: native SOURCE-SHA256.json, Hall files at 8dddec0, N/inputs/hall-motion-v2 byte identity,
single nlohmann/json header identity from the `-MD` dependency file.
Resources: core-7 lease, nice 15, one thread, no GPU, -j1, RLIMIT_AS 3 GiB, ~20 min run (allocation sweep dominates), < 1 GiB RSS.

Risks: native `census_sha` (inc:605-628) vs Hall `canonical_census_sha256` (`native_arena_motion.cpp:73-75`) never cross-checked on
native003 real data; duplicate JSON header ODR risk; Hall `final_approved_` stays set if native `try_publish` fails, so Hall publish
must only run inside `try_publish`; the 16-iteration midpoint loop can refuse.

## 2. Short coupled pilot (1-2 ns, one process, diagnostic checkpoints only)

Start state: the 0.749788 ns supported start (first C19 knot, Hall `dc_seed` time). The 26.4787 ns electron-only checkpoint is not usable:
legacy voltage-driven electrostatic state on a different mesh (conflicts with C19 current drive, A17), no importer for its field into
Hall E/B with conductor/exterior banks (A09), heating history in another engine, and past the 12.428 ns expansion comparison (A16).
C19 current is only 6.0-6.8 kA over 0.75-1.25 ns, so a pilot from 0.75 ns tests the machinery, not breakdown. A physically interesting
window needs a Hall-only vacuum pre-run (`prepare_vacuum` / `preflight_vacuum_final`, `owner.hpp:101-123`) to a chosen T_on, then a
fresh native genesis bound to that Hall state, with unreleased source quotes carried as unresolved scope (practical only after item H).

Must be built first:
| ID | Kind | Missing today | Build |
|---|---|---|---|
| H1 | DONE in native v5 (artifacts/native-v5-20261030) | A11 bounded history | Native v5 in-memory anchors every K steps; motion rows store endpoint hashes instead of full JSON. Disk restart can wait. |
| H2 | blocking | Native birth API | Native v5 source scope crediting source-v4 Trial leaves (exact debit `source_owner.hpp:183`) to native H2O and electron packets via one exact ExactParticleCount equality, with physical birth time and remaining flight; an open common predicate (V3 `accepts_closed_step` requires an unchanged source digest, contract:95). |
| H3 | blocking | Collection/escape | Native v5 sink records; Hall motion v3 with first-hit truncation (`kinetic_boundary.cpp:63`); collected charge free -> conductor bank in a new Hall field version, plus wall heat. |
| H4 | DONE: coverage v1 + v2 (artifacts/parent-coverage-v1-20261030, -v2-20261030; >10 keV provisional RBEB) | Parent-channel coverage | Qualified kernels for the onset gap and > 8 keV, or a versioned refusal scope recording unresolved counts (no silent skipping). |
| H5 | invalidating | One thermal owner | Hall thermal supplies sigma(T) and Joule; source v4 consumes Hall T and returns surface loads (new Hall version + source binding; no private Joule). |
| H6 | blocking | Arrival/count/RNG/ID owner | Counter-based RNG keyed by (run seed, step, patch, candidate index) and a global ID allocator. |
| H7 | blocking | 100,000-carrier native cap (`arena.cpp:33`) | Versioned packet merge/weighting policy conserving charge, momentum and energy, or the pilot stops at the cap. |
| H | practical | Per-step JSON hashing (~0.4-0.7 s/step) | Binary canonical encodings hashed once per immutable object (Hall `state_sha256` dumps the field ~10x per step, `field_model.hpp:65-66`; native `frame_prefix`, inc:232-240). |
| - | binding | Source v4 field view | Build `FinalFieldView{v3::TrialFrame, axial::Mesh, E, B}` (`source_owner.hpp:41,47`) from Hall's actual midpoint arrays. |

Source owner v5 (one version, after Stage 1 proves the joint interfaces): close A08 defects G1 (handoff zeroes representable current
loads, charge impulses and continuity for tiny admitted births) and G2 (subnormal expected-count scalar admitted with up to ~2x error;
admit through the finite log or refuse with scope), typed scoped refusals for stock drain and the v2 energy-closure path, H5 (consume
Hall temperature, no private Joule), the FinalFieldView binding to Hall midpoint arrays, and a decision on two-limb export for the
SS304 melt proposal. Evidence: artifacts/audit-reproducer-ports-20261030/source-v4.

Per step: first exposure (qualified parent channel only) -> source v4 prepare (emission + desorption births, exact debit = native credit)
-> Hall motion with first-hit sinks, ions as moving carriers (screening only from the Hall EM solve) -> field and heat proposal ->
second exposure -> seal -> preflight every owner -> allocation-free publication native -> chemistry -> Hall -> source.

Per-step ledgers (coordinator, canonical ID order): charge (free net/gross, conductor bank, exterior bank, collected, escaped, max absolute
Gauss, max divB); parent H2O exact inventory (source stock + native H2O + H2O+); H and O nuclei; energy (Hall field delta, plasma work,
Joule, outer loss/supply, circuit work, kinetic delta, thermal delta, wall heat, source-v4 neutral kinetic/internal/reference energy,
reaction defects, collected/escaped kinetic energy).

Resolution requirement: the purpose of the GPU-resident path is to make fine meshes affordable. Any screening or breakdown statement
from the pilot needs a refinement series on that path (200, 100, 50 and 25 um cells over the same window, more if the GPU memory
budget allows: ~200 B per edge, ~60M edges on 12 GB, ~110M on 24 GB), with the Debye length and the cell size reported side by side
per step. The 200 um mesh alone qualifies only the machinery.
Particle-count requirement: the GPU-resident path should also make larger macroparticle counts affordable (target 1e6-1e7 carriers
on 12-24 GB), reducing statistical noise in the EEDF tail that drives ionization and in the screening charge density. Pilot results
need a macroparticle series (e.g. 1e4, 1e5, 1e6 per species) alongside the mesh series (audit A14: macroparticle noise converged
independently of exact-count arithmetic). Enablers: a configurable, budget-checked carrier cap in native v5 (replacing the fixed 1e5,
item H7) and per-step sealing cost that stays cheap in carrier count (incremental block digests, GPU-computed block hashes with host
verification of a sample and the root).

Recorded limitations of the pilot: no disk restart (A05/A11); ionization only, ions do not react (A13); 1 keV seam (A14); emission phases,
unit accommodation, clean 4.3 eV surface (A15); SS304 solid <= 1500 K, refuses beyond, no melt (A16); C19-only drive (A17); 200 um mesh
does not resolve the Debye length above ~1e16 m^-3 at a few eV, so a refinement series is required before any screening claim; no
ion-neutral collisions or wall secondary emission; 0.75-2.75 ns covers only the low-current window.

## 3. Parallel and GPU architecture

Replicated mesh + particle decomposition with one coordinator rank. Rank 0 holds every authoritative owner and performs every commit,
hash, ledger and checkpoint; workers receive their ID-sorted census slice and the midpoint E/B and return endpoints, paths and
current/charge contributions. Contiguous ID blocks balanced by particle count: any N (including 3 or 5), no slabs or migration; for
N >= 4 rank 0 pushes no particles. Mesh operators and midpoint fields in per-node MPI-3 shared windows. Spatial decomposition is deferred
until the mesh exceeds node memory or the population far exceeds the native cap. Stage 1 runs serially first; the parallel versions must
reproduce Stage 1 byte for byte at N = 1. Parallel motion needs a new Hall motion provider producing the identical v2 record and a native
candidate004 (or reviewed re-pin) with an unchanged v2 recognizer.

Per fixed-point iteration: broadcast midpoint E/B; each worker hashes its copy and returns the digest (must equal the motion
`trial_field_sha256`); Gatherv + k-way merge in (edge/node, particle ID) order = the serial sum order; coordinator field solve. One
collective agreement before any swap; no MPI between swaps. Chemistry SSA, source and thermal stay serial on the coordinator.

Field-solve backends (site config; identical solver contract):
| Backend | Use |
|---|---|
| cpu-pcg | Existing sealed Jacobi PCG on the coordinator; default for the 45-cell mesh. |
| gpu1 | One GPU: deterministic FP64 CG on an explicitly assembled Hall CSR (row-per-thread SpMV, fixed-tree dots, no atomics); optional FP32 cuDSS factorization of the sigma-frozen fixed-dt operator as preconditioner, refactored when iterations grow. |
| hypre-cpu | Root-assembled CSR row-scatter into HypreParMatrix + PCG (Jacobi or AMS), reusing `mfem_pic_distributed.hpp:17-71` with static hypre 2.33.0; large N or multi-node. |
| multi-GPU | hypre rebuilt with CUDA (ParCSR PCG on device); cuDSS 0.6.0 multi-GPU/multi-node only as a cached preconditioner. |

Crossover: at ~4e4 edges (45-cell mesh) CPU PCG takes ~1-2 ms, below GPU transfer and launch overhead, so the GPU does not help there;
the GPU wins from about 5e5 edges (cells <= ~60 um). One GPU stops being enough above ~200 B/edge of device memory (~60M edges on 12 GB,
~110M on 24 GB) or when the solve exceeds the push at large N. Kernels select `--gpu-architecture` from the actual device at runtime
(the experimental NVRTC kernel hard-codes compute_86). RTX PRO 4000 Blackwell is sm_120: pinned cuDSS 0.6.0 carries sm_120 SASS and
PTX; CUDA 12.8 needs driver >= 570.

Rank-count independence: bitwise identical for every N with cpu-pcg and gpu1 (pure per-particle push with the same binary and
-ffp-contract=off, ordered merges, order-free exact-count chemistry, coordinator SSA, keyed birth RNG, ledgers and seals over ID ranges).
hypre-cpu and multi-GPU are not bitwise across N or GPU count; their receipt is the coordinator's serial recomputation of the independent
residual (<= max(1e-30, 2e-13 ||b||)), absolute Gauss, divB and field/Joule/Gauss-increment bounds (`field_model.hpp:121-127`), all
N-independent. Validation: N in {1,2,3,4,5,8} with cpu-pcg or gpu1 over 20 steps must be byte-identical (excluding the execution
envelope); hypre/multi-GPU pass guards for each N and match the first-step field within 1e-12 relative.

Every receipt and checkpoint carries an unsealed execution envelope (N, CPU list, coordinator/workers, backend and library hashes, GPU
count/name/UUID/capability/driver, rank-to-GPU map, site-config digest, image digest). Diagnostic checkpoints are written by the
coordinator from canonical owner states, so they are independent of N; a later restart may change N or backend after re-admission.

Launcher: one build plus a per-site JSON config (mpirun/srun, `--cpus 4-7` / `0,2,5`, one rank per listed physical core checked against
`thread_siblings_list`, optional frozen root and per-core leases, GPU identity rules, environment, expected kernel/driver). Workstation
profile keeps the frozen root, nice 15, per-core leases, GPU match by UUID, and uses only cores 6-7 (0-5 are allocated to Spheral).

## 4. Cost model per step (estimates)

dt <= 10 fs (>= 1e5 steps/ns); push ~3 us per carrier per iteration x 2 iterations; coordinator floor after item H ~15 ms on the
45-cell mesh; workers = N for N <= 3, N-1 for N >= 4.

| | N=1 | N=2 | N=4 | N=5 | N=8 | N=64 |
|---|---|---|---|---|---|---|
| P=1e4, 45-cell mesh (ms/step) | 75 | 45 | 35 | 30 | 24 | ~16 |
| P=1e5 (ms/step) | ~650 | ~330 | ~230 | ~180 | ~115 | ~45 |
| 25 um mesh, gpu1, P=1e4 (ms/step) | | | | | ~100 | |

75 ms ~ 2.1 h/ns; 24 ms ~ 0.7 h/ns; 115 ms ~ 3.2 h/ns; 25 um ~ 3 h/ns. With today's sealed code a step costs ~0.4-0.7 s (JSON
hashing) and the run stops at <= 818 steps. The electron-only reference (55 steps/s x 5.3 fs ~ 57 min/ns) used a different solver,
model and dt and is not comparable.

## 5. Portable deployment

Status 2026-10-30: paid cloud compute is out of scope. Near-term targets are the workstation (GPU-resident step on the
RTX 3080 Ti) and the 64-core RTX PRO 4000 Blackwell station; further owned machines and one or more clusters are expected later.
Portability therefore remains a design requirement (any N, per-site config and admission, container image of the frozen
environment, hypre-cpu and multi-GPU backends, multi-node launch), built when that hardware is available. Cloud shapes below are
reference only.


Does not delay the workstation pilot; restart (A05/A11) becomes required only for preemptible cloud runs. Container image (Apptainer for
Slurm, OCI) from `arch-base-20261028` root-files.json plus pinned MFEM/hypre/METIS/GKlib, cuDSS 0.6.0 and CUDA 12.8, image digest in
every receipt; the fixed-kernel guard becomes per-site admission (kernel, driver, GPU set, CPU model/flags, MPI fabric); each new site
passes Stage 1 and the N-parity suite before science runs. Cloud spending is minimized: own hardware first; Spot CPU VMs (at most one
low-cost GPU) for overflow only; 10-15 min checkpoint cadence, atomic checkpoints with SHA256 manifests in object storage, budget cap and
auto-shutdown. Cross-machine determinism: same binary, -ffp-contract=off, no fast-math, no -march=native, frozen libm; pin BLAS code
paths or keep BLAS out of sealed paths; elsewhere use the field-guard receipt plus ledger identities.

GCP shapes (A100/H100 excluded on cost; public specs, not verified here; see also 5b):
| Role | Shape | Usable ranks | Field backend | Notes |
|---|---|---|---|---|
| Portability smoke test | g2-standard-4 (1x L4 sm_89 24 GB, 4 vCPU, 16 GiB), Spot | N <= 2 (SMT siblings refused) | gpu1 (FP64 CG on device, optional FP32 cuDSS preconditioner) | Qualifies image, site admission, launcher, runtime architecture selection, repeat-solve determinism and cross-machine parity; not for timing. |
| Many ranks, push dominated | h3 / c3 / c2d HPC sizes, Spot, compact placement, Tier_1 networking | 64+, multi-node | cpu-pcg (45-cell mesh) or hypre-cpu (refined / multi-node) | The 45-cell field solve (1-2 ms) does not justify a GPU; cores and network set scaling. |
| Many ranks + low-cost GPU | g2-standard-48/96 (4/8x L4) | up to 24/48 | gpu1 on one L4; extra L4s only with a multi-GPU backend | Middle ground when a refined mesh fits one L4. |

L4 smoke test: Stage 1 binary at N = 1 and 2; repeat-solve determinism; a 20-step N-parity window compared with the workstation
(CPU-side seals and ledgers bitwise, field through the guard receipt); a few hundred pilot steps to measure rank-GPU communication.
Deterministic GPU kernels use `-fmad=false`, row-per-thread SpMV and fixed-tree reductions so FP64 CG can match across sm_86, sm_89
and sm_120 (to be shown, not assumed); runs using cuDSS internals compare across architectures only through the field guard receipt.

## 5b. Cost constraints (supersede any high-end-GPU option above)

- Backends cpu-pcg, hypre-cpu and gpu1 are equal in status and pass one solver contract (independent FP64 residual
  <= max(1e-30, 2e-13 ||b||), absolute Gauss, divB, energy and Gauss-increment bounds). The GPU is used only after a measurement shows
  the refined-mesh field solve exceeds the CPU push. Multi-GPU is deferred to the multi-GPU station.
- Hall edge operator after Jacobi scaling: kappa ~ 1.0007 (200 um) to 1.04 (25 um) in vacuum; the SS304 conduction term
  (sigma dt / 2 eps0 ~ 790) sits on the diagonal and is removed by the scaling; kappa O(1-10) (estimate). FP64 Jacobi PCG converges in
  ~10-30 iterations and sigma(T) changes each step, so no factorization for this operator; on low-FP64 GPUs FP64 CG is bandwidth bound.
- Poisson-type solves (Hall DC seed, legacy comparisons): kappa ~ (L/h)^2 ~ 2.6e6 on the legacy 40 mm / 25 um mesh, kappa u_fp32 ~ 0.15.
  Use an FP32 cuDSS factorization as a preconditioner inside FP64 CG (valid for SPD; never a correctness risk), cached while the
  matrix is unchanged. Earlier FP64 cuDSS tuning showed no end-to-end win, so FP32 must earn its place by measured time per step.
- FP64-accuracy receipt per solve: factorization precision and backend version; refactorization key (sha256 of FP64 CSR and dt); CG
  iterations, residual history, final FP64 residual recomputed independently on the CPU in canonical order against the same 2e-13 guard;
  Lanczos kappa estimate and forward-error bound; downstream Gauss, divB and energy guards; repeat-solve bit identity on one GPU.
- Measured electron-only phase split (pic-profile/phase-breakdown-summary.json, 210 steps, ~29.5k particles, 74 ms/step baseline,
  55 ms/step cpu-v2): field gather + Boris 42-55%; GPU Poisson incl. transfers 19.5-27%; CPU field assembly + ledgers 12%;
  drifts + continuity 9-11%; residual checks 2.5%; other 2.4%.
- CPU priorities (new versions; parity with sealed versions at N = 1 is the gate): (1) Hall push: structure-of-arrays carriers,
  precomputed cell index, fused gather/push/deposit, no per-particle heap allocation, same double-word deposit arithmetic
  (`native_arena_motion.cpp:104-132`); (2) Hall field apply and PCG: preassembled scaled CSR and fixed workspace instead of three
  allocated vectors per apply via std::function (`field_model.hpp:82-83`, `foil_em.hpp:56-63`); (3) binary canonical hashing;
  (4) native snapshot structural sharing instead of 4-6 deep copies per step (`inc:841,1073,1350,1430`); (5) allocations per phase as a
  gated metric; (6) `cudaDeviceScheduleBlockingSync` so GPU waits do not spin the coordinator core. Stage 1 records per-phase timers in
  the same groups as the electron-only profile.
- Recalibrated estimate with items 1-4: 1e4 carriers, 45-cell mesh, cpu-pcg: ~55 ms/step at N = 1 (~1.5 h/ns), ~20 ms/step at N = 8,
  coordinator floor ~12-15 ms/step.
- Cloud: own hardware first (workstation, then the 64-core station); cloud only for overflow or site qualification on Spot CPU VMs
  (8-32 physical cores, compact placement) and at most one L4 (g2-standard-4/8). Controls: billing budget alerts at 50/90/100% with a
  Pub/Sub-triggered stop of all VMs at 100%, per-VM `--max-run-duration`, power-off on job exit or failure, launcher wall-time cap, no
  standing GPU reservations, storage lifecycle rules; 1-2 h checkpointed runs with rank-independent atomic checkpoints.

## 5c. GPU-resident step (target on any machine with a GPU)

The "field solve not worth a GPU" conclusion above holds only while particles stay on the CPU. The particle work is the largest phase
(field gather + Boris 42-55% of an electron-only step, 23-40 ms for ~29.5k particles) and the earlier GPU gather microfixture took
~0.65 ms median for the same particle count (synthetic self-field; pic-profiling-resume.json). FP64 push cost (~100-200 operations per
particle) is well under 1 ms for 1e5 particles even at a 1/64 FP64 rate; the work is bandwidth bound.

- Device-resident: carriers (structure of arrays), midpoint E/B, the assembled scaled Hall edge operator (CSR), gather, Boris push,
  path/first-hit tracking, current and charge deposit, and the field CG. Iterations of the exact-midpoint fixed point stay on device.
- Host (coordinator): exact-count chemistry, source accounting and births, sealing, ledgers, preflight and publication. Per step only
  digests and compact carrier records cross PCIe (~1-10 MB).
- Determinism: deposits by sort on (edge/node, carrier ID) plus fixed-order segmented reduction, no atomics; CG with fixed-tree dot
  products and row-per-thread SpMV; compile with `-fmad=false`, IEEE division and square root. The push uses only + - x / and sqrt, so
  bitwise equality with the sealed CPU push (same double-word deposit arithmetic) is the acceptance test at every N and on every GPU
  architecture (sm_86 / sm_89 / sm_120 selected at runtime).
- Solver options: own deterministic Jacobi CG (default, bitwise reproducible); library solvers (AmgX, hypre-CUDA incl. AMS, PETSc,
  Ginkgo, cuDSS FP32 preconditioner) allowed where run-to-run bit identity is not required, accepted through the field guard receipt.
- The earlier GPU force candidate (pic-cuda-force-evolving-v3, NVRTC gather) is input material only: its two open review findings
  (stream ordering, async error teardown) must be closed, and it hard-codes compute_86.
- Placement in the gate order: after S1 (the CPU path is the reference). Expected effect: the push phase drops by more than an order of
  magnitude, the coordinator floor (~12-15 ms after hashing and snapshot rework) dominates, and fewer CPU cores are needed.

## 6. Gate order

1. S1: serial Stage 1 on the workstation.
2. S1b: new Hall versions (binary hashing, solver backend interface, distributed motion assembler) and native v5 (anchors, births,
   sinks, binary prefix), byte-identical to S1 at N = 1.
3. H2-H7.
4. GPU-resident step (5c) with bitwise parity to the CPU push, then N-parity (1,2,3,4,5,8) with cpu-pcg and the GPU path.
5. Workstation pilot (no restart).
6. Restart (A05/A11), cloud site admission, cloud runs.

## 7. Launch interface: portable, with core and GPU counts chosen per run

Requirement: on any owned machine or cluster, the core count, the exact cores and the GPUs used by a run are launch options,
changeable from run to run without editing code or wrapper scripts.

Command (one tool, same on every machine; version 1 is in tools/pf-launch with its README and site files):

    pf-launch [--site SITE.json] [--cores N | --cpus LIST] [--gpus none|auto|N|LIST] [--ranks-per-gpu K]
              [--cores-available] [--wait] [--dry-run] [--receipt PATH] -- PROGRAM ARGS...

- `--cores N` picks N free physical cores; `--cpus 4-7` / `0,2,5` names them; default `auto` = every physical core not reserved
  or leased. SMT siblings are never used as separate ranks (checked against `/sys/devices/system/cpu/*/topology/thread_siblings_list`).
- `--gpus none|auto|count:N|0,1|GPU-UUID,...`: default `auto` = every GPU not leased; `none` forces CPU backends; a bare integer is a
  GPU index (`--gpus 0` = GPU 0), a count is `count:N`. Ranks map to GPUs by
  `--ranks-per-gpu` (default: one solver rank per GPU); each rank gets its own `CUDA_VISIBLE_DEVICES`.
- One MPI rank per selected physical core, pinned (taskset / `--cpu-bind` / Slurm binding), one thread per rank, nice level from the
  site file. N may be any count, including 1, 2, 3 or 5.
- Inside a Slurm allocation the tool reads the allocation (SLURM_* variables, `srun`) instead of picking cores itself.
- Changing counts for a long run = checkpoint, then relaunch with new flags (needs the rank-independent restart after the pilot).

Site file (one per machine, `sites/<hostname>.json`; generated by `pf-launch --detect` and then edited):
- detected topology (physical cores, SMT siblings, NUMA), GPU inventory (index, UUID, name, compute capability, memory, driver);
- reserved cores and an optional reservation probe (on the workstation: physical cores 0-5 are allocated to the Spheral work and
  6-7 plus the GPU to Power-Flow, as of 2026-10-30);
- leases: per-core and per-GPU lock files (on the workstation core 6 = `pic-profile/testing-core.lock`, core 7 =
  `pic-profile/testing-core7.lock`, other cores `testing-core<K>.lock`, GPUs `testing-gpu<UUID>.lock`), held for the whole run;
- launcher (`mpirun` or `srun`), execution environment (frozen root via run-frozen on the workstation, container image elsewhere),
  nice level, environment variables;
- identity admission policy: strict on the workstation (kernel, driver, GPU UUIDs, frozen root digest must match), recorded-only or
  per-site records elsewhere.

Every run writes an execution envelope receipt: host, site-file digest, selected cores and GPUs, rank count, rank-to-core and
rank-to-GPU maps, GPU name/UUID/capability/driver, launcher argv, environment digest, root or image digest, start/end times, exit codes.
`--dry-run` prints the resolved plan and refuses if requested resources are reserved, leased or absent.
