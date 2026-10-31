Coordinator v6 -- one accepted coupled step with surface births, first-hit sinks and one thermal owner -- 2026-10-31

Owners (all read in place, sealed sources unchanged):
  native     candidate006 (native v6, identity aa4d8f03..., core archive libpic_native_v6.a 1c252595...): birth scope,
             first-hit sink records, binary flows motion record (no JSON motion route), anchors
  Hall       corrective v3 (runtimes/hall-corrective-v3-20261030): kinetic sinks with conductor/exterior banks, face
             newborns, ThermalV3 (the only thermal owner), CPU/GPU motion providers, cpu-pcg/gpu1 field backends
  source     heat/source owner v6 (runtimes/pic-electrode-heat-source-v6-20261030), one owner per Hall foil slice,
             quantized neutral release, loads to ThermalV3 through hall_source_binding.hpp
  arrival    arrival owner v1.1 (runtimes/pic-arrival-owner-v1.1-20261030): electron streams (n = SR(N/w0)) and
             quantized neutral streams (n = SR(D/w0)), counter-based RNG, native ID staging
  chemistry  parent coverage v2 (runtimes/pic-water-parent-coverage-v2-20261030): Song 2021 <= 1 keV, outer3 BEB
             1-10 keV, provisional RBEB above 10 keV (model-derived, tagged), refusal limit; composite channel
             registered at genesis (replaces the v3 ParentSubflow)
  predicate  native v6 accepts_closed_or_birth_step over the unchanged V3 common contract

Layout
  source/identity_v6.hpp     identities (species catalog, frame), physical setup (geometry, 8 foil-slice source owners
                             with 4 oxide-water patches each, 64 arrival streams), helpers
  source/joint_step_v6.hpp   observe, prepare_v6, can_publish_v6, try_publish_v6 (the step)
  source/genesis_v6.hpp      joint genesis with A09 admission (provisional Hall/source/arrival owners, native
                             create_full_time with configured caps, actual-census Hall bootstrap, real owners)
  tests/driver_v6.cpp        modes run | checks | vacuum; per-step ledger, digests, midplane frames, timings
  tests/alloc_hook.cpp       global allocation counter (zero allocations in the swap window are measured)
  run_v6.py                  build | run, receipted, every compile/link/run inside arch-base-20261028/run-frozen,
                             admission: cores 6-7 with their leases held, GPU lock files for GPU options
  tools/                     ledger_view.py, compare_runs.py (identity / CPU-GPU), render_pilot.py (PNG/GIF/tables)

One step (prepare_v6)
  observe all owners (live digests; the source seal is the coordinator-held accepted seal, valid only while every
  source owner and the arrival owner still have their recorded digests) -> common clock -> native outer
  (maximum_records 4096 >= 2M+6) -> first exposure (coverage v2) -> native birth scope (B_s) -> source field: the
  accepted start-of-step Hall E/B -> per slice: epoch and desorption quotes; per patch: emission exposure
  (source v6 prepare_patch: thermal-field spectrum at the Hall surface temperature and the one-sided field) ->
  arrival plan and staging at B_s -> electron births (kernel phase space, face point) and their Nottingham
  loads -> source trials for slices with leaves or loads (quantized leaves of exactly w0; exact debit n w0 =
  arrival debit) -> neutral release times bound -> source transition seal -> native birth credit, fork ->
  census and newborns -> Hall v3 closed proposal (newborns, source loads, sinks, banks, ThermalV3, provider and
  backend from the configuration) -> native admission of the binary flows record -> second exposure -> seal.
  can_publish_v6: native and chemistry can_publish, exposure/record binding, fresh observation, midpoint digests,
  native v6 predicate with the arrival evidence, source preflights against the live accepted Hall arrays, arrival
  commit checks replayed, Hall preflight_final last.
  try_publish_v6: native -> chemistry -> Hall -> source owners -> coordinator source seal (no allocation in this
  window, measured), then the arrival ledger commit (it allocates its record; its checks were replayed before) and
  the record of the live owner digests.

Decisions and stated approximations
  1. Source coupling is explicit: emission exposures and the source trials bind the accepted start-of-step Hall
     arrays (FinalFieldView of the live accepted E/B at preflight), not the step midpoint. One Hall proposal per
     step; no outer source/field fixed point.
  2. Electron emission is face-born. Count: source v6 exposure (kernel spectrum) -> arrival v1.1 electron stream.
     Kinematics: surface_emission::sample at the surface step (retarding field) or at the Schottky-Nordheim
     matching plane (extracting field); the packet starts on the foil face (Hall v3 face newborn, macro-charge
     from that foil's conductor bank). The sub-cell matching displacement and the source-layer image/return
     dynamics are not represented (the matching field work is recorded per step as unrepresented energy). The
     source-layer kernel (Active layers, PIC handoffs) is not used: Hall v3 refuses charged interior newborns.
     Nottingham energy enters ThermalV3 as a source-transfer load through the slice source trial.
  3. Neutral desorption: oxide-water Temkin/Eyring patches (coverage 0.33) on all four faces of every foil slice;
     quantized release with w0 per stream. A slice owner without leaf and load in a step has no trial (its stock,
     layers and energies are unchanged; the quantized remainder stays in the stock).
  4. Native foil planes: the lower-foil outer face and the left x faces are moved by 1-2 ulp (to 1.8e-3 m and
     2.0e-3 m) so that conductor-bank deposits at those faces put no rounding weight on vacuum nodes (see open
     issue 1). Every owner uses the same options.
  5. Emission kernel domain: the thermal-field kernel evaluates F <= 0 or 1 <= F <= 5e9 V/m. A patch whose
     one-sided start-of-step field lies in (0, 1) V/m (a sign change of the surface field) gets no electron
     decision in that step; the exposure is recorded as unresolved and bounded by the kernel values at F = 0 and
     F = 1 V/m (ledger columns arr_unresolved_emission_*). Pilot attempt 1 (build q2) refused at step 3849 on this.
  6. Qualification fixture (qual): work function 2.0 eV and electron w0 2e-23, neutral w0 2, 160 authored
     electron/H2O+ pairs (5 eV - 20 keV) around four 1e15 water packets and 10 pairs whose electrons hit the lower
     foil in the first step. It exercises every path and is not a physical configuration. Pilot fixture (pilot):
     4.3 eV, w0 1e6 for both kinds, empty genesis.

Build and run
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v6.py build LABEL
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v6.py run LABEL run NAME --fixture qual --steps 20
  flock GPU locks + any-core.sh ... run_v6.py run LABEL run NAME --fixture pilot --field gpu --push gpu --parity-every K
  run_v6.py run LABEL checks NAME           T1-T4, T6 and spot checks
  run_v6.py run LABEL vacuum NAME --steps N [--field gpu]
  Pilot (SimQueue): simq submit -p powerflow -c 2 --gpus 1 --mem 10G -t 10:00:00 -n pf-pilot-1 -C <this dir>
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 -- /usr/bin/python -I run_v6.py run q3 run
    <out>/run --fixture pilot --end-time 2.749788055492e-9 --field gpu --push gpu --parity-every 1000
    --frame-every 1000 --progress-every 1000 --checkpoint-every 50000 --stop-file <out>/STOP
  Builds of record: q3 = q4 (driver eb1936f1..., decision 5; q4 only reads the water binding through the neutral
  workspace path); q1 = q2 (abcc26cb...) before decision 5, used by pilot
  attempt 1 and the first qualification set.

Outputs
  qualification: /home/treefor/projects/.workspace/outputs/pic-water-common-commit-v6-20261031/qualification/
  pilot:         /home/treefor/projects/.workspace/outputs/pf-pilot-1-20261031/

Open issues
  1. Hall v3 bank residue: admission compares the double-word magnitude |hi|+|lo| of an accumulated bank value with
     the summed per-deposit residue bound; normalization of the sum can exceed that bound by about one ulp, so
     repeated rounding residue at one vacuum node refuses a later step ("conductor charge outside declared
     material/short"). Avoided at the foil faces by the plane choice above; the short/axial planes at
     z = -depth/2 cannot be made exact on this mesh, so repeated charged sinks there remain exposed.
  2. Source v6 Field::matches recomputes the whole-mesh digest on every call (about 10 ms per call on the 200 um
     mesh); a step with births in all 8 slices spends about 160 ms there. Steps without source trials avoid it.
  3. Emission exposures (32 patches) cost about 39 ms per step, the Hall host field preparation about 28-31 ms;
     the GPU field and push backends change the step cost by less than 2 % at this mesh and carrier count.
  4. The emission model is face-born (decision 2); the matching work is not represented.
