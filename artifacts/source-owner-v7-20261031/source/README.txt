Heat/source owner v7 -- electrode surface physics for 304 stainless foils with native oxide and adsorbed water
2026-10-31

Contract: API-v7.txt. Data, tags, citations, bounded and open items: DATA-v7.txt. Integration (not performed
here): INTEGRATION.txt. Sealed inputs (owner v6.1 and v6, heat/source v3 kernel, Hall heat-v2 source kernel,
arrival v1/v1.1 headers, native candidate003, hall-gpu-resident-v1 common layer) are read in place and hashed by
run_source_v7.py in every receipt.

What v7 is
  Owner v6.1 (versioned copy, namespace electrode_heat_v7) plus:
  * trial interval bound as a sealed configuration parameter (default 1e-14 s), with a first-order accuracy bound
    on emission/desorption exposures (typed refusal, receipt in the trial);
  * surface state: work function of oxidized SS304 with coverage-dependent water shift from the exact stock;
  * field-enhancement classes (exposures of the sealed planar SN thermal-field kernel at beta F), phase map;
  * emission births at the macroscopic matching plane (sub-cell displacement, image tail restored) and
    analytic sub-cell return in retarding fields;
  * wall responses replacing full absorption: Furman-Pivi SEE/backscatter (stainless), ion neutralization with
    potential/kinetic emission and reflection, neutral sticking/accommodation, ESD/ISD of adsorbed water through
    the approved quantized release (n = SR(D/w0)), optional photo hook; exact integer ledgers;
  * one wall kernel compiled for CPU and NVRTC, bit-identical, Philox4x32-10 keyed by event.
  With every option at its benchmark value v7 is v6.1 (digests, decisions, checkpoint bytes).

Build and run (receipted; every compile and run inside arch-base-20261028/run-frozen)
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v7.py build LABEL [normal|ubsan|asan]
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v7.py run LABEL [SANITIZER_LABEL ...]
  flock .../testing-gpu.lock flock .../testing-gpu-GPU-a7b9b458-....lock core-lease-20261030/any-core.sh --wait \
      /usr/bin/python -I run_source_v7.py gpu LABEL
  The lent cores 0, 1, 3 are admitted as well (lease held, taskset, nice 15, no RETURN-CORES-0-3 flag); the
  wrapper used for this qualification is /home/treefor/projects/.workspace/outputs/source-owner-v7-20261031/
  lent-cores.sh.
  Outputs: /home/treefor/projects/.workspace/outputs/source-owner-v7-20261031/LABEL/

Results of record (/home/treefor/projects/.workspace/outputs/source-owner-v7-20261031/)
  builds b3 (normal), b3u (UBSan), b3a (ASan): build/build-receipt.json, all passed (78 sealed input files hashed;
  the only v6 owner header in the binaries is compat-v6's alias header)
  run b3/run-111717 (normal + UBSan + ASan, lent cores 1/3): all passed, identical case counts in every variant
    test_surface_v7        33/33  (normal 48.7 s, UBSan 56.5 s, ASan 62.6 s)
    test_source_v7         15/15
    test_source_v6_on_v7    9/9   (the sealed v6 qualification against v7 at benchmark settings)
  run b3/run-111621 (normal, cores 6/7): same results.
  gpu b3/gpu-111501 (RTX 3080 Ti, NVRTC sm_86 cubin, GPU leases): test_wall_gpu_v7 5/5
  Key numbers (normal run):
    port        A02 probes: v7 verdicts and observations equal v6.1 (5/5); every digest, trial content and checkpoint
                byte of v7 at benchmark settings equals v6.1.
    interval    1e-14 default refuses 2e-14 with the v6.1 text; configured 1e-13: a 2^-45 s (28.4 fs) step publishes
                (birth, continuation, trial, commit); error-target policy bound 1.11e-14 s from a 0.9e11 /s rate:
                2e-14 refused (kind interval-bound), 1e-14 published with its receipt. Explicit-exposure order 0.998;
                error at the bound 5.0e-4 for eps 1e-3.
    emission    kernel v(f) within 0.0024 of Forbes (max 0.00236); 0 K current vs Murphy-Good FN within 4.9 %
                (f <= 0.35), equal to the neglected G(E) curvature within 1.4 %; MG temperature factor within 0.77 %;
                Richardson at F = 0 within 1e-6; Richardson-Schottky within 0.12 %; exact Fermi-Dirac law within 1e-8.
    hot         extended evaluate bit-identical to the sealed kernel <= 3000 K; 3200-5000 K: Fermi-Dirac law 1e-8,
                Richardson-Schottky 20 % bound, quadrature converged to 4e-15; validity map 256 points (32 refused
                by |F| > 5e9 V/m, no barrier collapse inside the domain).
    see         delta_max 2.051 at 292.5 eV (Table I: 2.05 at 292); 300 eV split 0.092/0.740/1.219 vs measured
                0.12/0.75/1.17; as-received SEY table 50-1950 eV max deviation 13.7 %, RMS 6.0 %; sampled yields,
                spectra (KS, gamma marginals per n, rediffused law), cosine law all within statistics.
    ions        TJS H on Fe within 0.25 % of the recommended table; channel statistics within 5 sigma.
    ledger      1e6 mixed events (2.7 s CPU): every per-event identity exact; int128 sums close; 1,246,906 births,
                1462 stuck, 48 energy-capped, 47,222 yield-clamped (oblique electrons), 375,025 outside-data flags
                (energies beyond the cited ranges, recorded).
    partition   bit-identical outputs and sums over 1/2/5/8 partitions, contiguous/strided/reversed.
    gpu         1e6 events: every SW7Out byte equal to the CPU; fixed-tree integer sums equal; 128- and 256-thread
                blocks and a repeat identical; kernel 1.26 s (FP64-limited consumer GPU) vs 2.7 s CPU.
    owner       class exposures bit-equal to J(T, beta F, phi) a_k A dt/qe; matching-plane births: K = E + Et - phi
                + F L exactly, ledger residual 0; retarding sub-cell return closes with net heat -1e-28 J; wall
                responses: 1000 hits -> 1223 births, stimulated hazard 0.221 releases 1.98e14 molecules as 198 leaves
                of w0 = 1e12 (thermal alone 0.24 molecules), stock exact, theta 0.80 -> 0.64, phi 3.479 -> 3.576 eV;
                phase-owner view: molten patch at 4000 K emits 1.73e9 electrons per 10 fs through the liquid class
                (phi 4.30 eV, thermionic-Schottky); water above 3000 K refused (surface-domain); 2e11-molecule
                residue written off at 2950 K with its binding-energy load; stale phase report refused; metal vapor
                recorded, no evaporation load.
    registry    165 parameters: measured 49, recommended 68, model_derived 28, scenario 19, benchmark 1
                (b3/run-111717/normal-test_surface_v7/parameters-v7.json).

Findings
  * The sealed v3 advance_source compares (end - t) with dt_s; when fl(t + dt) rounds up, end - t exceeds dt_s and
    a continuation is refused ("active source duration outside common interval"). It does not occur at 1e-14 in the
    fixtures but did at 4e-14 from t = 5e-9. Steps for which t + dt is exact in binary64 (dt a multiple of
    ulp(t)) avoid it; v7k keeps the sealed check (identical decisions at the v6.1 setting). INTEGRATION.txt D5.
  * The Furman-Pivi stainless angular parameters (from TiN/Al) exceed the model's own limits at oblique incidence
    (delta'_ts > M beyond ~18 deg, delta_e + delta_r -> 1 near 43 deg); clamped and counted per event.
  * The FP backscatter plateau (0.81) is ~3x the measured Fe backscatter coefficient above the fitted range; v7
    blends to the measured value between 1.1 and 3.93 keV.
  * No measured liquid steel work function, H2O+/H3O+ electron emission, ISD of water or beta distribution of
    rolled foil exists in the literature found; these are bounded scenario parameters (DATA-v7.txt).

Development scratch (not of record): /home/treefor/projects/.workspace/outputs/source-owner-v7-20261031/scratch,
research notes and downloads: .../research. Documentation (README, API, DATA, INTEGRATION) was completed after the
run of record; every code file hash in the receipts is unchanged.
