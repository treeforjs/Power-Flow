Complete saved PIC-state comparison

Scope
This isolated tool compares baseline-cpu6-run1 and cpu-candidate-cpu6-run1 against seed-current. It reads saved files, reconstructs the final RHS on CPU, and never solves a field system, invokes CUDA, signals a process or changes a source/run output. The completed CPU candidate and reviewed DOF patch are not edited. Python uses only its standard library. The C++ probes use the frozen compiler, headers and MFEM libraries. All execution uses CPU 6, nice 15 and one thread/job; builds and heavy checks acquire the shared testing-core.lock.

Entry point
compare_saved_pic.py BASELINE_DIRECTORY CANDIDATE_DIRECTORY SEED_DIRECTORY --output REPORT.json
Run with taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 ...
The Python CLI checks affinity/nice/thread settings and takes the lease internally. Do not wrap it in another flock on the same file. build_schema.sh and build_guard_probe.sh take the lease for compilation; invoke those scripts through the same CPU/nice/thread prefix. No dependencies are installed.

Exit codes: 0 = complete comparison passed; 1 = finite, well-formed data exceed an explicit comparison gate; 2 = invalid/nonfinite/missing/inconsistent data or a retained guard failure. JSON uses allow_nan=False. A failed comparison never reports passed. The fixture harness checks all three exit codes.

What is compared
- A C++ offsetof/sizeof probe binds the decoder to the frozen little-endian PIC2 ABI: 168-byte options, 200-byte diagnostics and 80-byte particle records. All 23 option components and all 25 diagnostic doubles are decoded explicitly. Padding is not physical state and is not compared. Unsupported ABI sizes are rejected.
- Magic, exact file size, particle resource guard, valid boolean bytes, unique/in-range identities, positive mass/weight, finite represented charge and kinetic energy, domain/conductor bounds, and finite metadata/CSV/JSON values are checked before numerical norms. Counts, identities, next ID, clocks, charge, mass, weight and all options are exact between paired final states.
- Reservoir spans, integrals and weighted-temperature arrays, checkpoint/meta clocks and voltage, seed hash/identity, run completion, solve count, scope and progress metadata are checked. Reservoirs are exact across runs; array norms are recorded.
- Every raw diagnostic and every seed-to-final public-ledger increment is compared. The checkpoint stores the raw ledger; diagnostics() derives kinetic energy and particle-charge defect. The comparator reconstructs both using saved particles and the original summation/ledger order before matching history.
- All history rows receive schema/finite/shape/clock checks. Only matching step numbers are compared. Endpoint increments are compared separately. The baseline has an extra wall-clock progress sample at step 403837; no candidate state at that step is invented or interpolated.
- Final frame_27.csv is checked for exact grid coordinates/order, time, solid markers and full row count. Potential and both electric components receive absolute, L2, RMS and global-scale relative norms. Both density components are independently reconstructed from each checkpoint, including the strict |z|<0.5 mm numerical slab and the full-depth average. Density integrals are checked against represented particle counts.
- Density differences caused by matched particles moving across occupancy bin/slab boundaries are reported explicitly rather than rejected by a near-zero pointwise relative threshold. This is occupancy bin analysis, not a reconstruction of crossing events. For these real replay outputs, all density cells are identical and no particle changes a bin or slab membership.

Retained residual guard
residual_guard_probe.cpp reconstructs the Cartesian Q1 mesh, diffusion matrix, charge deposition, essential conductors and constrained RHS using frozen MFEM and the production assembly order. It does not perform a linear solve or link the cuDSS wrapper. Persisted residuals must remain finite/nonnegative and satisfy the unchanged CUDA acceptance formula:
    residual_norm <= max(1e-30, 2e-13 * rhs_norm)
The final RHS norms agree at 8.227406016442554e-5, giving a 1.645481203288511e-17 C/m guard. Saved baseline/candidate residuals are 4.9391429701803693e-20 and 4.946397827285137e-20 C/m; the seed also passes. No production tolerances or guards are altered.

Tolerance rationale
These are analysis gates for a short replay, not changes to the runtime or claims of convergence.
- Particle position: the prior checkpoint regression limit of 1e-12 m, or 4e-8 of one 25 um cell. Actual maximum is 3.469446951953614e-18 m; L2 is 6.504208955270767e-17 m.
- Proper velocity: the prior global component-scale relative limit of 1e-10. Actual maximum difference is 5.287583917379379e-7 m/s, with global-scale relative Linf 9.931744773116695e-15 and L2 1.840346983681717e-5 m/s.
- Fields: 1e-10 times each component global amplitude, with a 256-epsilon cancellation floor based on potential amplitude (divided by cell size for derivatives). The full field-error budget over this 1 ps window corresponds to at most about 8.83e-4 m/s electron impulse, below the 5.32e-3 m/s velocity comparison budget. This is a conservative consistency budget; it is not a bound on physical-model or discretization error.
- Charge/energy ledgers: 1024 epsilon times the relevant full charge or energy scale. Electrode/Gauss comparisons use electrode charge, while particle/continuity/charge-balance comparisons use particle/source/sink charge. Kinetic and source/sink energy use particle-energy scales independently of the much larger field energy. A cancellation residual near zero is never used as its own relative denominator.
- Reservoir values and options are exact. Density reconstruction uses 256 epsilon times the larger of the component scale and the density equivalent of one physical electron in the sampling volume. All density differences are also checked through particle identity, weights, positions and bin assignments.
- Clocks use floating-point roundoff on the global time with a 1e-23 s floor; paired final local clocks and common-step history clocks are exact.

Actual field/density effect
Potential: maximum 9.167706593871117e-10 V; L2 1.198485706178077e-7 V; relative L2 6.740288763557199e-15.
Ex: maximum 2.9103830456733704e-6 V/m; L2 1.4702135477350596e-4 V/m; relative L2 1.3969903660292323e-13.
Ey: maximum 3.0566006898880005e-6 V/m; L2 1.721026622059537e-4 V/m; relative L2 4.398728535745734e-14.
Both density differences: exactly zero at all 129600 cell centers.
Raw field-energy difference: 3.2959746043559335e-17 J. Gauss-defect difference: 6.418921553411495e-22 C, against an electrode charge L1 scale around 5.63e-8 C. Public kinetic-energy difference: 1.5146129380243427e-28 J. Charge source/sink increments, represented weights and reservoir arrays match exactly.

Seed accounting
29537 particles initially; 29613 finally. Next-ID increase is 192, 29515 seed IDs survive, 22 seed IDs are absent, and 98 newly allocated IDs are present at the end. Assuming the inspected source remains the only allocator, these counts imply 116 total removals, including newly emitted particles removed before the final snapshot. The comparator does not infer their individual crossing times, paths or energies.

Fixtures and limitations
Synthetic fixtures cover every floating-point checkpoint field with NaN and infinity, plus malformed counts/IDs/booleans/mass/weights, truncation/trailing bytes, nonfinite CSV/history/metadata/JSON, finite field/density/trajectory mismatches and residual-guard failure. A small field roundoff perturbation is accepted. A separate conservative rebin fixture permits a 4000 cm^-3 occupancy change caused by a roughly 4e-16 m particle position shift across a grid boundary. Fixture states test file/consistency/numerical behavior and do not assert a physical Poisson solution.

The saved checkpoint has no nodal potential, residual error vector, per-step maximum-relative statistic, trajectory/crossing log or complete intermediate history. The comparator cannot recreate those absent records. CSV checks concern final cell centers only. RHS guard reconstruction checks the persisted norm, not a newly evaluated residual vector. The 1 ps replay is an integration regression check; grid, timestep, physical and long-run validation remain separate.

Reviewable artifacts
comparison.json: all actual norms, values, tolerances, checks, history coverage, reconstructed guards and input SHA-256 manifest.
fixture-tests.json: fixture outcomes and CLI exit-code tests.
boundary-density-fixture-result.json: the accepted conservative rebin example.
abi-schema.json, schema_probe.cpp, residual_guard_probe.cpp: ABI and guard reconstruction sources.
compare_saved_pic.py and test_comparator.py: comparator and fixture harness.
