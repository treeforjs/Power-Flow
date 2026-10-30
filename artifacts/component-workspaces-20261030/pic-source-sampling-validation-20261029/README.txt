CPU-only native source sampling validation

Scope: source emission sampling only, M=64,32,16,8, 32 deterministic paired
seeds, 100 ps of common prescribed source events from the previously audited
checkpoint at absolute heating time 26.2999463372 ns. No coupled solver or GPU
work. No changes to the active run or shared repositories.

Reproduce with:
  bash /home/treefor/projects/runtimes/pic-source-sampling-validation-20261029/run.sh

The script serially extracts and verifies the unchanged native slices, builds
one executable with the installed frozen compiler, runs the bounded validation,
and analyzes the output. Every utility is pinned to CPU6 at nice15. The script
acquires the shared pic-profile/testing-core.lock lease before extraction,
build, replay or analysis, and exits75 without starting if that lease is busy.
Do not wrap
the script itself in another nice15: nice increments are additive. No new
dependencies, MPI, MFEM or GPU libraries are required. There is a 50 second
validation wall cap; an exception/nonzero exit means incomplete. Results are
written only to the owned source-sampling-validation output directory and are
replaced on reproduction. run-summary.json confirms source_only=true and
coupled_convergence=false. analysis.json records the final statistical guards.

The final results use the corrected CPU6 allocation. Earlier development runs
used CPU0 under the initial assignment; those results are superseded.

Sources:
  extract_native.py: byte-identical source extraction and input fingerprints
  native_slice.hpp: generated native routines and storage-only class host
  validate.cpp: event capture, common replay, moments and physical checks
  analyze.py: uncertainty, sparse-patch count-tail screen, report, provenance

Event binary format, events.bin (native little-endian host, no padding):
  uint64 magic = 0x5352434556543031, uint64 event_count
  Repeated per event:
    double absolute_time_s, double span_s
    uint64 native_seed, uint64 source_step
    int32 face (0 back,1 gap,2 left,3 right)
    int32 lower (0 upper,1 lower), uint64 patch_count
    double temperatures_K[patch_count]
    double mean_current_densities_A_m2[patch_count]

Events preserve the arguments actually passed to the native source sampler.
For a patch, current density times span recovers the integrated-current input
up to division/multiplication roundoff. Every count replays identical T/J/span
arguments. The thermal and reservoir clock is prescribed at10fs; it does not
reproduce the live adaptive trajectory clock. A separate5fs source quadrature
check reports total supply, flush counts, emitted charge and pending charge.
The refined event stream is never substituted into a count comparison.

Group numbers in CSVs: 0 back upper,1 back lower,2 gap upper,3 gap lower,
4 left upper,5 left lower,6 right upper,7 right lower,8 all sources.
Energy totals use electron volts times represented physical electron count.
The source_energy_eV column uses exact native relativistic kinetic energy;
nonrel_source_energy_eV checks the exact Gamma(2,kT) effusive expectation.
Slab charge is emitted source in |z|<0.5mm, not transported active density.
Patch_squared_error is summed per event before adding events; missed source
probability mass is averaged equally over events. Neither is a claim that
total charge is lost: selected patches carry the conserved event total.

Native source source/math equivalence is documented in native-equivalence.json.
The copied methods, constants, types, helper math, distributions and source
reservoir arithmetic are unchanged. Options are copied from the checkpoint.
IDs and cumulative ledger offsets are reset per event and do not affect source
sampling. Whole-executable/cross-toolchain binary equivalence is not claimed.
provenance.json verifies input and frozen-source before/after hashes, matching
sealed native-source digests, executable/code/result hashes and CPU-only
dynamic linkage. The active output directory is never opened by these scripts.

Read source-sampling-report.txt for results, confidence intervals, limits and
proposed observable gates for a later independently authorized coupled test.
