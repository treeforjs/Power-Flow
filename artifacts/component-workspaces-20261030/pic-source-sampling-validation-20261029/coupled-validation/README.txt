Isolated64/32/16 PIC validation harness -- GPU EXECUTION DISABLED

Build and run the CPU-only preflight:
  bash /home/treefor/projects/runtimes/pic-source-sampling-validation-20261029/coupled-validation/build-and-preflight.sh

Every compiler/preflight runs serially on CPU6 at nice15 with one thread.
The wrapper holds pic-profile/testing-core.lock for the whole job and exits75
without work if busy. It has a180s cap, inside the proposed30min CPU utility
budget. Do not wrap it in another nice15 increment. No dependencies are added.

Results are fresh attempt-NNN directories under:
  <workspace>/outputs/pic-profile/source-sampling-validation/coupled-preflight
latest.json identifies the most recent successful attempt. Each contains
report.json, three identical restart copies, native/source checks, negative
fixtures, metric/budget checks, build/link logs, a disabled39-case campaign
plan and a seal of the result/code hashes. Attempts remain distinct.

Execution fences

coupled_trial --run, unknown activation flags, GPU environment requests and
enabled configuration requests are rejected. There is no runtime override.
The compiled run_trial implementation starts with an unconditional execution
fence before MPI initialization, GapPic construction, output creation or
field/trajectory work. gpu_disabled.hpp supplies a second unconditional
throwing backend fence. Neither executable links CUDA, cuDSS or NVML.
The preflight performs no MPI initialization, real GapPic construction,
field solve, trajectory advance or GPU probe/context creation. The native
checkpoint loader is deliberately NOT called: it solves fields immediately.

Enabling actual coupled work would require a separately reviewed SOURCE AND
BUILD change plus an assigned resource window. Flipping a config/environment
variable cannot do it. The compiled disabled backend is not a qualified GPU
backend, and no runnable coupled trace/ceiling or memory measurement exists.

Implemented components

prepare.py verifies frozen source against sealed-files.json and creates only
owned copies. It extracts the native thermal/field/voltage definitions and
the complete timestep/source-reservoir statements verbatim. The candidate
native PIC differs by one passive absorption observer and the GPU-stub
include. Header/implementation reverse checks prove those are the only
changes; serialized options/diagnostics/particle ABI stays168/200/80 bytes.
The original MPI/cuDSS residual-guard file stays byte-identical.

trial.cpp contains the disabled coupled loop:
  - exact copied restart state and fresh tag id>=copied next_particle_id;
  - count64/32/16 and declared seed ranges/unsigned offsets;
  -64 reference capture,32/16 source-and-step replay, or native adaptive mode;
  - replay rejection at a stricter native candidate stability bound;
  - native reservoir clock/threshold/reversal and per-event/cumulative charge;
  - pre-advance injection timestamps and SHA-bound trace identity;
  - source, full/slab particles, fresh sinks, native ledgers and space-charge
    potential/field diagnostics at0,5,...,100ps;
  - wall/step/RSS/output/duration caps and at most two trial checkpoints;
  - incomplete/invalid classifications without accepting partial physical time.
Actual GPU memory remains zero in this build; later memory/step-ceiling
qualification is explicitly absent. Frozen guard values are never loosened.

records.hpp defines bounded binary input/event/trace codecs, exact restart
checks, resource controls and output-path containment. Trace continuity allows
only floating-point clock-roundoff; replay still requires the candidate's
actual accepted time to equal the recorded time. A source-only or manufactured
trace cannot be used as a verified coupled reference. The existing source-only
event file records source integration interval ends; it is used ONLY as an
emission-math fixture, not imported as a coupled step/injection schedule.

source_preflight.cpp executes unchanged native emission routines in a storage
host, separately from the native MFEM class. It checks744 identical event
arguments at64/32/16 and three deterministic paired seeds. Source prefix checks
concern positions/velocities/scaled weights. Ephemeral per-event storage IDs
are not a transported population or a restart-continuation test.

metrics.py reads CVFRAME1 diagnostics and calculates physical-weight full/slab
and five-region counts/densities, fresh fractions, N_eff, kinetic means/totals/
quantiles, CDF distances, unnormalized coarse-density L1 and space-charge RMS
differences with the proposal's floors. Native hit4 is an axial subset of
escape; it is added to hit2 once for fresh total escape, never counted twice.
Raw particles retain3D positions/proper velocities for further diagnostics.

Confidence kernels enforce six paired seed triplets before acceptance, signed
t bounds, fixed absolute margins for preclassified small signals, deterministic
seed-level bootstrap distance bounds, completeness of both stages/all primary
metrics/three contrasts, fresh20% region gates, and unresolved distribution
statistics. Empty or partial input remains inconclusive. Bootstrap coverage at
six seeds is approximate; no physical confidence conclusion is established.
An intake bundle must retain all run/clock/identity/metric provenance. It is
not an independent execution attestation or permission to deploy.

analyze.py accepts an already exported observable bundle and writes a new
owned result; it acquires the analysis lease itself. It never launches a
simulation. For heavy use:
  taskset -c6 nice -n15 python3 /home/treefor/projects/runtimes/pic-source-sampling-validation-20261029/coupled-validation/analyze.py INPUT_BUNDLE NEW_OWNED_OUTPUT
The data-intake records are documented by metrics.evaluate_bundle:
  runs: stage,count,seed,complete,native_guards_passed,accepted_duration_ps
  scalar_records: stage,endpoint,metric,count,seed,value,absolute_margin?
  distance_records: stage,endpoint,metric,candidate,reference,seed,value
  information_records: stage,count,seed,partition,region,near_zero,fresh_fraction
  data_status: manufactured or owner_supplied
Scalar names come from frame_observables/primary_names. Distances come from
paired_distances.50-100ps values use trapezoid over the fixed physical times;
time samples never add independent seeds. Near-zero margins must be assigned
by the proposal's reference-only rule and recorded; they are not fitted to
candidate differences. A pass means a candidate for scientific review only.

budget.py is pure campaign ACCOUNTING, not a queue or resource reservation.
It rejects duplicate attempts, invalid count/seed/stage, a40th attempt,
stage budget overruns and exposures beyond726,000 steps/32,670 wall seconds.
The self-test exercises the39-slot envelope in memory and launches zero jobs.

Binary diagnostics

All formats use the frozen native little-endian ABI and explicit integer widths.
CVTRACE1 magic=0x4356545241434531:
  uint64 magic; uint32 coupled_verified,seed; double initial_s,end_s;
  uint64 step_count;64 ASCII SHA256 identity bytes.
  Per step: uint64 step;double pre_advance_time,dt;uint64 event_count<=8.
  Per event: double injection_time,span;uint64 native_seed,step;
  uint32 electrode(Lower0/Upper1),surface(Gap0/Back1/Left2/Right3);
  uint64 K;double T[K],J[K]. K<=4096. Full input/config identity is bound
  by campaign-plan.json and its seal; the verified flag alone is not proof.

CVFRAME1 magic=0x43564652414d4531:
  uint64 magic;double absolute_time;uint64 particles,fresh_start;
  native GapPicDiagnostics(200B);native GapParticle[particles](80B each);
  double fresh_sink_charge[5],fresh_sink_energy[5],source_electrons,pending;
  uint64 native_cell_count;per cell double(phi_sc,Ex_sc,Ey_sc,vacuum_volume).
Cell order is native j*360+i. The field correction is the height-averaged
Poisson contribution at native cell centers; it is not a new3D field solve.
The1mm slab and25mm full population overlap and must not be added together.

Checks and limits

CPU tests exercise restart/hash/ABI/ID/positive-state/reservoir integrity;
identical initial particle bytes across counts; gyro/plasma/source/input-clock/
empty-to-emission bounds with manufactured conditions; source charge/prefixes;
trace and frame roundtrips; fresh sink classes; weighted metrics; confidence
and cap controls; and negative corrupted-input/path/activation/guard checks.
No manufactured trace, B bound, spectrum or confidence result is promoted to
coupled physics evidence. Full evolution, selected field gradients, wall/RAM/
GPU cost and scientific convergence require later authorized execution.

The compiler reports inherited misleading-indentation warnings in verbatim
native snippets. Those statements remain unchanged; owned driver warnings
are fixed. No simulation failure is hidden by changing a native source line.

The active60ns output, its sealed executable, shared repositories, timestep/
residual/conservation guards and queue/core reservations remain untouched.
