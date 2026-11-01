Source v7.2 prepared-state successor
See API-v7.2.txt for the publication contract and limits.

Inherited v7.1 description (final v7.1 artifact committed 47899e2):
Source v7.1: explicit electron mass, isolated successor, 2026-11-01

Purpose
The sealed native chemistry v3.1.1 blob constructs electron mass as
(mass_eV * 1.602176634e-19) / (299792458 * 299792458). Its binary64 result
is 9.1093837138687491e-31 kg. Source v7 default is 9.1093837139e-31 kg.
This successor binds the actual native operand without replacing global
physical constants, rescaling heating caches, or changing existing defaults.

API and integration
surface_emission::Config appends electron_mass_kg and electron_mass_identity.
Both omitted: legacy behavior and configuration/state/checkpoint bytes remain
unchanged, subject to the recorded default differential qualification.
Nonlegacy mass must be finite, positive, and have a bounded nonempty identity.
This is a declaration mechanism, not physical qualification of arbitrary masses.

Use adapters/native_source_v7_mass.hpp:
  auto registry = full_time_chemistry_registry(native_configuration, options);
  auto states = full_time_chemistry_states(registry);
  cfg = native_source_v7::electron_mass_configuration(cfg, blob_path,
          expected_blob_sha256, registry, states);
The adapter independently loads the explicitly pinned chemistry v3.1.1 blob,
validates every species and arena-state mass/charge/nucleus/reference operand,
recomputes the native registry hash, and records native configuration/policy/grid
provenance in the source mass-policy identity. Electron lookup uses pfc-v3/e,
not an array position. Supported blob SHA256:
70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8

The numerical owner namespace remains electrode_heat_v7 for a separately built
replacement. The const source_patches(), source_patch(j), and source_settings(j) accessors
expose the actual configured support case for adapters. Indexed accessors return
independent value snapshots (no mutable alias or borrowed-reference lifetime),
remain valid after owner destruction, and throw std::out_of_range without
changing accepted state.

All consuming translation units must be rebuilt against this
runtime's private vendor/heat-v2 and vendor/v3 headers and private emission
objects. Do not link old emission objects or mix old and new Config/Packet
layouts in one program. Public owner headers remain under
native/mfem_electrode_heat_source_v7. Shared frozen sources are untouched.

Mass participates in free-electron supply and tunneling action, launch and
matching velocities, layer q/m dynamics, magnetic admission guards, canonical
handoff and kinetic ledgers. Persistent packets carry their declared mass.
Nondefault mass layers use the version-2 layer codec; default layers retain the
exact old codec. Owner configuration binds mass and provenance; restore checks
packet mass against the configured patch. Existing neutral patch mass remains
explicitly configurable and is unchanged by this patch.

Performance and remaining gates
Legacy mass keeps the existing optimized emission/cache path. Nonlegacy mass
uses the scalar reference kernel before any legacy cache. No speedup or full
run-time claim is made; the new mode may be slower. Before production adoption,
extend and qualify the immutable emission-cache key and all dependent operands
for declared mass. The hot-temperature law carries the same declared mass.

Chemical stock/free-gas reference offsets are NOT implemented here. Stock and
emitted reference offsets must cancel explicitly in a later count-bound adapter;
no chemical-reference difference may become artificial electrode heat, kinetic
energy, or work function. Wall response with nonlegacy mass is refused because
its legacy species/reference catalogue has not been mapped. SS304 vapor masses,
metal species/references, phase coupling and full native8 genesis qualification
remain separate integration gates. This patch alone does not qualify plasma
breakdown, a populated coupled pilot, or 60 ns evolution.

Qualification
run_source_v7.py admits only live one-core SimQueue jobs with explicit memory,
lease/cgroup/affinity checks and at least 5 GiB physical memory available. Normal
Power-Flow owner cores or explicitly borrowed raikou-mfem cores 4/5 are allowed.
Compilation and numerical tests execute through arch-base-20261028/run-frozen.
All compiled dependencies, archives, binaries, source .inc files, and the frozen
root manifest are receipted; system header paths resolve inside the frozen root.
No GPU changes or GPU qualification are claimed.

Example (scheduler arbitrates borrowing):
 simq submit -p powerflow -c 1 --mem 6G -t 00:30:00 --borrow-from raikou-mfem \
   --max-borrow-cores 1 --on-preempt freeze -C <this-runtime> -- \
   /usr/bin/python -I run_source_v7.py build LABEL
 simq submit -p powerflow -c 1 --mem 6G -t 00:20:00 --borrow-from raikou-mfem \
   --max-borrow-cores 1 --on-preempt freeze -C <this-runtime> -- \
   /usr/bin/python -I run_source_v7.py run LABEL

The actual native registry fixture uses sealed native v8.4 and chemistry v3.1.1.
It tests provenance/operand refusals and actual source births/kinetic energy.
Owner tests use that exact derived configuration, exercise handoff and active
layer checkpoint/thermal-clock binding, and refuse the legacy wall model.
A separate binary built from the original kernels is compared byte-for-byte
against default successor law, packet, layer, refusal and energy output across
56 temperature/field/work-function combinations. Inherited v6/v7 owner and
surface suites include the one-million-event wall ledger check.

c1: canceled after concrete compile failures (missing private material headers,
wrong blob member spelling and JSON alias); preserved logs. c2: all 110 checks
passed, plus equal original/successor default binary streams. c2 system-header
provenance used host paths and is therefore preliminary evidence only. Final
build/run labels and source/dependency seals are reported separately; do not
promote from c2 alone. No accepted physical state, deployed configuration,
paused simulation, or shared frozen dependency was altered.
