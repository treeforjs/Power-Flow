Hall corrective owner v3 -- 2026-10-30 (pilot items H3 and H5)

Entry point: source/root_api.hpp, namespace hall_corrective_v3 (GPU backends:
source/gpu/gpu_backends.hpp). Versioned copy of Hall corrective v2
(runtimes/hall-corrective-v2-20261030, unchanged and read in place). Sealed inputs used unchanged:
Hall commit 8dddec0 (axial/foil field headers, kinetic_boundary.{hpp,cpp}, heat-v2 bindings),
hall-gpu-resident-v1 common/ and push/, native candidate003 (coordinator v5 link) and native
candidate006 (native v6, read in place by native-v6-link/). source/v2_compat.hpp maps
hall_corrective_v2 -> hall_corrective_v3 so unchanged v2 clients (coordinator v5) compile against v3.

1. Kinetic sinks (H3)
1.1 Model (sink_core.hpp, shared by the CPU provider and the NVRTC kernels from the same bytes).
  The represented step of a carrier is the straight segment from its start x to x + D, with D the
  full displacement of the sealed work-kick push over its flight F (dt for residents, end - tau
  for newborns). A carrier is truncated when the sealed push reports a field-surface crossing or
  first_kinetic_hit(x, D) of native-finite-target-kinetic-policy-v1 reports a hit (tie order
  x y lower upper axial, later exact tie wins; the face-aware variant for a newborn that leaves a
  foil face). The carrier ends at the hit point x + f D with the hit coordinate set to the plane
  (the native v6 arithmetic, no snapping); flight F f, discarded flight F - F f. Its path current
  is the straight path deposit from the field start to the hit point; the incident kinetic energy
  is K0 + q E_path . (hit - start), E_path the path gather of that deposit (so the motion work
  guard closes with the truncated deposit), and the incident proper velocity has that energy
  along D. A neutral carrier drifts straight with D = F u / gamma and keeps u; a neutral first
  hit is a sink record without charge (v2 refused it).
1.2 Banks. In the same field step the macro-charge w q of every charged sink moves from the free
  charge to the conductor bank (surfaces 0 lower foil, 1 upper foil, 3 short bridge) or to the
  exterior bank (2 outer x/y faces, 4 axial planes) by the Words CIC deposit at the hit point
  (census_charge.hpp bank_transfer; order: sinks then newborns, ascending native ID). The
  endpoint continuity guard of Model::prepare checks after + bank - before - div(J) dt; absolute
  Gauss is admitted with the bank supports. CIC weight that rounding puts on a node next to the
  bank support (a plane's cell coordinate a few ulps off its node plane) stays in the bank as
  residue, admitted only up to 1e-12 of the macro-charge per deposit and bounded per node in
  ChargeLedger::residue_C (part of the ledger record and checkpoint).
1.3 Energy. Collected kinetic energy fl(w KE(incident u)) on a foil is wall heat applied by the
  thermal owner (SurfaceLoad producer "hall-v3/kinetic-sink-wall-heat", one per sink, sequence
  (revision + 1) 2^20 + i); short-bridge and escaped energies are recorded (no thermal owner).
  StepEnergy closes field + plasma + Joule + wall heat + external loads (residual against the
  heat roundoff bound).
1.4 Newborns (native v6 birth model). StepInputs{newborns, loads, birth_scope}; a newborn is a
  census row with birth time tau; it moves only over [tau, end]. Its start is a strictly interior
  vacuum point or a point on exactly one foil face with the velocity leaving the metal (edges,
  corners, starts in metal and velocity into metal are refused). A charged face newborn takes its
  macro-charge from that foil's conductor bank at its start in the same step. A charged interior
  newborn is refused (it needs the source path current and charge impulses, not in this version).
  birth_scope = true (native step with a birth scope, arrival committed every step) makes the
  record the flows record even with no newborn and no sink; the physical state is unchanged.
1.5 Records. SinkRecord (motion_receipt.hpp): first the native v6 FullTimeSinkRecord fields with
  the same values (native_id, kind 1 collected_at_conductor / 2 escaped_exterior, surface, axis,
  plane_m, fraction, segment_m = D, hit_point_m, incident_proper_velocity_m_s, flight_s,
  discarded_flight_s), then the Hall side (packet, state, conductor, field-frame bank point,
  charge, mass, weight, kinetic_J, weighted_kinetic_J, guards). NewbornRecord: native_id,
  birth_time_s, flight_s = end - tau, origin (foil face or interior), start, weight, charge,
  birth kinetic energy.
  Binary motion record: a receipt without newborns, sinks and birth scope is v2-identical (plain
  v3 schema canonical-full-time-total-field-motion-v3-binary, bytes identical to Hall v2);
  otherwise the native v6 flows record (API-006.txt section 4): tag
  canonical-full-time-total-field-motion-v3-flows-binary, coverage
  all_native_carriers_births_and_first_hit_sinks, the plain fields, then "sinks" (u64 count; per
  sink u64 native_id, kind, surface, axis; f64 plane, fraction, segment[3], hit_point[3],
  incident_u[3], flight, discarded) and "births" (u64 count; per newborn u64 native_id, f64
  birth_time, flight). native-v6-link/hall_v3_flows_adapter.hpp maps a receipt to
  raikou::mfem_backend::FullTimeMotionRecordV3 and binds it with bind_total_field_motion_v3 (the
  native v6 adapter bind_hall_total_field_motion_v3 maps only the plain record of a v2 receipt).
1.6 Ledgers (sink_ledger.hpp; Owner::kinetic_ledgers(), cumulative, checkpointed). SinkLedger per
  surface (collected/lower-foil, collected/upper-foil, escaped/outer-xy, collected/short-bridge,
  escaped/axial): exact count, exact weight, exact +e and -e counts (ExactParticleCount), exact
  sum of fl(w KE), compensated charge and energy words. BirthLedger per origin (lower-foil face,
  upper-foil face, interior). Per surface the values equal native's FullTimeSinkLedger surface
  totals exactly (native-v6-link test).
1.7 GPU. push_bodies_v3.hpp / push_kernels_v3.cu run sink_core.hpp per carrier (newborn flights,
  face starts, neutral drifts); GpuMotionProvider "gpu-push-v3-sinks" builds the same receipts.

2. One thermal owner (H5)
  ThermalV3 (thermal_v3.hpp) is the only owner of T, sigma(T) and Joule heat: the v2 per-slice
  heat-v2 Couplers, unchanged and bit-identical when no loads are given.
  - ThermalDeclaration: source patches declared heat-only on Hall slices (bound into the Hall
    material and configuration identities); declared slices evolve bit-identically with no loads.
  - External surface loads: StepInputs.loads / Owner::prepare_vacuum(dt, loads), SurfaceLoad
    {body, slice, channel, heat-v2 Transfer}. Each load is applied exactly once (heat-v2
    producer/energy-owner watermarks; a replay, a duplicate or an energy owner shared by two
    producers is refused). Joule channels or owners and Hall-reserved producers are refused
    (external-joule, reserved-producer); volume loads and unnamed channels are refused.
  - Admission before the engine (typed ThermalRefusal {kind, quantity, identity, observed, bound}):
    nonfinite-load, load-exceeds-material-range, net-heating-exceeds-capacity,
    net-cooling-exceeds-capacity; engine exceptions on a step with loads are rethrown typed
    (thermal-engine); closure check load-accounting. The owner is unchanged after a refusal.
  - Source owner v6 (runtimes/pic-electrode-heat-source-v6-20261030) reads the Hall slice
    temperature through a const view, returns its loads (emission/collection layer exchanges,
    desorption binding and release energy) and has no engine and no Joule term; its binding
    audits that no second owner reports Joule heat or a thermal change.

3. Identity changes relative to v2
  api "hall-compensated-field-owner-v3"; configuration record "hall-v3-field-configuration-v1"
  (adds the sink/newborn/exterior-support rule texts and the thermal declaration); charge ledger
  record "hall-v3-charge-ledger-v1" (adds residue_C); checkpoint adds kinetic_ledgers,
  kinetic_ledgers_sha256 and bank_residue_C. Field, heat and motion arithmetic of a step without
  newborns, sinks, birth scope and loads is unchanged; the motion record is the v2 record. The
  coordinator v5 identity mapping v2 -> v3 (111 digests) is one-to-one with no non-digest
  difference and only the api text changed.

4. Layout
  source/sink_core.hpp, sink_ledger.hpp, thermal_v3.hpp     sinks, ledgers, thermal owner
  source/field_contract.hpp, field_model.hpp, census_charge.hpp, owner.hpp, root_api.hpp,
         checkpoint.hpp, motion_receipt.hpp, motion_provider.cpp, encoding.hpp, solver_backend.hpp
  source/gpu/                                                GPU field and motion providers (v3)
  source/v2_compat.hpp                                       v2 namespace alias for v2 clients
  tests/test_hall_v3.cpp, test_hall_v3_gpu.cpp, sink_fixtures.hpp, sink_dump.hpp, ...
  run_hall_v3.py                    build | cpu | gpu | multistep (receipted, frozen launcher, leases)
  run_coordinator_parity.py         coordinator v5 T1-T7 on Hall v3 vs the v2 reference run
  native-v6-link/                   flows adapter, native v6 admission test, run_native_v6_link.py

5. Results of record (/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/)
  b7/cpu-211226       components 30/30 and UBSan 30/30 (0 runtime errors): all five surfaces,
                      side/bottom faces, node planes, corners, grazing, first of two surfaces,
                      near miss, ion impact, neutral first hit, newborns on faces and inside,
                      newborn sinks within end - tau, 400 random carriers near the lower foil
                      (19 sink fixtures, 152 sinks; absolute Gauss, charge-bank partition,
                      energy closure incl. wall heat, native v6 record conformance, publication);
                      partial-path current = truncated deposit; omitted/halved bank transfer and
                      collected charge on the exterior bank refused; newborn refusals; external
                      loads exactly once; double Joule; typed refusals; declaration; empty
                      birth scope; checkpoint with banks and ledgers.
                      v2 parity: 14 steps of the v2 physical-state scenarios (z = 0 and z = 1 mm
                      pairs, w = 1e9, node planes, 1e3 pairs, vacuum) bit-identical; motion records
                      equal Hall v2 and the sealed 8dddec0 provider on 24 fixture entries.
  b7/gpu-211241       RTX 3080 Ti (sm_86), two fresh processes: GPU sink records, truncated
                      paths, current, endpoint charge, free-charge words and closed proposals
                      bitwise equal to the CPU provider on 21 sink fixtures (286 sink records),
                      repeat in-process identical, v2 fixture set (12) identical, GPU push + GPU
                      field closed loop over 6 steps with births and collection at both foils
                      (Gauss, energy, preflight, publication); teardown drained, not quarantined.
  b7/multistep-211257 Hall closed loop (CPU), 20 cathode electrons (2 MeV, w 1e6) born on the
                      lower-foil face every 50 steps, 2732 steps until the census is empty:
                      18 collected at the upper foil (wall heat 5.767e-6 J = their incident
                      energy), 2 escaped at the outer y face (6.41e-7 J); every step: absolute
                      Gauss (max 0.044 of bound), total charge change <= 3.7e-28 C, free -> bank
                      partition exact to 1e-12, energy residual <= 2.9e-15 J (heat roundoff bound
                      6.6e-12 J), native v6 record conformance, preflight and publication;
                      71 ms/step; per-step ledger step-ledger.jsonl.
  coordinator-parity/run-p3-221714 (also run-p3-211730)
                      unchanged coordinator v5 sources on Hall v3 (v2_compat.hpp): T1 101, T2 28,
                      T3 35, T4 10, T5 37, T6 6, T7 9 checks, all passed; T1 byte-identical in two
                      fresh processes; Hall field words, charge ledger and 8 heat slices
                      bit-identical to the Hall v2 reference run for P1, P2 and both z = 1 mm
                      Stage 1 populations; 111-digest identity mapping, no non-digest difference.
  native-v6-link/n2/run-211210
                      15/15 checks, native v6 (aa4d8f03...) compiled from its sealed sources:
                      plain record admitted; sinks at all five surfaces admitted with native sink
                      records equal to the Hall records and native per-surface totals equal to the
                      Hall SinkLedger exactly; native first-hit port equals every Hall sink;
                      upper-foil face births, sink-prone births (sinks within end - tau),
                      interior births, births + sinks in one step admitted; five tampered records
                      refused typed (sink_outside_step, sink_not_crossing, sink_record,
                      motion_record x2); 6052-step run with 103 electrons born on the lower-foil
                      face, all collected at the upper foil, native and Hall ledgers equal on
                      every step (2.6 ms/step); arrival v1.1 quantized neutral leaves through the
                      native birth scope over 40 steps (99 leaves of exactly w0, native stream
                      debit = n w0 exactly).

6. Open
  - A charged newborn at an interior point is refused: its charge comes from the source path
    (emission path current and charge impulses), which this version does not carry.
  - Neutral sink kinetic energy on a foil is applied as wall heat (full accommodation); no
    reflection or partial accommodation model.
  - The bank residue (rounding-level CIC weight next to a bank plane) stays in the bank; it is
    bounded per node and part of the ledger, not redistributed.
  - The native v6 JSON route (v2 JSON record) cannot carry births or sinks; flows go through the
    binary adapter only.
  - Cross-architecture GPU bit identity (sm_89, sm_120) not shown.
