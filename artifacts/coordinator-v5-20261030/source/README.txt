Coordinator v5 — Stage 1 flow on Hall corrective v2, CPU and GPU backends — 2026-10-30

Coordinator v5 = the coordinator v4 flow (runtimes/pic-water-common-commit-v4-20261030) linked against Hall corrective v2
(runtimes/hall-corrective-v2-20261030) instead of the 8dddec0 corrective owner and the sealed Hall motion TU. Two native
builds of the same sources:
  n3  native v5 (candidate005, source identity cca72822..., libpic_native_v5.a 71572911...)  main line
  c2  native candidate003 (c045c42d..., libpic_native_candidate003.a d2568af9...)            Stage 1 parity
(the native choice is a build option of run_v5.py; -DCOORDINATOR_NATIVE_V5 selects the v5 genesis configuration
overload and history status). V3 native_parent_subflow.hpp and common_commit_contract.hpp are included unchanged.
Fixture: back at the plan position z = 0 (Hall node plane) for P1 and P2; T1 additionally runs both populations at the
Stage 1 position z = 1 mm for the bit-for-bit comparison with the Stage 1 receipts.

Changes relative to v4
 - Hall identities from Hall v2 (binary, cached; motion mesh and material digests unchanged values, no per-step JSON).
 - Field-solve backend (cpu-pcg | gpu1-factored-cg) chosen by configuration and bound into the Hall configuration and
   therefore into the coordinator configuration identity; motion provider (cpu | gpu) chosen by configuration and
   recorded in the execution envelope and audit only (records are bit-identical).
 - prepare_closed_parent_joint_v5 takes the motion provider; PhaseTimesV5 carries the Hall proposal split (motion,
   census digest, motion v3 digest, v2 JSON record seal, deposit support, field prepare, solve, independent residual,
   field seal, heat). The trial audit can keep only the motion seals (multi-step runs) instead of the full record.
 - Coordinator record tags whose encoding did not change keep their v4 text (species, frame, configuration, inactive
   source, common clock, genesis clock); the trial audit and genesis record schemas are -v5.
 - Carrier caps are configuration: Hall CarrierLimits, native v5 FullTimeConfiguration, parent chemistry carrier limit
   (above 1e5 only). Genesis requires the native/Hall/chemistry caps to agree.
 - Test modes added: gpu-step (acceptance b), multistep (joint, per-step ledger and per-phase timings), hall-loop
   (Hall v2 closed loop without native/chemistry, for carrier counts or memory where the joint step is not run),
   sweep-point (diagnostic). T3 attributes every crash and every injected failure that was swallowed into a
   publication to the component of its innermost frame.

Layout: source/{identity,genesis,closed_parent_joint_step}_v5.hpp; tests/{fixture_v5.hpp,test_stage1_v5.cpp,
alloc_hook.cpp,declaration_v5.cpp}; run_v5.py (build [003|005] | cpu | gpu | multistep | hall-loop | points);
compare_stage1.py (Stage 1 comparison and identity mapping); run_series.sh; series_summary.py.
Outputs: /home/treefor/projects/.workspace/outputs/pic-water-common-commit-v5-20261030/ (RESULTS.txt there).
