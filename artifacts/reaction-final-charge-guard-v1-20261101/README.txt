Final post-chemistry charge guard v1

SimQueue616 rebuilt this guard/test and the Hall motion implementation with
ASan/UBSan. It rejects the actual native recombination counterexample from609,
accepts unchanged charge with changed velocity, and permits an extra neutral.
All accepted native/RNG roots remain unchanged after the discarded candidate.
The native archive and external libraries are pinned, not rebuilt by616.

require_consistent(mesh,options,motion_endpoints,final_native_census,cap)
compares actual final native CIC charge with Hall motion endpoint charge using
Hall's existing same_net_charge rounding bound without increasing it. It throws
before publication if they differ. Its receipt binds the actual final census,
including changed kinematics, and reports maximum admitted per-node difference.
It neither repairs fields nor supplies reaction current or missing field work.

This is a root-side numerical predicate, not a publication capability. Both
inputs must come from live authenticated Hall/native candidates. The coordinator
must retain/preflight those exact candidates after this check and before any
no-fail publication. Wiring into the new common coordinator remains required;
Hall preflight_final alone does not read actual post-stage2 native carriers.

The counterexample uses16 manufactured carriers withweight1e15 only to expose
the mismatch. Its charge amplitudes are not physical prediction or convergence
evidence. See artifacts/reaction-local-charge-audit-v1-20261101 for the separate
actual native reproduction and explicit diagnostic-target-energy scope.
