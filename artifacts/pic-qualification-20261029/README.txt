PIC qualification follow-up snapshot, 2026-10-29.

The full run remains explicitly paused at its preserved ~26.48 ns state. Nothing here is a production deployment or restart authorization.

Sources are first-party experimental work; frozen CUDA/cuDSS/NVRTC/MFEM toolchains and third-party binaries are not redistributed. Sources include the corrected v2 CUDA gather/Boris lifetime handling, bounded actual-RHS observer, passive native event observer and exact-native-operand analyzer, and copied-state particle tracker. Original absolute source paths and SHA256 values are recorded in manifest.json. Rebuild into a fresh prefix with pinned dependencies.

Fixed-midpoint GPU qualification passed for46726 captured particles, including shared cuDSS-primary context lifecycle and zero field/velocity differences. This does not establish an integrated timestep speedup. Median GPU kernel0.307776ms is contained within the reference-inclusive22.868514ms force pipeline, which includes20.1632045ms fresh independent CPU reference. Actual evolving CUDA-force integration remains separate.

The RHS observer passed full saved-state comparison over1ps and captured16 half-solves in8 populated advances with the unchanged two-refinement FP64 solver/residual contract; the matrices match the audited native CSR exactly. Raw captured vectors and fixed-stage fixture are retained so later solver investigations can be reproduced.

The event analyzer corrected native binary64 endpoint/prefix arithmetic without changing input logs, physical fluxes or tolerances. Its800 real-stream checks passed; unknown inherited source provenance and incomplete incidence-angle coverage remain unknown.

The trajectory replay is a separate short copied-state diagnostic and must not be confused with accepted progress of the paused full run. All conditional voltage, mirrored anode temperatures, coarse prescribed B, height-averaged2D space charge, absent ions/gas/collisions/phase changes and unresolved convergence limits remain. CHICAGO publications are reference studies, not validation of this model.
