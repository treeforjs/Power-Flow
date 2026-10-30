# Native electrode heat/source binding v2

Use [API-v2.txt](API-v2.txt) for the integration contract. This separate module adds a certified external-field-work source port, persisted surface origins, typed source-layer checkpoint payloads, direct sealed SS304 material-v2 callbacks and one conditional oxide-water inventory law. Generic v1 remains unchanged.

The water model is a conditional ideal oxide surrogate. Its literature provenance and explicit scenario choices are in [water-scenario-v2.json](water-scenario-v2.json). Gas H2O keeps reference0; binding energy belongs to the negative adsorbed reservoir. The electron source is copied exactly from its owner's usable-api-v1 kernel. Material inputs are copied from the sealed v2 owner; [copied-inputs.json](copied-inputs.json) and [binding_seals.hpp](binding_seals.hpp) bind the actual compiled inputs.

Build with the frozen native toolchain, MFEM4.10 and MPI, using one job. Run electrode-heat-feedback-v2-tests REPORT_DIR on CPU7, nice15, one thread under testing-core7.lock with CUDA_VISIBLE_DEVICES empty. CMake target electrode-heat-feedback-v2 is the reusable interface; the test executable is component qualification only. Parent owns integration, field/current deposition, the single electron inventory and common atomic acceptance.
