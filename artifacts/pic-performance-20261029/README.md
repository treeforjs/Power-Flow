# PIC performance and paused-run recovery snapshot — 2026-10-29

This captures research source, guarded benchmark tools, recovered thermal and native magnetic inputs, the latest durable accepted checkpoint, and the corrected progress animation. It is an archive of the workstation experiment, with paths recorded in the original receipts; it is not a qualified new production build.

The full simulation is manually paused near 26.48 ns. Do not resume it automatically. Its supported history starts at 0.749788 ns and the archived checkpoint is at 26.47868472199812 ns. The final GIF frame is a reconstruction of that checkpoint with no time evolution. Restore/verify the exclusive CPU partition before any future full-run resume.

`source/sealed-gpu` is the exact source of the prior qualified executable; the other source variants are isolated experiments. CPU cache changes pass short checkpoint regression and serial/MPI tests, but graph capture, alternative orderings and gather offload are not promoted. CUDA Graph first-RHS bitwise qualification failed; nondefault ND128/256 failed the unchanged residual guard; AMD was slower. Experimental status can change in later receipts.

`inputs/native-fields-domain40mm.tar.gz` contains the launch-manifest-bound native field input files. The thermal cache, binary thermal history, and conditional midpoint pulse are separately retained. Frozen third-party runtimes/libraries are identified by the original receipts; they are not copied into this Git snapshot. Rebuild with those pinned versions and rerun checkpoint regression before using a new executable.

`paused-checkpoint/restart.meta` preserves every source reservoir and references the unchanged checkpoint binary. `history.csv` is trimmed to the checkpoint clock for restart; `history-complete-through-pause.csv` preserves the later in-memory progress samples as historical evidence. The checkpoint's serialized options/diagnostics/particles use the frozen PIC2 ABI.

The render uses the native unscaled B, fixed absolute |B| contours every 2 T, direction-only normalized arrows, a white background, and electron density in cm^-3 in the |z|<0.5 mm numerical slab. It does not calibrate the optical slit. Source-level model limits remain: conditional voltage, mirrored anode temperature, coarse prescribed B, height-averaged 2D self-space-charge with 3D particle motion, no ion/gas/collision feedback or phase changes, and unresolved sheath/grid/material convergence.

The CHICAGO literature review distinguishes published immediate-water-ionization closures from finite-rate hydrogen chemistry; neither establishes channel-resolved H2O fragmentation for this run. Published sources and a staged benchmark plan are in `notes/chicago-pic-review`.
