# CPU candidate audit

The candidate removes measured duplicate work and allocation/connectivity overhead while retaining MFEM basis evaluation, matrix products, and each sample's actual inverse Jacobian. It does not change physics, tolerances, source particle counts, dependencies, or checkpoint data structures.

The measured baseline profile has SparseMatrix::AddMult at 7.17% exclusive CPU, drift at 6.08%, native coordinate binning at 3.71%, malloc at 3.07%, and Table::GetRow at 2.18%. These are measured costs, not a claim that any one dominates the full run. GPU driver polling also occupies substantial CPU time.

Changes:

- Reuse the drift split-time vector across particles. Clear and append the same endpoint/crossing values, then run the same sort and deposition order.
- Cache identical Fields::set and Voltage::at times. All interval search, roundoff lookup adjustment, clamping, and interpolation arithmetic remain unchanged on a cache miss.
- Reuse the native field element/reference point when its electric and magnetic callbacks sample the same position. Field values remain freshly gathered.
- Cache checked scalar Q1 connectivity and use stack scratch for values and gradients. The guard requires a two-dimensional scalar H1 quadrilateral Q1 space, one shared finite element, positive in-range connectivity, and identity DOF transformations. Unsupported or updated spaces use GridFunction methods. CalcShape, CalcDShape, dot/matrix products, and actual physical inverse Jacobians are retained.
- Compute full_matrix*trial once, copy it for reaction charge, and use the original action for field energy.
- Retain the current distributed residual norm, reuse it in the caller, and compute its norm once for the unchanged relative statistic and residual guard. The serial solve still computes its own residual.

Validation uses the frozen MFEM 4.10 installation and frozen launcher, logical CPU 0, nice 15, one compiler job. The original foundation tests pass. Added Q1 tests match MFEM bit for bit on nonuniform rectangles, distorted quadrilaterals, corners/interiors, and Q2, triangle, and stale-connectivity fallback. The 12-step serial PIC trajectory, potential, diagnostic, and edge-transfer trace is byte-identical to the baseline (455,537 bytes). Native electric/magnetic values, interval transitions, repeated times, and backward lookups are byte-identical (31,386 bytes).

Five alternating two-million-sample benchmarks give median Q1-gradient times of 0.217253138 s baseline and 0.139595281 s candidate (1.56x for that operation); native E/B and voltage gathering is 0.602059740 s baseline and 0.323792750 s candidate (1.86x for that operation). Checksums match. These are targeted benchmarks on CPU 0 under competing machine load, not end-to-end deployment speed estimates.

The distributed residual reuse branch was compiled but not executed here because the frozen distributed branch selects cuDSS. GPU operations and deployment were not performed. The deployment owner must replay the merged candidate against the production checkpoint and confirm this branch before rollout. cudss_poisson.hpp is unchanged.

The stable candidate is /home/treefor/projects/runtimes/pic-cpu-candidate-20261029/source. Apply cpu-hotpaths.patch from the source root; source-manifest.json binds the changed source files. handoff.json is the machine-readable review entry point.
