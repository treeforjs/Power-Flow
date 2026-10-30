# Matrix structure and trace evidence

The single-threaded CPU reconstruction used the same frozen MFEM mesh, Q1 diffusion integrator, permittivity, coordinate fields and foil-aligned essential-node test as the sealed source. Prescribed boundary values and particles alter the RHS, not this stiffness/elimination pattern. It did not factor the matrix or load/cuDSS-execute a solve.

| Structural quantity | CPU result |
| --- | ---: |
| Rows | 130,321 |
| Stored nonzero upper entries | 633,713 |
| Symmetric nonzero entries including both triangles | 1,137,105 |
| Essential/constrained nodes | 3,618 |
| Isolated diagonal rows | 3,618 |
| Connected components | 3,619 |
| Largest connected block | 126,703 rows |
| Maximum off-diagonal graph degree | 8 |
| Original matrix-index bandwidth | 362 |

The raw CSR array sizes are 521,288 bytes for offsets, 2,534,852 bytes for columns and 5,069,704 bytes for FP64 values. These match each of the two original matrix-upload sizes in the existing baseline trace. This is a strong dimensional/structural cross-check; exact native byte/hash parity remains for the parent's diagnostic export. File hashes are in the manifest. No matrix permutation or numerical operation is changed by this audit.

The source leaves ordering at its documented default, customized METIS nested dissection, and ND_NLEVELS at documented512. This is source/documentation evidence, not a queried live handle. Existing traces do not reveal the actual permutation, elimination tree or factor fill. The new optional audit queries effective config, LU_NNZ and memory estimates after factorization, and can export/validate the row permutation. A query status/byte count must be checked before interpreting its value; a failed query does not mean zero fill. Exact tree-level scheduling is not exposed by this header audit.

Baseline forward/backward sparse kernels use grids of3704/68574 blocks. The two TRSV variants use17 blocks with128 threads each on an80-SM GPU. The former have many launch blocks but can still be dependency- or memory-limited; the trace does not distinguish those causes. The latter cannot fill all80 SMs with17 simultaneous blocks, but together account for only about8% of total kernel duration. No achieved occupancy, bandwidth, FP64 efficiency or critical-path timing is established. Reordering/fill comparisons must connect new config/factor data with actual solve traces and unprofiled receipts.

Constrained singleton rows are2.78% of all rows. Removing them could simplify the matrix, but needs RHS/solution packing and reconstruction under the original outer residual check. It is not the first implementation choice while an ordering knob can be tested with much less code. Likewise, the unit-bias matrix is a separate initial construction; the420 later solves arise from the two dependent half steps of210 advances, so bundling those RHSs is not a valid steady-state shortcut.

The main source diff adds only a small pre-analysis configuration parser, optional diagnostics and nested setup timers. The `Solve` method remains byte-identical to the sealed source. No graph implementation, lower precision, reduced refinement, relaxed guard or solver skip is present.
