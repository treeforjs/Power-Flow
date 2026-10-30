The optimized CUDA force candidate passes the fixed-seed1ps regression under the established policy. The independently reviewed streaming result covers all210 accepted steps in all three mode pairs, in addition to the three complete164-check final saved-state comparisons, shadow per-kick checks and clean runtime/production-hold receipts. No further comparison run was needed.

The streaming parser is SHA256 `364070ea2f3c2b960a40e6fa79d6b748c2f598542a1cbbd9787615bade0da322`, matching the parent review and source-review publication. Its reference comparator, ABI schema and policy lock hashes match the original established sources. The34 manufactured/corruption cases have the expected outcomes; source inspection confirms that malformed/partial/mutated data cannot become a passing completion result. Tests and the real comparison were not rerun by this reviewer.

The parser bounds metadata, step/particle/source dimensions and bytes before allocating, retains only one accepted record per mode, rejects nonfinite/invalid physical values and duplicate IDs, and checks the exact frozen particle/diagnostic ABI. It validates contiguous step/time/dt sequences from the pinned seed, ordered IDs/counts/q/m/weight, all four source reservoirs, native metadata, force rows, clean shutdown, completed endpoints and exact EOF. It rechecks all required final saved-state input hashes and observes file identity/stat changes. The success receipt contains full streamed hashes for all three1,457,665,600-byte traces. This review verified46 non-trace manifest hashes and current trace sizes without redundantly rereading the full streams.

Policy comparison is consistent with the established position/global-velocity/potential budgets,1024-epsilon ledger scales including the seed, and matching charge-scale treatment of edge transfers. Public `diagnostics()` already contains particle kinetic energy; it is not added a second time. The source compares all three pairs, all arrays and diagnostics on every record, and reports failed numeric comparisons separately from invalid/incomplete input. There are422,836,254 checked scalar values and zero failures across steps437412–437621.

For optimized versus off/native, the largest observed differences are:

| Quantity | Maximum absolute difference |
| --- | ---: |
| Position |6.94e-18m |
| Proper velocity |6.20e-7m/s |
| Nodal potential |3.90e-9V |
| Horizontal charge transfer |3.82e-34C |
| Vertical charge transfer |5.87e-34C |

All differences satisfy the original budgets. The largest fraction of any comparison limit across all pairs is0.46365, for a Gauss diagnostic in shadow versus optimized. Particle/source identity and all recorded step/dt clocks match exactly. This is numerical regression agreement, not a claim of bitwise identical evolving states. The earlier fixed-midpoint device output was bitwise equal; the evolving field solves can vary within the unchanged regression limits.

The observed process timing in this single instrumented1ps experiment is13.079084s native/off,14.282559s shadow, and9.225693s optimized. Thus optimized used **29.46% less wall time**, or an observed **1.418x speed factor**, with passing scoped parity. Shadow was slower than off. These process times include initial setup, ordinary simulation work and approximately1.46GB of full-state trace/file observation per mode; whole-force sums retain initial setup/export costs. No repetitions, reversed run order, uninstrumented full-run timing or general production speed claim are established. The result is a measured improvement for this exact bounded trial, not an extrapolation to60ns.

The residual limitation remains explicit: per-step RHS norms and original guard limits were not exported. The parser compares residual diagnostics using the original guard limits recorded in the passed final saved-state reports, as the established comparator does. It does not independently reconstruct or prove each step's RHS guard. The original runtime residual guards were unchanged and executed; final-state guards were separately reconstructed by the established comparator. This limitation cannot justify skipping solves or accepting a future warm-start solver on a biased RHS threshold. Future solver work still needs charge/self-field accuracy observables.

Other limits remain: the late1ps window does not cross a real native-frame boundary; CPU boundary/epoch tests and earlier live shared-context tests do not imply general trajectory or hardware refresh coverage. Unrecorded midpoint E/B, raw source samples, next_id/crossing-event logs and independent electric-work ledgers are not inferred. No physical/grid/time convergence, deployment approval or production resume follows from this regression result. The optimized mode remains an isolated candidate, with accepted v2 and all original scientific guards unchanged.

The parent's profile summary is independently cross-checked in `profile-crosscheck.json`: mean advance time46.419990ms off ->28.115523ms optimized, a39.43% reduction inside the advance timer; mean whole force20.454194ms ->2.885525ms, including initial setup/export cost. The advance timer excludes the later full-state observation, so that reduction is distinct from the29.46% instrumented process-wall reduction. All modes retain210 advances,422 cuDSS solves and2 factor setups. This is a force-offload result, not a new winning cuDSS configuration.

The main chat has the parser review, complete comparison receipt, metadata/hash verification and this final scoped disposition. Production remains the original stopped/frozen process. No GPU execution, solver change, comparison rerun, production signal, source edit or tolerance change was performed by this reviewer.
