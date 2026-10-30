GPU particle push for the corrected Hall motion path — 2026-10-30

source/ holds the shared GPU layer (common/: device session by UUID, one ordered stream, pinned staging, NVRTC compilation with
runtime architecture selection and fixed precision options, drained teardown, deterministic scan/radix sort/segment primitives),
the push backend (push/) and tests. run1/ holds SUMMARY.txt, parity and timing summaries and the CPU/GPU receipts (raw array
dumps omitted; workspace paths written as /home/treefor/sim-workspace/powerflow).

Parity with the sealed CPU motion (Hall commit 8dddec0), bit for bit on every fixture: endpoints, paths, per-edge current,
initial and endpoint node charge (double-word net, x87-equivalent gross), the 7 budget values, refusal text and order.
Fixtures: field-port pairs at weights 1e9/1e11/1e12 at the DC-seed and converged midpoints (also swapped into the
prepare_closed_proposal fixed-point loop: same iteration counts, midpoint E/B and candidate field hash); synthetic 1e3-1e6
carriers (1e6 against a statement-by-statement copy of the sealed loop, since the sealed path refuses above its 1e5 cap);
boundary populations 1e3-1e6 (12,593 bisected trajectories at 1e6); a carrier in metal (same refusal). 45 fresh processes and
5 in-process repeats each were identical; all teardowns clean.

Timing (ms, GPU total including uploads, downloads and host sums): 1e4 carriers 1.32 vs 15.0 CPU arithmetic (48.2 sealed);
1e5 7.6 vs 151.6 (533.8 sealed, incl. 331 ms census hash and JSON seal); 1e6 66.7 vs 1557. Peak device memory 1.5 GiB at 1e6.
All GPU runs ran inside run-frozen (its sealed root carries libcuda 615.78.08 identical to the host's) with NVRTC 12.8 on sm_86.

Open: census ordering, duplicate-packet check, census hash, JSON seal and the 80-bit kinetic/work sums stay on the host; no
on-device acceptance of endpoints yet (needed for a resident multi-step loop); boundary bisection runs 49 trajectories per hitting
carrier in one thread (135 ms vs 23 ms at 1e6 with boundaries); sm_89 and sm_120 parity not yet shown.
