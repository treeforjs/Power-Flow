Chemistry v3.2.2 audit preparation -- NOT a production binding
===========================================================
This isolated successor starts from a coherent local snapshot of stable v3.2.
It applies only the qualified v3.1.1 audit corrections and the reviewed swarm
completeness/publication/sealing fixes. It does not adopt active v3.2.1's new
threshold/attachment/recombination estimates or convert unknown channels into
sampled reactions to bypass strict refusal.

Production and physical-model qualification are false. Existing v3.2 model
estimates remain inherited assumptions; this work supplies no new approval for
them or for the observed multi-seed scientific comparison criterion. The new
policy must refuse all reachable unbounded regions by default. Diagnostic test
fixtures explicitly admit unknowns where needed to test mechanics.

reference/v3.2 is an immutable before/after verified snapshot, including old
generated assets and historical documentation/seal. None of those old files is
a successor qualification receipt. data/sources and data/extracted point into
that local snapshot; licensed numerical source data must stay local.
data/tables is generated separately with the corrected builder and strict policy.
No active predecessor directories or jobs were changed.

The original qualified patches and their exact base manifests are in adoption/.
The only additional integration changes are successor path relocation, the
sealed audit policy/CPU tests, GPU coverage regression hunks (while keeping
v3.2's exp-parity checks), and bounded preparation launch/receipt tools.
run_checks.py is deliberately disabled for this preparation version. It retains
the adopted complete-stage driver for later review, but cannot launch the full
swarm matrix. Use run_prep.py, which submits one-core frozen SimQueue jobs:
  python3 run_prep.py regenerate --borrow
  python3 run_prep.py cpu --borrow
  python3 run_prep.py sanitized --borrow
  python3 run_prep.py gpu-build --borrow
  python3 run_prep.py gpu
Each successful prerequisite is checked before dependent work. The GPU command
requires an explicit shared GPU allocation. These are functional checks, not
performance measurements or full million-event/component qualification.

The inherited swarm contract remains pinned to the historical v3.2 blob and
scientific specification. A regenerated blob must have a separately reviewed
contract before any future swarm comparison; no automatic pin refresh occurs.
The sealer requires a complete bound qualification receipt and will refuse
these preparation receipts. Even a future technical seal reports
production_qualified=false and physics_model_validated=false.

No expensive swarm or coupled simulation is started by this preparation.
