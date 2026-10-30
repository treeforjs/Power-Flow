Offline occupancy diagnostic and conservative EEDF/H2O/ion instrumentation plan

Start with DIAGNOSTICS-PLAN.txt. diagnostic-contract.json gives the proposed units, metadata/event fields and explicit missing inputs.

snapshot_occupancy.py reads the already validated frozen PIC2 checkpoint schema through the preserved v1 decoder. It computes weighted instantaneous energy/angular OCCUPANCY only. It does not estimate crossing events, time averages, fragmentation rates, neutral density or ions. It rejects non-electron charge/mass and nonfinite input. It preserves energy overflow and zero-speed angular weight.

Run on CPU6 at nice15 with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1. The tool acquires outputs/pic-profile/testing-core.lock internally. It writes only to a new requested output directory.

taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES= python3 /home/treefor/sim-workspace/powerflow/outputs/pic-profile/eedf-diagnostics-plan/snapshot_occupancy.py /home/treefor/sim-workspace/powerflow/outputs/pic-profile/cpu-candidate-cpu6-run1 --out /home/treefor/sim-workspace/powerflow/outputs/pic-profile/eedf-diagnostics-plan/another-occupancy --grid 9

snapshot-audit.json binds the three demo inputs and records count/energy/angular normalization. snapshots/{seed,baseline-final,candidate-final}/ contains occupancy.json, energy_occupancy.csv and angular_occupancy.csv. The two final inputs represent the same physical endpoint and must not be pooled as extra samples/time. fixture-tests.json records18 passing finite/unit/slab/normalization/overflow/angular/empty-state tests. The original decoder source hash is recorded with every snapshot.

Current model remains collisionless electron transport with a height-averaged2D self-field, imposed axial/native fields and finite target sinks. No physical/grid/timestep convergence or reaction prediction is claimed. References and source-evidence files distinguish checked sources from proposed diagnostics. No production or shared candidate source was changed, and no GPU work or process signals were performed.
