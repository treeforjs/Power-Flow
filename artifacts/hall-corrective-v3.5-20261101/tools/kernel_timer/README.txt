Kernel timer (measurement tool, not part of the Hall build or its qualification)

A stand-in libcuda.so.1 that forwards every driver entry point to the real driver library and times each
cuLaunchKernel with a pair of events recorded on the launch stream. At every cuStreamSynchronize/cuCtxSynchronize
the completed pairs are read (cuEventElapsedTime) and accumulated per kernel name (cuFuncGetName); at exit the table
(kernel, launches, total device ms, threads launched) is written to $HALL_KTIMER_OUT as JSON.

Use: build.sh OUTDIR (writes OUTDIR/libcuda.so.1), then run a test binary inside the frozen launcher with
LD_LIBRARY_PATH=OUTDIR set for the binary only (run-frozen env LD_LIBRARY_PATH=OUTDIR HALL_KTIMER_OUT=... BINARY ...).
The Session's dlopen("libcuda.so.1") then resolves to the stand-in. The event pairs add host overhead per launch and do
not change device results; timings from such a run are for the kernel breakdown only, not for host-cost tables.

Scripts:
  measure.sh BUILD OUTDIR        test_v34 bench, lane configuration, 200 um, 1e5 and 1e6 carriers, 2 and 5 steps
                                 (OUTDIR/run-N-S; expects the stand-in in OUTDIR/ktimer)
  measure_modes.sh BUILD         the same for modes v33, lane (50 um) and owner, 200 and 50 um (gpu-kernel-measure/m2/
                                 run-MODE-T-N-S in the v3.4 outputs)
  device_table.py DIR OUT_JSON   timing table (host CPU and wall from the bench tables of record, device time per step
                                 from the kernel-timer runs; steady step = (5-step - 2-step) / 3) and the per-kernel
                                 resource split (upper bounds of FP64, INT32 and DRAM use from ../ptx_count.py counts)
  fp64_vs_bandwidth.py DIR OUT   the first form of the split (lane, 200 um only); superseded by device_table.py
