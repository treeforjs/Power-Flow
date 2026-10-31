Heat/source owner v6 workspace -- 2026-10-30

Contract: API-v6.txt (owner v6, Hall binding, quantized neutral release, step protocol, typed
refusals). Arrival policy v1.1: runtimes/pic-arrival-owner-v1.1-20261030/API-v1.1.txt.
Hall side: runtimes/hall-corrective-v3-20261030/README.txt.

Build and run (receipted; every compile and run inside arch-base-20261028/run-frozen):
  core-lease-20261030/both-cores.sh (or any-core.sh --wait) /usr/bin/python -I run_source_v6.py build LABEL [ubsan]
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v6.py run LABEL [UBSAN_LABEL]
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v6.py growth LABEL T:STEPS[,T:STEPS...] W0
Sealed inputs (source v5, source kernel v3, heat-v2 headers, arrival v1 headers, native candidate003
library) are hashed before and after every build; they are only read.

Results of record (/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/source-v6/)
  s4/run-212341   9/9 normal and 9/9 UBSan (0 runtime errors, same outcomes):
    binding           v6 binds the declared heat-only Hall slice; another slice or an undeclared
                      Hall slice is refused
    temperature       after Hall heats the slice (patch mean 300.166 K), the desorption quote and
                      the emission exposure use the Hall temperature bit for bit
    exactly-once      emission layer and desorption loads enter Hall once: thermal delta = Joule +
                      loads within the heat roundoff bound, the source reports no Joule, the
                      single-owner audit passes, replaying the loads is refused
    double-joule      a sealed v5 private engine driven by a Joule coefficient is refused by the
                      audit (typed double-joule); v6 passes; an external Joule load is refused by
                      Hall (typed external-joule)
    typed-refusals    1e200 stock and ~5e21 release: refused typed by Hall before any engine call
                      (load-exceeds-material-range), Hall unchanged; 1e20 release: typed
                      net-cooling-exceeds-capacity; the sealed v5 engine on the 5e21 case: untyped
                      engine error (the 1e200 case aborts a v5 process and is not run)
    checkpoint        v6 checkpoint restores inventory and layers against the same Hall clock and
                      is refused after the Hall clock advanced
    protocol          a source trial preflights before a Hall publication and is stale after it
    quantized-release 6 Hall steps, w0 = 1e6, about 1e3 leaves per step: every leaf exactly w0,
                      exact source debit = arrival credit = packet sum, stock + released = initial
                      stock exactly, a non-w0 leaf is refused
    unbiasedness      2e5 keyed decisions at x = 1e-3, 0.3, 1.7, 12.5: |z| <= 1.8
  growth (test_source_v6 TRACE OUT growth T STEPS W0): one Hall-slice oxide-water patch (lower-foil
                  inner face, 5 mm x 6.25 mm, coverage 0.33, 1.34e14 molecules), source v6 quantized
                  mode with arrival v1.1, the one thermal owner a heat-v2 slice Coupler at T receiving
                  the desorption loads each step, dt = 10 fs; the v1 floor rule evaluated on the same
                  quotes for comparison.
    s3/growth-205937 (w0 = 1e6)
      T       steps  mean D/step  carriers  E[carriers]  z      v1 floor rule  ms/step
      300 K   1e5    5.24         0         0.52         -0.72  100000         29.4
      400 K   1e5    7.86e3       830       785.8        +1.58  100000         30.6
      600 K   1e4    1.32e7       132305    132301.5     +0.08  132305         29.9
    s4/growth-212659 (w0 = 1e5, 400 K, 1e5 steps): 8086 carriers (E 7858.2, z = +2.68; a 40-seed
      sweep of the same decision gives mean z 0.07, var 0.71, see API-v1.1.txt).
    Every run: released molecules = carriers x w0 exactly, stock + released = initial stock exactly.
    The s3 binary differs from s4 only in Hall v3 sources outside the growth path (birth-scope flag).
  Default w0_neutral = 1e6 and its basis: API-v1.1.txt section 3 of the arrival owner v1.1.
