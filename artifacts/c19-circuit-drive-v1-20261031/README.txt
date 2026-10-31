c19-circuit-drive-v1 -- circuit drive for the Hall port and validation-data inventory (2026-10-31)

Closes the construction part of pilot limitation L5 (drive was the prescribed C19 current only) and drafts the data
scope of audit items A17 and A19. Nothing is integrated into Hall or the coordinator (INTEGRATION.txt says how).
Tags used everywhere: measured / measured_published_display / measured_derived / published_estimate /
machine_data_supplied (values supplied directly for this work, 2026-10-31) / model_derived / recommended.

1. Naming
 "C19" is dissertation Figure C.19: the published current of Mykonos shot 17133 (repository data/dissertation/
 shot_17133_current.csv, sha256 2c54d2db...). The driver is the Mykonos linear transformer driver (LTD) at Sandia;
 there is no separate "C19" machine in the documents.

2. Machine description found, and what is missing
 Documents (dissertation, hash 1d7b1f9c..., the source of every repository recovery; calibration/external_sources.json):
   five 3-m cavities in series, 36 bricks each, 40 nF capacitors, 200 kV multi-stage spark-gap switches,
   ferromagnetic cores, coaxial water transmission line; cavity n+1 triggered 6.6 ns after cavity n (line transit);
   "up to 500 kV" cumulative; 1-MA class, 100 ns rise into a matched load (pp. 52-53). Brick = two capacitors charged
   +V/-V in series through the switch (p. 28, MAIZE, same brick type). Measured voltage: V_bank - L0 dI/dt curves in
   Figs. C.19 (17133), C.11 (15210), 4.7 (15778) and A.1-A.17; V_bank measured at the water-line output, L0 = total
   inductance from there through the target (value not published), dI/dt from calibrated B-dots (pp. 59-60).
   Not in the documents: charge voltage of any shot, stored energy, brick/cavity/line inductances, line impedances,
   switch model or jitter, core material/saturation, L0.
 Supplied machine data (machine_data_supplied, 2026-10-31): 5 stages x 36 bricks x 2 x 40 nF (confirmed: 20 nF/brick,
   720 nF/stage, 144 nF erected); +/-80 kV (160 kV/brick; 256 J/brick, 9.216 kJ/stage, 46.08 kJ; 800 kV open circuit);
   LTD inductance 36 nH total for all 5 stages; feed from the dI/dt monitor to the target entrance ~40 nH (36 nH
   feed-post value recalled); target ~25 nH (20-30); ~140 kV/stage under load, ~2/3 of it at the monitor; FWHM
   160 ns into 0.5 ohm; ~220 kV peak at the target entrance (validation only).
 Derived (from the supplied values): L_stage = 7.2 nH; L_brick <= 36 x 7.2 = 259 nH; sqrt(L/C) = sqrt(36 nH/144 nF)
   = 0.500 ohm (equals the 0.5 ohm rating); short-circuit quarter period (pi/2) sqrt(LC) = 113.1 ns.
   Split of the 36 nH (recommended, typical LTD scale, not from documents; the lumped model uses only the total):
   L_brick = 216 nH (6.0 nH per stage) + 1.2 nH per stage cavity stray/line.
 Consistency with documents: 40 nF per capacitor and two per brick agree with the dissertation and with Lamppa et
   al. SAND2021-12691 p. 15 (external public report, the "Lamppa" lead of HANDOFF.md; fetched, sha256 baabb5ff...).
   Charge voltage: no document gives the charge of 17133 or 15210, so +/-80 kV cannot be checked per shot;
   "up to 500 kV" (dissertation) and "100 kV peak output per cavity, 1 MA into 0.5 ohm" (SAND2021-12691) are maximum
   ratings consistent with +/-80 kV being below maximum. The dissertation calls V_bank both "the charge voltage" and
   "measured at the water-line output" (conflict recorded). The 36 nH feed-post value was not found in the repository
   or in SAND2021-12691 (which gives load-region inductances, see 6).

3. Circuit as built (src/c19_circuit.hpp; CPU, header only, a few scalars per step)
   C (144 nF, v_C(0) = 800 kV) -- R_s(t) -- L_ltd (36 nH) -- node A (dI/dt monitor = water-line output)
   node A -- [core branch: magnetizing 1/L_m and loss G_c to return]
   node A -- L_feed -- R_feed -- Hall port (target entrance; feed plane of the Hall mesh)
   R_s(t) = R_inf + R_sw exp(-(t - t_trigger)/tau_sw) (switch resistive phase / stage-timing onset, phenomenological);
   open switch before t_trigger. Implicit midpoint (trapezoidal) step: 5x5 midpoint system in 80-bit long double,
   fixed pivot order. Port relation per step I_mid = a - Y V_mid. Energy E = C v^2/2 + L_ltd i_s^2/2 + L_m i_m^2/2 +
   L_feed i_f^2/2; per-step identity dE + dt(R_s i_s^2 + G_c v_A^2 + R_f i_f^2) + dt i_f V_port = 0 (exact in exact
   arithmetic; ledger with error-free sums and a stated rounding bound).
   Modes: PrescribedCurrentDrive (regression: Hall's C19 parse and interpolation, bit for bit) and LumpedLtdCircuit +
   PortCoupler (coupled). Reduced Hall vacuum load for calibration: src/reduced_hall_load.hpp.
   Xyce: approved as optional; not built (runtimes/xyce-7.10.0-20261031/STATUS.txt: tarball verified, Trilinos 14.4,
   SuiteSparse and BLAS/LAPACK absent from the frozen root, multi-hour build; stopped per the scope decision).

4. Hall vacuum load (computed with the sealed Hall operators; tools/hall_port_load.cpp)
   DC on the Hall mesh of record (45x45 cells of 200 um, 9 mm box, 4 axial foil cells + 1 air cell each end):
   R_dc = 0.0359454 ohm at 300 K (two 200 um SS304 foils; Hall pilot 216.106 V / 6012.07 A), L_dc = 11.045 nH.
   Box study (same operators, 200 um cells): 9 mm 11.045, 12.6 mm 13.799, 27 mm 15.856, 34.6 mm 16.041, 45.4 mm 16.154 nH.
   Dynamic: a 1e13 A/s ramp from the DC state gives a line response (Z ~ 148 ohm for 0.17 ns), overshoot to 317 kV and
   slow ring-down toward L dI/dt = 110 kV (211 kV after 0.4 ns): the port is a shorted 25 mm two-strip line.
   Reduced load: L_H + R_H(T), R_H from the Hall SS304 law (calibration/hall_ss304_law_table.csv, generated from the
   sealed law), adiabatic uniform-J heating; sigma held at 1500 K beyond the law range (reached at 67-84 ns).

5. Calibration (calibration/calibrate.py, fit_results.json, calibration_currents.png, port_voltage.png)
   Fixed: supplied C, V0, L_ltd (and L_feed where stated); Hall L_H. Fitted (model_derived): R_inf, R_sw, tau_sw, G_c
   (+1/L_m in *_magnetizing variants), one plotted-time offset per shot; time-weighted least squares over every
   published vertex of both shots; 16-start multistart. Peak/rise errors = model/measured - 1.
   variant (load beyond the LTD)                  RMS/peak 17133,15210   peak 17133,15210   10-90 rise   V_target peak
   A  feed 40 nH + Hall 9 mm box (11.0 nH)        12.1 %, 10.5 %         -7.9 %, -10.9 %     +36, +37 %   101 kV (-54 %)
   A3 feed 40 nH + Hall 45 mm box (16.2 nH)       12.8 %, 11.2 %         -7.9 %, -10.9 %     +42, +43 %   139 kV (-37 %)
   B  feed 40 nH + 25 nH target                   13.7 %, 12.2 %         -7.8 %, -10.7 %     +55, +56 %   180 kV (-18 %)
   D  feed 36 nH + Hall 9 mm box                  11.5 %, 10.0 %         -7.8 %, -10.8 %     +31, +32 %   106 kV (-52 %)
   C  feed fitted (-> 0.1 nH bound) + 9 mm box     4.4 %,  3.6 %         -2.9 %,  -6.1 %     -4.3, -3.9 % 146 kV (-34 %)
   C3 feed fitted (-> 0.1 nH bound) + 45 mm box    5.5 %,  4.6 %         -4.2 %,  -7.3 %     -2.4, -2.0 % 204 kV (-7 %)
   (C2 with a magnetizing branch: 2.9 %/2.8 % RMS, rise +1.5/+2.0 %, but L_m = 29 nH and the 0.5 ohm FWHM falls to 110 ns.)
   Parameters: A R_inf 0.285 ohm, switch phase and core loss fitted to ~0; C3 R_inf 0.434 ohm, R_sw 146 ohm, tau_sw
   4.5 ns, G_c ~0; offsets 17133 -0.2 (C3) to -15.5 ns (A2), 15210 +0.1 (A2) to +15.3 ns (C3). Rating check (0.5 ohm, circuit only):
   FWHM 166 ns (A, +4 %), 173 ns (C3, +8 %). Formal 1-sigma values in fit_results.json are unreliable wherever a
   parameter sits at a bound or is degenerate; the spread across variants is the meaningful uncertainty:
   R_inf 0.26-0.44 ohm (formal +/-0.007 in C3); onset R_sw 96-146 ohm with tau_sw 4.5-6.2 ns (formal +/-0.7-1 ns; R_sw
   itself unconstrained, +/-100-900 ohm) when the convex early rise is fitted, ~0 otherwise; core loss G_c 0 in the
   preferred sets, up to 0.69 S where it substitutes for missing inductance (B); plotted-time offsets +/-2 to 28 ns
   formal, -0.2 to -15.5 ns (17133) and +0.1 to +15.3 ns (15210) across variants.
   Lumped switch closure: at t = 0+ the erected voltage divides inductively (A: ~100 kV on the port; C3: a 217 kV,
   ~0.5 ns spike while R_s is large). This is a property of the lumped onset, not a physical prediction; start a
   coupled run after the first 5 ns of the circuit clock or replace the onset law once switch data exist.
   Prescribed C19 current into the Hall 9 mm port: L_H dI/dt + R I peaks at 150 kV (-32 %; max dI/dt 1.28e13 A/s);
   16.5 nH would give 220 kV, close to the converged-box 16.15 nH.
   Accuracy targets (peak and rise within 10 %, 220 kV within 15-20 %): met by C3 only; not met with the supplied
   40 nH feed (rise 31-56 % slow) in any variant.

6. Inductance finding (needs an answer before a production coupled run)
   With 144 nF and 800 kV, the measured rise (10-90 % ~71 ns, peak 0.82/0.85 MA at ~103-106 ns after the 10 % point)
   needs about 47-52 nH in total, i.e. 11-16 nH beyond the 36 nH LTD. 36 + 40 + target (11-25 nH) = 87-101 nH gives a
   quarter period of 176-189 ns. SAND2021-12691 p. 15 states that a 20-25 nH load platform draws ~850 kA from Mykonos
   and 35-40 nH hardware draws 650-750 kA; the shots (0.82-0.85 MA) match the first. Either the ~40 nH feed applies to
   other hardware or it includes the target, or the 36 nH LTD value includes part of it.

7. Coupling, ledger and regression (tests/; logs in outputs/c19-circuit-drive-v1-20261031/tests/)
   test_circuit (sets A, C, C3): ledger closure full pulse (dt 100 ps, 10 ps) and 2e6 Hall-size steps (10 fs):
   per-step residual/bound <= 0.12, cumulative closure 2e-13 to 2e-15 of the 46 kJ store; midpoint order 2.00/2.00
   against the analytic RLC; reduced-load coupling <= 3 evaluations per step; open switch exact; bitwise determinism.
   test_hall_port: R1 2,000,133 drive queries (knots, 1e6 random times, 1e6 Hall step midpoints with Hall's time
   recurrence) bit-identical to hall_corrective_v3::CurrentDrive (v3.1 header); R2 2,000 Hall field steps (sealed
   axial operator) driven by the Hall drive vs the circuit module's prescribed mode: every receipt scalar and final
   E, B, E_low, B_low identical, source ledger = sum circuit_work exactly. H1 circuit + Hall field step, vacuum,
   2,000 steps from a 40 ns handover (I ~ 0.33-0.35 MA): 2.0 evaluations per step, max |F|/|I| <= 4.1e-14, join
   (circuit delivered - Hall circuit_work) <= 1.7e-14 J of 0.15-0.20 J delivered, circuit closure 0.9-4.1e-10 J of
   0.6-1.3 J flux (within the summed rounding bound, set by the 4.5e4 J store held in doubles), Hall balance within its
   admitted bound. H2 with an iterate-dependent gap conduction current switched on mid-run (breakdown stand-in), 1,000
   steps: 5.75 evaluations per step (max 8, set by the conduction fixed point), join <= 2.6e-15 J, all within bounds.
   All 11 checks pass for each of the three parameter sets, 33 in total (tests/RESULTS.txt).

8. Validation inventory and A19 plan: validation/INVENTORY.txt (+ .json, 24 items, hashes) and validation/A19-PLAN.txt.
   Recovered here from the dissertation (vector, same tick calibration as the repository current recoveries):
   resistive voltage V_bank - L0 dI/dt for 17133 (338 vertices) and 15210 (352), Normalized APD 2 (17133, 330),
   2x Normalized SVS (15210, 1058): validation/data/, recovery/recover_c19_c11_signals.py.
   15210: resistive voltage collapses between ~52 and 60 ns plotted (breakdown by the dissertation's definition);
   17133: no clean collapse (plateau ~370-383 kV to ~88 ns, then decline to 0 at ~168 ns; V/I ~1 ohm early, so the
   curve carries an inductive residual of the unpublished L0); APD 2 onset 60-64 ns plotted.

9. Data to supply (exact)
   a. The load inductance beyond the 36 nH LTD for the planar-MITL hardware of 17133/15210: feed (monitor to target
      entrance) and target separately, or the measured L0 of the V_bank - L0 dI/dt curves; whether the 36 nH LTD value
      includes the water line and output section; the source of the 36 nH feed-post value.
   b. Charge voltage of shots 17133 and 15210 (confirm +/-80 kV) and the trigger time of the published time origin
      (or the delay between switch trigger and current start); SAND2021-12691 gives ~10 ns trigger jitter.
   c. Water-line output voltage V_bank waveforms (raw) and calibrations, raw current waveforms, B-dot location and
      probe averaging for "Average Current".
   d. Per-brick inductance/resistance or switch resistive-phase data, core magnetizing inductance and saturation
      volt-seconds per cavity, line segment impedances and transit times (would replace R_sw/tau_sw/G_c fits).
   e. The shot, waveform and method behind the ~220 kV target-entrance value.
   f. A17: RGA model/sensitivities/cracking matrix/date and total pressure; the 12 unrecovered VUV optics inputs;
      interferometry (PDV-laser) n_e data; ICCD frames and gates.

10. Layout and commands (build and tests inside arch-base-20261028/run-frozen; outputs never in /tmp)
   src/c19_circuit.hpp, src/reduced_hall_load.hpp      library
   tools/c19_sim.cpp                                    circuit + reduced Hall load (calibration engine)
   tools/hall_port_load.cpp                             Hall DC inductance / ramp probe (sealed axial operators)
   tools/ss304_law_table.cpp                            Hall SS304 law table
   tools/build_tests.sh, tools/run_tests.sh, tools/calibrate_and_test.sh, tools/frozen-py, tools/lent-core.sh
   tests/test_circuit.cpp, tests/test_hall_port.cpp
   calibration/ (calibrate.py, fit_results.json, parameter files A_supplied / C_traces / C3_traces_45mm, figures,
                 external_sources.json, hall_ss304_law_table.csv)
   recovery/ (signal recovery from the dissertation PDF), validation/ (inventory, plan, recovered data)
   Build:  core-lease-20261030/any-core.sh --wait arch-base-20261028/run-frozen sh tools/build_tests.sh
   Tests:  tools/lent-core.sh 3 sh tools/run_tests.sh   (or any-core.sh --wait run-frozen sh tools/run_tests.sh)
   Binaries and logs: /home/treefor/projects/.workspace/outputs/c19-circuit-drive-v1-20261031/{build,tests,...}
   Read-only inputs: Hall v3.1 source/field_contract.hpp; sealed hall-finite-foil-em axial headers; sealed heat-v2
   material-v2 SS304 law; MFEM 4.10 static libraries (link only); Power-Flow data/dissertation (unchanged).
