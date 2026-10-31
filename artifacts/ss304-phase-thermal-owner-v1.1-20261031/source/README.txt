SS304 phase-change thermal owner v1.1 -- 2026-10-31

v1.1 extends the approximate parts of v1 (runtimes/ss304-phase-thermal-owner-v1-20261031, committed, unchanged and
read in place by the regression test). Every extension is an option; with the Config defaults (all options off)
v1.1 reproduces v1 bit for bit (tests/test_regression.cpp links both). Config::enable_v1_1(c) switches every
extension on. Provenance and tags: DATA-SOURCES.json (v1.1) and DATA-SOURCES-v1.json (v1 entries).

1. What changed (per item)
  1. Liquid transport from pulse-heating data (LawOptions::caloric = 1). Liquid enthalpy, resistivity and density
     from Pichler, Simonds, Sowards, Pottlacher, J. Mater. Sci. 55 (2020) 4081 (ohmic wire pulse heating of NIST
     SRM 1155a, a Cr-Ni-Mo austenitic steel; measured 1708-2900 K; 2.9 % / 4.5 % / 3 % at k = 2):
       H - H(298) = -335 + 0.847 T kJ/kg; rho_e = 1.154 + 1.893e-4 T uOhm m; D = 8065 - 0.661 T kg/m3;
     thermal conductivity k = L0 T / rho_e with the Sommerfeld Lorenz number L0 = 2.443e-8 W Ohm/K^2; linear in T
     across the 1672-1727 K mushy zone from the Ho & Chu / Mills solid values; the linear fits are extrapolated
     above 2900 K (stated). The solid law is the v1 law (304 data).
     Validation against the pulse-heating dataset (law and owner, test_law 7a, test_v11 pulse): the owner, heated by
     a uniform Joule pulse like the experiment, reproduces H(T) and rho_e(T) at 1800-2800 K to 1e-15 (the law is
     the data); the v1 law deviates by up to 1.8 % (H), 4.3 % (rho_e) and 3.5 % (density), all within or near the
     data uncertainty. WF conductivity at 1800 K: 29.4 W/(m K) versus Mills 29 (10 %).
  2. Composition and activities. Six inventories Fe, Cr, Ni, Mn, Si, C (Mills Table 1 typical 304: 70.16 / 19.0 /
     9.3 / 1.0 / 0.5 / 0.04 mass %); vapor pressures of Mn, Si (pure liquid) and C (graphite) from NIST-JANAF like Fe,
     Cr, Ni (1 bar at the JANAF boiling points of Fe, Cr, Ni, Mn, Si within 2e-5 decades). Activities (activity = 1):
     regular solution with Omega_FeNi = RT ln 0.66 at 1873 K (gamma_Ni infinite dilution 0.66, Zellars et al. 1960,
     Knudsen measurements); the other pairs 0 (Gilby & St. Pierre 1969: negligible interactions in iron-rich Fe-Cr-Ni;
     Roy & Hultgren 1965: Fe-Mn nearly ideal). Preferential depletion: each surface cell evaporates with its own
     composition, so Mn (and Cr) deplete in the surface layer; melt flow carries the composition of the donor cell.
     Mn carries 91 % of the evaporated atoms at 1700 K, 73 % at 2000 K, 16 % at 3000 K (nominal 304, vacuum).
  3. Gas-phase enthalpy (gas = 1). Energy removed per evaporated atom = eps0_i + E_int,i(T) + kinetic, with E_int
     tabulated so that the implied latent heat equals the JANAF vaporization enthalpy at every table row (relative
     deviation 3e-16; v1: up to 10.6 % at 5000 K). E_int carries the gas electronic excitation and the
     element-versus-alloy caloric difference; it is a new ledger column (evaporation internal energy) and, in the
     quantized route, the internal energy each released macro neutral must carry.
  4. Knudsen layer and radiation. knudsen = 1: Knight (1979) jump conditions (as reproduced in arXiv:2501.18781,
     App. D): the Mach number at which the Knudsen-layer pressure equals the patch back pressure (coupling input
     StepInput::back_pressure_Pa, the native neutral/plasma pressure at the surface; Config::back_pressure_Pa for
     standalone runs) sets the net flux and the recoil pressure; vacuum gives the sonic values (T 0.669, rho 0.308,
     p 0.206, net flux 0.816 of the one-sided flux, recoil 0.550 p_sat); at or above saturation no evaporation.
     radiation = 1: eps sigma (T^4 - T_env^4) on exposed faces, eps = 0.37 (Cezairliyan & Miiller 1980, 304L,
     abstract value), new ledger column.
  5. Melt dynamics (Config::flow). Lubrication flow of the liquid layer along x on the two large faces of every
     slice: q = delta^3/(3 mu) (f_x - dp/dx) + delta^2/(2 mu) dsigma/dT dT_s/dx, with f = J x B from Hall
     (StepInput::em_JB, per EM cell), p = recoil (Knudsen layer) - sigma x curvature of the displaced surface -
     normal Lorentz load x delta/2, sigma(T) and dsigma/dT from Mills Eq. 26/27 (S content configurable), mu(T)
     from Mills Eq. 28; delta = liquid thickness from mass and density (mean of the two columns, at most half the
     foil). Mass and energy move between surface cells as exact antisymmetric transfers (donor composition and
     specific energy), so both ledgers stay exact. Outputs per patch: displacement (mean, max), melt depth, melt
     speed, recoil pressure, largest melt Reynolds number, surface composition. Surface displacement per face column
     = mass gained through the face / (rho0 A) + half the thermal expansion of the column (density law). Typed
     refusals: geometry-departure (|displacement| > Config::geometry_refusal_fraction x EM cell thickness, 0.5 in
     enable_v1_1), melt-flow-step (a sub-step would move more than flow_max_fraction of a donor cell),
     melt-flow-reynolds (optional, Config::max_reynolds), phase-explosion (cell T >= 0.9 T_c; T_c default 8000 K,
     the lowest cited estimate, which puts the threshold at 7200 K, above the 6000 K law domain).
  6. Engineering. Block-Merkle state digest (state_digest.h; native v7 forms layout: 32-record chunks, 1024-record
     leaves, native-v5 node rule) over the cell words (E, six atom words, flags) and the face-column words, computed
     on the device by the GPU engine (1.4 ms for 262,500 cells) and checked against a host recomputation
     (Owner::verify_state_leaves). GpuOptions::session runs the engine on a caller-provided hall_gpu::Session (its
     stream), e.g. Hall's.
  7. Architectures. The NVRTC program compiles for sm_86, sm_89 and sm_120 (tests/nvrtc_archs.cpp, cubin digests in
     the run output). Bit identity is shown on sm_86 only; the test for the RTX PRO 4000 Blackwell station is in
     Sec. 4.

2. Law domain: 298.2-6000 K with the v1.1 law (JANAF gas tables end at 6000 K); 298.2-5000 K with the v1 caloric law.

3. Tests and results of record (RESULTS.txt in /home/treefor/projects/.workspace/outputs/ss304-phase-thermal-owner-v1.1-20261031/)
  test_regression  v1.1 (options off) = v1 bit for bit: evaporation/loads/deposits (600 steps), quantized releases,
                   Stefan sub-steps, through a surface-cell merge
  test_owner       the v1 suite on v1.1 (v1 mode): 49 checks
  test_law         v1 law checks plus pulse-heating comparison, boiling points, activities, internal-energy channel,
                   Knight jump values and limits, Mills viscosity/surface tension, Mn-bearing HKL table
  test_v11         pulse-heating owner validation; Knudsen limits (vacuum, 0.6 p_sat, above p_sat); Poiseuille layer
                   (q = delta^3 f/(3 mu), exact to 5e-16 in the first step) and Marangoni layer (q = delta^2 tau/(2 mu),
                   2e-14); capillary leveling (discrete eigen-rate to 1e-7, lubrication-theory rate within 0.2 %);
                   radiation; thermal-expansion displacement; Merkle; phase-explosion and geometry refusals; ledger
                   closure with every option (3000 steps, melt flow, Mn depletion, Knudsen, radiation, internal
                   energy)
  test_gpu         GPU = CPU bit for bit: v1 mode, every v1.1 option (300 steps), quantized releases with internal
                   energy, capillary flow with 251 sub-steps per step, device Merkle digest = host digest,
                   caller-provided session; full Hall geometry 262,500 cells: 22.6 ms per step GPU, 147 ms CPU
  nvrtc_archs      sm_86, sm_89, sm_120 compile
  Sanitizers: the CPU programs and the GPU parity program with -fsanitize=address,undefined (label c2s).

4. Test for the RTX PRO 4000 Blackwell station (sm_120): build with run_owner_v1_1.py build LABEL, then
     flock <testing-gpu.lock> flock <testing-gpu-UUID.lock> <core lease> /usr/bin/python -I run_owner_v1_1.py gpu \
       LABEL parity test_gpu {OUT} {NVRTC} <station GPU UUID>
     (GPU_UUID in run_owner_v1_1.py names this station's device and must be set for the station), and
     run_owner_v1_1.py run LABEL v11 test_v11 {OUT}. Pass criteria: test_gpu all checks (GPU = CPU on the station),
     and the CPU digests ledger11_state_sha256 (test_v11) and the regression results equal those of RESULTS.txt here;
     together they give GPU(sm_120) = CPU(station) = CPU(here) = GPU(sm_86).

5. Layout: as v1 plus source/state_digest.h, tests/test_v11.cpp, tests/test_regression.cpp, tests/nvrtc_archs.cpp,
   data/janaf/{Mn,Si,C}-*.txt, run_owner_v1_1.py.

6. Open and approximate
  O1  The liquid data are for SRM 1155a (316L-type, 2.2 % Mo, 12.5 % Ni), not 304; above 2900 K the linear fits are
      extrapolated (rho_e 2.29 uOhm m and k 64 W/(m K) at 6000 K, not bounded by data).
  O2  Activities: only the Fe-Ni binary is measured-based; Si and C coefficients not acquired (gamma = 1, an upper
      bound of fluxes that are negligible); partial excess enthalpies are not in the per-atom energy; Mo, P, S not
      tracked; no species diffusion in the melt (the depletion depth below the surface-cell thickness is not resolved).
  O3  Knudsen layer: single-species monatomic jump conditions applied to the mixture with one net fraction; the back
      pressure must come from the coupling; no condensation from the background (returning atoms only as deposits).
  O4  Emissivity: one abstract value (0.37) for solid and liquid; range ambiguity noted in DATA-SOURCES.json.
  O5  Melt flow: lubrication (Stokes) model along x only (no axial z flow between slices, no flow on the 0.2 mm edge
      faces); inertia neglected (the Reynolds number is reported; with molten 100 um layers and Lorentz forces it can
      exceed 1e5, outside the model); viscous heating and the Lorentz work on the flow are not in the energy ledger
      (the Lorentz work is not debited from Hall); fully molten columns are split half and half between their faces;
      merged cells keep reference geometry.
  O6  Phase explosion: the 0.9 T_c onset is the stated convention (not verified from a source); T_c estimates
      8000-10000 K put it above the law domain.
  O7  The GPU step is 3.4x slower than v1 (more inventories and the flow passes); per-step uploads are synchronous;
      super-cycling (accumulated Joule energy) remains the cheap option.
  O8  Cross-architecture bit identity shown on sm_86 only (sm_89/sm_120 compile only).
