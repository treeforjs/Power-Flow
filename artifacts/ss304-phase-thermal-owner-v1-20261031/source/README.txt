SS304 phase-change thermal owner v1 -- 2026-10-31

Purpose: one owner of foil metal mass and condensed-phase energy for the Power-Flow coupled PIC water-breakdown
simulation, valid from 298.2 K through melting into strong vaporization (law domain 298.2-5000 K), with exact mass
and energy ledgers, a deterministic CPU reference and NVRTC GPU kernels that reproduce it bit for bit. Closes the
thermal part of pf-pilot-1 LIMITATIONS L4 and audit A16 once integrated (INTEGRATION.txt). Not integrated here.

1. Material law (source/phase_law.h; provenance and tags of every value: DATA-SOURCES.json)
  Caloric (J/kg, zero at 300 K solid, fixed reference volume rho_ref = 8000 kg/m3 as the sealed heat-v2 law):
    298.2-1672 K  solid   Cp = 472 + 0.136 T - 2.82e6/T^2 (Mills 2004 Eq. 10, <5%)                    recommended
    1672-1727 K   mushy   L_f = 160 J/(K kg) x 1727 K = 276,320 J/kg released linearly in T (Mills Eq. 13; producer
                          melting range of TP 304, 1399-1454 C); Cp blended to 800                 recommended / model
    1727-5000 K   liquid  Cp = 800 J/(kg K) (Mills Eq. 12, +-50, measured basis to 1900 K; constant above, supported
                          by the JANAF constant liquid Cp of Fe/Cr/Ni, Kopp-rule value 796)       recommended / model
  Transport: Ho & Chu 1977 Table 2 (k, rho_e; 5% solid, 15% liquid; >1500 K extrapolated by the source) to 1800 K;
    above 1800 K k = 28.99 + 0.01 (T - 1800) W/(m K) (Mills Table 3(b) slope 1800-1900 K, linear continuation) and
    rho_e = L T / k with the constant Lorenz ratio anchored at the 1800 K row (2.4432e-8 W Ohm/K^2). The electrical
    conductivity per Hall EM cell (Owner::sigma_em) is the volume average of 1/rho_e over its thermal sub-cells.
  Density (diagnostic only; the owner keeps the fixed reference volume): Mills Table 3, solid 7921 - 0.5 (T - 298.2),
    liquid 6908 - 0.86 (T - 1723) kg/m3.
  Bit identity: h, apparent Cp, k, rho_e, sigma, liquid fraction and 8000 h are bit-identical to the sealed
    ss304v2::Model evaluated_solid_1500 on 298.2-1500 K (479,962 points) and to nominal_1800 on 298.2-1800 K
    (600,013 points), 0 differences. The inverse T(h) (Newton / quadratic / linear per phase, then nearest-double
    polish) is within 1 ulp of the sealed bisection and never further from h.
  Vapor: pure-liquid vapor pressures of Fe, Cr, Ni from NIST-JANAF gas and liquid tables (Gibbs energy at the table
    rows, ln p linear in 1/T between rows; liquid continued with its tabulated constant Cp beyond the table end:
    Fe > 4000 K, Cr > 3500 K, Ni > 3600 K). Raoult activity a_i = gamma_i x_i (gamma_i = 1; Mills: Fe-Cr, Fe-Ni nearly
    ideal) with the surface cell's own atom fractions (composition Fe/Cr/Ni 71.7/19.0/9.3 mass %, Mills Table 1 with
    Mn, Si, C lumped into Fe). Cross-checks: Fe and Ni liquid within 1.2% and 2.4% of Alcock et al. 1984 over their fit
    ranges; above the boiling points equal to -logKf(liquid) within 0.0005 decades; Cr consistent with Alcock's solid
    equation plus the JANAF fusion Gibbs energy.
  Evaporation: Hertz-Knudsen-Langmuir j_i = (1 - beta_r) alpha_e a_i p_i / sqrt(2 pi m_i k T), alpha_e = 1,
    recondensation beta_r = 0.18 (Anisimov 1968; set 0 when the gas owner returns atoms explicitly). Energy removed
    per atom = eps0_i + kinetic, eps0_i constant (Fe 4.254, Cr 4.032, Ni 4.607 eV; anchored to the JANAF vaporization
    enthalpy at 3000 K), kinetic 2 kT (expected mode) or the actual packet energy (quantized mode). The implied latent
    heat is within 3% of JANAF over 2000-3500 K, within 6% over 1700-4000 K, 10.6% (Cr) at 5000 K.
    Nominal 304 at 3000 / 4000 / 5000 K: 33 / 893 / 5860 kg/m2/s, evaporative cooling 2.0e8 / 5.1e9 / 3.0e10 W/m2,
    recession 5.7e-3 / 0.18 / 1.4 m/s (test_law hkl table).
  Domain limit 5000 K: about half the published critical temperatures of steel/304L (9500-10000 K).

2. Owner (source/phase_owner.hpp; engines source/engine.hpp, source/gpu_engine.{hpp,cpp}; kernels
   source/phase_kernels.h, arithmetic source/det_math.h)
  Geometry: per slice (body, axial cell k) a 2-D cross-section, nx = nx_em x rx columns across the 5 mm width and
    graded rows through the 0.2 mm thickness (EM row edges included); fixed reference volume; no axial conduction
    (as Hall v3). Recommended rows: graded_rows(0.2 mm, 0.25 um, 1.3, 10 um) = 42 rows.
  State per cell: energy (128-bit fixed point, 2^-96 J) and Fe/Cr/Ni atoms (128-bit, 2^-32 atom); T = h^-1(E/m).
  Step (explicit, sub-steps of equal length below a state-independent stability bound that includes the
    evaporation stiffness at 5000 K): harmonic-mean face conduction, each face energy rounded once and applied with
    opposite signs to the two cells; Hall Joule energy per EM cell split exactly over its sub-cells (sigma V or V);
    surface loads, metal deposits, quantized releases or expected evaporation on surface elements; every step total
    divided over the sub-steps by exact integer division; prescribed-temperature faces for benchmarks.
  Surface recession: surface cells below 25% of their reference mass are merged into the inward neighbour at the
    end of the step (exact transfer); exposure is a pure function of the void flags (rebuilt after restart).
  Ledgers: E = E0 + Joule + loads + deposits + boundary - evaporation(reference + kinetic) - releases, and
    N_i = N_i,0 + deposited - evaporated - released, as exact 128-bit integer identities, checked inside every
    prepare (std::logic_error on violation) and by Owner::ledger_closes().
  Protocol: quote(dt) -> prepare(StepInput) -> preflight -> commit (Trial; accepted state untouched until commit;
    a commit invalidates every other trial); producer watermarks for exactly-once loads and deposits; save/load.
  Refusals (PhaseRefusal kinds): temperature-above-domain (5000 K), temperature-below-domain (298.2 K),
    surface-depletion-within-step, column-depleted, release-exceeds-inventory, release-binding, load-identity,
    off-surface, nonfinite-input, invalid-input, invalid-step, substep-budget, fixed-point-range. None at 1500 K.
  Outputs: per patch (slice, face, fraction range) area, mean/max surface T, liquid fraction, melt depth (contiguous
    liquid layer along the inward face normal: through the thickness for the 5 mm faces, along the width for the
    0.2 mm edge faces), HKL flux and expected atoms per element, cumulative evaporated/deposited atoms, recession;
    sigma per EM cell; quotes for the quantized release route n = SR(D/w0).
  GPU: the same functions as NVRTC kernels (one thread per cell / EM cell / element / patch / slice; no atomics;
    integer sums), compiled by the sealed hall-gpu-resident-v1 session with --fmad=false --prec-div=true
    --prec-sqrt=true --ftz=false; transcendental functions are ports of fdlibm exp/log in basic IEEE operations
    (det_math.h), so host and device agree bit for bit.

3. Tests and results of record (run_owner_v1.py; RESULTS.txt in
   /home/treefor/projects/.workspace/outputs/ss304-phase-thermal-owner-v1-20261031/)
  test_law    bit identity with the sealed law; inverse; continuity; Mills Table 3 spot values; JANAF transcription;
              vapor cross-checks (Alcock 1984, JANAF logKf); latent-heat model; deterministic exp/log (<= 1 ulp from
              glibc); fixed-point round trip and exact division; water-desorption cap analysis
  test_owner  uniform heating (exact enthalpy solution through melting to 5000 K, typed refusal above), constant-flux
              conduction (erfc solution, 2nd order), 2-D mode decay (exact discrete eigenvalue to 4e-13 K), two-phase
              Stefan problem (front within 0.06%), refusals and rollback, checkpoint restart (bit-identical
              continuation), quantized releases with binding checks, super-cycling, Hall geometry, Bennett 2023
              published surface-heating case, surface-cell merges, ledger closure through melting and strong
              vaporization (4500 K)
  test_gpu    CPU versus GPU bit identity of every accepted state, step record, patch report, quote and sigma:
              evaporation + loads + deposits (300 steps to 4978 K), quantized releases (150 steps), Stefan with prescribed temperature and 24
              sub-steps per step, a surface-cell merge (27,835 steps), the full Hall foil geometry (262,500 cells:
              GPU 6.7 ms per step, CPU 93 ms per step)
  Sanitizers: the same three programs built with -fsanitize=address,undefined (label b2s).

4. Layout
  source/        det_math.h, phase_law.h, phase_kernels.h (host + NVRTC); phase_law_host.hpp, janaf_data.h
                 (generated by tools/gen_janaf.py from data/janaf), engine.hpp (CPU reference), phase_owner.hpp,
                 gpu_engine.{hpp,cpp}
  tests/         test_law.cpp, test_owner.cpp, test_gpu.cpp, check.hpp
  data/          janaf/*.txt (NIST-JANAF tables as fetched), metals-vapor-pressure.txt (Alcock transcription)
  DATA-SOURCES.json, INTEGRATION.txt, run_owner_v1.py

5. Open and bounded items
  O1  Liquid transport above 1900 K is a continuation (k linear, rho_e constant Lorenz ratio): no data bound; Joule
      heating above melting depends on it.
  O2  Activity coefficients are 1 (no data acquired); Mn (about 1 mass %, more volatile than Fe), Si, C are not
      tracked; surface-cell composition changes by selective evaporation without species diffusion in the melt.
  O3  Constant per-atom vapor reference energy: gas electronic excitation (JANAF gas Cp above 5/2 R) is not carried;
      latent-heat deviation bounded above (O3 numbers in Sec. 1).
  O4  Solid alloy evaporation uses the pure-liquid standard state (overestimate up to x1.8 Fe, x5.1 Cr, x2.3 Ni at
      1000-1672 K; fluxes there are below 1e-3 kg/m2/s).
  O5  No hydrodynamics, pressure or expansion geometry: the melt layer stays in place (no j x B flow, ejection or
      surface tension), density is diagnostic, superheated liquid above the 1-bar boiling points is kept condensed,
      phase explosion is not modelled.
  O6  Knudsen layer as a fixed recondensation fraction; no plasma/vapor back pressure; condensation only through
      explicit deposits. Thermal radiation is not included (sigma T^4 = 1.45e7 W/m2 at 4000 K, below the HKL
      evaporative cooling of 5.1e9 W/m2 there; it can be added as a load).
  O7  Per-step device-side state digest (Merkle, as Hall v3.1) not built: state_sha256 downloads the state.
  O8  Electron emission and water desorption laws of the source owner end at 3000 K (desorption cap justified in
      INTEGRATION.txt 2.2; emission open).
  O9  GPU engine opens its own session (no shared Hall session yet); per-step transfers are synchronous.
  O10 Merged cells keep the reference geometry (surface position error up to one cell); side walls toward voids
      are adiabatic.
  O11 No axial conduction between slices (as Hall v3).
  O12 The published comparison is Bennett et al. 2023 Eq. (1) (constant sigma, no conduction; the owner reproduces it
      to 5e-9 of the exact cell integral and gives solidus at 114.5 T versus the published 104 T with constant
      c_v). No resolved pulsed-heating melt/vaporization measurement for 304 was acquired.
  O13 Bit identity shown on one device (RTX 3080 Ti, sm_86) only.
