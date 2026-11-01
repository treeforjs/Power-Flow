SS304 phase-change thermal owner v1.2 -- 2026-10-31

v1.2 continues the extension of the approximate parts of the owner: it closes the energy ledger of the melt (Lorentz
work, surface-force work, dissipation), replaces the lubrication melt model where inertia matters, keeps ejected melt
in the ledgers, and makes the owner affordable next to the Hall field step (sub-cycling). v1.1
(runtimes/ss304-phase-thermal-owner-v1.1-20261031, committed fe235c5, unchanged and read in place by the regression
test) is reproduced bit for bit when the v1.2 options are off. Config::enable_v1_2(c) switches every v1.2 option on
over the v1.1 options. Provenance: DATA-SOURCES.json (v1.2), DATA-SOURCES-v1.1.json, DATA-SOURCES-v1.json. The v1.1
text is README-v1.1.txt; the integration contract is INTEGRATION.txt (v1.2) over INTEGRATION-v1.1.txt and
INTEGRATION-v1.txt.

1. What changed (per item of the v1.2 scope)

1.1 Energy ledger completeness
  The owner's energy is now thermal E plus melt kinetic energy K (one 128-bit word per face column). Exact closure,
  checked every step:
    dE + dK = Joule + loads + deposits + boundary + W_L + W_S
              - evaporation (reference + kinetic + internal) - radiation - releases - ejected,
  with new exact columns
    W_L  Lorentz work on the melt (taken from the field): sum over columns of tau A [h u_mid . f_t
         - (1/2) f_n,in h^2 div(u)] (tangential J x B on the depth-averaged velocity, plus the work of the inward
         normal load on the normal motion of the layer, w(y) = -y div u);
    W_S  work of the surface forces (recoil and capillary pressure gradient, Marangoni shear): tau A [u_mid .
         (3/2 tau_s - h grad p_s) + |tau_s|^2 h / (4 mu)];
    Q    dissipation booked as heat (internal, kinetic -> thermal), the closure of each column's kinetic-energy
         balance K_new - K_old = W_L + W_S + K_in - K_out - K_ejected - Q. In a steady layer Q equals the viscous
         dissipation of the Poiseuille-Couette profile, 3 mu u^2/h + tau_s^2 h/(4 mu) per area (test_v12 work:
         interior columns to 4e-5 (Poiseuille) and 1.4e-3 (Marangoni)); it also carries the dissipation of jumps
         and wall stagnation of the shallow-layer scheme.
  Lubrication model with Config::book_work: per interface the Stokes dissipation G q + tau_m u_s (exactly W_L +
  W_S, all dissipated at once) is booked as heat in the two surface cells (exact half split); W_L = (f_x - dp_n/dx)
  q dx dz tau with p_n the normal-load part of the layer pressure (test_v12 lubwork: f q dx dz dt to 2e-13, dE = W_L
  + W_S exactly). With book_work off the lubrication model is v1.1 (no work columns).
  Exchange term for the field solver: W_L (StepRecord::lorentz, per patch patch_lorentz) is energy the melt received
  from the electromagnetic field; Hall must debit its field-energy ledger by W_L (INTEGRATION.txt 1.3). W_S is not a
  field term (vapor momentum and surface energy). No other energy enters or leaves the owner unbooked.

1.2 Inertia: depth-averaged melt layer (Config::flow_model = 2)
  Layer of a face column: the contiguous fully liquid cells from the surface cell inward inside that face's half of
  the column; mass m, mass-based thickness h, density, viscosity at the layer's mass-mean temperature. Equations on
  each large face, along x (width) and z (the height through the slices k of one body: flow between slices):
    dm/dt + div(m u) = 0
    d(m u)/dt + div(m u u + P) = A [h f_t + 3/2 tau_s - h grad p_s] - 3 mu u A / h
  P = rho g h^2/2 with rho g = -min(f_n, 0) the inward normal Lorentz load (J x B from Hall); p_s = recoil (Knight
  Knudsen layer with back pressure) - sigma kappa (curvature of the displacement field along x and z); tau_s =
  dsigma/dT grad T_s; bed friction 3 mu u/h - tau_s/2 from the Poiseuille-Couette profile (no slip at the solid,
  shear tau_s at the surface), which reproduces the lubrication fluxes h^3 f/(3 mu) and h^2 tau_s/(2 mu) in the
  steady limit. Scheme: Rusanov interface fluxes with hydrostatic reconstruction (bed = displacement - h; a flat
  surface over a varying melt depth stays at rest), semi-implicit friction, explicit sub-steps limited by
  (|u| + a) tau <= cfl x column width per interface (a^2 = g h + sigma h k_max^2/rho, cfl 0.25); a sub-step
  estimate from the last step's largest signal speed, and on a CFL refusal inside the step the step is redone with
  twice the sub-steps (deterministic, up to Config::max_cfl_retries). Walls: mirror states (x edges, z ends) unless
  open with ejection. Transfers: a donor gives the fraction dm/m_layer of every layer cell's atoms and energy (the
  interface words are the sums of the cells' words); the receiver distributes by mass share (remainder to the surface
  cell); dissipation is distributed the same way. Solidification or evaporation of layer liquid takes its share of
  the momentum; new melt joins at rest. The lubrication model stays selectable (flow_model = 1, the v1.1 default).
  Benchmarks (test_v12, inviscid layer, g = 1000 m/s2 from J_z B_x):
    dam break (Stoker 1957, h_L/h_R = 1.4/0.6 mm, t = 40 ms): middle depth 0.953 mm vs 0.954 mm (0.10 % at 400
      columns, 0.17 % at 200), shock position within one column of the Rankine-Hugoniot shock, L1 depth error 8.6 %
      (200) and 5.0 % (400) of the disturbance (first-order scheme, ratio 1.73);
    moving bore (Rankine-Hugoniot speed 1.528 m/s into still water): measured 1.540 m/s over 25 ms (within one
      column of the predicted position), upstream state kept to 0.09 % (depth) and 0.08 % (velocity);
    the same dam break along z (200 slices of one body) = along x, bit for bit;
    flat surface over a melt depth varying 0.2-1.0 mm: |u| < 1.2e-16 m/s after 200 steps;
    low Reynolds number (Re 0.5 and 0.16): steady Poiseuille velocity f h^2/(3 mu) to 5e-7 and Marangoni velocity
      tau h/(2 mu) to 1.4e-7; depth-averaged flux = the v1.1 lubrication flux of the same configuration to 5e-7.
  On the v1.1 full-options foil configuration the layer reaches Re 2.4 (v1.2, 3000 windows of 10 ps) where the
  lubrication model reported 7.1e5: the Stokes flux had no inertia limit; with inertia the layer accelerates as f t/rho
  over these times.

1.3 Melt ejection (Config::eject)
  Open edges (Config::open_x, x edges of each face, default open; Config::open_z, z ends of a body, default wall):
  transmissive outflow (zero-gradient state, outflow only). Departure: liquid of a column whose outward displacement
  exceeds Config::eject_fraction x EM cell thickness (0.25; the geometry refusal stays at 0.5) leaves at the rate
  (eta - limit) rho0 A per sub-step (at most half the admitted donor fraction). Ejected atoms, thermal energy and
  kinetic energy leave the owner through an explicit reservoir: ledger columns eject_N (per element), eject_E
  (thermal + kinetic) and eject_K (kinetic part); per patch in StepRecord::patch_eject_N/E/K and cumulative in
  PatchReport::cumulative_ejected_*; this is the output for the native owner (future droplet species). Tests: edge
  outflow = rho u h dz dt of the steady layer (3.5e-11), per-patch ejecta add up to the ledger; departure keeps the
  outward displacement at 0.029 of the EM cell with the limit at 0.02 (one sub-step of inflow above the limit);
  exact ledgers in both.

1.4 Cost: thermal sub-cycling, fused transfers, device residency
  Owner::accumulate(power, J_B, dt_field) or Owner::accumulate_device(d_power, d_JB, dt_field) (device pointers on
  the engine's session, stream order, no synchronization) add each field step's Joule energy as exact 128-bit words
  (from_double(P dt) per EM cell, the same words a per-step call books) and (J x B) dt; StepInput::use_accumulated
  makes one owner step consume the window (dt = window length, force = window mean of J x B). The window Joule word
  equals the sum of the per-field-step words exactly (test_v12 subcycle). Accuracy bound reported per window:
  StepRecord::sigma_change_max (largest relative change of an EM cell's conductivity over the window; the field used
  the start value for K steps) and window_heating_K (window Joule / (m c_p,min), the largest temperature change that
  can occur inside a window). SS304 foil, K = 4 and 16 against a per-field-step owner: temperatures agree to 5e-5 and
  7e-5 K (bounds 0.25 and 0.79 K), sigma change 1.6e-5 and 5.2e-5.
  GPU engine: the surface physics (boundary exchange, loads, deposits, releases, HKL evaporation with the Knudsen
  layer, radiation) moved from the cell kernel to one thread per surface element (full warps for the double-precision
  transcendental work, which in the cell kernel ran in nearly every warp for one or two edge cells); a two-level slice
  reduction (rows, then slices) instead of one serial thread per slice; fused kernels (one zeroing pass, one
  start-of-step pass per EM cell for force, start conductivity and Joule partition, the flag copy inside the cell step);
  one staging area and one synchronization per step; the start-of-step thermodynamic pass skipped when the scratch
  holds the accepted slot; a light patch-sum kernel in the step (full report only on request); the conductivity of the
  new state returned by the step; the slice sums of the accepted state cached by the owner (no evaluation before each
  step); release inputs uploaded only in quantized mode. All of it keeps every word bit-identical (exact integer sums;
  same expressions). Optional device-time profile per kernel group (GpuOptions::profile).
  Cost on the Hall foil (262,500 cells, RTX 3080 Ti, every v1.2 option, field step 10 fs): owner step 3.0 ms (device
  2.4 ms; v1.1 engine 22.6 ms per step, v1.1 options on the v1.2 engine 3.2 ms); per field step 3.0 ms (K = 1),
  0.75 ms (K = 4), 0.38 ms (K = 8), 0.19 ms (K = 16), 0.093 ms (K = 32); accumulation 1.5 us per field step. The
  <= 2 ms target is met from K = 2; at 10 fs field steps the stability sub-step (7.5 ps) allows one owner sub-step per
  window up to K = 746.

1.5 Sources (DATA-SOURCES.json)
  Phase explosion: X. Xu, Appl. Surf. Sci. 197-198 (2002) 61-66 (read; cached copy and digest recorded): nucleation
  negligible below 0.9 T_c, significant at about 0.9 T_c, measured surface temperature during phase explosion about
  0.9 T_c. The convention is now sourced (v1.1 O6). T_c default unchanged (8000 K, the SESAME 4272 value as stated by
  Mattsson & Desjarlais 2004, the lowest estimate cited). Emissivity 0.37 (Cezairliyan & Miiller 1980): two index
  copies of the abstract give it for 1700-1900 K; used at all temperatures and therefore extrapolated outside that
  range.

2. Selectability and identity
  Config defaults: flow_model 1, book_work off, eject off, diagnostics off: every state word, record, ledger word,
  Merkle root (cells and face-column sections keep their v1.1 tags) and patch report equals v1.1 (test_regression:
  v1-mode evaporation, quantized releases, all v1.1 options for 600 steps, through a surface-cell merge). The layer
  state (momentum, layer mass, K word) is a third Merkle section ("ss304-phase-v1.2/layer-flow"), digested when
  flow_model = 2; checkpoints carry it and refuse to be written with a non-empty accumulation window.

3. Tests and results of record (RESULTS.txt in /home/treefor/projects/.workspace/outputs/ss304-phase-thermal-owner-v1.2-20261031/)
  test_regression  v1.2 (options off) = v1.1 bit for bit (state, ledger, records, roots, reports)
  test_owner       the v1 suite on v1.2 (v1 mode)
  test_law         the v1/v1.1 law checks (law unchanged)
  test_v11         the v1.1 suites on v1.2 (v1.1 options)
  test_v12         dam break, bore, z-flow, well-balanced rest state, low-Re limits, work and dissipation, lubrication
                   work booking, edge and departure ejection, sub-cycling, checkpoint, full v1.2 ledger (3000 windows,
                   every option, K = 4, exact every step, device-style Merkle leaves = host)
  test_gpu         GPU = CPU bit for bit for every v1.1 suite and for v1.2 (inertial layer along x and z, ejection,
                   sub-cycling with host and with device accumulation, dam break with CFL halvings): 19 checks; timing
                   on the Hall foil (Sec. 3 of RESULTS.txt)
  nvrtc_archs      sm_86, sm_89, sm_120 compile
  Sanitizers: CPU programs and the GPU parity program with -fsanitize=address,undefined, no report.
  Labels of record: e3 and e3s (RESULTS.txt).

4. Layout: as v1.1 plus source/gpu_program.h (NVRTC main program shared with tests/nvrtc_archs.cpp), tests/test_v12.cpp,
   run_owner_v1_2.py; tests/test_regression.cpp now links the v1.1 owner.

5. Open and approximate
  A1  Depth-averaged layer: one velocity per column with an assumed Poiseuille-Couette profile (no profile
      evolution, no vertical inertia, no turbulence model); first-order Rusanov scheme (diffusive at jumps: 5 % L1
      error at 400 columns in the dam break); the layer is the fully liquid cells of the face's half column (mushy
      cells do not flow; a fully molten column is split between its faces); no flow on the 0.2 mm edge faces (liquid
      reaching an x edge leaves through ejection or meets a wall).
  A2  Normal Lorentz load: only the inward part acts (a pressing load, g >= 0); an outward load (pull-off,
      Rayleigh-Taylor-like) is not represented. The Lorentz work uses the time-centred velocity and the column
      divergence for the normal part; the field solver does not yet include the back-EMF u x B in its Ohm's law.
  A3  Surface-force work W_S comes from the vapor momentum flux and the surface free energy, which are not owner
      states; it is a ledger column, not a conversion between two owned energies.
  A4  Ejection: the reservoir has mass, composition, thermal and kinetic energy but no droplet size, direction
      distribution or flight; the departure limit is a geometric criterion (fixed Hall mesh), not a breakup model;
      crater depth (inward displacement) still ends in the geometry refusal at 0.5 EM cell.
  A5  Sub-cycling: the field solver sees the conductivity of the window start for K field steps (bound:
      sigma_change_max); the Joule energy is exact but applied uniformly over the window inside the owner (bound:
      window_heating_K); the melt sees the window-mean J x B.
  A6  Phase explosion: refused at 0.9 T_c without the nucleation time lag (about 5 ns in Xu 2002); T_c itself is
      uncertain by 2000 K across the estimates.
  A7  Emissivity: one value, measured range 1700-1900 K, extrapolated to 300-6000 K.
  A8  v1.1 items O1-O3 (liquid data for SRM 1155a, activities, Knudsen-layer mixture) unchanged; cross-architecture
      bit identity shown on sm_86 only (sm_89/sm_120 compile only).
