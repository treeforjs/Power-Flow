# Current-driven foil, return conductor and vacuum source

The source in `native/mfem_edge_source` replaces the earlier metal-only box
heating case with a finite-width electromagnetic field solve. Current crowding
and its Joule heating follow from that solve; no edge-temperature multiplier or
imposed edge profile is present. This is a new MFEM magnetoquasistatic driver,
not the unchanged stock `joule.cpp` kernel and not a full MHD calculation.

The current production record is
`runs/mfem_edge_current_source_v3/source_run.json`. Its `top_surface.json` and
`bottom_surface.json` feed the neutral transport calculation. All run products
are ignored by Git and remain in the shared workspace.

## Geometry and formulation

The current direction is invariant z, with a 10 mm length used once to convert
two-dimensional per-length powers and energies to W and J. The upper foil spans
x = +/-2.5 mm and y = +2.5 to +2.7 mm. Its identical return conductor spans
y = -2.7 to -2.5 mm. The gap between facing surfaces is 5 mm. The baseline
exterior box extends to +/-10 mm in both coordinates, with A_z = 0 on its outer
boundary. This boundary truncates the field and its inductance; a larger-domain
case checks its influence.

With B = (dA/dy, -dA/dx, 0) and E_z = U_i - dA/dt in conductor i, the weak
backward-Euler system is

```
H A_new = M_sigma A_old / dt + F_foil U_foil + F_return U_return
H = K_(1/mu0) + M_sigma / dt
J_i = sigma_i (U_i - (A_new - A_old)/dt)
integral(foil J) = I(t), integral(return J) = -I(t).
```

Both voltage gradients U_i are unknowns. Three MFEM linear solves and a 2 by 2
Schur complement enforce the two integrated current constraints each step.
Conductivity is exactly zero in the vacuum; no small numerical vacuum
conductivity is introduced. H1 bilinear finite elements represent A_z on the
conductor-aligned rectangular mesh. J is evaluated from Ohm's law, not from an
elementwise second derivative of A_z.

Each conductor has its own H1 thermal problem with insulated boundaries:

```
(M_Cv / dt + K_k) deltaT_new = M_Cv deltaT_old / dt + load(J^2/sigma).
```

The squared current is integrated by quadrature before assembling the source.
It is not the square of cell-average current. Vacuum has no thermal degrees of
freedom. The top output is the upper foil back face at y = +2.7 mm, normal +y;
the bottom output is its gap-facing surface at y = +2.5 mm, normal -y. Both are
actual material-side boundary traces averaged over each width patch and the
invariant z length. Their CSV values are temperature increments; the surface
adapter adds the 300 K initial temperature once.

## Drive and conditional initial state

The input is the recovered dissertation Fig. C.19 shot 17133 trace,
`data/dissertation/shot_17133_current.csv`, with the adjacent recovery JSON as
provenance. It is digitized published data, not raw oscilloscope data. Its first
supported point is 0.749788055492 ns at 6012.07136047 A. The calculation starts
there and stops at 60 ns, inside the known support. Piecewise-linear
interpolation includes the original trace knots and exact 40/60 ns snapshots;
intervals are divided to remain at or below 0.5 ns. The code does not pad an
unknown trigger interval, retime the trace, or extrapolate outside support.

At the first supported point, the model assumes uniform DC current in each
conductor, solves the corresponding static magnetic field, and sets T = 300 K.
This is an explicit conditional initial state. The magnetic and thermal
prehistory is unknown. The published current clock is assumed comparable with
the published ALEGRA clock for comparison, but that registration is unverified.
The runner currently requires exactly the 300 K baseline to keep its native
absolute snapshots and surface metadata consistent.

## Numerical records and energy checks

The build uses installed MFEM 4.10, release commit
`d964264cdb9a13e94a201b6c236c7721e0c8765f`. The driver uses one MPI rank and
single-rank Hypre AMG-preconditioned conjugate gradients. It is a serial spatial
solve, despite linking MPI-enabled MFEM. The library installation is reused
without modifying MFEM's source or its build.

The run record binds configuration, current input, binary, native source and
output hashes. The provenance CSV hash is compared both as raw bytes and with
CRLF changed to LF. That newline-only normalization is explicitly recorded and
does not change the input numeric values.

For every time interval, the electromagnetic ledger checks

```
dt Lz sum_i(U_i I_i) = delta magnetic_energy
                    + dt total_Joule_power
                    + 0.5 Lz deltaA^T K deltaA.
```

The final term is backward-Euler magnetic damping, retained separately from
physical Joule heating. The thermal ledger checks that the increase of heat
content equals integrated Joule power. Initial magnetic energy is retained; it
is not charged again as drive work. The whole-loop voltage is
`(U_foil - U_return) Lz`, whereas each foil voltage uses its own gradient.

At 60 ns the baseline target is 500734.775996 A. The maximum reconstructed
current error over the run is 2.48e-8 A. The maximum absolute per-step EM and
thermal ledger errors are 6.14e-12 J and 8.00e-14 J. The final reported linear
residual is `||H x - b||_2 / ||b||_2 = 8.36e-14`, maximized over the three
electromagnetic solves that step. Solver convergence is checked for the static
initial state and for thermal solves as well, but their residuals are not
included in this column. The initial diagnostic row uses zero placeholders for
step-only quantities and must not be interpreted as a measured zero residual.

`foil_snapshots.csv` contains upper-foil cell averages at exactly 40 and 60 ns.
Its depth coordinate increases into the upper foil from the gap-facing surface.
The JSON sidecar declares SI units, sampling, geometry, current provenance and
the CSV hash for `tools/compare_mfem_alegra.py`.

## Scope and physical limits

The values sigma = 3e6 S/m, Cv = 2.8e6 J/(m^3 K) and k = 16 W/(m K) are
constant illustrative coefficients. They are not a qualified SS304 property
model or parameters fitted to the ALEGRA curves. The calculation omits material
motion, phase change, latent heat, conductivity feedback, radiation and Hall
terms. It does not simulate a plasma or infer a measured wall temperature.

At 60 ns the baseline outer-edge-band mean is about 2625 K and its central-band
mean about 570 K. The corresponding mean Joule-power-density ratio is 7.81.
Those means hide larger corner temperatures: the largest gap-face patch is
about 9811 K, the back-face patch about 8449 K, and the largest foil cell average
about 9089 K. Even at 40 ns the largest cell average is about 3354 K. These
temperatures exceed the solid constant-property model's useful regime. Values
are preserved rather than clipped, but high-temperature neutral release and
ALEGRA differences cannot be interpreted as qualified predictions.

The rectangular sharp corners create mesh-sensitive field and heating maxima.
Band diagnostics use the outer 10% of width on each side and the central 20%,
integrated through the full thickness. Native band integrals classify
quadrature points by these fixed physical cuts, so their effective band area
has a small quadrature discretization error. Broad band means, corner maxima,
and a lineout at a specified physical depth are distinct measurements. Mesh and
outer-boundary sensitivities are required before claiming edge convergence.

## Completed sensitivity checks

The baseline, doubled mesh and expanded exterior runs all completed on the same
bounded current schedule, with their records in
`runs/mfem_edge_current_source_v3`, `runs/mfem_edge_current_mesh_refined` and
`runs/mfem_edge_current_domain_extended`. The reproducible comparison is
`runs/mfem_edge_current_sensitivity/sensitivity.json`, with conservatively
projected lineouts in its adjacent CSV. No coefficients or input timing were
adjusted between cases.

The table below uses exact width overlap of exported cell averages for band
means, rather than the native quadrature-point cuts. It is a piecewise-constant
reconstruction, so its baseline band values differ slightly from the native
values above. Both edge bands and the center band retain fixed physical widths.

| At 60 ns | Baseline 64 x 8 | Mesh 128 x 16 | Exterior +/-20 mm |
|---|---:|---:|---:|
| Edge/center mean q | 7.886 | 7.579 | 7.870 |
| Edge-band mean T (K) | 2650.5 | 2559.7 | 2658.8 |
| Center-band mean T (K) | 569.8 | 572.5 | 572.9 |
| Maximum cell-average T (K) | 9089 | 10562 | 9318 |
| Maximum gap-face patch T (K) | 9811 | 10988 | 10084 |

For the finer mesh, gap/back-face temperature-increment lineouts change by
2.46%/4.26% RMS after averaging onto the baseline width bins. At a fixed depth
of 12.5 micrometers, reconstructed temperature increment, J and q lineouts
change by 3.45%, 0.934% and 2.53% RMS. For the larger exterior the corresponding
surface changes are 4.30%/3.83%, and the interior changes are 3.87%, 2.41% and
3.27%. RMS percentages are normalized by the baseline RMS of the same quantity;
for temperature this uses the increment above 300 K. Interior lineouts are
linearly reconstructed between cell-average depth centers, with no depth
extrapolation.

The larger exterior changes final magnetic energy from 677.8 to 720.9 J and
the inferred foil voltage gradient from 2.792 to 2.956 MV/m. Its effect on the
temperature asymmetry between gap and back faces is measurable even while the
broad band q ratio changes little. The refined mesh increases sharp-corner
maxima while modestly changing broad bands. These checks support the presence
of field-derived edge crowding, but neither corner nor boundary convergence is
established. The largest current reconstruction error among these runs is
5.78e-8 A, with per-step EM/thermal ledger errors below 1.23e-11/1.64e-12 J.

## Reproduction

From this repository in Windows PowerShell, with Ubuntu-24.04 and MFEM installed:

```powershell
python tools/run_mfem_edge_source.py --wsl-distro Ubuntu-24.04 --output-dir runs/new_edge_source
python tools/run_mfem_edge_source.py --wsl-distro Ubuntu-24.04 --config configs/mfem_edge_current_mesh_refined.json --output-dir runs/new_edge_refined
python tools/run_mfem_edge_source.py --wsl-distro Ubuntu-24.04 --config configs/mfem_edge_current_domain_extended.json --output-dir runs/new_edge_domain
python tools/summarize_mfem_edge_sensitivity.py --baseline runs/new_edge_source --case runs/new_edge_refined --case runs/new_edge_domain --output-dir runs/new_edge_sensitivity
```

The runner requires a new output directory and preserves unsuccessful runs.
The refined configuration doubles all mesh counts, including 128 width and 16
thickness cells in each conductor. The domain case changes only the exterior
half-extent from 10 to 20 mm and its margin counts from 16 to 38, keeping the
exterior cell spacing approximately constant; this is a boundary/mesh
sensitivity rather than a pure isolated boundary change. Both retain the same
trace, timestep cap, material assumptions and 60 ns endpoint.
