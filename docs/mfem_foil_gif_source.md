# Native MFEM foil source for the provisional transport animation

`native/mfem_foil_source/foil_driver.cpp` constructs a rectangular conductor and
uses the **unchanged** external MFEM 4.10 `MagneticDiffusionEOperator` from
`miniapps/electromagnetics/joule_solver.cpp`. It links the installed MFEM library;
the source tree and installation are not modified or rebuilt.

This is a manufactured/provisional SI calculation, not a measured MITL shot or
qualified SS304 prediction. The [configuration](../configs/mfem_foil_gif_source.json)
declares a 5-mm width, 200-µm thickness, 10-mm current length, constant conductivity
3e6 S/m, volumetric heat capacity 2.8e6 J/(m³ K), thermal conductivity 16 W/(m K),
and mu0. A prescribed terminal-potential difference ramps linearly to 150 V in
100 ns. These values are explicit illustrative assumptions rather than fitted
material or drive data.

The 3D mesh has 32×4×2 hexahedra with coordinates x across width, y through
thickness, and z along current. x and y are centered about zero; z runs 0–10 mm.
The two broad material faces are y=±100 µm. Voltage is essential at the z ends;
the induced tangential electric-field boundary is homogeneous. All thermal
boundaries impose zero normal heat flux. The stock operator's linear thermal
field starts at deltaT=0; the declared initial 300 K is added exactly once by the
existing SI surface adapter. Fixed material coefficients make this a consistent
temperature-rise formulation. No phase change, motion, Hall induction,
temperature-dependent property feedback, exterior magnetic-domain qualification,
or measured-current drive is claimed.

```powershell
python tools/run_mfem_foil_source.py --wsl-distro Ubuntu-24.04 --output-dir runs/new_mfem_foil_source
```

The output directory must be new. Native build/runtime stdout and stderr,
commands and exits, configuration, source and executable hashes are retained.
The current source release reference is official v4.10, commit
`d964264cdb9a13e94a201b6c236c7721e0c8765f`; the installed binary identity is
independent and no source/build attestation is inferred.

At every snapshot, mapped MFEM boundary-face quadrature evaluates the
**one-sided material L2 trace** using `Elem1` and its face-to-element integration
point. Area-weighted averages are accumulated into aligned x/z patches and
reduced over MPI. Every patch's summed quadrature area must match dx*dz to
relative 1e-10. This is boundary trace evaluation, not relabeling volume averages.

Files `top.csv` and `bottom.csv` contain
`time_s,z_index,x_index,temperature_value`, with temperature values in kelvin
increments. `top_surface.json` and `bottom_surface.json` contain the corresponding
hash-bound MFEM surface configurations, including exact native snapshot times,
partitions, face coordinates/normals and the explicit initial-temperature
convention. They load directly through
`mhdlab.mfem_surface.load_mfem_surface_samples` and preserve `source_kind=mfem`.
That origin label records an actual MFEM execution; it does not turn its
manufactured parameters into physical qualification.

`diagnostics.csv` retains projected Joule power, independently assembled direct
electric power, thermal and magnetic energies, current at both ends, maximum
linear iterations/residuals, raw temperature extrema and accumulated projected
heat. The runtime rejects nonfinite diagnostics, missing snapshots and nonpositive
absolute temperature. It preserves numerical differences and does not clip fields.
Closing thermal increase against projected source does not establish the complete
electromagnetic boundary-energy balance or mesh/time convergence.

The accepted 32×4×2/40-step/two-rank run ended at 2 µs with absolute temperature
300–771.514459 K. End-current mismatch was at most 1.84e-8 A, maximum projected
versus direct Joule-power difference was 0.003465 relative, and thermal/projected-
source energy residual was below1.5e-13 J. Its explicit status remains provisional.
Surface release and vacuum transport are separate downstream calculations; no
desorption, chemistry or fragment production occurs in this native source driver.
