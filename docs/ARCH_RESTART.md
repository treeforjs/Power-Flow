# Arch restart checkpoint

Work stopped on request before restarting the computer. No new
high-resolution native simulation has run. The selected restart mesh is
**2048 x 256 cells per foil**, superseding the earlier 8192/16384 requests.

## First run after the Linux environment is ready

Install/configure MFEM 4.10 with MPI/Hypre/Metis, CMake, Ninja and the Python
requirements. `MFEM_DIR` may point to its installed CMake configuration.
The new native snapshot-export changes have not yet been compiled in Linux.
Do not mistake the saved configuration for a completed simulation.

```bash
mkdir -p runs
python3 tools/benchmark_mfem_edge_source.py --config configs/mfem_edge_current_arch_2048x256.json --output-dir runs/arch_2048_preflight --mode preflight
/usr/bin/time -v -o runs/arch_2048_walltime.txt python3 tools/run_mfem_edge_source.py --config configs/mfem_edge_current_arch_2048x256.json --output-dir runs/arch_2048_heating --ranks 8
python3 tools/render_mfem_foil_temperature.py --source-run runs/arch_2048_heating --output-dir runs/arch_2048_foil_temperature
```

The runner also records configure/build/native-command wall times separately.
The driver distributes every electromagnetic and thermal linear system over
the requested MPI ranks. Assembly, field evaluation and output retain rank-zero
ownership. Each rank owns a disjoint contiguous matrix-row block; a completed
`parallel_execution.json` checks equal system/solve participation on every rank.
The numerical equations, source trace, solver tolerances and current checks
are unchanged. Single-rank results reproduced the earlier driver exactly in
two differential controls; eight-rank field, face and diagnostic comparisons
passed a declared 1e-8 column-scaled limit with separate energy checks.
Use the frozen environment launcher for builds and runs on this installation.
The case retains the shot 17133 trace, 60 ns cap, four material-face exports,
0.25 ns maximum step and 66 selected actual field snapshots. Check free RAM
and disk before the full run: field CSV output alone may occupy several GB.
No startup automation or new Windows/WSL configuration was installed.

The existing independent H/H2/H2O/OH-pool GIF remains a visualization control.
The physical spectral fit instead uses adjustable initial H2O monolayers as
its primary parameter, with fragments generated through conserved chemistry.

## Saved temperature views and unfinished work

- Existing density overview/edge GIFs remain in the previously pushed artifacts.
- New saved-data foil views are in `artifacts/cathode_temperature_saved_20261028`.
  The interior GIF contains ONLY the two recorded 40/60 ns refined fields;
  the separate surface GIF contains actual recorded boundary frames.
- The gas-temperature renderer is implemented; its 10 focused tests passed.
  It shows supported, centered species-resolved translational kinetic T,
  not bulk-flow energy, equilibrium temperature or molecular internal T.
- The unfinished producer is preserved under
  `review_only/temperature_continuation_20261028`. It is not production-ready.
  Fix corner histogram ownership using original unshifted points; advance
  x/y only for invariant-current observation; add finite-overflow guards;
  then test drift invariance, 3V moments, tiny weights and statistical masks.
- Preserve equal physical axes and both overview and edge close-ups for
  density, foil temperature and gas kinetic temperature.

The shown upper foil corresponds to the dissertation cathode geometry.
Its top/back and bottom/gap-facing faces are surfaces of that SAME electrode.
The lower/anode return is outside the near-foil gas movie. See
`electrode_polarity_mapping.md` for evidence and electrical-sign limitations.

## Questions and requested physics for the next session

1. **Thermal expansion of the solid:** add this to the Arch work plan. Couple
   the temperature history to transient solid mechanics with thermal strain,
   deformation and updated emitting-face position/normal/area/velocity. MFEM
   supplies elasticity capabilities, but the present heating run has fixed
   geometry. Resolve SS304 expansion/elastic/plastic properties, supports,
   inertial response and Lorentz forces. Melting/vaporization require phase/EOS
   and flow treatment beyond a solid thermal-strain model. This is requested
   work, not an implemented or enabled feature.
2. **Expansion into the stated background pressure `5e-5`:** confirm its units
   and gas temperature before converting to number density. Determine whether
   residual-gas collisions, plume self-collisions, backscatter, sticking or
   re-adsorption affect the desorption/advection maps as material leaves the
   surface. Compare collision times/mean free paths with the actual plume
   scales and 0-60 ns window before retaining collisionless flight.
3. **Dissertation atomic/molecular weight spectra:** locate and recover the
   relevant mass/RGA spectra, their diagnostic response, pressure conditions,
   background subtraction and fragmentation pattern. Assess whether they
   constrain residual-gas composition or source chemistry; do not equate mass
   peaks directly with neutral branch fractions without that response.
4. Keep initial water monolayer count adjustable through the complete source,
   transport, CR and detector forward model. Retain optical-throughput and
   monolayer-definition uncertainty; do not claim a unique fitted coverage
   while those scales and physical collision rates remain unresolved.

MFEM references: https://mfem.org/examples/ and
https://mfem.org/electromagnetics/. No new simulation was launched after the
instruction to defer until Arch and then stop.
