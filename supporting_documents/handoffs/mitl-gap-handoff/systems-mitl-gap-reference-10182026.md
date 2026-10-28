# MITL kinetic transport and state-exchange handoff

Date: 18 October 2026

Status: independently executable reference components delivered for a native
integration. Quark production source, its integration branch and its runners
have not been modified or used. Integration remains with the Quark owner.

## Delivered behavior

- Three-component charged/neutral particle motion, including relativistic and
  nonrelativistic Boris reference updates and orbit/work diagnostics.
- Species and internal-state catalogs with source-defined polyatomic assignments,
  common chemical energy references and strict conservation/data-domain checks.
- Conservative two-reactant to N-product kinematics with supplied energy sharing
  and scattering directions. Internal-state changes affect available recoil.
- A homogeneous-cell collision operator coupled to prescribed-field streaming;
  equal-weight pair hazards, changing states, product creation and energy ledgers.
- Flux-weighted thermal launches with explicit internal states, unit surface
  normals and wall-frame velocity distributions.
- Rejected reference steps preserve particle state and restore random state.
- Standalone C++20 orbit/recoil candidates, analytic verification executables,
  and a manufactured coupled example with Make and CMake harnesses.

The production requirement covers every retained molecular neutral, cation and
anion, not only H2. The current mechanism inventory contains 12 molecular species.
These are H2O, OH, H2, O2, H2O+, OH+, H2+, H3+, O2+, H3O+, OH- and O2-.
The generic catalog supports their different level-label schemes; physical
state-to-state data and a converged molecular basis are not supplied by the
manufactured example.

## Source and package layout

The reference source is pinned to Power-Flow commit
`05a3d2f3ec2c3d536eb47c1522f9b51c2e940f8e` on local branch
`research/flash-publication-foundation`.
The exact packaged bytes are pinned by `reference-source.sha256`; the commit
identifier is informative and does not replace that manifest.

- `reference/`: self-contained component source, tests, example and run guide.
- `native/`: independently executable C++20 header candidates, tests and example;
  these do not configure the Quark source tree.
- `evidence/`: local verification logs and manufactured reference outputs.
- `electrothermal-port-map.md`: private reviewed source/component map for the
  existing energy, population and boundary implementation.
- `reference-source.sha256`: source/test snapshot identity.
- `native-source.sha256`: native candidate and standalone harness identity.
- `verification.json`: environment versions, commands and observed outcomes.
- `MANIFEST.sha256`: all delivered file hashes, excluding this manifest itself.

The extracted reference package has a minimal package initializer so it does
not import the old reduced-MHD runner. The four component modules and example
are exact copies of the pinned implementation. No private implementation was
copied into those components.

## Local reproduction

From `reference/`, in a separate Python environment:

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python tools/run_gap_reference.py --output-dir ../evidence/recreated --plot
```

These commands run the independent reference package. They are not Quark
runner submissions, native Quark builds or production qualification commands.

From `native/`, the standalone C++ components can be checked with:

```bash
make check BUILD_DIR=/tmp/mitl-gap-reference-build
```

The native example advances prescribed-field orbits with an explicitly scheduled
manufactured excitation. It checks recoil and internal-energy exchange without
claiming a physical event rate or spatial collision calculation.

The reference example uses a manufactured molecule, manufactured coefficients
and an explicitly selected isotropic COM scattering assumption. Source hashes,
input fields, catalog identity, random seed, events and energy residuals are
recorded with the result. It is not a water-plasma or dissertation reconstruction.

## Observed verification

- Power-Flow full project suite: 210 passed, 2 skipped.
- Extracted Python reference suite: 135 passed.
- Native CMake/CTest harness: all 3 targets passed; orbit and reaction executables
  each cover 9 verification groups, plus the coupled example.
- Python manufactured case: 25 excitations and 25 deexcitations; maximum step
  energy-balance residual about 4.7e-15 eV.
- Native manufactured case: one scheduled excitation; maximum step balance
  residual about 3.1e-16 eV.

These residuals describe numerical conservation in the manufactured cases.
They do not establish the accuracy of water fragmentation or emission physics.
The logs and exact environment are included in `evidence/` and `verification.json`.

## Native integration contract

1. Retain one species/state/energy registry shared by source, transport and
   diagnostics. Every molecule uses a source-appropriate rovibrational basis.
2. Preserve one event extent for particle changes, electron redistribution,
   internal energy, recoil, radiation and boundary exchange. Do not count the
   same formation/excitation energy as both stored species energy and heat.
3. Return trial updates and residuals without mutating accepted host state.
   Acceptance and rollback cover particles, populations, fields, random state,
   time integration history and caches together.
   A reaction adapter must change state identity, mass, charge, momenta and
   particle count together, including ionization and dissociation products.
4. Reuse continuity's actual flux/current in electrical work and boundary
   energy accounting. Keep particle, momentum, translation, internal, chemical
   and photon transfer moments distinct.
5. Couple current to the field owner for Hall feedback. Lorentz trajectories
   alone in prescribed fields are not a predictive Hall gap calculation.
6. Add spatial collision cells, material intersections and conservative
   surface inventory/energy exchange before a foil-width calculation.

The private port map identifies existing components and their current validity
domains. Port the equations and behavioral contracts into the native interfaces;
do not substitute the existing bounded nitrogen data for water rates. Preserve
the listed source and licensing provenance for any direct source reuse.

## Explicit remaining work

There is no native Quark integration or self-consistent field/circuit solver in
this packet. The coupled reference operator uses uniform prescribed fields, one
homogeneous collision volume, equal particle weights and nonrelativistic reaction
kinematics with a speed-validity bound. Spatial mesh deposition, walls, radiation
evolution, the full physical water mechanism and measured-spectra comparison
remain unimplemented. The standalone orbit update supports relativistic motion;
that does not extend the reaction kernels to high-energy electron collisions.

No Quark files are claimed by this packet and no production queue work is requested.
