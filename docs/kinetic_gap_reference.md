# Executable kinetic transport and state-exchange reference

This is the first executable component set for the gap model. The integration
target is a native solver port delivered through a handoff. Verification runs
in the current workspace; it requires no Quark runner, external plasma executable,
or FLASH installation.

## Implemented components

| Component | Behavior |
| --- | --- |
| `mhdlab/kinetic_orbits.py` | Three-position/three-momentum Boris orbits for neutrals, electrons and ions; relativistic or nonrelativistic motion; explicit gyro and field-work diagnostics |
| `mhdlab/state_kinetics.py` | Immutable species and internal-state records, arbitrary polyatomic quantum labels, common-reference energies, exact state lookup, conservative channel checks and bounded cross-section interpolation |
| `mhdlab/reaction_kinematics.py` | Deterministic 2-to-N reaction momenta from supplied Jacobi energy shares and directions; state-energy changes, recoil and conservation ledgers |
| `mhdlab/kinetic_gap.py` | Prescribed-field streaming coupled to homogeneous-cell stochastic state-changing collisions; outward flux-weighted thermal emission with supplied internal states |

The state catalog can represent molecular neutrals, cations and anions with
source-defined rovibrational assignments. Each reaction names explicit reactant
and product state IDs. There is no fallback from a missing state to a bulk species
or a prescribed rotational temperature. This capability does not supply the still
missing physical water rates or demonstrate convergence of all molecular bases.

The reaction constructor computes available translational energy from incoming
relative kinetic energy and the difference in common-reference state energies.
It conserves nuclei, charge and the declared additive nonrelativistic masses.
Explicit ordered Jacobi modes then determine recoil. A mode-energy fraction is
not the kinetic-energy fraction of an individual outgoing electron or fragment.
Actual differential cross sections, state branching, electron energy sharing and
angular correlations remain supplying-kernel responsibilities.

## Run the reference case

With normal project requirements installed, from the repository root:

```bash
python tools/run_gap_reference.py --plot
python -m pytest tests/test_kinetic_orbits.py tests/test_state_kinetics.py \
  tests/test_reaction_kinematics.py tests/test_kinetic_gap.py -q
```

The case writes `runs/kinetic_gap_reference/history.csv`, `summary.json` and,
with `--plot`, `overview.png`. The summary pins the source-file hashes, catalog,
seed, prescribed fields and parameters. It reports energy residuals and event
counts, without declaring physical qualification.

The example uses a **manufactured molecule M**, manufactured excitation rates,
and an explicitly chosen isotropic COM scattering assumption. It is an executable
coupling check, not a water mechanism, measured shot, or synthetic VUV result.
It includes an electron, a neutral and an ion to expose their different dynamics.

## Numerical and physical scope

- The collision operator is one homogeneous cell with equal macro-particle
  weights W. Each unordered pair's event hazard is `W sigma(g) g / V`, with
  cross sections evaluated at relative COM kinetic energy. Rates are rebuilt
  after every stochastic event. Unknown supported-energy data raise an error;
  zero below a state threshold follows the energy constraint.
- Streaming then collisions is a first-order operator split. The collision
  interval freezes the cell's spatial representation. A spatial collision mesh,
  wall intersections, finite inventory injection and surface energy feedback
  still require implementation. Cell-local reaction placement is not a
  charge-conserving subcell PIC deposition algorithm.
- The coupled operator is nonrelativistic, with a maximum speed cap of 0.05 c
  or a stricter requested value. The standalone orbit kernel supports
  relativistic motion; a relativistic pusher does not make the reaction kernel
  relativistic. A high-voltage production calculation must provide valid
  reaction kinematics over its actual electron-energy range.
- Uniform E/B are prescribed. Gyro motion and species currents are evaluated;
  self-consistent electromagnetic evolution, Hall feedback, circuit coupling,
  charge-conserving current deposition and sheaths remain separate work.
  No extra Hall force is added to Lorentz trajectories.
- Gyro-angle-controlled subcycling resolves rotation up to a configured bound
  and rejects work exceeding the substep limit. Convergence in orbit accuracy,
  split timestep and collision-cell size is still required for a physical case.
- Thermal emission samples an outward Rayleigh normal velocity and Gaussian
  tangential components, giving mean wall-frame kinetic energy `2 kB T`.
  The caller supplies internal states and must debit surface inventory and
  energy. Fragment products follow reaction recoil, not this launch law.
- Radiative-transition records are validated inputs, but photon-event evolution,
  absorption, full CR kinetics and detector synthesis are not integrated here.

Reaction, field-work and total-energy diagnostics retain separate ownership.
The step balance uses kinetic changes and integer state-population changes rather
than subtracting large rounded formation-energy totals. Trial operations do not
mutate input particle states and restore the supplied random-generator state on
failure. Differential callbacks must be free of external side effects. A native
transactional host must extend this acceptance contract to populations, fields,
caches, random state and time across the complete coupled proposal.

## Verification evidence

Checks include ballistic neutrals, electron/ion gyromotion, magnetic energy
preservation, common crossed-field drift, relativistic speed/energy, thermal
emission energy and mass scaling, two-body recoil, three-body breakup/ionization,
moving COM and rotation invariance, state-dependent thresholds, strict data domains,
equal-weight scaling and analytic event survival/pair counting. The manufactured
reaction data are explicitly separated from physical water data.

Next implementation work is the native field/cell interface and conservative
surface/boundary events, followed by supplied water state-to-state data and
electromagnetic current feedback. Preserve the energy-ownership and acceptance
contracts when porting existing transport/source components.
