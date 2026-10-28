# Water-monolayer fit, collection optics and CR starting point

**initial_h2o_monolayers is the primary adjustable spectral-regression
parameter.** Each candidate binds a common H2O inventory over actual face
areas and recomputes depletion. H/H2/OH fragments must come from conserved
chemistry. The earlier four-independent-pool movie is a separate control.

Start with `replay/docs/target_collection_and_cr.md` for equations, assumptions,
numerical results and replay commands. Main figures:

- `optics/edge_collection_overlap.png`: collection by current geometric region.
- `origin/source_face_collection.png`: collection by original emitting face.
- `optics/collection_acceptance.png`: finite-pupil acceptance and foil shadow.
- `cr/neutral_cr_radiative_control.png`: sourced H2 radiative branching control.

The optics use the user's approximately19-inch target distance,8-inch mirror
spacing and50um by1inch slit height along foil thickness. The few-inch
first-mirror-to-slit separation is straight-line geometry, not the last
optical leg. Final-leg length, focal lengths/curvature, mirror tilts and camera
mapping remain explicit conditional parameters. Three focus-compatible relay
scenarios are included; none is claimed to be the recovered optical prescription.

The nominal first pupil intercepts1.7304e-4 of isotropic point-source photons
before losses and slit clipping. The spatial overlap coefficients are evaluated
at unit photon rate per historical neutral, not a physical CR emissivity.
Geometric edge positions and side-face origins are kept separate. In the
reference57-60ns scenario, collected side-face-origin proxy fractions are
7.3-8.0%; geometric near/far edge profiles overlap substantially. These are
conditional diagnostics, not measured light fractions or fitted monolayers.

The finite-pupil kernel has accepted-volume units m3 and includes1/(4pi),
source volume, the accepted current-direction interval, both mirror pupils,
slit and spectrometer cone exactly once. `mhdlab.detector_emission` maps
supplied CR upper-state densities and explicit throughput to expected photons
per detector row and spectral bin. It supplies no excitation model or ADU gain.
Actual mirror/window/grating response, the final camera mapping and calibrated
gain are still required for absolute iCMOS counts. Throughput and coverage
can be confounded; the inference must retain that uncertainty. The number of
H2O molecules per monolayer is a separate explicit area convention.

Raikou native-validation accepted the1673-state/28466-branch H2 radiative
graph. A15-branch transient control passes population, inventory and photon
energy checks. Its initial excited density is authored. Full physical CR
still needs state-resolved collision/formation/loss data, electron histories,
nascent internal states and transfer/opacity treatment. The external Raikou
repository was not changed; binary and source identities remain separate.

The SS304 EOS path includes a tested bounded P/u table reader and T(rho,u)
inverse, recovered4273 model parameters and a generation request. No numeric
4273/electronic-free-energy table was acquired, no physical SS304 extrapolation
was installed, and the MFEM material model was not changed. ETHOS-informed
multi-fidelity extension is a next generation step after those inputs exist.

Validation:34 focused checks pass across optics geometry, photon projection,
EOS tables and water-source binding. Tagged-face densities reconstruct the
saved maps within4.9e-17 relative L1; tagged photon profiles sum to the total
within5.2e-16. Pupil refinement changes total reference coefficients by1.4-1.5%
and row profiles by1.4-1.6%; this is sensitivity, not demonstrated convergence.
The collection and source-face plots were visually inspected.

## Replay

Use Python>=3.11 with numpy, scipy, matplotlib, pillow, h5py, PyYAML and imageio.
From `replay`, with fresh output directories:

```powershell
python tools/compute_target_collection.py --maps ../inputs/neutral_maps.npz --metadata ../inputs/metadata.json --output-dir runs/replayed_collection
python tools/neutral_cr_probe.py --raikou-root C:/Users/tsmit/Documents/Raikou --executable C:/Users/tsmit/Documents/Raikou/build-v1-install/bin/raikou-native.exe --output-dir runs/replayed_cr_control
```

The source-face profiles are saved for replay/plotting. Recreating their tagged
flight observation needs the original205MB source_packets.npz, left in the
shared project's `runs/neutral_desorption_current_17133_perimeter_v2`, plus its
metadata. No external native runtime is bundled. Sources, cases and input data
are retained with SHA256 identities in `manifest.json`. Results are provisional;
no experimental monolayer, neutral temperature or calibrated detector fit is
claimed.
