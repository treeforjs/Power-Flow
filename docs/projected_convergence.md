# Qualification of projected spectral features

`mhdlab.projected_convergence.qualify_projected_convergence` checks independent
mesh, time, spectral and state-basis refinements of two hypotheses on common
projected features. Each coarse case must change exactly one resolution count
relative to its hypothesis's common finest case. Run identities must be unique;
features, units, gates, excitation, collection and calibration share a declared
comparison contract. The report retains all input evidence.

For each feature and axis, require
`abs(coarse - finest) / abs(finest) < 0.01`. A zero-to-zero change is zero;
nonzero change to a zero finest value has undefined relative change and fails.
There is no arbitrary normalization floor to conceal a changing feature.
Choose meaningful common features before examining contrasts; do not compare
different bins or silently integrate over missing experimental support.

Successive differences are sensitivity indicators. They do not establish a
numerical uncertainty bound. To qualify a proposed contrast, both model owners
must supply justified per-feature uncertainty bounds and their method/provenance.
Their summed bound must be strictly below 10% of the absolute contrast for every
declared claim feature, in addition to passing all four refinement axes. A zero
contrast cannot qualify discrimination. Declared bounds and contracts still
require scientific review; the utility never asserts physical discrimination.

```powershell
python tools/qualify_projected_convergence.py configs/projected_convergence_manufactured.json --output runs/projected_convergence_manufactured/report.json
```

This fixture exercises the checks with authored values and an explicit synthetic
exact-solution bound. It performs no FLASH, transport or population evolution.
Use a new output file for each request. Physical forward-model convergence
evidence is still absent.
