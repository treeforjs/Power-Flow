Emission integration qualification v2

This separate package links immutable usable-api-v1. It does not change that
API, run the GPU, resume production, or admit a native physical foil model.
The sealed80 source checks and4 microcases remain unchanged. All compilation
and numerical runs require the exclusive testing-core7.lock, CPU7, nice15,
one job and arch-base-20261028/run-frozen.

The32-stage pilot uses one common10fs clock, one electron inventory, a32-cell
1D scalar Poisson solve from resident charge and an elliptic vector-potential
solve from resident current. External boundary drivers and the accepted state
determine field stages. The caloric model evolves1490K solid electrodes using
sampled complete wall heat and an explicitly authored Joule scenario. The
vendored SS304 v1 header/scenario are byte-identical to their material-owner
source at packaging; fixed reference density8000kg/m3 and Cp-as-Cv are explicit
conditional assumptions. This is neither native heat conduction nor a Maxwell
energy proof. It uses start-of-interval frozen fields; a converged common
physical midpoint corrector, work-conjugate field/current coupling and native
one-sided actual surface fields still require parent qualification.

An accepted-state trial is pure. Field/thermal/source reservoirs, unique IDs,
packet phase-space, full LayerState, ledger and field arrays are serialized
fieldwise. An active-layer checkpoint at epoch16 replays byte-exactly through
epoch32. A deliberately mismatched mesh epoch rejects and leaves its actual
rollback base byte-identical. A count reservoir emits a weighted packet when
it reaches one physical electron; its pending count is checkpointed. Fixed
caller-owned quantiles exercise the interface, not Monte Carlo convergence.
Only the chosen sampled wall heat is debited, never the mean diagnostic too.

Nine matching/substep variants retain common dt10fs. Nine further variants
evolve200stages from10MV/m with25T tangential B and zero normal B. They resolve
extracting-side magnetic recollection. Nominal magnetic applicability remains
0.2 throughout. The earlier0.3 manufactured microcase is not physical default.

Returning PIC particles retain their originating LayerState after handoff.
At inward crossing of its calculated handoff distance, inward_layer_adapter.hpp
adds |Uimage| to kinetic energy, restores negative image-reference energy,
preserves transverse proper momentum, and marks the particle Active. No charge
change, displacement current, field work or wall heat is added by that local
reference transfer. advance_layer then resolves the remaining same interval.
Failing to perform this inverse of the outward tail restoration missed magnetic
recollection in build-and-run-v3.log. The failed log is preserved. The v4 run
includes the fix. This is an adapter behavior outside immutable v1 and must be
implemented in native transport before any relevant integration is admitted.
Physical trajectories and quantum bridges still own q dx/current and q E dx
once. A zero-duration bridge is integrated over accepted dt, not divided by0.
There is still one inventory; a packet cannot be both ordinary and Active.

The local field removes the physical single-electron discrete Poisson self
response only when that electron was in the solved accepted inventory. It
does not remove collective space charge or all macro-packet charge. This is
a fixture convention, not a validation of the native conductor-image gather.

Nominal magnetic scan:16legacy finite-offset probe locations,3temperatures
(300,773.15,1490K),3normal-energy quantiles, fixed transverse/azimuth quantiles.
All144phase points are retained:98admitted,46rejected by the0.2 criterion.
Failures occur on extracting upper faces, with signed F between0.213175 and
97.724679MV/m. These counts are not physical probabilities or a missing-current
bound; the probe fields are not admitted one-sided surface fields. Do not
discard or resample failures. The guard is a conditional energy-scale estimate,
not a rigorous20% current error guarantee.

canonical_barrier_diagnostic.cpp derives a purely tangential-B alternative
from the minimal-coupling Hamiltonian, retaining the canonical-momentum cross
term: Ueff/e = phi -(F+B*v_z)*x - c/x + e*B^2*x^2/(2m).
It resolves the first finite forbidden barrier and the subsequent allowed
pocket, with a stated parabolic continuation above the local saddle. It
recovers sealed B=0 SN action and tests quadrature refinement. Three of18
chosen phase cases have no local allowed pocket. A local pocket does not
imply terminal escape; subsequent magnetic return belongs to transport.
This diagnostic is not a finite-T source law: it lacks a magnetized solid
supply, a joint energy/angle sampler and arbitrary oblique-B scattering.
It cannot replace failures of the existing guard. No E*=sqrt(E^2-c^2 B^2)
substitution is made; the general magnetic-dominated range is outside that
electric-dominated estimate's real domain.

Primary derivation: Lebedynskyi et al., Problems of Atomic Science and
Technology98(4),62-66(2015), https://arxiv.org/abs/1702.02713 Eqs4,20,28.
The implementation expands the full squared canonical momentum rather than
using a phase-independent scalar barrier shift. The stable unmagnetized law
and its other primary references remain documented in usable-api-v1/MODEL.txt.

Run reproduce.sh inside the required frozen/locked/affined environment. The
second argument is the original probes.csv. Evidence is component-level;
only the parent can admit native midpoint/Maxwell/heat/mechanics integration.
The optional third argument selects the immutable API source root. The local
sealed library is reused; a source-only published checkout builds its v1
library in a build directory before linking the qualification package.
