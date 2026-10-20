# MITL source convergence, heating and recovered evidence

Date: 18 October 2026

This packet pins Power-Flow commit `8eae133fe48baae40511bc1e89cf33e29a69fb96` on
`research/flash-publication-foundation`. It extends the previous boundary/material
packet at `/workspace/handoffs/mitl-boundary-extension-20261018T1826.tar.gz`.
Local commits separate conditional-time quadrature, source convergence,
recoverable foil markers, ramp heating, host requirements and source assessment.
They are included in `convergence-heating.patch`. The extracted subset contains
41 byte-identical public files plus four named packaging wrappers.

## Delivered calculations

The manufactured planar source now has deterministic integrated observables:
particle and energy budgets, position means/variances/covariance with the
within-patch spread, outgoing flux-Maxwellian velocity moments and optional
conditional emission-time moments. Actual width/current patch areas enter
before integration. The uniform-width hypothesis preserves each current,
state, time and temperature mixture; it changes width positions only.
Exact-state masses remain supplied inputs. Expected outgoing translation is
`2*kB*T` per particle; sampled energy and its difference from that expectation
remain distinct at the packet interface.

The source refinement CLI varies width, current, original snapshot count and
midpoint subdivisions in separate families. It retains finite-reference
comparisons and actual three-solution orders. Reversal, unresolved roundoff or
nondecreasing differences produce warnings and null orders. The first snapshot
family keeps subdivisions per knot interval fixed, so its changes include the
shortening integration step. A separate follow-up represents every original
piecewise-linear history on exactly the same 512-interval integration grid.
It isolates snapshot interpolation more cleanly while retaining finite time
and spatial errors.

| Follow-up family | Actual three numerical levels | Release-number order |
| --- | --- | ---: |
| Width | 32 / 64 / 128 bins | 2.0003 |
| Current direction | 8 / 16 / 32 bins | 1.9933 |
| Original snapshots, common 512-step integration | 16 / 32 / 64 intervals | 1.9893 |

The smaller current-direction correction to width second moment has order
1.8794 and remains less settled. Comparisons of the highest study level with
the finite finest release number are about 0.0134%, 0.00753% and 0.0140% in
these respective families. These are finite numerical differences, not
certified error bounds or uncertainty in the experiment. Five historical
cases were reused without repeating depletion; only six new cases were
integrated. New compressed arrays were reopened and hash-verified. Historical
source-array digests have no saved array bytes and cannot independently be
recomputed; matching reports/code/configuration and regenerated input hashes
qualify that reuse explicitly.

Conditional-time quadrature uses open Gauss-Legendre probability nodes and
actual log-seconds. Centered variance avoids raw-moment cancellation. Log
moments, representation flags and log-based comparisons expose differences
even when positive SI moments round to zero. Quadrature changes and the timing
kernel's numerical qualification remain separate; agreement of rounded
estimates does not establish exact accuracy.

The linear-ramp half-space reference supplies point and finite-volume B, J
and Joule power, incoming Poynting energy, magnetic storage and accumulated
Joule heat. It integrates mean J squared independently of mean J. Its wall
adiabatic rise is finite, unlike the instantaneous-step reference. The
manufactured example closes global energy to roughly 5.43e-16 relative
residual and exports canonical hash-associated profiles. This is analytic
reference generation and interchange verification, not a numerical FLASH
convergence result. Constant properties, stationary material, an indefinitely
continued ramp and adiabatic point temperatures are explicit limits.

## Recovered and newly acquired material evidence

The partial Wilkes recovery retains ten identifiable published Figure 1
absolute-resistivity markers for specimen 2, 0.203 mm Type 304 foil: six
heating and four cooling. An affine calibration handles native scan skew;
open/filled morphology distinguishes source branches from drawn lines.
Three merged clusters remain masked. Figure 2 is percent deviation, so it
does not establish the 0.305 mm specimen's absolute resistivity. Readout
envelopes are display reconstruction limits, not physical error bars. Sparse
coordinates spanning about 298.82–778.16 K do not establish a continuous
curve, a full cooling branch or a runtime property law. The source-pixel QA
overlay and original PDFs remain local and are identified by hash.

The new official Mills et al. 2004 full text supplies evaluated austenitic
mass-specific Cp/enthalpy and composition-dependent density leads. The printed
Cp equals the derivative of the printed solid enthalpy, but the rounded
enthalpy constant does not give exact zero at the nominal reference; that
offset is retained. Illustrative 304 composition, family-table melting
assumptions and liquid fusion contributions remain separate. Mass-specific
Cp is neither Cv nor an EOS. No composition, density, phase treatment,
high-temperature adapter or runtime default is selected here. The earlier
Flynn table remains bounded to its separate 253.15–423.15 K specimen domain.

## Owner integration requirements

`source/docs/boundary_transaction_requirements.md` specifies missing live-host
inventory, nuclei/charge, mass, momentum and energy transactions. It separates
actual weighted launch energy, kinetic rate barriers, independently referenced
adsorption/internal energies, departing mass and support work. Moving-wall
Galilean identities and a two-packet exact-arithmetic example close the stated
accounts. Atomic accept/rollback, rejected-step replay and restart requirements
are handed to the integrating owner. No live transaction or moving-wall
implementation is supplied.

No Quark source, runner, build or production job is included or used. Boltzmann,
collisional-radiative, state-resolved collision/radiative data and kinetic gap
development remain outside this packet. No private implementation is copied.
No remote push, merge, publication or external message occurred.

## Verification and next work

The new scopes passed **194 targeted checks**. The independent extracted subset
passed **193**, with one optional source-PDF integration skipped because that
external file was deliberately unavailable. This packaging check exercises new
import, dependency and data boundaries; the unchanged whole-project suite was
not repeated. The first extraction run passed 192 and caught an unclosed CSV
handle in one recovery test. After a context-manager correction, only that
failed test was rerun and passed; the 193 count is unique coverage across those
two commands. `verification.json` and `evidence/` retain exact commands, exits,
captured or explicitly reported logs, hashes, runtime versions, resolved findings
and limitations. Original targeted runs without persisted process logs are
identified as such; combined captures are not described as separate stderr.

Use `source/README.md` for independently runnable public examples. The advanced
follow-up script in `evidence/` preserves its actual recorded workspace paths
and historical-cache dependencies. Checksums identify source, generated
reports/arrays, delivery files and the archive.

Next upstream work is joint source time/position/thermal-energy statistics.
For symmetric edge heating, Cov(time,x) can vanish even when release varies
strongly across width; Cov(time,x²), conditional source moments and their
matched-width controls are more useful. A separate guarded enthalpy-input
contract can then qualify the new material evidence without equating Cp with
Cv or inventing the actual foil's phase/composition. Actual FLASH execution
still requires installed licensed source and the configured workstation.
