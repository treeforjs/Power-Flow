# Independent native-v3 thermodynamic review

Reviewed on 2026-10-29 by the material/thermodynamics agent. This is a source review; the new v2 adapter and tests were not executed by this agent, following the central CPU7 test-owner instruction.

## Implemented adapter boundary

- `mhdlab/solid_foil_caloric_v2.py` selects only the sealed thermal-v2 q1 `nominal_1800` law, validates source, transitive headers, library, and frozen runtime hashes, and uses the owner's C ABI.
- `native/mfem_thermoelastic_v3/ss304_heat_material_v2.hpp` exposes `electrode_heat::conditional_ss304_v2(...)` and `ss304_v2_source_seal_sha256()`. It calls the sealed owner's h, capacity, transport, and phase functions. Its stable enthalpy increment integrates the exact sealed capacity piecewise.
- `tests/test_solid_foil_caloric_v2.py` covers the cold negative reference energy, inverse, independent capacity integral, cold transport, phase support, density/mass convention, provenance receipt, and rejection paths.
- Exact source seal: `1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d`.
- The admitted thermal interval is 298.2–1800 K. Energy remains zero at 300 K, so energy below 300 K is negative. This does not relax the mechanical 773.15 K cap or the smooth-segment restrictions.

## Energy sign and work accounting

Let `W` denote the summed midpoint stress work over actual mechanical substeps and Gauss points. The mechanics use `Delta K = -W` for free or stationary homogeneous constrained boundaries. Let `Delta Uel` denote the change of the conditional elastic internal energy, with `Uel = Phi - T*dPhi/dT`.

The exchange packet defines the caloric debit as `debit = Delta Uel - W`, and sends `qnew = -debit` to native heat. Native heat advances `Delta H = sum(q)` within its stated nodal/global solve bound. Therefore

`Delta H + Delta Uel + Delta K = sum(q - qnew) + (Delta K + W) + native residual`.

The sign is correct and stress work is counted once. Substituting recoverable strain energy `Phi` for `Uel` would change this thermodynamic model and is not equivalent. The code retains the conditional Cp-as-Cv caloric approximation and correctly avoids claiming a physical compressible EOS.

The local heat law belongs to native mass-lumped nodal enthalpy. Averaging nodal h to mechanical cells and inverting creates a spatially homogenized mechanical temperature; this does not establish a cellwise native conduction balance. The current receipts disclose this distinction.

## Acceptance-resolution recommendation

The reviewed corrector initially used `2e-30 + 5e-10*sum(abs(q)+abs(qnew))` as one scalar binding tolerance. This scales with mechanical exchange only. Near zero exchange, it can be smaller than the energy resolution of the native nodal temperatures and caloric inverse even when heat conduction is nonzero.

For each native cell i, use a reported local bound based on its actual reference volume and admitted temperature/capacity:

`r_i = 512*eps*V_i*max(C0 at old/new corners)*max(abs(T at old/new corners)) + 5e-10*(abs(q_i)+abs(qnew_i))`.

Retain the independent native/global residual tests. Require each `abs(qnew_i-q_i) <= r_i` and `abs(sum(qnew-q)) <= sum(r_i)`. Where the caloric projection has a larger local inverse bound, include that measured/model bound, with its axial fractions summed once. Its present formula is `256*eps*V*C0*abs(T) + RTOL*abs(H)`. Do not assign the whole global bound independently to every cell. Report local residuals and local bounds separately. This qualifies an exchange below heat-state resolution, not relative convergence of arbitrarily small exchange.

## Restart and identity review

The latest reviewed load path creates `native.fresh()`, loads into that staged instance, verifies combined state identity, and swaps only after validation. This resolves the earlier concern about a rejected state-digest check occurring after accepted native state mutation.

Further recommendations sent to the root owner:

1. Bind an explicit proof-screen scenario/hash to the joint configuration. The reviewed configuration used a fixed descriptive string while accepting an arbitrary callback, which permits a changed proof screen to retain the same restart identity.
2. Compare native reference density, caloric reference density, and mechanical density at the stress-free reference during construction. A reconstruction at exactly 300 K cannot identify density mismatch because all base reference energies equal zero.
3. Before a staged restart is swapped, also validate the staged mechanical state, proof screen, and homogeneous constrained u/v at its projected temperatures. Shape/finiteness and byte hashes alone do not establish physical admissibility of a consistently re-sealed checkpoint.

Suggested central checks are a near-zero-exchange/conduction fixture, positive and negative work fixtures, independent total-energy reconstruction, rejected staged restart preserving the original digest, changed proof-screen rejection, and density mismatch rejection at exactly 300 K.

## Stable native increment

For solid endpoint temperatures a and b, the exact sealed-model increment is

`Delta h_specific = (b-a)*(472 + 0.068*(a+b) - 2820000/(a*b))`.

The derivative equals `472 + 0.136*T - 2820000/T^2`. On the mushy interval 1672–1727 K, integrate the affine sensible capacity plus `latent_fusion/55` once. Above 1727 K, capacity is 800 J/(kg K). Multiply by caller reference density once to obtain J/m3; the heat solver then applies reference volume once. The formula accepts negative increments and avoids cancellation of large reference energies for a one-ULP temperature increment.
