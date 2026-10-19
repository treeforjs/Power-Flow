# Conservative ledgers for conditioned surface emission

`mhdlab.window_source_ledger.window_source_ledger(source,
include_components=True)` reduces a validated `TemkinPlanarSourceWindows`
record. It retains the original selected intervals, surface geometry, pools,
exact state identifiers, temperatures and source labels. It computes new
bookkeeping from those source windows. Canonical validation replays the declared
stored release law; no replacement source history, packet realization or gap
evolution is generated.

Five extensive quantities retain independent qualifications:

| Quantity | Meaning |
| --- | --- |
| `released_number` | Released parent count |
| `translational_energy_j` | Expected outgoing translational energy |
| `internal_energy_j` | Supplied outgoing internal excitation cost |
| `binding_energy_j` | Supplied adsorbate-to-gas removal cost |
| `total_wall_energy_debit_j` | Declared total source wall debit |

The removal cost retains the interchange name `binding_energy_j`. Its supplied
energy reference is separate from the kinetic activation barrier. Adding the
activation barrier to the wall debit would count an additional cost that this
source model has not declared. Outgoing translation uses the outward-flux
Maxwellian expectation `2*kB*T`; internal and removal costs are separate terms.
An outgoing internal cost is not a resolved molecular level distribution,
daughter yield or collision/radiation calculation.

The report contains resolved and matched-uniform-width ledgers. Both include a
total, width, current, state, reservoir and selection marginals, together with
current/state/width patch groups. Multiple reservoirs with the same exact state
identifier contribute to one state group while their reservoir and component
associations remain available. Exact zero inventories are explicit source skips.
Known geometry/state groups with zero release remain present.
Component-free output is a budget summary. Its separate marginals do not
replace the original joint time/state/temperature source distribution. Preserve
the canonical source and component associations when transporting particles or
predicting spectra.

Represented extensive values are summed as exact rational numbers. An
unqualified contributor makes that quantity's complete group total unavailable;
its reason, logarithm, zero status and source association remain in an
unqualified-contributor catalog. A separately named qualified subtotal may still
be reported. It is not the full group total. No unavailable positive count or
energy becomes a structural zero, and one unavailable channel does not erase a
separately qualified channel.
The catalog covers the five consumed extensive fields. Optional component
manifests also retain the original areal-quantity qualification diagnostics.

Source qualification and floating-point projection are distinct. Exact rational
sums can remain meaningful when a displayed SI value underflows, overflows or
rounds. Projection flags and exact fractions retain that distinction. An
`available` source total does not by itself establish that its float projection
is suitable for a host transaction. The report always sets
`host_budget_accepted` to false.

Each density divides by the group's full declared physical support once:

| Group | Area denominator |
| --- | --- |
| Total, state, reservoir, selection | Full original face `W*L` |
| Width | Original width-bin span times full current length |
| Current | Full width times original current-bin span |
| Current/state/width patch | Original patch area |

These are window-integrated areal densities, rather than fluxes per second.
Multiple selected times, reservoirs and zero-release subpatches do not multiply
or shrink the area denominator. The full-support convention remains the same
when a source emits from only part of the face. Exact coordinate-product areas
and original rounded scalar patch areas retain their signed differences.

The matched-width control allocates each original component's number and every
energy channel by the exact target width fraction `dx/W`. It preserves the
original time/current/state/temperature/cost mixture. Factored current/state
mixtures avoid materializing a component-by-target-width product. This
redistribution is a source comparison; it does not evolve a uniformly heated
surface or change the original depletion.

Signed residuals expose the declared wall debit minus translation, internal and
removal terms, and translation minus the represented count-weighted `2*kB*T`
expectation when their inputs are available. Floating-point source ledgers can
have nonzero residuals. The reducer retains those residuals without normalizing
weights, changing temperatures or repairing energy. Source parameter uncertainty,
kinetic numerical qualification and area/projection rounding remain separate.

Export ledgers from a preserved fixed-window study:

```sh
python tools/export_window_source_ledger.py --study saved_window_study --output-dir new_ledgers
```

The command checks the exact original report/configuration, source/evidence
bytes, selected producer code pins and canonical source association. It
reconstructs the original global-window plan and conditioned source manifest,
then evaluates the new ledger. Saved moment estimates remain input evidence;
they are not recalculated. Copied original files retain their original names
under `study/`. Per-case JSON contains exact group fractions and qualifications;
`ledger.csv` exposes numeric projections and explicit flags/reasons.

A saved source-window packet export can supply the same source association:

```sh
python tools/export_packet_source_ledger.py --packet-export saved_packets --output-dir new_packet_ledgers
```

This command checks the packet producer pins, exact input bytes, packet schema,
dimensions and canonical source selection. It calculates new deterministic
source ledgers and preserves the original files under `packet_export/`. It does
not regenerate the saved realization or recalculate source moments. Its
`ledger.json` retains the original local selections; blank global-window columns
in its CSV mean a composite source selection rather than one inferred photon
gate. Expected source translation remains separate from sampled packet energy.
Original state-mass declarations are retained as provenance and do not alter
the new budget calculation.

Existing destinations and unsafe input names are refused. Staging and final
input/code rehashes assume a coordinated single writer. Selected numerical code
pins do not include every eagerly imported legacy module; the portable
verification subset labels its empty package-initializer wrapper explicitly.

Status labels describe the supplied provenance. Manufactured cases remain
manufactured; grouping them by state does not assign physical water populations
or a monolayer convention. These ledgers do not establish absolute VUV
calibration, photon-gate timing, recorded spectral discrimination, wall feedback,
or an accepted host transaction. Separately owned molecular kinetics and
collision/radiative models consume the source association through their owner
interfaces.
