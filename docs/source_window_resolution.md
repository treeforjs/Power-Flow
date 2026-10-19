# Fixed source-emission windows across saved histories

`mhdlab.source_window_plan.plan_source_windows(time_edges_s, windows,
relative_time_tolerance=...)` maps `SourceTimeWindow(start_s, stop_s)` records
in physical source time to `SourceWindowSelection` offsets in the **original**
frozen-temperature intervals. It reads coordinates, rather than evolving a new
inventory. Source-emission windows do not specify photon recording gates:
material launched earlier can emit during a later recording gate.

Intersections use exact rational arithmetic on the supplied binary64 values.
The result retains the original window index for every local selection, exact
intersection endpoints and duration, projected offsets, canonical parent/local
durations, and signed rounding records. Positive intersections must remain
positive and inside their original interval. Windows must be disjoint and fully
supported by the history; exact adjacency is allowed and input order is kept.
No missing support, collapsed sliver or overlap is silently clipped or repaired.
Accepted endpoint rounding can move cuts or close an exact tiny gap between
long windows; all such motion remains in the projection records.

The caller supplies a relative coordinate tolerance. Local endpoint and duration
errors are measured against each positive intersection duration. A parent-end
duration discrepancy is qualified when that endpoint participates in the
selection; unrelated parent-duration rounding is still recorded. Rounded
absolute endpoints have separate addition and precision records. The immutable
plan replays its input and derived associations during validation and manifest
creation. These are coordinate representation checks on declared times, rather
than bounds on source release, kinetic time sensitivity, physical timing
uncertainty or diagnostic synchronization.

Compare preserved canonical HDF sources with:

```sh
python tools/compare_windowed_source_resolution.py --config study.json --output-dir new_results
```

The configuration has these exact fields:

```json
{
  "schema": "mhdlab-windowed-source-resolution-config-v1",
  "parameter_status": "manufactured",
  "conditioning_provenance": "Declared source-emission restriction only.",
  "comparison_provenance": "Common manufactured drive; preserved source arrays.",
  "relative_time_tolerance": 1e-12,
  "log_mass_tolerance": 1e-8,
  "windows": [{"start_s": 4.91e-8, "stop_s": 4.97e-8}],
  "orders": [8, 16, 32],
  "cases": [{
    "case_id": "saved_source",
    "source_file": "saved.h5",
    "source_sha256": "REPLACE_WITH_EXACT_LOWERCASE_SHA256",
    "evidence_file": "saved_association.json",
    "evidence_sha256": "REPLACE_WITH_EXACT_LOWERCASE_SHA256",
    "resolution": {
      "width_bins": 16, "current_bins": 4,
      "snapshot_intervals": 8, "subdivisions_per_interval": 64
    }
  }],
  "comparisons": []
}
```

Replace the example paths, hashes, resolution and time support with the actual
preserved inputs. Relative paths resolve from the configuration's directory.
Every case needs a hash-associated evidence JSON with schema
`mhdlab-preserved-source-association-v1`, `parameter_status`, exact canonical
`parameter_provenance`, `source_sha256`, the same `resolution`, a nonempty
`preserved_inputs_sha256` digest map, and `source_evolution_repeated: false`.
Additional recovery and qualification details may be retained. The command
checks this association and the canonical HDF, without independently proving
the truth of an evidence declaration or retrieving its ancestral files. The
recovery record must therefore retain its own hash replay and provenance.
Manufactured, sourced, reconstructed and unvalidated labels stay explicit.

Optional comparisons specify `estimate_case_id`, `reference_case_id` and
`family`. For `width`, `current`, `snapshot` and `integration`, exactly the
declared count increases toward the reference while the others stay fixed.
For `snapshot_common_intervals`, finer snapshots have fewer subdivisions so
the total source interval count stays fixed; the actual source time edges
must be identical. This separates snapshot interpolation changes from a change
in source integration intervals. Counts alone do not prove a common physical
temperature drive. The supplied comparison provenance must establish that
association.

The preserved common-grid transformation record, when present, must cross-check
the original snapshot count, its equivalent subdivisions and the actual
integrator's common grid with one subdivision. Its explicit no-smoothing and
no-additional-information flags remain false. Only the two qualified varying
original/equivalent counts are omitted from metadata equality for the
`snapshot_common_intervals` family; all original physical attributes and other
transformation details remain compared. Resampled temperatures do not become
additional independent analytic snapshots.

Compared sources must share pool parameters, parameter provenance, surface
metadata, physical domains and unchanged spatial coordinates. Each pool's
original areal inventory must be uniform across the face and identical between
cases. General spatial inventory remapping is unsupported. Source dimensions
must agree with declared width/current counts and
`snapshot_intervals * subdivisions_per_interval`. Spatial refinement comparisons
also require identical source time edges.
For current refinement, the canonical geometry-derived
`current_patch_intervals_m` list may change with the validated z partition;
other physical metadata remains compared.

The command conditions each requested physical window separately and saves its
joint time/position/thermal moments for every requested quadrature order.
`case_CASE.json` keeps the exact plan, original component associations and
quadrature changes. `report.json` reports estimate-minus-reference source
differences at each fixed order. Quadrature changes are refined-minus-coarse.
Time differences use each result's retained interval origin, window start and
projected local mean before subtraction. Raw mean offsets with different
anchors are never compared directly. Nonzero changes that underflow their
display retain exact represented-input fractions and explicit flags.

These are finite resolution differences. They are neither certified bounds nor
experimental uncertainty; a small quadrature change does not prove accuracy.
The code infers no convergence order. Outputs copy the exact input HDF/evidence
and configuration bytes and bind input/code/output hashes. Existing destinations
and indirect input files are refused. Staging and rehash checks assume a
coordinated single writer, rather than a concurrent transaction protocol.
The `source_sha256` map pins the selected numerical modules listed in
`CODE_PATHS`, rather than every eagerly imported legacy module. The portable
verification subset uses an explicitly identified empty package initializer to
isolate those upstream modules. It retains byte-identical numerical files.

No packet generation, source re-evolution, host budget acceptance, wall cooling,
gap trajectory, molecular chemistry, radiation or actual FLASH run occurs.
The matched-width control changes launch width only, retaining the complete
conditional time/current/state/temperature mixture. Spectral discrimination,
physical H2O coverage and absolute VUV normalization require their own qualified
physical inputs and observation mapping.
