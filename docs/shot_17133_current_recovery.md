# Shot 17133 published current input

[The CSV](../data/dissertation/shot_17133_current.csv) and
[provenance](../data/dissertation/shot_17133_current.json) recover the red
“Average Current [MA]” vector curve in dissertation Fig. C.19, printed p. 165 /
PDF p. 190. The title and caption identify Mykonos shot 17133 and a 200 um SS304
foil. This is the shot named as the drive for the legacy heating profiles in
Figs. 4.28–4.30. The recovered display supplies same-shot current context;
the original simulation drive file and its clock registration remain missing.

The recovery retains 132 vertices over approximately 0.749788–196.75 ns. The largest retained current is 0.820547 MA at 119.25 ns; 0 vertices lie inside the later legend's bounding rectangle. The display bounds are approximately +/-0.327196 ns and +/-2738.32 A. The independent gate rectangle spans approximately 26.0393–29.0394 ns. These are published display quantities, with physical trigger and simulation-clock alignment unresolved.

The shared ordinate is “Signal Height”; the current's units come from its
own MA legend. Other curves on that axis are separate quantities. The recovery
retains native PDF vertices, their page coordinates, signed current, and SI
columns `time_s` and `current_A`. It inserts no endpoints or prehistory and
applies no smoothing, baseline correction or extrapolation. These vertices
are correlated published display geometry with unknown original filtering,
rather than acquisition samples.

The sidecar pins the PDF, decoded page content, current path, tick coordinates,
generator and CSV. Its display envelopes combine half the plotted stroke width
with the largest tick-fit residual. They are readout bounds, not diagnostic
calibration errors, physical uncertainty or independent likelihood variances.
The VUV rectangle is retained separately. Its plotted location establishes
neither breakdown nor a physical trigger or a registered simulation clock.

Recreate the recovery from the hash-matched dissertation:

```sh
python tools/recover_dissertation_drive_17133.py --pdf /path/to/smtrevor_1.pdf --qa-plot runs/recovery/shot_17133_current_qa.png
```

Inspect this shot with an explicit coordinate choice:

```sh
python tools/inspect_current_drive.py --csv data/dissertation/shot_17133_current.csv --provenance data/dissertation/shot_17133_current.json --preserve-source-axis --json-out runs/current_17133/summary.json --plot runs/current_17133/current.png
```

`--preserve-source-axis` retains the published clock with zero offset. It does
not align the current to the 40/60 ns heating profiles or to a hardware trigger.
An explicitly justified later offset must remain a separate declaration.
`DriveTrace` interpolates only inside the retained finite support. Its optional
integrals describe that interpolant; an integral of I² is not deposited heat.

This recovery supplies no circuit voltage, measured dI/dt, electromagnetic
normalization, complete current-density cross-section, Joule-heating closure or
FLASH result. The [shot 15210 input](../data/dissertation/shot_15210_current.json)
remains an independent shot.
