# Published simulation-current insets

The upper-left panels of dissertation Figs. 4.28 and 4.29 provide finite
current traces displayed with the legacy “Allegra” simulations. Section 4.5
identifies Mykonos shot 17133 as their drive. These reconstructed raster
traces are distinct from the experimental current display in Fig. C.19;
the original simulation input file and clock registration remain unavailable.

| Figure, printed / PDF page | Data and provenance | Retained samples | Displayed time support [ns] |
| --- | --- | --- | --- |
| 4.28, 94 / 119 | [CSV](../data/dissertation/figure_4_28_simulation_drive.csv), [metadata](../data/dissertation/figure_4_28_simulation_drive.json) | 178 | −2.5 to 39.642857 |
| 4.29, 96 / 121 | [CSV](../data/dissertation/figure_4_29_simulation_drive.csv), [metadata](../data/dissertation/figure_4_29_simulation_drive.json) | 262 | −2.5 to 59.642857 |

The literal horizontal axis is **Time [ns]** and the vertical axis is
**Current [MA]**. Five time ticks and six current ticks per panel supply
independent linear calibrations. The zero-time tick lies to the right of
the left frame. Eleven negative-time columns per figure remain unchanged;
no exact zero-time sample, pulse onset or pre-pulse baseline is inserted.
SI columns apply only the declared ns-to-s and MA-to-A conversions.

Each native column retains its exact RGB `(178, 34, 34)` run, median ordinate
and run endpoints. One contiguous run of 2–8 pixels is accepted within the
inspected panel. All 388 candidate columns per figure are catalogued:
178 retained and 210 excluded in Fig. 4.28, and 262 retained and 126 excluded
in Fig. 4.29. A single-pixel terminal cap in Fig. 4.28 is excluded and recorded.
The recovered current maxima are approximately 0.291613 and 0.483268 MA,
respectively, within the retained display segments.

Both curves reach the first interior column of the left frame. Their earlier
history is unknown. The 40/60 ns caption times label the simulation snapshots;
they are not recovered waveform endpoints. No point is padded to those times,
and neither finite segment supplies the later pulse tail. Recovery applies
no smoothing, resampling, interpolation, baseline subtraction or current repair.

Time readout bounds extend each native column by two pixels, approximately
±0.476190 ns. Current bounds contain the complete matched run, extend it by
two pixels on each side and add the maximum current tick-fit residual.
Early negative lower bounds remain signed. These correlated display envelopes
are neither physical uncertainty nor independent errors for a likelihood.
An envelope extending beyond the first or last retained column does not
recover additional source support.

Reproduce the two datasets and their PNG/PDF QA figures into new directories:

```sh
python tools/recover_dissertation_simulation_drive.py --pdf /path/to/smtrevor_1.pdf --output-dir runs/simulation_drive/data --qa-dir runs/simulation_drive/qa
```

The tool checks the PDF, decoded page contents, native RGB pixels and tick
geometry. It refuses to overwrite recovered datasets. Sidecars bind the CSV,
executed generator, source geometry and QA hashes. The source PDF remains
outside this repository.

These curves do not establish the original drive normalization, experimental
trigger alignment, a VUV gate, current-density closure or deposited energy.
They are not FLASH results or validated physical heating inputs.

## Comparison on literal plotted coordinates

The companion verification handoff retains a standalone comparison script,
six exact CSV/JSON inputs, a difference table and PNG/PDF figures. It compares
piecewise-linear displayed curves only on their shared finite support. It
performs no time shift, amplitude fit, time warp or physical registration.
For the central interpolants, extrema are evaluated at every union knot and
overlap endpoint using exact rational arithmetic on the CSV decimal cells.

| Difference | Central overlap [ns] | Largest absolute central difference [kA] | Supported outer-envelope queries |
| --- | --- | --- | --- |
| Fig. 4.28 minus Fig. 4.29 | −2.5 to 39.642857 | 1.52107 | 292 |
| Fig. 4.28 minus C.19 | 0.749788 to 39.642857 | 16.2914 | 355 |
| Fig. 4.29 minus C.19 | 0.749788 to 59.642857 | 16.2914 | 502 |

Both legacy central interpolants lie below the C.19 central interpolant by
approximately 1.88118–16.2914 kA over their respective overlaps; the largest
absolute difference occurs near 33.74983 ns. This comparison describes the
published coordinates, not the actual simulation input or a calibrated current
error. The two legacy insets differ by at most about 1.52107 kA centrally.

For a query to be supported, each entire horizontal readout window must lie
inside its curve's retained support. The outer current interval contains the
minimum lower and maximum upper display envelopes across that window,
including internal knots. This is conservative under a monotone retiming of
piecewise-linear vertex displays within their estimated time/current boxes.
Edge queries without full support are flagged rather than clipped or extrapolated.

None of the supported finite queries excludes zero difference. This is
non-exclusion by an outer set; it proves neither continuous agreement nor the
existence of one jointly consistent retiming, and it does not register physical
clocks. Pixel/axis correlations remain unknown. The display bounds supply no
statistical likelihood, diagnostic uncertainty or justified drive repair.

The comparison report SHA-256 is
`a6b71d79349bb9021fb67f5101daaf01f9f03089ad04b418c04eb466576ddfba`.
