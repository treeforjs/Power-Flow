# Current-drive inputs for FLASH heating

Inspect the recovered shot with an explicit choice of time coordinate:

```bash
python tools/inspect_current_drive.py --preserve-source-axis \
  --json-out runs/current_inspection/summary.json \
  --plot runs/current_inspection/current.png
```

Alternatively, `--time-offset-s` records a chosen translation. The inspection
checks the CSV against its source hash and reports the finite support, signed
values, legend-covered vertices and independent VUV gate geometry. It never
identifies breakdown or establishes trigger alignment from that rectangle.

`mhdlab.drive_trace` supplies a bounded, explicitly sourced SI current history.
It is independent of the legacy trace loader and does not fit a circuit or infer
missing voltage. Interpolation outside the published or measured support raises
an error rather than padding, clamping or assigning zero current.

The first recovered input is
[`shot_15210_current.csv`](../data/dissertation/shot_15210_current.csv), with
[`provenance`](../data/dissertation/shot_15210_current.json). It contains 118
vertices from the native vector current curve in dissertation Fig. C.11, printed
p. 157 / PDF p. 182. Its supported times are approximately 2.250–198.250 ns and
the plotted peak is approximately 0.848 MA. Encoded portions hidden by the legend
are flagged. It is a recovered published curve, not the original oscilloscope
waveform. Early signed values are retained.

```python
from mhdlab.drive_trace import load_current_csv

drive = load_current_csv(
    "data/dissertation/shot_15210_current.csv",
    time_column="time_s", current_column="current_A",
    time_scale_to_s=1.0, current_scale_to_a=1.0,
    diagnostic_time_offset_s=0.0,
    source_kind="reconstructed", shot_id="15210",
    source_reference="data/dissertation/shot_15210_current.json",
    figure_reference="Dissertation Fig. C.11, printed p.157 / PDF p.182",
)
current = drive.current_at(55e-9)
integral_i2 = drive.integrate_current_squared_a2_s(10e-9, 100e-9)
```

The zero offset in this example preserves the published plot coordinate. It
does not establish the hardware trigger origin or align the waveform with FLASH
initial conditions. Supply any later timing correction explicitly and retain
the supporting calibration. A FLASH drive callback must query only the known
time support. Initializing earlier requires additional measured evidence or an
explicitly labeled model assumption; it cannot silently acquire a zero-current
prehistory from this contract.

The methods integrate I and I² exactly for the piecewise-linear interpolant,
including partial segments and sign changes. I² is a current-drive diagnostic;
it is not deposited Joule energy without the actual current distribution and
resistive material state.

`slope_a_s(time, knot_policy="left"|"right"|"raise")` differentiates the
interpolation. The dissertation reports calibrated B-dot measurements for dI/dt
and propagation timing corrections. A derivative of the recovered plotted
current is not that original measurement. Do not use it to claim a recovered
resistive voltage or measured circuit inductance.

Shot 15210 belongs to the VUV campaign. The legacy heating profiles in
Figs. 4.28–4.30 used drive shot 17133. These inputs cannot be combined and labeled
a same-shot reproduction. Their role is separate source/thermal benchmarking
until a consistent heating drive and material setup are available.
