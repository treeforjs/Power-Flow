# Electrode identity of the current foil GIF

The movie shows **one geometric upper foil**, at y=+2.5 to +2.7 mm. That is the
foil corresponding to the **cathode in the dissertation's ALEGRA comparison
geometry**. The native source itself currently labels it only as `upper foil`;
it does not solve or record the cross-gap electrode-potential polarity.

The movie's `top` and `bottom` labels are two broad faces of this same foil.
They are not two different electrodes:

| Native source label | Physical surface | Intended dissertation association |
| --- | --- | --- |
| Upper foil | x=-2.5 to +2.5 mm, y=+2.5 to +2.7 mm | Cathode comparison foil |
| `top` face | y=+2.7 mm, outward +y | Back face of that upper foil, away from the AK gap |
| `bottom` face | y=+2.5 mm, outward -y | Gap-facing face of that upper foil |
| `left` and `right` faces | x=-2.5/+2.5 mm, through the same 200 um thickness | Thickness sides of that upper foil |
| Lower return foil | y=-2.7 to -2.5 mm | Lower/anode counterpart in the comparison geometry |

The current near-foil GIF does not include the lower return electrode or its
desorption. Its normal-axis observation box surrounds only the upper foil.
Upper/lower here are simulation coordinates, not a new claim about laboratory
vertical orientation.

## Dissertation evidence

Source: Trevor Johannes Smith's 2025 dissertation, cached at
`C:/Users/tsmit/Documents/Power Flow/.cache/dissertation/smtrevor_1.pdf`, SHA256
`1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7`.
Complete relevant pages were extracted and visually inspected.

| Figure | Printed / one-based PDF page | Evidence |
| --- | --- | --- |
| 4.2 | 54 / 79 | Target photograph/CAD and 5 mm gap, 5 mm width, 25 mm height. This page does not independently assign electrical polarity to the current native axes. |
| 4.27 | 93 / 118 | The right-hand cross-section explicitly labels the upper red foil `Cathode` and the lower blue foil `Anode`. Current is drawn into the page in the cathode and out of the page in the anode. |
| 4.28 | 94 / 119 | Caption identifies current-density and temperature profiles as cathode profiles. Their plotted foil occupies approximately +2.5 to +2.7 mm from the AK-gap middle. |
| 4.29 | 96 / 121 | Same cathode identification and positive gap-normal foil coordinate for the 60 ns case. |
| 4.30 | 98 / 123 | Caption identifies the width profiles as being across the cathode, sampled 12.5 um inward from the AK gap. |

These geometric/figure associations support describing the shown foil as the
upper cathode-reference foil. This association does not come from the sign of
the native prescribed current alone.

## Native sign and polarity limitation

`native/mfem_edge_source/edge_driver.cpp` constrains the upper foil to +I along
its local z direction and the lower return to -I. The dissertation's displayed
cathode J_z profiles are negative, consistent with the into-page marker in
Fig. 4.27. A mapping between these z/current orientations has not been recorded.
The scalar current trace is used as a magnitude with opposite conductor signs;
it is not a cross-gap voltage measurement.

The solved conductor voltage gradients enforce those longitudinal current
constraints. They do not establish the electrostatic potential difference
between electrodes across the gap. The model also contains no cathode-specific
electron emission, sheath, space-charge or charged-particle heating closure.
Under the present stationary, fixed-property Joule model, a global current/field
sign reversal leaves J squared heating unchanged, which is why a thermal
comparison cannot independently settle electrical polarity.

Future metadata should distinguish the intended geometric role
(`cathode_by_dissertation_geometry`) from a polarity established by an explicit
circuit/voltage model, and should state the current-direction mapping. No
existing source records, fields or historical GIFs were relabeled or rerun in
this audit.
