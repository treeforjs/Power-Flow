# Target-to-iCMOS collection geometry evidence

The dissertation supports a conditional photon-budget and edge-overlap
calculation. It does not provide a unique numerical target-to-detector transfer
matrix or absolute ADU prediction. The recovered geometry is recorded in
[`target_optics_geometry.json`](../data/dissertation/target_optics_geometry.json).
The existing [calibration audit](absolute_spectral_calibration.md) and
[manufacturer follow-up](calibration_primary_source_followup.md) remain valid.

## Recovered experimental facts

Page numbers below are printed dissertation page / one-based PDF page.

| Item | Recovered description | Source |
| --- | --- | --- |
| Target | 5 mm gap, 5 mm electrode width, 25 mm height | 54 / 79, Fig.4.2 |
| VUV view | At the 3 o'clock position, looking across the gap; separate from the SVS collection | 56 / 81, Fig.4.3; 58 / 83, Fig.4.5 |
| Relay | Two concave mirrors described as 100 mm and diameter 1 inch, with Al+MgF2 coatings | 82 / 107 |
| Window | An added MgF2 window; reference [69] names Thorlabs WG61050 | 82 / 107; 178 / 203 |
| Limiting path | Text identifies the first-mirror/window path as the Mykonos limiting etendue | 82 / 107 |
| Spectrometer | McPherson 234/302, 0.2 m focal length, concave holographic grating, 1200 grooves/mm, Al+MgF2 | 78-79 / 103-104 |
| Camera | Photek iCMOS 160 with S20 intensifier and MgF2 window at the spectrometer focal plane | 79 / 104 |
| Band / sampling | Approximately 115-165 nm, centered at 140 nm; approximately 0.03 nm/pixel | 79 / 104 |
| Exposure | One frame integrated over a stated 3 ns gate | 79 / 104 |
| Reported camera QE | 19-23% between 115 and 180 nm | 79 / 104 |
| Reported grating reflectance | 75-85% above 115 nm; this is not an order-efficiency measurement | 79 / 104 |
| Lamp setup | 50 micrometre slit and quoted 0.037 mm2 sr etendue | 80,82 / 105,107 |

The last row belongs to the D2 lamp calibration. The user has now confirmed
the same slit width for the collection setup, but that does not transfer the
lamp etendue to the shot relay. The 16-inch chamber radius and mechanical adapters on
page82 are not recovered object distances or pupil locations. Figs.4.3 and4.5
are layout schematics; scaling their drawn rays would not recover metrology.

The phrase describing the relay mirrors does not explicitly label 100 mm as
focal length versus radius of curvature. The
[eSource spherical-mirror catalog](https://www.esourceoptics.com/vuv-uv-spherical-mirrors.html)
lists mirrors by radius of curvature and states f=ROC/2. Therefore, f=100 mm is
an explicit scenario until the installed part or drawing is identified. Its
catalog clear aperture and coating curves are generic product data, not proof
of the historical units' illuminated apertures or incidence angles.

## Human-confirmed collection dimensions

The user supplied the following geometry in this session:

> first collection mirror is approximately 19 inches away. Slit is 50 um by 1 inch height aligned parallel to the foil thickness.

These are trusted session inputs, separate from PDF recovery: approximate
target-to-first-mirror distance 0.4826 m, slit width 50 micrometres, slit height
25.4 mm, with the slit height along the foil-thickness/gap-normal direction y.
The conversion from inches is exact; the physical distance remains approximate
and no measurement tolerance has been supplied. The user subsequently gave an
approximately8-inch (203.2-mm) mirror separation and confirmed that the first
mirror's few-inch distance from the slit is a straight-line distance. The
final second-mirror-to-slit optical path, mirror100-mm focal-length/ROC
interpretation and relay magnification remain unknown.

## Manufacturer constraints that can inform scenarios

The [legacy 234/302 datasheet](https://www.mcphersoninc.com/pdf/234302.pdf) gives
f/4.5, 32-degree incidence, a 40x45 mm grating, adjustable slit width from
10 micrometres to 3 mm with settable height, and a 25 mm focal plane. Its
1200-groove/mm table lists 4 nm/mm dispersion. The quoted scanning resolution
uses particular wavelengths and 10-micrometre slits, so it does not establish
the experimental VUV line-spread function.

The [current iCMOS family page](https://www.photek.com/icmos-camera/) lists a
1920x1200 IMX174 sensor and 5.86-micrometre pixels. These current family values
do not establish the historical custom intensifier-to-sensor mapping. Combining
5.86 micrometres with the generic 4 nm/mm gives 0.02344 nm/pixel, different from
the dissertation's approximate 0.03. Neither number supplies the spatial row
scale or a justified resampling rule.

The dissertation references
[eSource normal-incidence reflectance curves](https://www.esourceoptics.com/vuvmirrorcurves.html)
and the [WG61050 product lead](https://www.thorlabs.com/item/WG61050).
The exact installed coating, lot, incidence and window transmission remain
unrecovered. No flat wavelength throughput or ADU-per-photon factor is adopted.

## Remaining missing inputs that affect the calculation

The human confirmation now supplies the first-mirror distance, approximately
8-inch mirror separation, and slit dimensions/orientation. The first mirror
is a few inches in straight-line distance from the slit; this does not give
the final second-mirror-to-slit optical leg. That leg,
mirror tilts, actual clear pupils and mounted window aperture/position remain
unrecovered. Magnification, spatial row scale, focus plane, point-spread
functions and spectral-spatial coupling also remain unknown. The slit-height
axis is confirmed, while exact line-of-sight registration to simulation x/z
is still a scenario assumption.

For a numerical conditional calculation, keep these as named scenario inputs.
Also retain separate wavelength-dependent mirror reflectances, added-window
transmission, useful-order grating efficiency, detector assembly response and
ADU conversion. Camera QE alone ends at expected photoelectrons; it does not
define the digital counts. An empirical full-train calibration can own some of
these factors already, so their ownership must prevent double multiplication.

Do not borrow the 175-micrometre spot radius, 450-micrometre cathode offset or
0.57 demagnification from page71 / PDF96. Those describe the SVS lens/fiber
system, not the VUV mirror relay.

## Useful side-on overlap scenario

An explicit scenario consistent with the qualitative VUV layout takes the
line of sight along simulation x (foil width), the human-confirmed slit height
along y (gap normal), and current direction z as the other transverse coordinate.
The line of sight remains a coordinate-registration assumption. In that scenario, near-edge,
far-edge and interior volumes can map to the same gap rows; a 2D density-map
pixel is not a separate detector pixel.

Opaque foil occlusion must be included. Far-side side plumes at y within the
foil's thickness can be hidden behind solid material, while near-side plumes
are visible. Finite-pupil rays from broad-face plumes may also clip the foil.
Edge-origin labels and current edge-position masks answer different questions;
preserve both when attributing collected light.

Optical-depth defocus can also mix gap rows. A superseded conditional example
uses the human first-mirror distance u=0.4826 m, assumes f=0.1 m and a symmetric
unit-magnification two-mirror relay, and focuses at the foil-width midplane.
The thin-mirror equation gives intermediate conjugate q=fu/(u-f)=0.126137 m,
mirror separation 2q=0.252274 m, and second-mirror-to-slit distance u. The latter
distances and magnification are derived scenario choices. Its 252-mm mirror
separation conflicts with the later approximately203-mm human estimate, so
it is not the current optical prescription. The current parameter family fixes
the confirmed distances and derives focus-compatible effective focal lengths
for several explicitly assumed final optical legs; these do not identify the
installed mirrors' focal lengths, tilts or curvature.

With the full nominal diameter25.4 mm, the object cone has numerical aperture
about0.02631. A width edge at depth offset2.5 mm then has a paraxial object-space
blur radius about65.8 micrometres. The earlier0.318-mm illustration assumed an
object distance of100 mm and is superseded by the19-inch input. At unit
magnification, the nominal f/4.5 spectrometer cone is wider and does not limit
this on-axis cone; actual stops, tilts and vignetting still require tracing.

The on-axis circular first pupil subtends about0.00217449 sr, or1.7304e-4 of
an isotropic point emitter's photons, before window/mirror/grating losses and
slit clipping. This is geometric collection, not photons reaching the iCMOS or
ADU. These are scenario estimates, not recovered spatial resolution or serial
throughput. A finite-pupil depth kernel is preferable to assigning one
depth-independent Gaussian blur.

Photon transport should begin with CR upper-state emissivity, not neutral
number density alone. For isotropic spontaneous emission, a cell contributes
A_ul*n_u*V*dt photons times its accepted solid-angle fraction and spectral-bin
fraction, followed by explicitly owned losses. A unit-emissivity or unit-radiance
response map can be computed before an absolute source strength is known and
later contracted with the CR solution without adding a fitted slab depth.
