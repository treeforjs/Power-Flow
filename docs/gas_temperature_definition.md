# Neutral kinetic translational temperature

The neutral temperature diagnostic is computed separately for H, H2, H2O and
OH from the particles present in each gas cell at the stated physical time.
For particle-number weights w and one species mass m:

```math
u_i = \frac{\sum_p w_p v_{p,i}}{\sum_p w_p},
\qquad
C_{ij} = \frac{\sum_p w_p(v_{p,i}-u_i)(v_{p,j}-u_j)}{\sum_p w_p},
```

```math
\Theta_{ij}=\frac{m}{k_B}C_{ij},
\qquad
T_{\mathrm{tr}}=\frac{\mathrm{tr}(\Theta)}{3}.
```

The species bulk velocity is removed before forming the covariance. Mean drift
energy is not counted as thermal energy. The diagonal components of Theta
preserve directional anisotropy; its trace is the scalar displayed in the GIF.
Use a centered second pass over particle velocities to avoid cancellation in
`<v²> - |u|²` for narrow directed distributions. Do not pool different species
using one mass or one bulk velocity.

This is an empirical second-moment diagnostic of collisionless particles. It
does not establish local equilibrium, collisional heating/cooling, an EOS
temperature, electron temperature, or rotational, vibrational or electronic
state temperatures. Ballistic velocity selection and mixtures of streams can
change the local diagnostic without collisions. Foil temperature remains the
separate material-temperature field from MFEM.

## Statistical support

```math
N_{\mathrm{eff}}=\frac{(\sum_p w_p)^2}{\sum_p w_p^2}.
```

The default display requires effective sample size at least 50 and expected
particle number at least one in the cell. Normalize weights by their cell's
largest weight before computing moments and effective sample size; this avoids
underflow of squared weights without changing either result. Unsupported,
empty and solid cells retain NaN temperature. They are not assigned zero kelvin
or filled by interpolation. Supported zero temperature is distinct from NaN.

The covariance above uses the total number weight as its denominator; it does
not silently apply a finite-sample unbiasedness correction. Effective sample
size is a sampling diagnostic, not a physical particle count or an uncertainty
model. Fine spatial cells can lose statistical support even when the source
mesh is refined.

## Flux-weighted source check

An outward thermal flux launch has Gaussian tangential velocities and a
Rayleigh-distributed outward normal velocity. With sigma² = kB*Twall/m, its
normal mean is sigma*sqrt(pi/2), tangential variances are sigma², and centered
normal variance is (2-pi/2)*sigma². Therefore the launch-ensemble scalar is

```math
\frac{T_{\mathrm{tr,launch}}}{T_{\mathrm{wall}}}
=\frac{4-\pi/2}{3}\simeq0.8097346.
```

This is a source-ensemble check, conditioned on the actual birth temperature.
It is not a requirement on every spatial cell after flight. The uncentered
second moment would give 4*Twall/3 and would include the outward drift.

## Saved maps and rendering

`tools/render_neutral_temperature_gif.py` reads HDF5 datasets `time_s`,
`x_edges_m`, `y_edges_m`, UTF8 `species` ordered H,H2,H2O,OH, and
`foil_bounds_m=[xmin,xmax,ymin,ymax]`. `translational_temperature_k`,
`effective_sample_size` and `number_density_m3` each have shape
`(time,4,normal_bin,width_bin)`. Temperature allows unsupported NaNs; finite
values must be nonnegative and satisfy the declared support. Infinite values
are rejected. Optional bulk-velocity and component-temperature arrays need not
be loaded by the renderer.

HDF attributes require `status` and a nonempty JSON object encoded as
`definition_metadata_json`. The JSON sidecar supplies matching status,
`model_description`, `time_registration`, `invariant_current_length_m`, a
nonempty `assumptions` or `source_hypotheses` list, and the support thresholds
`minimum_effective_sample_size` and `minimum_expected_cell_particles`.
The thresholds default to 50 and 1; the shorter `min_Neff` and
`min_expected_cell_particles` aliases are also accepted. The invariant length
converts density and cell area into expected cell particle number once.

```sh
python tools/render_neutral_temperature_gif.py --input runs/case/kinetic_moments.h5 --metadata runs/case/metadata.json --output-dir runs/case/new_temperature_gif --fps 12 --dpi 110
```

The four panels share one fixed linear kelvin scale over all times. Its
default minimum is zero and its maximum is the largest supported input value;
temperatures below 300 K are preserved. An explicitly requested
`--temperature-max-k` changes the color range only, with over-range bins
recorded. Axes use equal physical length scales. Optional
`--view-bounds-mm XMIN XMAX YMIN YMAX` crops only the displayed axes.

HDF frames are streamed rather than loading the whole temperature cube.
Every supplied time is encoded once, without interpolation or smoothing.
The output contains the GIF, first/middle/last previews decoded from that GIF,
and a report of input/source hashes, timestamps, support, display limits and
encoding checks. A successful render does not qualify the underlying source,
temperature estimator or physical model.
