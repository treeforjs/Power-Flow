"""Build the validation-data inventory (INVENTORY.json and INVENTORY.txt) for the Power-Flow coupled PIC
water-breakdown model. Every referenced local file is hashed at build time. Tags:
  measured                     direct experimental record (none of the raw records is in the repository)
  measured_published_display   experimental quantity recovered from a published figure (display geometry)
  measured_derived             experimental construction from several measured signals (e.g. V_bank - L0 dI/dt)
  published_estimate           author-derived estimate from measured spectra by a model fit (e.g. Yacora by eye)
  machine_data_supplied        values supplied directly for this work (2026-10-31), approximate as stated
  model_derived                simulation output or fit (code-to-code reference only)
  recommended                  a value or threshold proposed here
Usage: python -I build_inventory.py OUTDIR
"""
from __future__ import annotations

import hashlib, json, sys
from pathlib import Path

PF = Path("/home/treefor/projects/Power-Flow")
D = PF / "data/dissertation"
HERE = Path(__file__).resolve().parent
WS = Path("/home/treefor/projects/.workspace/outputs")
PDF_SHA = "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"


def sha(p: Path):
    try:
        return hashlib.sha256(p.read_bytes()).hexdigest()
    except OSError:
        return None


def files(*paths):
    return [{"path": str(p), "sha256": sha(p)} for p in paths]


LOS_VUV = ("VUV synthetic diagnostic: (1) state-resolved emissivity of H Lyman-alpha (121.6 nm) and the H2 ~161 nm Lyman band "
           "from the simulated neutral/electron state (CR model with sourced rates; not the whole-level Yacora assignment); "
           "(2) integrate along the line of sight across the 5 mm foil width (the published reduction assumes a 5 mm "
           "emitting depth, Eq. 4.8, p. 92) and over the slit (50 um wide x 25.4 mm tall, slit height along the gap normal, "
           "supplied geometry); (3) gate: 3 ns rectangle, timing registered to the shot clock; (4) spatial PSF/magnification "
           "and the instrumental line-spread function (not recovered; 0.03 nm/pixel is sampling only); (5) apply the "
           "published reduction (transition density per cm^3 with 5 mm depth) to the synthetic radiance, or compare "
           "radiance after the absolute response is recovered; (6) keep the 3D collision volume separate from the "
           "cm^-3 display/slit sampling (A17).")

ITEMS = [
    # ---------------- electrical ----------------
    {"id": "E1", "category": "electrical", "quantity": "load current (published 'Average Current [MA]')", "shot": 17133,
     "tag": "measured_published_display", "availability": "recovered in repository (vector polyline, 132 vertices, 0.75-196.75 ns plotted)",
     "source": {"figure": "C.19", "printed_page": 165, "pdf_page": 190, "files": files(D / "shot_17133_current.csv", D / "shot_17133_current.json")},
     "uncertainty": "display envelope +/-0.327 ns, +/-2.74 kA (stroke + tick fit); B-dot calibration, filtering, averaging over probes and "
                    "trigger origin unpublished; physical uncertainty unknown",
     "simulated_observable": "port current I(t) of the circuit-coupled run (dI/dt-monitor node = feed current), shifted by the calibration "
                             "offset (t_plot = t_circuit + offset_17133)",
     "synthetic_steps": "none beyond the clock offset; the B-dot location (feed start) is the circuit node A",
     "metrics": "A19 M-E1", "gaps": ["trigger/clock origin", "B-dot calibration and probe averaging", "raw scope samples"]},
    {"id": "E2", "category": "electrical", "quantity": "load current", "shot": 15210, "tag": "measured_published_display",
     "availability": "recovered in repository (118 vertices, 20 behind the later legend flagged)",
     "source": {"figure": "C.11", "printed_page": 157, "pdf_page": 182, "files": files(D / "shot_15210_current.csv", D / "shot_15210_current.json")},
     "uncertainty": "display envelope +/-0.327 ns, +/-2.0 kA; same unpublished items as E1; foil thickness caption conflict 200/500 um",
     "simulated_observable": "port current I(t), offset_15210", "synthetic_steps": "clock offset", "metrics": "A19 M-E1",
     "gaps": ["as E1", "foil thickness of the shot (200 vs 500 um captions)"]},
    {"id": "E3", "category": "electrical", "quantity": "resistive voltage V_bank - L0 dI/dt [MV]", "shot": 17133,
     "tag": "measured_derived", "availability": "recovered here (vector, 338 vertices, 0.25-199.75 ns plotted)",
     "source": {"figure": "C.19 (green)", "printed_page": 165, "construction": "dissertation Eq. 4.1-4.2, pp. 59-60",
                "files": files(HERE / "data/fig_C_19_shot_17133_resistive_voltage.csv", HERE / "data/signal_recovery_manifest.json")},
     "uncertainty": "display +/-0.327 ns, +/-2.7 kV; L0 value, V_bank divider calibration and timing corrections unpublished. Shape: rises "
                    "with I to ~370 kV by 56 ns, plateau ~370-383 kV to ~88 ns, falls to ~0 at ~168 ns; V/I ~1 ohm early, far above the "
                    "200 um foil resistance (0.036 ohm at 300 K), so the curve carries an inductive residual",
     "simulated_observable": "V_A(t) - L0 dI/dt with V_A the circuit node-A voltage (water-line output) and L0 the experimenters' inductance "
                             "(treated as a nuisance parameter); equivalently the resistive part of the Hall port voltage plus the feed "
                             "inductance mismatch",
     "synthetic_steps": "clock offset; L0 fit or supplied; same differentiation of the simulated current as the experimenters' B-dot (no extra filtering)",
     "metrics": "A19 M-E3 (timing of features, not amplitude)", "gaps": ["L0 value and its calibration-shot basis", "V_bank monitor calibration and location"]},
    {"id": "E4", "category": "electrical", "quantity": "resistive voltage V_bank - L0 dI/dt [MV]", "shot": 15210,
     "tag": "measured_derived", "availability": "recovered here (vector, 352 vertices)",
     "source": {"figure": "C.11 (green)", "printed_page": 157, "files": files(HERE / "data/fig_C_11_shot_15210_resistive_voltage.csv")},
     "uncertainty": "display +/-0.327 ns, +/-2.0 kV; L0 unpublished. Shape: rises to ~94 kV at 32 ns, collapses to ~0 between ~52 and 60 ns "
                    "(breakdown by the dissertation definition), regrows after 70 ns to ~216 kV at 196 ns",
     "simulated_observable": "as E3; breakdown time = resistive-voltage collapse of the simulated load", "synthetic_steps": "as E3",
     "metrics": "A19 M-T1", "gaps": ["as E3"]},
    {"id": "E5", "category": "timing", "quantity": "Normalized APD 2 (anode-surface light, 2 o'clock view)", "shot": 17133,
     "tag": "measured_published_display", "availability": "recovered here (dashed vector, 330 vertices; 67 behind the legend flagged)",
     "source": {"figure": "C.19 (orange dashed)", "files": files(HERE / "data/fig_C_19_shot_17133_apd2_normalized.csv")},
     "uncertainty": "normalized; onset ~60-64 ns plotted (0.01 at 60 ns, 0.115 at 64 ns); 31 m fiber delay and propagation correction unpublished",
     "simulated_observable": "time history of emission from the anode-surface region in the APD view",
     "synthetic_steps": "collimated view of the anode surface (1-inch mirror), broadband visible response (unknown spectral response): use onset time only",
     "metrics": "A19 M-T2", "gaps": ["APD spectral response, view footprint, fiber delay"]},
    {"id": "E6", "category": "timing", "quantity": "2x Normalized SVS intensity (streaked visible, side-on)", "shot": 15210,
     "tag": "measured_published_display", "availability": "recovered here (vector, 1058 vertices, 65.6-170.2 ns)",
     "source": {"figure": "C.11 (blue)", "files": files(HERE / "data/fig_C_11_shot_15210_svs_2x_normalized.csv")},
     "uncertainty": "normalized; rises from 0.009 at 68 ns; collection 450 um off the cathode surface, 175 um spot (setup Fig. 4.13 is for 15750-15802)",
     "simulated_observable": "visible emission near the cathode surface in the SVS spot", "synthetic_steps": "spot-volume integration, timing only",
     "metrics": "A19 M-T2", "gaps": ["SVS absolute response for 15210", "streak time base calibration"]},
    {"id": "E7", "category": "electrical", "quantity": "current, resistive voltage, APD, SVS for shot 15778 (Fig. 4.7) and current/resistive "
     "voltage for the 17 VUV shots (Figs. A.1-A.17)", "shot": "15778; 15156...15802", "tag": "measured_derived",
     "availability": "published only, not digitized (A-figures are raster)", "source": {"figure": "4.7, A.1-A.17", "printed_pages": "61, 102-118"},
     "uncertainty": "unknown until recovered", "simulated_observable": "as E1-E4", "synthetic_steps": "as E1-E4", "metrics": "A19 M-T1",
     "gaps": ["digitization of Fig. 4.7 and A.1-A.17 electrical panels (breakdown times of every VUV shot)"]},
    {"id": "E8", "category": "electrical", "quantity": "target-entrance AK voltage ~220 kV peak (0 at the short, ~110 kV at the simulated midplane)",
     "shot": "Mykonos planar target (shot not specified)", "tag": "machine_data_supplied",
     "availability": "supplied value (also recorded in artifacts/pre-integration-audit-20261030/reviews/root-data/field-dissertation-comparison/comparison.json)",
     "source": {"files": files(PF / "artifacts/pre-integration-audit-20261030/reviews/root-data/field-dissertation-comparison/comparison.json")},
     "uncertainty": "approximate; use as validation target with +/-15-20 % (scope note 2026-10-31); time waveform not supplied",
     "simulated_observable": "Hall terminal_V (feed-plane gap voltage) peak", "synthetic_steps": "none (direct port quantity)",
     "metrics": "A19 M-E2", "gaps": ["waveform, shot, measurement method of the 220 kV value"]},
    {"id": "E9", "category": "electrical", "quantity": "machine rating and load-voltage division: FWHM 160 ns into 0.5 ohm; ~140 kV/stage "
     "(~700 kV) under load; ~2/3 reaches the dI/dt monitor", "shot": "machine", "tag": "machine_data_supplied",
     "availability": "supplied values", "source": {}, "uncertainty": "approximate (stated)",
     "simulated_observable": "circuit-only: 0.5 ohm load FWHM; LTD terminal voltage v_C - R_s i_s; node-A / LTD-terminal ratio",
     "synthetic_steps": "none", "metrics": "A19 M-E4", "gaps": []},
    {"id": "E10", "category": "timing", "quantity": "published gate rectangles: VUV 26.04-29.04 ns (17133), VUV 52.05-55.05 ns (15210), ICCD 334/340 "
     "windows (15210, not recovered)", "shot": "17133, 15210", "tag": "measured_published_display",
     "availability": "VUV rectangles recovered in repository sidecars; ICCD rectangles not recovered",
     "source": {"files": files(D / "shot_17133_current.json", D / "shot_15210_current.json")},
     "uncertainty": "+/-0.327 ns display; electronics edge definition and trigger reference unpublished; Table 4.1 gives 55 ns 'After' for 15210",
     "simulated_observable": "gate window for every synthetic optical diagnostic", "synthetic_steps": "rectangular 3 ns gate on the shot clock",
     "metrics": "used by A19 M-S*", "gaps": ["gate electronics definition", "ICCD gate coordinates"]},
    # ---------------- spectroscopy ----------------
    {"id": "S1", "category": "spectra/temperature/density", "quantity": "Table 4.1 Yacora estimates: T, n(H2), n(H), gate time and before/after "
     "breakdown for shots 15750, 15210, 15765, 15795, 15778, 15802", "shot": "6 shots", "tag": "published_estimate",
     "availability": "transcribed in repository (recovery_manifest.yaml published_estimates_table_4_1)",
     "source": {"table": "4.1", "printed_page": 89, "files": files(D / "recovery_manifest.yaml")},
     "uncertainty": "fits by eye (p. 90), no covariance; whole-level population assignment; density ranges conflict with the abstract",
     "simulated_observable": "n(H), n(H2) and their temperatures in the VUV view at the gate time", "synthetic_steps": LOS_VUV,
     "metrics": "A19 M-S2", "gaps": ["Yacora decks and electron inputs", "fit covariance", "per-shot gate time on the shot clock"]},
    {"id": "S2", "category": "spectra", "quantity": "H Lyman-alpha and H2 161 nm transition-density gap profiles (central lineouts)",
     "shot": 15750, "tag": "measured_published_display", "availability": "recovered in repository (106 samples)",
     "source": {"figure": "4.20 (= A.13)", "printed_page": 84, "files": files(D / "figure_4_20_profiles.csv", D / "figure_4_20_profiles.json")},
     "uncertainty": "digitization +/-0.0262 dex, +/-0.018 mm; published physical bands not recovered; derived ordinate (transition density with 5 mm depth)",
     "simulated_observable": "gap profile of the synthetic transition density", "synthetic_steps": LOS_VUV, "metrics": "A19 M-S1",
     "gaps": ["published uncertainty bands", "calibrated spectral image", "lineout averaging regions"]},
    {"id": "S3", "category": "spectra", "quantity": "VUV wavelength lineouts (top/middle/bottom)", "shot": "15750, 15765",
     "tag": "measured_published_display", "availability": "recovered in repository (157 and 183 samples)",
     "source": {"figure": "4.20, 4.21", "files": files(D / "figure_4_20_wavelength_lineouts.csv", D / "figure_4_21_wavelength_lineouts.csv")},
     "uncertainty": "display envelopes; LSF unknown; broad strokes at sharp features omitted",
     "simulated_observable": "synthetic spectrum 115-165 nm in the three bands", "synthetic_steps": LOS_VUV + " Convolve with a measured LSF.",
     "metrics": "A19 M-S1 (line ratios: Lyman-alpha to 161 nm band; band shapes)", "gaps": ["instrument LSF", "band collection kernels"]},
    {"id": "S4", "category": "spectra", "quantity": "VUV profiles/lineouts for 15778, 15210, 15795, 15802 (Figs. 4.22-4.25) and 15802 Yacora "
     "comparison (Fig. 4.26)", "shot": "15778, 15210, 15795, 15802", "tag": "measured_published_display",
     "availability": "published only, not digitized", "source": {"figure": "4.22-4.26", "printed_pages": "86-91"},
     "uncertainty": "unknown until recovered", "simulated_observable": "as S2/S3", "synthetic_steps": LOS_VUV, "metrics": "A19 M-S1",
     "gaps": ["digitization of Figs. 4.22-4.26 (15210 is the only VUV shot with a recovered current: Fig. 4.23 first)"]},
    {"id": "S5", "category": "calibration", "quantity": "D2 lamp calibration display and camera gain response", "shot": "instrument",
     "tag": "measured_published_display", "availability": "recovered in repository",
     "source": {"figure": "4.18, 4.19", "files": files(D / "figure_4_18_calibration.csv", D / "figure_4_19_gain.csv",
                                                       D / "calibration_primary_source_manifest.json")},
     "uncertainty": "calibration described as approximate; lamp certificate and camera serial response not recovered",
     "simulated_observable": "absolute radiance -> counts mapping for S2-S4", "synthetic_steps": "instrument response chain", "metrics": "input to M-S1",
     "gaps": ["lamp certificate", "per-shot gain", "absolute response plane"]},
    {"id": "S6", "category": "spectra/temperature", "quantity": "visible SVS spectrum lineout (330 ns) and SVS temperatures with uncertainties "
     "(text pp. 66-75: e.g. 6600+/-1650 K before breakdown and 26200+/-6600 K peak for 15080; 18500+/-3500 K and 6200+/-1300 K for 15166; "
     "4700+/-800 K and 10500+/-3000 K for 15174)", "shot": "15765 (spectrum); 15080, 15166, 15174 (temperatures)",
     "tag": "measured_published_display (spectrum); published_estimate (graybody temperatures, ~25 % error stated)",
     "availability": "spectrum recovered in repository (136 samples); temperatures only in the text",
     "source": {"figure": "4.15; 4.10-4.12", "files": files(D / "observed_spectrum_4_15.csv", D / "observed_spectrum_4_15.json")},
     "uncertainty": "arbitrary intensity; temperature errors as published (~25 %)",
     "simulated_observable": "surface/near-surface plasma temperature in the SVS spot; line identities (O II 464-466 nm, H-alpha 656 nm)",
     "synthetic_steps": "spot-volume emission, graybody + line fit applied identically to the synthetic spectrum",
     "metrics": "A19 M-S3", "gaps": ["SVS absolute response", "15080/15166/15174 are 25-100 um foils (not the 200 um target of the model)"]},
    {"id": "S7", "category": "optics", "quantity": "VUV relay geometry: first mirror ~19 in (0.483 m) from the target, mirror separation ~8 in, "
     "f = 0.1 m mirrors, slit 50 um x 1 in parallel to the foil thickness", "shot": "VUV campaign", "tag": "machine_data_supplied",
     "availability": "recorded in repository (target_optics_geometry.json human_confirmed_geometry)",
     "source": {"files": files(D / "target_optics_geometry.json")}, "uncertainty": "approximate distances; 12 inputs listed as unrecovered",
     "simulated_observable": "LOS registration for S1-S4", "synthetic_steps": "ray-cell projection", "metrics": "input",
     "gaps": ["second-mirror-to-slit distance", "magnification", "LOS angle", "spatial PSF", "LSF", "grating efficiency", "ADU per photon"]},
    # ---------------- field / thermal (legacy simulation; code-to-code) ----------------
    {"id": "M1", "category": "magnetic gap / current depth / temperature (code-to-code)", "quantity": "Allegra B_x gap profiles, J_z depth profiles, "
     "temperature depth profiles at 40 and 60 ns; drive insets", "shot": "17133 drive", "tag": "model_derived",
     "availability": "partially recovered in repository",
     "source": {"figure": "4.28, 4.29", "files": files(D / "figure_4_28_magnetic_gap.csv", D / "figure_4_29_magnetic_gap.csv",
                                                       D / "figure_4_28_current_depth.csv", D / "figure_4_29_current_depth.csv",
                                                       D / "figure_4_28_temperature.csv", D / "figure_4_29_temperature.csv",
                                                       D / "figure_4_28_simulation_drive.csv", D / "figure_4_29_simulation_drive.csv")},
     "uncertainty": "legend units (B_x [100 T], J_z [MA/cm^2]) unverified; clock registration and original drive unknown; digitization bounds only",
     "simulated_observable": "Hall B_x across the gap at the foil middle and edge lineouts; J_z and T vs depth", "synthetic_steps": "same lineout positions; "
     "frame rotation as in comparison.json", "metrics": "A19 M-F1 (informational)",
     "gaps": ["Allegra decks, mesh, normalization", "edge lineout offset"]},
    {"id": "M2", "category": "magnetic / current / temperature across width (code-to-code)", "quantity": "Allegra B_x, J_z and T across the foil width "
     "at 12.5 um depth, 40 ns", "shot": "17133 drive", "tag": "model_derived", "availability": "partially recovered in repository",
     "source": {"figure": "4.30", "files": files(D / "figure_4_30_electromagnetic.csv", D / "figure_4_30_heating.csv")},
     "uncertainty": "as M1; T legend eV interpreted literally", "simulated_observable": "Hall fields and thermal owner T at 12.5 um depth",
     "synthetic_steps": "depth interpolation; Hall 200 um cells do not resolve 12.5 um (refinement series)", "metrics": "A19 M-F1",
     "gaps": ["as M1"]},
    {"id": "M3", "category": "neutral desorption maps (code-to-code)", "quantity": "H, H2, H2O, OH desorption and collisionless flight maps "
     "0.75-60 ns, broad faces and four-face perimeter", "shot": "17133 drive", "tag": "model_derived",
     "availability": "repository artifacts", "source": {"files": files(PF / "configs/neutral_desorption_map_17133_provisional.json",
                                                                   PF / "artifacts/neutral_desorption_20261028/manifest.json",
                                                                   PF / "artifacts/neutral_desorption_perimeter_20261028/manifest.json")},
     "uncertainty": "illustrative kinetics (0.8 eV, 1e12 s^-1, independent pools), constant properties; not measured",
     "simulated_observable": "native neutral census by face and width position", "synthetic_steps": "same binning", "metrics": "A19 M-N1",
     "gaps": ["no measured desorption map exists"]},
    {"id": "N1", "category": "neutral composition", "quantity": "background residual-gas mass spectrum (50 mass channels)", "shot": "chamber background",
     "tag": "measured_published_display", "availability": "vector recovery in the workspace (not in the repository data tree)",
     "source": {"figure": "2.3", "printed_page": 19,
                "files": files(WS / "rga-channel-recovery/figure_2_3_rga_recovery.json", WS / "rga-channel-recovery/figure_2_3_rga_mass_channels.csv",
                               WS / "rga-channel-recovery/figure_2_3_rga_peak_channels.csv")},
     "uncertainty": "plotted partial pressures; supplied total-pressure context ~5e-5 Torr (approximate); peaks include ionizer fragments; no calibration",
     "simulated_observable": "initial background gas composition and pressure (not the per-shot desorbed inventory)",
     "synthetic_steps": "convert to neutral fractions with an RGA cracking/sensitivity model", "metrics": "A19 M-N2",
     "gaps": ["RGA model, electron energy, relative sensitivities, cracking matrix, measurement date, independent total gauge reading"]},
    {"id": "N2", "category": "electron density", "quantity": "fiber-coupled interferometry (Mykonos PDV laser) of gap electron density vs time",
     "shot": "first campaign", "tag": "measured", "availability": "mentioned (Fig. 4.4 caption, ref. [67]); no data in the repository",
     "source": {"figure": "4.4", "printed_page": 57}, "uncertainty": "unknown",
     "simulated_observable": "line-integrated n_e along the probe chord vs time", "synthetic_steps": "chord integration, phase-to-density conversion",
     "metrics": "A19 M-P1", "gaps": ["the interferometry data and chord geometry (would directly test the plasma density, absent otherwise)"]},
    {"id": "N3", "category": "imaging", "quantity": "gated ICCD images (cathode view, 2 o'clock / 4 o'clock)", "shot": "VUV campaigns",
     "tag": "measured", "availability": "not in the repository", "source": {"figure": "4.3, 4.5"}, "uncertainty": "unknown",
     "simulated_observable": "line-of-sight visible emission image", "synthetic_steps": "camera projection, gate", "metrics": "A19 M-T2",
     "gaps": ["images, gate coordinates, response"]},
    {"id": "T1", "category": "timing", "quantity": "breakdown definition: resistive-voltage collapse; visible onset supporting", "shot": "all",
     "tag": "measured_derived", "availability": "text (pp. 59-61); collapse visible in E4 (15210, ~52-60 ns plotted)",
     "source": {"printed_pages": "59-61", "files": files(D / "recovery_manifest.yaml")},
     "uncertainty": "per-shot; the ~50 ns of the conclusions is not a per-shot value", "simulated_observable": "time of simulated resistive-voltage collapse",
     "synthetic_steps": "same definition applied to the simulated load voltage", "metrics": "A19 M-T1", "gaps": ["trigger registration of the plotted clocks"]},
]

A17 = {
    "C19_binding": {
        "bound": "C19 CSV bytes (sha256 2c54d2db...) and parse (Hall CurrentDrive; reproduced bit for bit, tests R1/R2); unshifted plotted clock",
        "open": ["physical trigger/clock origin of the plotted time (calibration offsets are fits, model dependent: -2 to -16 ns for 17133)",
                 "B-dot location, probe averaging and calibration ('Average Current')",
                 "L0 of the V_bank - L0 dI/dt construction and the V_bank monitor calibration",
                 "machine inductance consistency: supplied 36 nH LTD + 40 nH feed + target (11-16 nH Hall, ~25 nH supplied) vs the "
                 "~47 nH total implied by the measured rise with 144 nF / 800 kV",
                 "charge voltage of shots 17133 and 15210 (documents give none; +/-80 kV supplied as machine setting)"]},
    "RGA_calibration": {"open": ["RGA model/electron energy/sensitivities/cracking matrix", "date and relation to the shots",
                                 "total-pressure gauge reading (5e-5 Torr supplied approximate)", "no per-shot inventory measurement"]},
    "slit_volume_separation": {"open": ["slit function and LOS registration (12 unrecovered optics inputs, target_optics_geometry.json)",
                                        "published densities assume a 5 mm emitting depth; keep the 3D model volume separate and apply the reduction to the synthetic radiance",
                                        "instrument LSF and absolute response"]},
}


def main(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    doc = {"schema": "c19-circuit-drive-v1-validation-inventory-v1", "dissertation_pdf_sha256": PDF_SHA,
           "tags": {k: v for k, v in [l.strip().split(None, 1) for l in __doc__.splitlines()[2:9]]},
           "items": ITEMS, "A17_gaps": A17}
    (out / "INVENTORY.json").write_text(json.dumps(doc, indent=1, ensure_ascii=False))
    lines = ["Validation-data inventory for the coupled PIC water-breakdown model (generated by build_inventory.py; machine-readable: INVENTORY.json)",
             "Tags: " + "; ".join(f"{k} = {v}" for k, v in doc["tags"].items()), ""]
    for it in ITEMS:
        lines.append(f"[{it['id']}] {it['quantity']}  (shot {it['shot']}; {it['category']})")
        lines.append(f"    tag: {it['tag']}   availability: {it['availability']}")
        src = it["source"]; fl = src.get("files", [])
        lines.append("    source: " + ", ".join(f"{k} {v}" for k, v in src.items() if k != "files") + ("" if not fl else ""))
        for f in fl:
            lines.append(f"      {f['path']}  sha256 {f['sha256'] or 'MISSING'}")
        lines.append(f"    uncertainty: {it['uncertainty']}")
        lines.append(f"    simulated observable: {it['simulated_observable']}")
        lines.append(f"    synthetic diagnostic: {it['synthetic_steps']}")
        lines.append(f"    metric: {it['metrics']}")
        if it["gaps"]:
            lines.append("    missing: " + "; ".join(it["gaps"]))
        lines.append("")
    lines.append("A17 (data scope) -- what is bound and what is missing")
    for k, v in A17.items():
        lines.append(f"  {k}:")
        if "bound" in v:
            lines.append(f"    bound: {v['bound']}")
        for o in v["open"]:
            lines.append(f"    missing: {o}")
    (out / "INVENTORY.txt").write_text("\n".join(lines) + "\n")
    missing = [f["path"] for it in ITEMS for f in it["source"].get("files", []) if not f["sha256"]]
    print("items", len(ITEMS), "missing files", missing)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
