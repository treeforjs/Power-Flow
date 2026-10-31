"""Calibrate the lumped Mykonos LTD circuit (src/c19_circuit.hpp via build/c19_sim) against the published currents of
Mykonos shots 17133 (Fig. C.19) and 15210 (Fig. C.11), with the reduced Hall vacuum load (L_H from the sealed Hall
operators, R_H(T) from the Hall SS304 law, adiabatic heating).

Fixed (supplied machine data, 2026-10-31): C = 144 nF, V0 = 800 kV, L_ltd = 36 nH, L_feed = 40 nH.
Fitted (model_derived): R_inf, R_sw, tau_sw (series resistance with a switch resistive phase), core branch
1/L_m and G_c at the dI/dt monitor node, and one plotted-time offset per shot (t_plot = t_circuit + offset; the
published time origin is not a trigger reference). Variants change the inductance assumptions (see VARIANTS).

Objective: time-weighted (trapezoid) least squares of the current over every retained published vertex of both shots.
The vertices are correlated display geometry, so the formal 1-sigma values from the Jacobian understate the
uncertainty; the variants give the systematic spread.

Usage (inside run-frozen): python -I calibrate.py C19_SIM LAW_TABLE OUTDIR SCRATCH
"""
from __future__ import annotations

import csv, hashlib, json, subprocess, sys, time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/home/treefor/projects/Power-Flow/data/dissertation")
SHOTS = {17133: DATA / "shot_17133_current.csv", 15210: DATA / "shot_15210_current.csv"}
SHA = {17133: "2c54d2dbaec6b6c8fcf2e45b5abd837aa72999f5fd14217f530277ffd5949076",
       15210: "7e330390e4c678008e1fae0eed983a0283a04c87cd01d3b6118bda99ee1e8cd8"}
L_H_HALL = 1.1045115428330026e-08          # Hall mesh of record: 45 x 45 cells, 9 mm transverse box
L_H_HALL_45MM = 1.6154162174912792e-08     # same operators and 200 um cells, 45.4 mm box (box study: 9/12.6/27/34.6/45.4 mm)

# name -> (fixed overrides, extra free parameter or None, use magnetizing branch)
# The requested minimal model is: series R (with a switch resistive phase) and a core LOSS branch (G_c); the
# magnetizing inductance is an optional extension (variants *_magnetizing).
VARIANTS = {
    "A_supplied_inductances_hall_target": ({"L_feed_H": 40e-9, "L_H": L_H_HALL}, None, False),
    "A2_supplied_inductances_hall_target_magnetizing": ({"L_feed_H": 40e-9, "L_H": L_H_HALL}, None, True),
    "B_supplied_inductances_25nH_target": ({"L_feed_H": 40e-9, "L_H": 25e-9}, None, False),
    "C_feed_inductance_free_hall_target": ({"L_H": L_H_HALL}, "L_feed_H", False),
    "C2_feed_inductance_free_hall_target_magnetizing": ({"L_H": L_H_HALL}, "L_feed_H", True),
    "D_feed_36nH_hall_target": ({"L_feed_H": 36e-9, "L_H": L_H_HALL}, None, False),
    "A3_supplied_inductances_hall_target_45mm_box": ({"L_feed_H": 40e-9, "L_H": L_H_HALL_45MM}, None, False),
    "C3_feed_inductance_free_hall_target_45mm_box": ({"L_H": L_H_HALL_45MM}, "L_feed_H", False),
}
# free parameters: name, initial, lower, upper, scale
FREE = [("R_inf_ohm", .30, 0., 3., .1), ("R_sw_ohm", 1.0, 0., 200., 1.), ("tau_sw_s", 10e-9, .3e-9, 80e-9, 1e-9),
        ("G_core_S", .05, 0., 5., .1)]
MAGNETIZING = ("inv_L_core_per_H", 1e6, 0., 2e8, 1e6)
OFFSETS = [("offset_17133_s", -5e-9, -60e-9, 40e-9, 1e-9), ("offset_15210_s", 5e-9, -60e-9, 40e-9, 1e-9)]
DT, T_END = 1e-10, 280e-9


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_shot(shot):
    p = SHOTS[shot]
    if sha(p) != SHA[shot]:
        raise ValueError(f"shot {shot} CSV changed")
    rows = list(csv.DictReader(p.open()))
    t = np.array([float(r["time_s"]) for r in rows]); I = np.array([float(r["current_A"]) for r in rows])
    w = np.empty_like(t); w[1:-1] = .5 * (t[2:] - t[:-2]); w[0] = .5 * (t[1] - t[0]); w[-1] = .5 * (t[-1] - t[-2])
    return t, I, w


class Model:
    def __init__(self, sim, law, scratch):
        self.sim, self.law, self.scratch = sim, law, Path(scratch)
        self.scratch.mkdir(parents=True, exist_ok=True)
        self.n = 0
        self.cache = {}

    def run(self, params: dict, tag="eval", every=1):
        key = (tuple(sorted(params.items())), every)
        if key in self.cache:
            return self.cache[key]
        pf = self.scratch / f"{tag}.params"; out = self.scratch / f"{tag}.csv"
        pf.write_text("".join(f"{k} {float(v)!r}\n" for k, v in params.items()))
        r = subprocess.run([self.sim, str(pf), self.law, str(out), str(every)], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-2000:])
        a = np.genfromtxt(out, delimiter=",", names=True)
        self.n += 1
        res = (a, r.stdout)
        if len(self.cache) > 64:
            self.cache.clear()
        self.cache[key] = res
        return res


def params_of(x, names, fixed):
    p = {"dt_s": DT, "t_end_s": T_END}
    for k, v in fixed.items():
        if k != "no_core":
            p[k] = v
    for (n, *_), v in zip(names, x):
        p[n] = v
    if fixed.get("no_core"):
        p["inv_L_core_per_H"] = 0.; p["G_core_S"] = 0.
    return p


def metrics(t, I):
    k = int(np.argmax(I)); Ip = I[k]
    def cross(level):
        j = int(np.argmax(I >= level)); return t[j - 1] + (level - I[j - 1]) / (I[j] - I[j - 1]) * (t[j] - t[j - 1])
    t10, t90 = cross(.1 * Ip), cross(.9 * Ip)
    return {"peak_A": float(Ip), "peak_time_s": float(t[k]), "t10_s": float(t10), "t90_s": float(t90),
            "rise_10_90_s": float(t90 - t10), "t10_to_peak_s": float(t[k] - t10)}


def fit_variant(model, name, fixed, free_L, magnetizing, shots):
    names = list(FREE) + ([MAGNETIZING] if magnetizing else [])
    if free_L:
        names = names + [(free_L, 20e-9, .1e-9, 80e-9, 1e-9)]
    allnames = names + OFFSETS
    x0 = np.array([n[1] for n in allnames]); lo = np.array([n[2] for n in allnames]); hi = np.array([n[3] for n in allnames])
    sc = np.array([n[4] for n in allnames])
    Iscale = 1e5

    def resid(z):
        x = z * sc
        a, _ = model.run(params_of(x[:len(names)], names, fixed))
        out = []
        for j, shot in enumerate((17133, 15210)):
            t, I, w = shots[shot]
            off = x[len(names) + j]
            Im = np.interp(t - off, a["t_s"], a["I_port_A"], left=0.)
            out.append(np.sqrt(w / 1e-9) * (Im - I) / Iscale)
        return np.concatenate(out)

    best = None
    import itertools
    starts = []
    for Rinf, Rsw, tau, Gc in itertools.product((.05, .3), (.5, 20.), (1.5e-9, 8e-9), (0., .5)):
        z = x0.copy()
        for k, n in enumerate(allnames):
            z[k] = {"R_inf_ohm": Rinf, "R_sw_ohm": Rsw, "tau_sw_s": tau, "G_core_S": Gc}.get(n[0], z[k])
        starts.append(z)
    for s in starts:
        s = np.clip(s, lo, hi)
        r = least_squares(resid, s / sc, bounds=(lo / sc, hi / sc), x_scale=1., diff_step=1e-4, max_nfev=300)
        if best is None or r.cost < best.cost:
            best = r
    x = best.x * sc
    J = best.jac; dof = max(1, J.shape[0] - J.shape[1]); s2 = 2 * best.cost / dof
    try:
        cov = np.linalg.pinv(J.T @ J) * s2; sig = np.sqrt(np.clip(np.diag(cov), 0, None)) * sc
    except Exception:
        sig = np.full_like(x, np.nan)
    pars = params_of(x[:len(names)], names, fixed)
    a, stdout = model.run(pars, tag=f"best_{name}")
    res = {"variant": name, "fixed": {k: v for k, v in fixed.items()}, "free_inductance": free_L,
           "parameters": {n[0]: float(v) for n, v in zip(allnames, x)},
           "formal_1sigma": {n[0]: float(v) for n, v in zip(allnames, sig)},
           "at_bounds": [n[0] for n, v in zip(allnames, x) if v <= n[2] * (1 + 1e-6) + 1e-30 or v >= n[3] * (1 - 1e-6)],
           "cost": float(best.cost), "nfev": int(best.nfev), "sim_stdout": stdout, "shots": {}}
    tm = a["t_s"]; Im = a["I_port_A"]; Vt = a["V_port_V"]; vA = a["V_port_V"] * 0 + a["v_A_V"]
    mm = metrics(tm, Im)
    for j, shot in enumerate((17133, 15210)):
        t, I, w = shots[shot]
        off = x[len(names) + j]
        Ipred = np.interp(t - off, tm, Im, left=0.)
        r_ = Ipred - I
        mmeas = metrics(t, I)
        rmsw = float(np.sqrt(np.sum(w * r_ ** 2) / np.sum(w)))
        res["shots"][shot] = {
            "offset_s": float(off), "n_vertices": int(len(t)),
            "rms_residual_A_time_weighted": rmsw, "rms_residual_over_peak": rmsw / mmeas["peak_A"],
            "max_abs_residual_A": float(np.max(np.abs(r_))), "max_abs_residual_time_s": float(t[int(np.argmax(np.abs(r_)))]),
            "measured": mmeas,
            "model": {"peak_A": mm["peak_A"], "peak_time_plotted_s": mm["peak_time_s"] + float(off),
                      "rise_10_90_s": mm["rise_10_90_s"], "t10_to_peak_s": mm["t10_to_peak_s"]},
            "peak_error_rel": mm["peak_A"] / mmeas["peak_A"] - 1,
            "rise_10_90_error_rel": mm["rise_10_90_s"] / mmeas["rise_10_90_s"] - 1,
            "t10_to_peak_error_rel": mm["t10_to_peak_s"] / mmeas["t10_to_peak_s"] - 1,
            "peak_time_error_s": mm["peak_time_s"] + float(off) - mmeas["peak_time_s"]}
        with (Path(OUT) / f"fit_{name}_shot_{shot}.csv").open("w", newline="") as h:
            wr = csv.writer(h); wr.writerow(["time_plotted_s", "I_measured_A", "I_model_A", "residual_A"])
            for ti, Ii, Ip in zip(t, I, Ipred):
                wr.writerow([repr(float(ti)), repr(float(Ii)), repr(float(Ip)), repr(float(Ip - Ii))])
    # Main-pulse peaks exclude the first 5 ns after closure: an ideal lumped switch closure divides the full erected
    # voltage inductively at t = 0+ (a sub-ns transient of the lumped model, reported separately).
    main = tm >= 5e-9
    k = int(np.argmax(np.where(main, Vt, -np.inf))); kA = int(np.argmax(np.where(main, vA, -np.inf)))
    res["voltages"] = {"V_target_entrance_peak_V": float(Vt[k]), "V_target_entrance_peak_time_circuit_s": float(tm[k]),
                       "V_target_entrance_peak_over_220kV": float(Vt[k] / 220e3 - 1),
                       "V_target_entrance_closure_transient_max_V": float(np.max(Vt[~main])) if np.any(~main) else None,
                       "V_monitor_peak_V": float(vA[kA]), "V_monitor_peak_time_circuit_s": float(tm[kA]),
                       "V_monitor_peak_over_supplied_two_thirds_of_700kV": float(vA[kA] / (700e3 * 2 / 3) - 1),
                       "ltd_terminal_voltage_peak_V": None}
    # LTD terminal voltage behind the LTD inductance: v_C - R_s(t) i_s ("applied" voltage; supplied: ~140 kV/stage, ~700 kV)
    Rs = pars.get("R_inf_ohm", 0.) + pars.get("R_sw_ohm", 0.) * np.exp(-tm / pars.get("tau_sw_s", 1e-9))
    emf = a["v_C_V"] - Rs * a["i_s_A"]
    ke = int(np.argmax(np.where(main, emf, -np.inf)))
    res["voltages"]["ltd_terminal_voltage_peak_V"] = float(emf[ke])
    res["voltages"]["ltd_terminal_voltage_at_peak_current_V"] = float(np.interp(tm[int(np.argmax(Im))], tm, emf))
    res["voltages"]["monitor_over_ltd_terminal_at_Vmon_peak"] = float(vA[kA] / np.interp(tm[kA], tm, emf))
    res["_trace"] = {"t": tm, "I": Im, "Vt": Vt}
    res["foil_T_final_K"] = float(a["T_foil_K"][-1])
    return res, pars, a


def rating_check(model, pars):
    p = dict(pars); p["resistive_load_ohm"] = .5; p["L_feed_H"] = 1e-12; p["t_end_s"] = 500e-9
    a, _ = model.run(p, tag="rating_0p5ohm")
    t, I = a["t_s"], a["I_port_A"]; k = int(np.argmax(I)); half = .5 * I[k]
    above = np.where(I >= half)[0]; i0, i1 = above[0], above[-1]
    tl = t[i0 - 1] + (half - I[i0 - 1]) / (I[i0] - I[i0 - 1]) * (t[i0] - t[i0 - 1])
    tr = t[i1] + (half - I[i1]) / (I[i1 + 1] - I[i1]) * (t[i1 + 1] - t[i1])
    return {"load": "0.5 ohm resistor at the LTD output (L_feed -> 1 pH)", "peak_A": float(I[k]), "peak_time_s": float(t[k]),
            "fwhm_s": float(tr - tl), "supplied_rating_fwhm_s": 160e-9, "fwhm_error_rel": float((tr - tl) / 160e-9 - 1),
            "peak_load_voltage_V": float(.5 * I[k])}


if __name__ == "__main__":
    sim, law, OUT, scratch = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    Path(OUT).mkdir(parents=True, exist_ok=True)
    shots = {s: load_shot(s) for s in SHOTS}
    model = Model(sim, law, scratch)
    t0 = time.time()
    results = {"schema": "c19-circuit-drive-v1-calibration-v1", "dt_s": DT, "t_end_s": T_END,
               "machine_data_supplied": {"C_F": 144e-9, "V0_V": 800e3, "L_ltd_H": 36e-9, "L_feed_H": 40e-9,
                                              "L_feed_alternative_H": 36e-9, "target_L_estimate_H": [20e-9, 30e-9]},
               "hall_vacuum_load": {"L_H": L_H_HALL, "R_H_300K_ohm": 0.03594536304816677,
                                    "source": "tools/hall_port_load.cpp dc mode (sealed Hall dc_initial, 45x45x(4+2) mesh)"},
               "inputs_sha256": {str(SHOTS[s]): SHA[s] for s in SHOTS}, "variants": {}}
    chosen = chosen_c = chosen_c3 = None
    for name, (fixed, free_L, magnetizing) in VARIANTS.items():
        r, pars, a = fit_variant(model, name, fixed, free_L, magnetizing, shots)
        r["rating_check"] = rating_check(model, pars)
        results["variants"][name] = r
        if name.startswith("A_"):
            chosen = pars
        if name.startswith("C_"):
            chosen_c = pars
        if name.startswith("C3_"):
            chosen_c3 = pars
        print(name, "cost", r["cost"], {s: (round(v["rms_residual_over_peak"], 4), round(v["peak_error_rel"], 4),
              round(v["rise_10_90_error_rel"], 4)) for s, v in r["shots"].items()}, r["voltages"]["V_target_entrance_peak_V"], flush=True)
    # Target-entrance voltage that the Hall vacuum load gives under today's prescribed C19 current:
    # V = L_H dI/dt + R_H(300 K..) I, dI/dt from 2-ns secants of the published polyline (diagnostic only).
    t, I, _ = shots[17133]
    tt = np.arange(t[0], t[-1] - 2e-9, .25e-9); didt = (np.interp(tt + 2e-9, t, I) - np.interp(tt, t, I)) / 2e-9
    Vp = L_H_HALL * didt + 0.0359454 * np.interp(tt + 1e-9, t, I)
    results["prescribed_C19_target_voltage_estimate"] = {
        "method": "L_H * (2-ns secant dI/dt of the published polyline) + R_H(300 K) * I; Hall L_H = 11.045 nH",
        "peak_V": float(Vp.max()), "peak_time_plotted_s": float(tt[int(np.argmax(Vp))] + 1e-9),
        "max_didt_A_per_s": float(didt.max()), "peak_over_220kV": float(Vp.max() / 220e3 - 1),
        "L_target_needed_for_220kV_at_max_didt_H": float((220e3 - 0.0359454 * np.interp(tt[int(np.argmax(didt))] + 1e-9, t, I)) / didt.max())}
    # figures: plot_results.py (re-simulates the stored fits)
    for v in results["variants"].values():
        v.pop("_trace", None)
    results["evaluations"] = model.n; results["wall_s"] = time.time() - t0
    Path(OUT, "fit_results.json").write_text(json.dumps(results, indent=1))
    for fname, vkey, pars, head in (
            ("c19_circuit_v1_parameters_A_supplied.txt", "A_supplied_inductances_hall_target", chosen,
             "variant A: supplied inductances (36 nH LTD, 40 nH feed), Hall target load; does not meet the 10 % rise-time target"),
            ("c19_circuit_v1_parameters_C_traces.txt", "C_feed_inductance_free_hall_target", chosen_c,
             "variant C: feed inductance fitted to the shot traces (model_derived), Hall target load, 9 mm box"),
            ("c19_circuit_v1_parameters_C3_traces_45mm.txt", "C3_feed_inductance_free_hall_target_45mm_box", chosen_c3,
             "variant C3: feed inductance fitted to the shot traces (model_derived), Hall target load with a 45.4 mm box")):
        keep = {k: v for k, v in pars.items() if k not in ("dt_s", "t_end_s")}
        offs = {k: v for k, v in results["variants"][vkey]["parameters"].items() if k.startswith("offset")}
        Path(OUT, fname).write_text(
            f"# c19-circuit-drive-v1 lumped Mykonos LTD parameters ({head})\n"
            "# supplied machine data: C_F V0_V L_ltd_H (and L_feed_H in variant A); Hall-computed: L_H; "
            "model_derived (fit): the rest\n# t_trigger_s = 0 on the circuit clock; plotted time = circuit time + offset_<shot>_s\n"
            + "".join(f"{k} {float(v)!r}\n" for k, v in keep.items())
            + "".join(f"# {k} {float(v)!r}\n" for k, v in offs.items()))
