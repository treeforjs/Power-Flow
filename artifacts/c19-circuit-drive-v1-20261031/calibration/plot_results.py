"""Figures of the calibration (re-simulated from fit_results.json with build/c19_sim, so they show the stored fits).
  calibration_currents.png  measured currents of 17133 and 15210 with variants A (supplied inductances) and C3
  port_voltage.png          target-entrance (Hall port) voltage of variants A, C and C3 against the supplied ~220 kV
Usage (inside run-frozen): python -I plot_results.py C19_SIM LAW_TABLE CALIBRATION_DIR SCRATCH
"""
from __future__ import annotations

import json, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate as C  # noqa: E402

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
SERIES = {"A_supplied_inductances_hall_target": ("#2a78d6", "A: supplied 40 nH feed, Hall 9 mm box"),
          "C_feed_inductance_free_hall_target": ("#1baf7a", "C: feed fitted, Hall 9 mm box"),
          "C3_feed_inductance_free_hall_target_45mm_box": ("#eb6834", "C3: feed fitted, Hall 45 mm box")}


def main(sim, law, cal, scratch):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    res = json.loads(Path(cal, "fit_results.json").read_text())
    model = C.Model(sim, law, scratch)
    runs = {}
    for name in SERIES:
        v = res["variants"][name]
        pars = {"dt_s": C.DT, "t_end_s": C.T_END}
        pars.update({k: val for k, val in v["fixed"].items()})
        pars.update({k: val for k, val in v["parameters"].items() if not k.startswith("offset")})
        a, _ = model.run(pars, tag="plot_" + name)
        runs[name] = (a, v["parameters"])
    shots = {s: C.load_shot(s) for s in C.SHOTS}
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    for ax, shot in zip(axes, (17133, 15210)):
        t, I, _ = shots[shot]
        ax.plot(t * 1e9, I / 1e6, color=INK, lw=2, label="measured (published polyline)")
        for name in ("A_supplied_inductances_hall_target", "C3_feed_inductance_free_hall_target_45mm_box"):
            a, p = runs[name]; col, lab = SERIES[name]
            ax.plot((a["t_s"] + p[f"offset_{shot}_s"]) * 1e9, a["I_port_A"] / 1e6, color=col, lw=1.6, label=lab)
        ax.set_xlim(0, 200); ax.set_xlabel("plotted time (ns)"); ax.set_ylabel("current (MA)")
        ax.set_title(f"Mykonos shot {shot}", color=INK, loc="left"); ax.grid(color=GRID, lw=.6)
    axes[0].legend(frameon=False, loc="lower right")
    fig.savefig(Path(cal) / "calibration_currents.png", dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 3.8), layout="constrained")
    for name, (col, lab) in SERIES.items():
        a, p = runs[name]
        ax.plot((a["t_s"] + p["offset_17133_s"]) * 1e9, a["V_port_V"] / 1e3, color=col, lw=1.6, label=lab)
    ax.axhline(220, color=MUTED, lw=1, ls="--"); ax.text(4, 226, "supplied value ~220 kV (validation target, +/-15-20 %)", color=MUTED)
    ax.set_xlim(0, 200); ax.set_ylim(-50, 320); ax.set_xlabel("plotted time, shot 17133 (ns)")
    ax.set_ylabel("target-entrance voltage (kV)")
    ax.set_title("Hall port voltage with the reduced Hall vacuum load", color=INK, loc="left")
    ax.grid(color=GRID, lw=.6); ax.legend(frameon=False, loc="upper right")
    fig.savefig(Path(cal) / "port_voltage.png", dpi=150); plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    main(*sys.argv[1:5])
