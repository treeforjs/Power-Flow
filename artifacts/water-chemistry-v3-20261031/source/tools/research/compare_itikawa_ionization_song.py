"""Compare LXCat Itikawa-database H2O ionization blocks against Song 2021 Table 12
values held in the frozen dataset (common energies only, no interpolation).

Usage: python3 -I compare_itikawa_ionization_song.py <lxcat_itikawa_txt> <frozen_dataset_json>
"""
import json
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lxcat_window_scan import blocks  # noqa: E402


def main():
    lx, frozen = sys.argv[1], sys.argv[2]
    d = json.load(open(frozen))
    ch = {c["id"]: c for c in d["channels"]}
    want = {
        "Total Ionization.": "H2O:ion_yield:total",
    }
    out = []
    for kind, meta, pts in blocks(lx):
        if kind != "IONIZATION":
            continue
        proc = meta.get("PROCESS", [""])[0]
        comment = " ".join(meta.get("COMMENT", []))
        if "Total Ionization" in comment:
            cid = "H2O:ion_yield:total"
        elif "H2O+" in proc:
            cid = "H2O:ion_yield:H2O+"
        elif "OH+" in proc:
            cid = "H2O:ion_yield:OH+"
        elif "O++" in proc:
            cid = "H2O:ion_yield:O_atomic_2plus"
        elif "O+" in proc:
            cid = "H2O:ion_yield:O+"
        elif "H2+" in proc:
            cid = "H2O:ion_yield:H2+"
        elif "H+" in proc:
            cid = "H2O:ion_yield:H+"
        else:
            continue
        c = ch[cid]
        scale = float(c["ordinate_scale_to_SI"])
        song = {float(e): float(v) * scale for e, v in zip(c["energy_eV"], c["values"])}
        n = 0
        mism = []
        for e, s in pts:
            e = float(e); s = float(s)
            if e in song:
                n += 1
                ref = song[e]
                if ref == 0 and s == 0:
                    continue
                rel = (s - ref) / ref if ref else float("inf")
                if abs(rel) > 1e-6:
                    mism.append((e, s, ref, round(rel, 4)))
        out.append({"lxcat_process": proc, "lxcat_comment": comment, "song_channel": cid,
                    "song_locator": c.get("source_locator"), "n_lxcat_points": len(pts),
                    "n_common_energies": n, "mismatches(E, lxcat_m2, song_m2, rel)": mism,
                    "lxcat_E_range": [float(pts[0][0]), float(pts[-1][0])] if pts else None})
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
