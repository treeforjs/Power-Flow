"""List LXCat cross-section points of IONIZATION blocks within an energy window.

Usage: python3 -I lxcat_window_scan.py <emin_eV> <emax_eV> <lxcat_txt> [<lxcat_txt> ...]
Prints the block header (type, PROCESS, PARAM) and every (E, sigma) pair with
emin <= E <= emax, verbatim as text from the file, plus the nearest points
bracketing the window.
"""
import sys


def blocks(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() in ("IONIZATION", "EXCITATION", "ATTACHMENT", "ELASTIC", "EFFECTIVE", "ROTATIONAL"):
            kind = lines[i].strip()
            j = i + 1
            meta = {}
            while j < len(lines) and not lines[j].startswith("-----"):
                if ":" in lines[j] and lines[j].split(":")[0] in ("PROCESS", "PARAM.", "COMMENT"):
                    meta.setdefault(lines[j].split(":")[0], []).append(lines[j].split(":", 1)[1].strip())
                j += 1
            k = j + 1
            pts = []
            while k < len(lines) and not lines[k].startswith("-----"):
                parts = lines[k].split()
                if len(parts) == 2:
                    pts.append((parts[0], parts[1]))
                k += 1
            yield kind, meta, pts
            i = k + 1
        else:
            i += 1


def main():
    emin, emax = float(sys.argv[1]), float(sys.argv[2])
    for path in sys.argv[3:]:
        for kind, meta, pts in blocks(path):
            if kind != "IONIZATION":
                continue
            inwin = [p for p in pts if emin <= float(p[0]) <= emax]
            below = [p for p in pts if float(p[0]) < emin]
            above = [p for p in pts if float(p[0]) > emax]
            print(path.split("/")[-2], "|", meta.get("PROCESS"), "|", meta.get("PARAM."))
            print("   last below:", below[-1] if below else None, " in window:", inwin, " first above:", above[0] if above else None)


if __name__ == "__main__":
    main()
