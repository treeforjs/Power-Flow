"""Build data/sources/ELECTRON-SOURCES.json for the electron + H2O collision and swarm sources.

Usage: python3 -I build_electron_sources_index.py <component_root> <manifest_swarm.json> <manifest_lit.json>

Computes sha256 of every file in each listed source directory, collects every extracted
table whose source_id is one of these sources, and merges the missing-items lists.
"""
import datetime
import glob
import hashlib
import json
import os
import sys

ROOT, MAN_SWARM, MAN_LIT = sys.argv[1:4]
SRC = os.path.join(ROOT, "data", "sources")
EXT = os.path.join(ROOT, "data", "extracted")
ms = json.load(open(MAN_SWARM, encoding="utf-8"))
ml = json.load(open(MAN_LIT, encoding="utf-8"))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read(d, name):
    p = os.path.join(SRC, d, name)
    return open(p, encoding="utf-8").read().strip() if os.path.exists(p) else None


def crossref_match(relpath):
    d = json.load(open(os.path.join(SRC, relpath), encoding="utf-8"))
    m = d["message"]
    if "items" in m:
        m = m["items"][0]
    au = ", ".join(a.get("family", "") for a in (m.get("author") or [])[:8])
    return (f"{m.get('DOI')}: '{(m.get('title') or [''])[0]}', {(m.get('container-title') or [''])[0]} "
            f"{m.get('volume') or ''}, {m.get('page') or m.get('article-number') or ''} "
            f"({(m.get('issued', {}).get('date-parts') or [[None]])[0][0]}); authors {au} [record: data/sources/{relpath}]")


LXCAT_PLATFORM = [
    "S. Pancheshnyi et al., Chem. Phys. 398, 148-153 (2012), doi:10.1016/j.chemphys.2011.04.020",
    "L. C. Pitchford et al., Plasma Process. Polym. 14, 1600098 (2017), doi:10.1002/ppap.201600098",
    "E. Carbone et al., Atoms 9, 16 (2021), doi:10.3390/atoms9010016"]
LXCAT_PLATFORM_CR = [crossref_match("crossref-electron-cs-20261031/pancheshnyi-2012-chemphys.json"),
                     crossref_match("crossref-electron-cs-20261031/pitchford-2017-ppp.json"),
                     crossref_match("crossref-electron-cs-20261031/carbone-2021-atoms.json")]
LXCAT_TERMS = ("LXCat terms of use (accepted on the site before download): databases remain property of their "
               "contributors, are not to be redistributed by third parties or used commercially; cite the database, "
               "www.lxcat.net and the retrieval date. Pages saved as lxcat_how_to_reference_page.html and "
               "lxcat_redistribution_page.html in each LXCat cross-section directory.")

S = []


def lx(db, dirname, citation, cr, url, summary):
    S.append({
        "source_id": dirname,
        "citation": citation + f" Retrieved from www.lxcat.net (nl.lxcat.net mirror) on 31 October 2026.",
        "doi": None,
        "crossref_verified": cr,
        "url": url,
        "contents_summary": summary + " " + LXCAT_TERMS,
        "lxcat_platform_references": LXCAT_PLATFORM,
        "lxcat_platform_references_crossref": LXCAT_PLATFORM_CR,
    })


lx("Itikawa", "lxcat-itikawa-h2o-cs-20261031",
   "Itikawa database, www.lxcat.net/Itikawa (contact Y. Itikawa); data group: Y. Itikawa and N. Mason, Cross Sections for "
   "Electron Collisions with Water Molecules, J. Phys. Chem. Ref. Data 34, 1-22 (2005), doi:10.1063/1.1799251.",
   {"value": True, "matched": crossref_match("crossref-electron-lit-20261031/itikawa-mason-2005.json"),
    "note": "The database itself has no DOI; the cited data-group reference was verified. Full text of Itikawa & Mason 2005 not obtained, so the LXCat blocks could not be mapped to the paper's table numbers."},
   "https://nl.lxcat.net/data/ (browse flow; db 13, species H2O, all 28 processes); cache https://nl.lxcat.net/cache/6ae60f2f19116/Download.zip",
   "28 e/H2O blocks: 3 dissociative attachment (H-, O-, OH-), elastic integral, momentum transfer, total scattering, rotational "
   "J0->0 (labelled excitation, 0 eV), J0->1/2/3, vibrational (010) and (100+001), OH(X) production (Harb et al. 2001), "
   "O(1S) production (Kedzierski et al. 1998), OH(A-X) emission (Beenakker et al. 1974), H Balmer alpha (Mohlmann and de Heer "
   "1979), H Balmer beta, O 777.4 nm and O 844.7 nm emission (Beenakker et al. 1974), H Lyman alpha and O 130.4 nm emission "
   "(Morgan and Mentall 1974, renormalized), partial ionization H2O+, H+, OH+, O+, H2+, O++ and total ionization (13.5-1000 eV). "
   "Not marked 'complete set' on LXCat.")
lx("Hayashi", "lxcat-hayashi-h2o-cs-20261031",
   "Hayashi database, www.lxcat.net/Hayashi; data group: M. Hayashi, 'Electron collision cross sections for molecules "
   "determined from beam and swarm data', in Swarm Studies and Inelastic Electron-Molecule Collisions, eds. L. C. Pitchford, "
   "B. V. McKoy, A. Chutjian, S. Trajmar (Springer, New York, 1987), pp. 167-187; digitized by S. Chowdhury, LAPLACE.",
   {"value": True, "matched": crossref_match("crossref-electron-cs-20261031/hayashi-1987-chapter.json"),
    "note": "Book chapter verified; database has no DOI."},
   "https://nl.lxcat.net/data/ (browse flow; db 11, species H2O, all 8 processes); cache https://nl.lxcat.net/cache/6ae60f38a4d73/Download.zip",
   "8 blocks: total attachment, elastic momentum transfer, excitation at 0.206, 0.459, 1.058, 8.445 and 14.052 eV, total "
   "ionization (13.76 eV). Swarm-derived set; LXCat warns that rotational cross sections are not included.")
lx("Morgan", "lxcat-morgan-h2o-cs-20261031",
   "Morgan database (W. L. Morgan, Kinema Research & Software), www.lxcat.net/Morgan; data-group comment refers to Morgan and "
   "Penetrante (1990) for the CAR rotational parameters; several blocks cite Hayashi.",
   {"value": True, "matched": crossref_match("crossref-electron-cs-20261031/morgan-penetrante-1990-cpc.json"),
    "note": "Only the 1990 paper named in the data-group comment could be verified; the database has no DOI and no per-process primary reference."},
   "https://nl.lxcat.net/data/?pre=cs_set (complete-set download, s[1:44]=16:142)",
   "10 blocks (complete set): 3 dissociative attachment, effective momentum transfer, vibrational V2 0.2 eV, V1 0.43 eV, V3 "
   "0.43 eV, dissociative excitation H+OH 7 eV, O*+H2 13 eV, ionization 13 eV. Data-group comment gives CAR parameters "
   "'1.78e-3 0.50 2' and warns that rotational cross sections are not included explicitly.")
lx("MuroranIT", "lxcat-muroranit-h2o-cs-20261031",
   "MuroranIT database (S. Kawaguchi, Muroran Institute of Technology), www.lxcat.net/MuroranIT; data group: S. Kawaguchi, "
   "K. Takahashi, K. Satoh, H. Itoh, 'Electron transport analysis in water vapor', Jpn. J. Appl. Phys. 55, 07LD03 (2016).",
   {"value": True, "matched": crossref_match("crossref-swarm-20261031/kawaguchi-2016.json"),
    "note": "Data-group reference verified; database has no DOI."},
   "https://nl.lxcat.net/data/?pre=cs_set (complete-set download, s[1:44]=836:11880)",
   "41 blocks (complete set, updated 2025-07-11): 3 dissociative attachment, elastic momentum transfer, rotational J=1/2/3, "
   "vibrational (010) and (100)+(001), 26 electronic-excitation blocks (a~3B1 7.14 eV ... H2O* 12.5 eV, including B~1A1 9.73 eV), "
   "partial ionization H2O+ (12.65 eV), H++OH, OH++H, O++H2, H2++O, O2++H2. Swarm-validated set (Monte Carlo vs measured "
   "swarm parameters).")
lx("Phelps", "lxcat-phelps-h2o-cs-20261031",
   "Phelps database (A. V. Phelps, JILA), www.lxcat.net/Phelps; data group 'H2O NOV 1983' (Davies, unpublished); ionization "
   "block 'Total ionization FROM SHUTTEN ET AL.' i.e. J. Schutten et al., J. Chem. Phys. 44, 3924 (1966).",
   {"value": True, "matched": crossref_match("crossref-electron-lit-20261031/schutten-1966.json"),
    "note": "Only the Schutten et al. reference named in the ionization block is verifiable; the set itself is unpublished (Davies)."},
   "https://nl.lxcat.net/data/ (browse flow; db 19, species H2O, all 10 processes); cache https://nl.lxcat.net/cache/6ae60f419bfec/Download.zip",
   "10 blocks: O- attachment, effective momentum transfer, vibrational 0.198 and 0.469 eV (Linder), energy-loss processes 1.1, "
   "6.3, 9, 12, 12.5 eV, total ionization (12.6-200 eV). LXCat notes (Oct 2024) that the set is no longer marked complete "
   "because rotational cross sections are not included.")
lx("QUANTEMOL", "lxcat-quantemol-h2o-cs-20261031",
   "QUANTEMOL database (Quantemol Ltd), www.lxcat.net/QUANTEMOL; BEB total ionization computed with Quantemol-N (UK R-matrix "
   "target models), how to reference: http://quantemol.com.",
   {"value": False, "matched": "No DOI or journal reference is given by the database for these blocks; nothing to verify."},
   "https://nl.lxcat.net/data/ (browse flow; db 21, species H2O, all 5 processes); cache https://nl.lxcat.net/cache/6ae60f678c789/Download.zip",
   "5 BEB total-ionization blocks, 13.67/13.72-4993.7 eV, 499 points each, differing in target model (CI vs Hartree-Fock, "
   "DZP vs 6-311G basis, C2/Cs/C2v point group; Koopmans binding energies listed in each COMMENT).")
lx("TRINITI", "lxcat-triniti-h2o-cs-20261031",
   "TRINITI database (EEDF package, A. P. Napartovich, N. A. Dyatko, I. V. Kochetov, A. G. Sukharev, TRINITI, Troitsk), "
   "www.lxcat.net/TRINITI; H2O blocks cite 'M. Yousfi et al., Cent. de Phys. Atom. de Toulouse, Univ. di Padova 1987'.",
   {"value": False, "matched": "The cited Yousfi et al. 1987 report has no DOI; not verifiable via Crossref."},
   "https://nl.lxcat.net/data/?pre=cs_set (complete-set download, s[1:44]=24:0)",
   "8 blocks (complete set): O- and H- attachment, effective momentum transfer, rotational (lumped, 0.04 eV), vibrational "
   "(010) and (101), dissociative excitation H+OH 7.1 eV, ionization 12.61 eV.")

S.append({
    "source_id": "song2021-jpcrd-reference",
    "citation": ("M.-Y. Song, H. Cho, G. P. Karwasz, V. Kokoouline, Y. Nakamura, J. Tennyson, A. Faure, N. J. Mason, Y. Itikawa, "
                 "Cross Sections for Electron Collisions with H2O, J. Phys. Chem. Ref. Data 50, 023103 (2021)."),
    "doi": "10.1063/5.0035315",
    "crossref_verified": {"value": True, "matched": crossref_match("crossref-electron-lit-20261031/song-2021.json")},
    "url": "https://discovery.ucl.ac.uk/id/eprint/10128988/1/5.0035315.pdf (download URL recorded by the frozen dataset)",
    "contents_summary": ("Pointer directory only: the PDF is read in place at "
                         "/home/treefor/sim-workspace/powerflow/outputs/h2o-collision-data-20261029/source-documents/song-2021.pdf "
                         "(sha256 d94ae0c73ab745a3e5b958d096f2c774b123b8257877e2d71911990e9b83b595, identical to the frozen "
                         "dataset's song2021 pdf_sha256) and is not copied. Tables 3, 4, 5, 6/7, 8/9 (A1B1), 10, 12, 13 are in "
                         "the frozen dataset /home/treefor/projects/runtimes/h2o-collision-data-20261029/data/collisions/"
                         "h2o-gas-20261029.json and are not re-transcribed. Transcribed here: Table 1 (TCS), Table 2 (elastic "
                         "DCS, 10 energies), Table 8 (six state energies), Table 11 (O(1S0) production), running-text numeric "
                         "statements (Secs. 2, 9, 10)."),
    "primary_references_verified": [
        "Table 2: " + crossref_match("crossref-electron-cs-20261031/matsui-2016-epjd.json"),
        "Table 8: " + crossref_match("crossref-electron-cs-20261031/thorn-2009-pmcphysb.json"),
        "Table 11: " + crossref_match("crossref-electron-cs-20261031/kedzierski-1998-jpb.json"),
        "Table 11 compilation: " + crossref_match("crossref-electron-cs-20261031/mcconkey-2008-physrep.json"),
        "Sec. 8 states: " + crossref_match("crossref-electron-cs-20261031/ralphs-2013-jpb.json"),
        "Table 1 low-energy basis: " + crossref_match("crossref-electron-cs-20261031/faure-2004-jpb.json"),
        "Table 10 (frozen dataset) basis: " + crossref_match("crossref-electron-cs-20261031/harb-2001-jcp.json"),
        "OH(A) emission: " + crossref_match("crossref-electron-cs-20261031/beenakker-1974-chemphys.json")],
})

lit_by_id = {s["source_id"]: s for s in ml["sources"]}
sw_by_id = {s["source_id"]: s for s in ms["sources"]}
b = lit_by_id["bolorizadeh-rudd1986-pra"]
S.append({
    "source_id": "bolorizadeh-rudd1986-pra",
    "citation": ("M. A. Bolorizadeh and M. E. Rudd, Angular and energy dependence of cross sections for ejection of electrons "
                 "from water vapor. I. 50-2000-eV electron impact, Phys. Rev. A 33, 882-887 (1986)."),
    "doi": "10.1103/PhysRevA.33.882",
    "crossref_verified": b["crossref_verified"],
    "url": read("bolorizadeh-rudd1986-pra", "source-url.txt"),
    "contents_summary": ("Author-institution open copy (University of Nebraska-Lincoln DigitalCommons, via Internet Archive "
                         "snapshot). Table III (printed p. 887): total ionization cross sections 50-2000 eV from integration "
                         "of measured ejected-electron SDCS, 15% uncertainty (printed p. 883); the paper reports them 15-35% "
                         "above Schutten et al. 1966."),
})
S.append({
    "source_id": "crossref-electron-lit-20261031",
    "citation": "Crossref REST API work records and bibliographic searches for the ionization literature (tasks 1b, 1d, 1e).",
    "doi": None,
    "crossref_verified": {"value": True, "matched": "The records are the Crossref responses; per-citation matches below.",
                          "per_citation": ml.get("crossref_verification")},
    "url": "https://api.crossref.org/works/<doi> and https://api.crossref.org/works?query.bibliographic=...",
    "contents_summary": ("Records for Itikawa & Mason 2005, Rao/Iga/Srivastava 1995, Straub et al. 1998, Lindsay & Mangan "
                         "(Landolt-Bornstein I/17C), Djuric/Cadez/Kurepa 1988, Schutten et al. 1966, Rieke & Prepejchal 1972, "
                         "Bolorizadeh & Rudd 1986, Mark & Egger 1976, Orient & Srivastava 1987, Hanel et al. 2002, Montenegro "
                         "et al. 2007, Scully et al. 2006, King & Price 2008, Song et al. 2021, Munoz et al. 2007; plus search "
                         "responses (search-*.json)."),
})
S.append({
    "source_id": "crossref-electron-cs-20261031",
    "citation": "Crossref REST API work records for references behind the LXCat cross-section sets and the Song 2021 tables transcribed here.",
    "doi": None,
    "crossref_verified": {"value": True, "matched": "The records are the Crossref responses. Two DOIs initially guessed for Faure et al. 2004 and Harb et al. 2001 resolved to unrelated papers; those responses were discarded and the correct DOIs (10.1088/0953-4075/37/4/007, 10.1063/1.1397327) were found by bibliographic search and fetched."},
    "url": "https://api.crossref.org/works/<doi> and https://api.crossref.org/works?query.bibliographic=...",
    "contents_summary": ("Records: LXCat platform papers (Pancheshnyi 2012, Pitchford 2017, Carbone 2021), Hayashi 1987 "
                         "chapter, Morgan & Penetrante 1990, Matsui et al. 2016, Thorn et al. 2009, Kedzierski et al. 1998, "
                         "McConkey et al. 2008, Ralphs et al. 2013, Faure et al. 2004, Harb et al. 2001, Beenakker et al. 1974; "
                         "search-*.json are the bibliographic queries."),
})
sw = sw_by_id["lxcat-swarm-laplace-h2o-20261031"]
S.append({
    "source_id": "lxcat-swarm-laplace-h2o-20261031",
    "citation": ("LAPLACE database (measurements after 1975; digitized by S. Chowdhury, LAPLACE), www.lxcat.net/LAPLACE, "
                 "swarm/transport data for e in H2O. Entries cite R. W. Crompton, J. A. Rees, R. L. Jory, Aust. J. Phys. 18, "
                 "541 (1965) and 'Hasegawa, H., Date, H., Shimozuma, M., Hokkaido University' (probably J. Phys. D 40, 2495 "
                 "(2007)). Retrieved from www.lxcat.net (nl.lxcat.net mirror) on 31 October 2026."),
    "doi": None,
    "crossref_verified": sw["crossref_verified"],
    "url": "https://nl.lxcat.net/data/ (swarm/transport browse flow), see request.txt",
    "contents_summary": sw.get("contents_summary"),
    "lxcat_platform_references": LXCAT_PLATFORM,
})
cs = sw_by_id["crossref-swarm-20261031"]
S.append({
    "source_id": "crossref-swarm-20261031",
    "citation": "Crossref REST API work records and searches for water-vapour swarm literature.",
    "doi": None,
    "crossref_verified": cs["crossref_verified"],
    "url": "https://api.crossref.org/works/<doi>",
    "contents_summary": cs.get("contents_summary"),
})

# per-source file hashes, retrieval time, tables
tables_by_source = {}
for p in sorted(glob.glob(os.path.join(EXT, "*.json"))):
    try:
        t = json.load(open(p, encoding="utf-8"))
    except Exception:
        continue
    tables_by_source.setdefault(t.get("source_id"), []).append((p, t))
for s in S:
    d = os.path.join(SRC, s["source_id"])
    s["directory"] = os.path.relpath(d, ROOT)
    rt = read(s["source_id"], "retrieved_utc.txt") or read(s["source_id"], "created_utc.txt")
    s["retrieved_utc"] = rt.split() if rt and len(rt.split()) > 1 else rt
    files = []
    for f in sorted(os.listdir(d)):
        fp = os.path.join(d, f)
        if os.path.isfile(fp):
            files.append({"path": os.path.relpath(fp, ROOT), "bytes": os.path.getsize(fp), "sha256": sha(fp)})
    s["raw_files"] = files
    tl = []
    for p, t in tables_by_source.get(s["source_id"], []):
        xs = [v for v in t.get("x", []) if isinstance(v, (int, float))]
        tl.append({
            "table_id": t["table_id"],
            "file": os.path.relpath(p, ROOT),
            "file_sha256": sha(p),
            "tag": t.get("tag"),
            "quantity": t.get("quantity"),
            "process": t.get("process"),
            "x_name": t.get("x_name"),
            "x_range": [min(xs), max(xs)] if xs else None,
            "n_points": len(t.get("x", [])),
            "threshold_eV": t.get("threshold_eV"),
            "locator": t.get("locator"),
        })
    s["tables"] = tl

missing = [
    {"item": "LXCat cross-section sets for H2O named IST-Lisbon, Budapest, Biagi or any database other than Itikawa, Hayashi, Morgan, MuroranIT, Phelps, QUANTEMOL, TRINITI",
     "reason": "On 2026-10-31 the legacy LXCat browse flow with all 25 cross-section databases selected listed e/H2O processes only in these 7 databases. IST-Lisbon has no H2O set; Budapest is a swarm database (no H2O entries there either, see swarm items)."},
    {"item": "New LXCat (LXCat 3) API or download endpoints",
     "reason": "www.lxcat.net and lxcat.net redirect to the legacy nl.lxcat.net site; no API endpoint found. All data were taken from the legacy site."},
    {"item": "Song 2021 tabulated cross sections for a~3B1, 3A2, 1A2, b~3A1, B~1A1 excitation, the six-state sum, total dissociation into neutrals, OH(A) emission, H(n) and O emission",
     "reason": "Song 2021 gives these only in figures (Figs. 14-16) or as single numbers in the text (recorded in song2021-text-numeric-statements). Only A1B1 (Table 9, frozen dataset), OH (Table 10, frozen dataset) and O(1S0) (Table 11, transcribed here) are tabulated. Tabulated electronic-state cross sections are available only as model-derived LXCat MuroranIT blocks; measured emission cross sections are available as LXCat Itikawa blocks. Figures were not digitized."},
    {"item": "Song 2021 counting ionization cross section",
     "reason": "Defined in Sec. 2 (printed p. 4) but not tabulated; Table 12 'Total' is described as gross total ionization (Fig. 17 caption, printed p. 18)."},
    {"item": "Song 2021 tabulated Djuric et al. 1988 and Schutten et al. 1966 total ionization",
     "reason": "Shown only in Fig. 17; no table. Figure not digitized."},
    {"item": "Elastic DCS beyond 130 deg and outside 2-100 eV from Song 2021",
     "reason": "Table 2 covers 2-100 eV and 10-130 deg only (Matsui et al. 2016)."},
] + [dict(m, found_by="literature search (tasks 1b/1d/1e)") for m in ml.get("missing", [])] \
  + [dict(m, found_by="swarm search (task 2)") for m in ms.get("missing", [])]

findings = {
    "H2O_plus_12p621_to_13p5_eV": {
        "measured_values_strictly_inside_window": [],
        "statement": ("No measured H2O+ (or total ionization) cross section strictly between 12.621 and 13.5 eV was obtained. "
                      "The candidate primary measurements (Rao et al. 1995, Straub et al. 1998, Djuric et al. 1988, Mark & Egger "
                      "1976, Lindsay & Mangan) could not be retrieved. The lowest measurement-based tabulated value is at the "
                      "window edge: 13.5 eV, 0.025 x 1e-16 cm2 = 2.5e-22 m2 (Song 2021 Table 12, printed p. 20, H2O+ and Total "
                      "columns; identical in LXCat Itikawa blocks 22 and 23, preceded there by a (13.5 eV, 0) point)."),
        "lxcat_points_inside_window_not_measured": {
            "lxcat-muroranit-h2o-36-ionization-h2op (model_derived, swarm-validated; threshold 12.65 eV)": [
                ["1.265000e+1", "0.000000e+0"], ["1.267700e+1", "2.834300e-24"], ["1.273500e+1", "9.601400e-24"],
                ["1.285300e+1", "2.605800e-23"], ["1.303200e+1", "6.130600e-23"], ["1.327400e+1", "1.377600e-22"],
                ["1.350000e+1", "2.500000e-22"]],
            "lxcat-morgan-h2o-10-ionization-h2op (model_derived; threshold 13 eV)": [
                ["1.300000e+1", "0.000000e+0"], ["1.323500e+1", "1.000000e-23"]],
            "threshold-only zero points": "Phelps total ionization (12.6 eV, 0) then (14.5 eV, 7.3e-22); TRINITI (12.61 eV, 0) then (20 eV, 5.0e-21); Hayashi, QUANTEMOL and Itikawa have no points below 13.5 eV.",
            "remark": "MuroranIT reaches 2.5e-22 m2 at 13.5 eV, the same value as the Song/Lindsay-Mangan recommended point there."},
        "reference_energies_from_Song2021_text": "experimental ionization threshold 12.61 eV; Quantemol vertical threshold 13.67 eV (Sec. 10, printed p. 18).",
    },
    "total_ionization_above_1keV": {
        "Bolorizadeh_Rudd_1986_TableIII_p887_measured": {"energy_eV": [1000, 1200, 1500, 2000], "printed_1e-20_m2": [0.78, 0.70, 0.58, 0.44], "uncertainty": "15%", "table": "lit-bolorizadeh-rudd1986-total-ionization"},
        "Munoz_2007_TableI_printed_p4_measured_frozen_dataset": {"energy_eV": [1000, 1500, 2000, 3000, 4000, 5000], "printed_a0^2": [2.54, 1.83, 1.45, 1.05, 0.845, 0.702], "a0^2_to_m2": 2.800285205390781e-21, "channel": "H2O:munoz2007:measured_ionization in h2o-gas-20261029.json (not re-transcribed)"},
        "Song2021_Table12_last_point": "1000 eV total 0.705 x 1e-16 cm2 (frozen dataset H2O:ion_yield:total)",
        "calculated_only_above_1keV": "LXCat QUANTEMOL BEB blocks extend to 4993.7 eV (calculated, not measured).",
        "not_obtained": "Schutten et al. 1966 (0.1-20 keV) values; Rieke & Prepejchal 1972 Bethe parameters M^2 and C for H2O. No measured total ionization above 5 keV was obtained, so nothing measured covers the 10 keV seam.",
        "Phelps_Schutten_block": "LXCat Phelps total ionization ('from Shutten et al.') ends at 200 eV.",
    },
    "itikawa_lxcat_vs_song_table12": ml.get("itikawa_lxcat_comparison"),
    "song2021_table11_unit_conflict": ("Table 11 header prints 10^-16 cm2 but the Sec. 9 text (printed p. 16) gives a maximum of "
                                       "1.5 x 10^-18 cm2 and the LXCat Itikawa O(1S) block has the same mantissas in units of "
                                       "1e-18 cm2; the extracted table uses 1e-22 m2 per printed unit and keeps the literal "
                                       "reading in a separate field."),
    "song2021_double_ionization_statement": ("Sec. 10 (printed p. 19) prints 0.7 x 10^-16 cm2 each for (H+ + OH+) and (H+ + O+ + H) "
                                             "at 200 eV, inconsistent with the same section's ~7% double-ionization fraction "
                                             "and with Table 12; recorded verbatim and flagged."),
    "song2021_table10_vs_lxcat_itikawa_OH": ("Song Table 10 (frozen dataset) equals the LXCat Itikawa OH(X) block (Harb et al. "
                                             "2001) multiplied by 4.6/6.7 = 0.687 at all 11 energies to within rounding, "
                                             "consistent with the H- renormalization stated in Sec. 9 (printed p. 15)."),
    "song2021_table12_total_definition": ("Table 12 'Total' equals the sum of the six partial-ion columns within rounding and is "
                                          "called 'gross total ionization' in the Fig. 17 caption; no counting cross section is "
                                          "tabulated."),
    "swarm_data": ("Only 3 measured pure-H2O swarm sets were obtained, all LXCat LAPLACE transcriptions: D_T/mu 64.0-188.4 Td "
                   "(Crompton et al. 1965), muN 59.7-800.0 Td and (alpha-eta)/N 152.5-2999.2 Td (Hasegawa et al.). No separate "
                   "alpha/N, eta/N, D_L N or measured W tables were obtained."),
}

out = {
    "schema": "electron-sources-index-v1",
    "description": ("Index of electron + H2O collision and water-vapour swarm sources retrieved for this component, with the "
                    "extracted tables derived from each. Tags: measured | recommended | calculated | model_derived. "
                    "Values in data/extracted are verbatim transcriptions; see each table's extraction_method and notes."),
    "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "frozen_dataset_reference": {
        "path": "/home/treefor/projects/runtimes/h2o-collision-data-20261029/data/collisions/h2o-gas-20261029.json",
        "note": "Song 2021 Tables 3, 4, 5, 6/7, 8/9, 10, 12, 13 and Munoz 2007 Tables I/II are transcribed there and are not repeated here."},
    "tag_policy": {
        "measured": "direct experimental data reproduced without model adjustment",
        "recommended": "evaluated/recommended values from a review (Itikawa database, Song 2021)",
        "calculated": "ab initio or BEB calculation (QUANTEMOL)",
        "model_derived": "cross-section sets adjusted to reproduce swarm data or assembled for Boltzmann solvers (Hayashi, Morgan, MuroranIT, Phelps, TRINITI)"},
    "sources": S,
    "findings": findings,
    "missing": missing,
    "tools": ["tools/research/lxcat_browse_fetch.sh", "tools/research/unzip_safe.py", "tools/research/parse_lxcat_cs.py",
              "tools/research/verify_lxcat_cs.py", "tools/research/build_song2021_tables.py",
              "tools/research/verify_song2021_tables.py", "tools/research/song_rotated_table_words.py",
              "tools/research/crossref_record_summary.py", "tools/research/lxcat_swarm_fetch.sh",
              "tools/research/parse_lxcat_swarm.py", "tools/research/compare_itikawa_ionization_song.py",
              "tools/research/lxcat_window_scan.py", "tools/research/build_electron_sources_index.py"],
}
dest = os.path.join(SRC, "ELECTRON-SOURCES.json")
with open(dest, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
    fh.write("\n")
print(dest, "sources:", len(S), "tables:", sum(len(s["tables"]) for s in S), "missing:", len(missing))
for s in S:
    print(f"  {s['source_id']}: files={len(s['raw_files'])} tables={len(s['tables'])} retrieved={s['retrieved_utc']} crossref={s['crossref_verified'].get('value')}")
