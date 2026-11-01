"""Assemble the literature-source manifest (tasks 1b, 1d, 1e) with file hashes.

Usage: python3 -I build_lit_manifest.py <component_root> <out_json>
Static descriptive content is defined below; sha256 values and retrieval
times are read from disk at run time.
"""
import hashlib
import json
import os
import sys


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def cr_match(root, name):
    m = json.load(open(os.path.join(root, "data/sources/crossref-electron-lit-20261031", name)))["message"]
    return "Crossref %s: '%s', %s %s, %s, %s; authors %s" % (
        m.get("DOI"), (m.get("title") or [""])[0], (m.get("container-title") or [""])[0],
        m.get("volume"), m.get("page") or m.get("article-number"),
        (m.get("issued", {}).get("date-parts") or [[None]])[0],
        ", ".join(a.get("family", "") for a in m.get("author", [])))


def main():
    root, out = sys.argv[1], sys.argv[2]
    src = os.path.join(root, "data/sources")
    br_dir = os.path.join(src, "bolorizadeh-rudd1986-pra")
    cr_dir = os.path.join(src, "crossref-electron-lit-20261031")
    rel = lambda p: os.path.relpath(p, root)

    sources = []
    sources.append({
        "source_id": "bolorizadeh-rudd1986-pra",
        "citation": "M. A. Bolorizadeh and M. E. Rudd, 'Angular and energy dependence of cross sections for ejection of electrons from water vapor. I. 50-2000-eV electron impact', Phys. Rev. A 33, 882-887 (1986)",
        "doi": "10.1103/PhysRevA.33.882",
        "crossref_verified": {"value": True, "matched": cr_match(root, "bolorizadeh-rudd-1986.json")},
        "url": open(os.path.join(br_dir, "source-url.txt")).read().strip(),
        "retrieved_utc": open(os.path.join(br_dir, "retrieved_utc.txt")).read().strip(),
        "raw_files": [{"path": rel(os.path.join(br_dir, f)), "sha256": sha(os.path.join(br_dir, f))}
                      for f in sorted(os.listdir(br_dir))],
        "contents_summary": "Open author-institution copy (UNL DigitalCommons, M. Eugene Rudd Publications no. 10) of the published article, 7-page PDF (cover page + 6 article pages). Table I: angle-integrated singly differential cross sections; Table II: orbital constants; Table III: total ionization cross sections 50-2000 eV (absolute, 15% uncertainty), obtained by integrating SDCS.",
        "tables": [{"table_id": "lit-bolorizadeh-rudd1986-total-ionization", "tag": "measured",
                    "quantity": "cross_section_m2", "x_range": [50, 2000], "n_points": 14,
                    "locator": "Table III, printed page 887"}],
    })
    cr_files = sorted(f for f in os.listdir(cr_dir))
    sources.append({
        "source_id": "crossref-electron-lit-20261031",
        "citation": "Crossref REST API work records (api.crossref.org/works/<doi>) and bibliographic-search responses for the literature cited in tasks 1b, 1d, 1e",
        "doi": None,
        "crossref_verified": {"value": True, "matched": "records themselves"},
        "url": "https://api.crossref.org/works/",
        "retrieved_utc": open(os.path.join(cr_dir, "retrieved_utc.txt")).read().split(),
        "raw_files": [{"path": rel(os.path.join(cr_dir, f)), "sha256": sha(os.path.join(cr_dir, f))} for f in cr_files],
        "contents_summary": "Direct work records: " + "; ".join(f for f in cr_files if f.endswith('.json') and not f.startswith('search-')) + ". Files prefixed search- are query.bibliographic responses used to locate DOIs (search-fiegele-2000.json: no water-specific Fiegele 2000 paper found; nearest hit Hanel et al. 2002).",
        "tables": [],
    })

    verification = {}
    for key, fname in [
        ("Itikawa & Mason 2005", "itikawa-mason-2005.json"),
        ("Rao, Iga & Srivastava 1995", "rao-iga-srivastava-1995.json"),
        ("Straub, Lindsay, Smith & Stebbings 1998", "straub-1998.json"),
        ("Lindsay & Mangan, Landolt-Boernstein I/17C, ch. 5.1", "lindsay-mangan-2003-lb-I17C.json"),
        ("Djuric, Cadez & Kurepa 1988", "djuric-1988.json"),
        ("Schutten, de Heer, Moustafa, Boerboom & Kistemaker 1966", "schutten-1966.json"),
        ("Rieke & Prepejchal 1972", "rieke-prepejchal-1972.json"),
        ("Bolorizadeh & Rudd 1986", "bolorizadeh-rudd-1986.json"),
        ("Mark & Egger 1976", "mark-egger-1976.json"),
        ("Orient & Srivastava 1987", "orient-srivastava-1987.json"),
        ("Hanel et al. 2002 (near-threshold appearance energies, H2O/D2O)", "hanel-2002.json"),
        ("Montenegro, Scully, Wyer, Senthil & Shah 2007", "montenegro-2007-jesrp.json"),
        ("Scully, Wyer, Senthil, Shah & Montenegro 2006", "scully-2006-pra.json"),
        ("King & Price 2008", "king-price-2008.json"),
        ("Song et al. 2021", "song-2021.json"),
        ("Munoz et al. 2007", "munoz-2007.json"),
    ]:
        verification[key] = {"crossref_verified": True, "matched": cr_match(root, fname),
                             "record": rel(os.path.join(cr_dir, fname))}

    manifest = {
        "sources": sources,
        "crossref_verification": verification,
        "itikawa_lxcat_comparison": {
            "full_text_obtained": False,
            "status": "Only the LXCat Itikawa-database transcription is available (data/sources/lxcat-itikawa-h2o-cs-20261031). The Itikawa & Mason (2005) full text could not be retrieved, so the LXCat blocks cannot be mapped to the paper's table numbers or checked against the printed tables.",
            "cross_check_performed": "LXCat Itikawa ionization blocks compared point by point (common energies, no interpolation) with Song et al. 2021 Table 12 (printed p. 20, from Lindsay & Mangan, based on Straub et al. 1998) as held in the frozen dataset h2o-gas-20261029.json; Song Table 12 O2+ entry at 100 eV re-read from the Song 2021 text layer.",
            "cross_check_result": {
                "H2O+ partial (31 pts, 13.5-1000 eV)": "identical to Song Table 12 H2O+ column at all 31 energies; LXCat additionally lists (13.5 eV, 0) before (13.5 eV, 2.5e-22 m2)",
                "Total ionization (31 pts, 13.5-1000 eV)": "identical to Song Table 12 Total column; same leading (13.5 eV, 0) point",
                "H+ (16.9-1000 eV)": "identical at the 27 common energies; LXCat adds threshold point (16.9 eV, 0)",
                "OH+ (18.116-1000 eV)": "identical at the 27 common energies; LXCat adds threshold point (18.116 eV, 0)",
                "O+ (19-1000 eV)": "identical at the 25 common energies; LXCat adds threshold point (19 eV, 0)",
                "H2+ (20.7-1000 eV)": "identical at the 24 common energies; LXCat adds threshold point (20.7 eV, 0)",
                "O2+ (80-1000 eV)": "one difference: 100 eV LXCat 1.9e-24 m2 vs Song Table 12 printed 1.94e-4 x 1e-16 cm2 = 1.94e-24 m2 (-2.1%); other 15 common energies identical. Which value the Itikawa & Mason table prints cannot be checked without the paper.",
            },
            "interpretation": "The LXCat Itikawa ionization data reproduce the Lindsay & Mangan / Straub et al. 1998 recommended values (as also tabulated by Song 2021 Table 12), with a recommended H2O+ onset grid starting at 13.5 eV. Non-ionization LXCat Itikawa blocks were not cross-checked here.",
            "tool": "tools/research/compare_itikawa_ionization_song.py",
        },
        "onset_12p621_to_13p5_eV_H2O_plus": {
            "measured_values_strictly_inside_(12.621,13.5)_obtained": [],
            "boundary_value_at_13.5_eV": {
                "value_m2": 2.5e-22,
                "printed": "0.025 x 1e-16 cm2",
                "where": ["Song et al. 2021 Table 12, printed p. 20 (H2O+ and Total columns; from Lindsay & Mangan, based on Straub et al. 1998, partials from D2O measurements)",
                          "LXCat Itikawa database H2O+ and total-ionization blocks (13.5 eV, 2.5e-22 m2), preceded by (13.5 eV, 0)"],
                "tag": "recommended (measurement-based)",
            },
            "lxcat_points_inside_window_for_reference": {
                "note": "Parent extraction covers these; listed only so that the literature search result is not misread as 'no numbers anywhere'. These are swarm-adjusted or assembled sets, not direct measurements.",
                "MuroranIT (Kawaguchi et al. 2016), threshold 12.65 eV": [[12.65, 0.0], [12.677, 2.8343e-24], [12.735, 9.6014e-24], [12.853, 2.6058e-23], [13.032, 6.1306e-23], [13.274, 1.3776e-22], [13.5, 2.5e-22]],
                "Morgan, threshold 13 eV": [[13.0, 0.0], [13.235, 1.0e-23]],
                "TRINITI, threshold 12.61 eV": "no points inside; (12.61, 0) then (20, 5.0e-21)",
                "Phelps (total ionization from Schutten et al.), threshold 12.6 eV": "no points inside; (12.6, 0) then (14.5, 7.3e-22)",
                "Hayashi (threshold 13.76 eV), QUANTEMOL (13.67/13.72 eV)": "no points inside",
            },
            "candidate_primary_sources_not_obtained": [
                "Rao, Iga & Srivastava 1995 (measured from threshold to 1 keV; table content unknown)",
                "Straub et al. 1998 (threshold to 1000 eV; table content below 13.5 eV unknown)",
                "Djuric, Cadez & Kurepa 1988 (total ionization from threshold; Song 2021 Sec. 10 says up to 150 eV)",
                "Mark & Egger 1976 (H2O+ from threshold to 170 eV)",
                "Hanel et al. 2002 (near-threshold ion-yield scans; abstract reports appearance energies, not absolute cross sections)",
                "Lindsay & Mangan, Landolt-Boernstein I/17C (recommended values; reproduced by Song 2021 Table 12 from 13.5 eV)",
            ],
        },
        "total_ionization_above_1keV": {
            "obtained": {
                "Bolorizadeh & Rudd 1986, Table III, printed p. 887 (primary, measured, 15%)": {
                    "x_eV": [1000, 1200, 1500, 2000], "y_1e-20_m2": [0.78, 0.70, 0.58, 0.44],
                    "extracted_table": "data/extracted/lit-bolorizadeh-rudd1986-total-ionization.json"},
                "Munoz et al. 2007, Table I, printed p. 4 (primary, measured; already in frozen dataset channel H2O:munoz2007:measured_ionization, not re-transcribed)": {
                    "x_eV": [1000, 1500, 2000, 3000, 4000, 5000], "y_a0^2": [2.54, 1.83, 1.45, 1.05, 0.845, 0.702],
                    "a0^2_to_m2": 2.800285205390781e-21,
                    "notes": "Relative ion-current measurement 50-5000 eV normalized to N2 ionization at 1000 eV assumed (3.05 +/- 0.19) a0^2 (refs 22-26 of the paper). The Munoz text does not quote numerical values of Schutten et al. or Rieke & Prepejchal (neither is in its reference list)."},
            },
            "lxcat_phelps_schutten_block": "LXCat Phelps H2O ionization block ('Total ionization FROM SHUTTEN ET AL.') ends at 200 eV: points 12.6, 14.5, 18, 20, 30, 50, 100, 150, 200 eV; no points above 1 keV.",
            "not_obtained": ["Schutten et al. 1966 (0.1-20 keV gross ionization)", "Rieke & Prepejchal 1972 (0.1-2.7 MeV; Bethe parameters M^2 and C)", "Djuric et al. 1988 (does not extend above 1 keV per Song 2021 Sec. 10: threshold to 150 eV)", "Montenegro/Scully/Shah water papers (closed; electron-impact dissociative double ionization 45-1500 eV, not total ionization)"],
        },
        "missing": [
            {"item": "Itikawa & Mason, J. Phys. Chem. Ref. Data 34, 1 (2005) full text", "reason": "Semantic Scholar lists a free-to-read (bronze) publisher PDF at pubs.aip.org/aip/jpr/article-pdf/34/1/1/8183723/1_1_online.pdf, but the host returns a Cloudflare bot challenge (HTTP 403) to scripted requests; OpenAlex reports closed; Open University repository record oro.open.ac.uk/18224 returns HTTP 403 live and its archived copy states 'No digital document available to download for this item'; no NIST srd.nist.gov/jpcrdreprint copy (HTTP 404); no Internet Archive PDF capture. Only the LXCat Itikawa transcription is available. A manual browser download from AIP may succeed."},
            {"item": "Rao, Iga & Srivastava, J. Geophys. Res. 100(E12), 26421 (1995)", "reason": "agupubs.onlinelibrary.wiley.com /doi/pdf, /doi/pdfdirect and /doi/epdf all return HTTP 403 (Cloudflare challenge); OpenAlex and Semantic Scholar report closed; no Internet Archive PDF capture."},
            {"item": "Straub, Lindsay, Smith & Stebbings, J. Chem. Phys. 108, 109 (1998)", "reason": "Paywalled (OpenAlex and Semantic Scholar: closed); no archived open copy found."},
            {"item": "Lindsay & Mangan, Landolt-Boernstein I/17C ch. 5.1 'Ionization' (doi 10.1007/10874891_2)", "reason": "Springer paywalled; Crossref record verified. Values are available via Song 2021 Table 12 (13.5-1000 eV), which the parent handles."},
            {"item": "Djuric, Cadez & Kurepa, Int. J. Mass Spectrom. Ion Processes 83, R7 (1988)", "reason": "Elsevier paywalled (closed in OpenAlex and Semantic Scholar)."},
            {"item": "Schutten, de Heer, Moustafa, Boerboom & Kistemaker, J. Chem. Phys. 44, 3924 (1966)", "reason": "AIP paywalled (closed); no archived copy; the only numeric derivative found is the LXCat Phelps block, which stops at 200 eV. No authoritative compilation with a table locator quoting its keV values was obtained."},
            {"item": "Rieke & Prepejchal, Phys. Rev. A 6, 1507 (1972): H2O Bethe parameters M^2 and C", "reason": "APS paywalled (closed); no open compilation quoting the H2O M^2 and C with a table locator was found (Inokuti 1979 ICPEAC report, OSTI 5971494, cites the paper; a search of its scanned text layer found no tabulated H2O parameters; Kim, Santos & Parente PRA 62, 052710 (2000) is closed). Values not recorded."},
            {"item": "Mark & Egger, Int. J. Mass Spectrom. Ion Phys. 20, 89 (1976)", "reason": "Elsevier paywalled (closed)."},
            {"item": "Orient & Srivastava, J. Phys. B 20, 3923 (1987)", "reason": "Listed as free-to-read (bronze) at iopscience.iop.org/article/10.1088/0022-3700/20/15/036/pdf, but IOP returns a Radware bot-manager redirect to scripted requests; no archived copy."},
            {"item": "Hanel et al., J. Chem. Phys. 116, 2456 (2002)", "reason": "AIP paywalled (closed). Abstract indicates appearance-energy determinations from near-threshold ion-yield scans rather than absolute cross sections."},
            {"item": "Montenegro et al., J. Electron Spectrosc. Relat. Phenom. 155, 81 (2007); Scully et al., Phys. Rev. A 73, 040701 (2006)", "reason": "Closed. Abstract of the 2006 paper: electron-impact dissociative double-ionization cross sections 45-1500 eV (ion pairs H+ + OH+, H+ + O+), not total or partial single ionization above 1 keV. No Montenegro electron-impact partial ionization of water above 1 keV with tabulated values was found."},
            {"item": "Fiegele et al. 2000 water near-threshold data", "reason": "Crossref search found no water-specific Fiegele 2000 paper (the 2000 J. Phys. B paper covers CF4, CHF3, CH4, C3H8); the related water work is Hanel et al. 2002 (closed)."},
            {"item": "Bull et al. 2012 (water near-threshold)", "reason": "No matching water paper located in the time available; not pursued further."},
        ],
    }
    with open(out, "w") as f:
        json.dump(manifest, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
