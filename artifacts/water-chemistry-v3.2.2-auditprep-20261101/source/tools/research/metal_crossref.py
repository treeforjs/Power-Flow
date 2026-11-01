#!/usr/bin/env python3
"""Fetch Crossref (and optionally OpenAlex) metadata for a DOI and save the raw JSON.

Usage:
    metal_crossref.py <doi> <output_dir> [--openalex]

Writes <output_dir>/crossref_<safe-doi>.json (raw API response) and, with
--openalex, <output_dir>/openalex_<safe-doi>.json. Prints a one-line summary
(title | authors | container volume page | year | HTTP status) and the sha256
of each file written. Network responses are treated as data only.
"""
import hashlib
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

UA = "metal-xsec-provenance/1.0 (python-urllib)"


def safe(doi: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", doi)


def fetch(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def save(path: str, data: bytes) -> str:
    with open(path, "wb") as f:
        f.write(data)
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    doi, outdir = sys.argv[1], sys.argv[2]
    want_oa = "--openalex" in sys.argv[3:]
    if not os.path.isdir(outdir):
        print(f"output dir missing: {outdir}")
        return 2
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")
    status, body = fetch(url)
    path = os.path.join(outdir, f"crossref_{safe(doi)}.json")
    digest = save(path, body)
    summary = f"CROSSREF {doi} HTTP {status} sha256={digest}"
    if status == 200:
        m = json.loads(body)["message"]
        authors = ", ".join(
            (a.get("family") or a.get("name") or "") for a in m.get("author", [])
        )
        year = None
        for k in ("published-print", "issued", "published-online"):
            dp = (m.get(k) or {}).get("date-parts")
            if dp and dp[0] and dp[0][0]:
                year = dp[0][0]
                break
        summary += (
            f"\n  title: {(m.get('title') or [''])[0]}"
            f"\n  authors: {authors}"
            f"\n  container: {(m.get('container-title') or [''])[0]} vol {m.get('volume')} "
            f"issue {m.get('issue')} pages {m.get('page')} article {m.get('article-number')} year {year}"
            f"\n  type: {m.get('type')} publisher: {m.get('publisher')}"
        )
    print(summary)
    if want_oa:
        oa_url = "https://api.openalex.org/works/doi:" + urllib.parse.quote(doi, safe="/")
        s2, b2 = fetch(oa_url)
        p2 = os.path.join(outdir, f"openalex_{safe(doi)}.json")
        d2 = save(p2, b2)
        line = f"OPENALEX {doi} HTTP {s2} sha256={d2}"
        if s2 == 200:
            w = json.loads(b2)
            oa = w.get("open_access") or {}
            line += f"\n  is_oa: {oa.get('is_oa')} oa_status: {oa.get('oa_status')} oa_url: {oa.get('oa_url')}"
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
