"""Bibliographic query against Crossref; prints top matches (no files written unless output path given).

Usage: python3 -I crossref_query.py "<bibliographic string>" [rows] [output_json_path]
"""
import json
import sys
import urllib.parse
import urllib.request
import urllib.error
import time


def _get(req):
    delay = 3.0
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < 5:
                time.sleep(delay)
                delay *= 2
                continue
            raise


def main():
    q = sys.argv[1]
    rows = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode({"query.bibliographic": q, "rows": rows})
    req = urllib.request.Request(url, headers={"User-Agent": "pic-water-chemistry-source-check/1.0"})
    raw = _get(req)
    if len(sys.argv) > 3:
        with open(sys.argv[3], "wb") as f:
            f.write(raw)
    for it in json.loads(raw)["message"]["items"]:
        a = it.get("author", [])
        first = a[0].get("family", "") if a else ""
        year = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        print(it.get("DOI"), "|", first, "|", year, "|", (it.get("container-title") or [""])[0], it.get("volume"), it.get("page"), "|", (it.get("title") or [""])[0][:120])


if __name__ == "__main__":
    main()
