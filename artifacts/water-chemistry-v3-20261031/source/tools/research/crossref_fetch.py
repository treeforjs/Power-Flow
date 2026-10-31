"""Fetch a Crossref work record by DOI and save the raw JSON.

Usage: python3 -I crossref_fetch.py <doi> <output_json_path>
Prints a short bibliographic summary of the saved record.
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
    doi, out = sys.argv[1], sys.argv[2]
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")
    req = urllib.request.Request(url, headers={"User-Agent": "pic-water-chemistry-source-check/1.0"})
    raw = _get(req)
    with open(out, "wb") as f:
        f.write(raw)
    msg = json.loads(raw)["message"]
    authors = "; ".join(
        (a.get("family", "") + ", " + a.get("given", "")) for a in msg.get("author", [])
    )
    year = None
    for key in ("published-print", "issued", "published-online"):
        if key in msg and msg[key].get("date-parts"):
            year = msg[key]["date-parts"][0][0]
            break
    print(json.dumps({
        "DOI": msg.get("DOI"),
        "title": (msg.get("title") or [""])[0],
        "authors": authors,
        "container": (msg.get("container-title") or [""])[0],
        "volume": msg.get("volume"),
        "issue": msg.get("issue"),
        "page": msg.get("page"),
        "year": year,
        "type": msg.get("type"),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
