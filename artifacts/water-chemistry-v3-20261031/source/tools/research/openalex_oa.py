"""Look up open-access locations for DOIs via the OpenAlex API (prints only; writes nothing).

Usage: python3 -I openalex_oa.py <doi> [<doi> ...]
"""
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "pic-water-chemistry-source-check/1.0"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt < 3:
                time.sleep(3 * (attempt + 1))
                continue
            raise


def main():
    for doi in sys.argv[1:]:
        d = get("https://api.openalex.org/works/doi:" + urllib.parse.quote(doi, safe="/"))
        if d is None:
            print(doi, "| not in OpenAlex")
            continue
        oa = d.get("open_access", {})
        print(doi, "|", (d.get("title") or "")[:90], "| is_oa:", oa.get("is_oa"), "| status:", oa.get("oa_status"), "| oa_url:", oa.get("oa_url"))
        for loc in d.get("locations", []):
            if loc.get("pdf_url") or loc.get("is_oa"):
                print("    loc:", loc.get("landing_page_url"), "| pdf:", loc.get("pdf_url"), "| oa:", loc.get("is_oa"), "| version:", loc.get("version"))
        time.sleep(0.3)


if __name__ == "__main__":
    main()
