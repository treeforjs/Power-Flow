"""Summarize Crossref work records and bibliographic search results.

Usage: python3 -I crossref_record_summary.py <json> [<json> ...]
Prints DOI, title, container, volume, page, issued date and first authors.
"""
import json
import sys


def line(it):
    au = "; ".join(f"{a.get('family', '')}" for a in (it.get("author") or [])[:6])
    return (f"{it.get('DOI')} | {(it.get('title') or [''])[0][:120]} | {(it.get('container-title') or [''])[0][:60]} | "
            f"vol {it.get('volume')} | p {it.get('page') or it.get('article-number')} | "
            f"{(it.get('issued', {}).get('date-parts') or [[None]])[0]} | {au}")


for p in sys.argv[1:]:
    d = json.load(open(p, encoding="utf-8"))
    m = d.get("message", {})
    print("==", p)
    if "items" in m:
        for it in m["items"]:
            print("  search:", line(it))
    else:
        print("  record:", line(m))
