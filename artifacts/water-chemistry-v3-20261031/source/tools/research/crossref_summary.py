"""Print a one-line bibliographic summary of Crossref work JSON files.

Usage: python3 -I crossref_summary.py <file.json> [...]
Handles both single-work responses (message is a work) and query responses
(message.items is a list of works).
"""
import json
import sys


def summarize(w):
    authors = "; ".join(
        f"{a.get('family', '')}, {a.get('given', '')}".strip(", ") for a in w.get("author", [])
    )
    year = None
    for k in ("published-print", "issued", "published-online"):
        dp = w.get(k, {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            year = dp[0][0]
            break
    return {
        "doi": w.get("DOI"),
        "title": (w.get("title") or [""])[0],
        "journal": (w.get("container-title") or [""])[0],
        "volume": w.get("volume"),
        "issue": w.get("issue"),
        "page": w.get("page") or w.get("article-number"),
        "year": year,
        "authors": authors,
    }


def main():
    for p in sys.argv[1:]:
        with open(p, encoding="utf-8") as fh:
            d = json.load(fh)
        msg = d.get("message", {})
        works = msg.get("items") if isinstance(msg.get("items"), list) else [msg]
        for w in works:
            print(p.rsplit("/", 1)[-1], json.dumps(summarize(w), ensure_ascii=False))


if __name__ == "__main__":
    main()
