"""Extract the text of a PDF to a text file (page markers included). Usage: python -I pdf_text.py PDF OUT_TXT"""
import sys, pymupdf
with pymupdf.open(sys.argv[1]) as d, open(sys.argv[2], "w") as o:
    for i, p in enumerate(d):
        o.write(f"\n=== page {i + 1}\n"); o.write(p.get_text())
