"""Extract rows of the ATcT v1.202 main table (HTML file argv[1]) whose formula cell matches one of argv[2] ('|'-separated, e.g. 'H2O  (g)')."""
import re, sys, json, html
text = open(sys.argv[1], encoding='utf-8', errors='replace').read()
want = set(sys.argv[2].split('|'))
rows = re.findall(r'<tr id="([^"]*)">(.*?)</tr>', text, flags=re.S)
found = {}
for rid, r in rows:
    def span(cls):
        m = re.search(r'<span class="%s">(.*?)</span>' % cls, r, flags=re.S)
        return html.unescape(re.sub(r'<[^>]+>', ' ', m.group(1))).split() if m else None
    fm = re.search(r'type="button">(.*?)</button>', r, flags=re.S)
    formula = html.unescape(fm.group(1)).strip() if fm else ''
    formula = re.sub(r'\s+', ' ', formula)
    if formula in want and formula not in found:
        name = span('Name')
        found[formula] = {'row_id': rid, 'name': ' '.join(name or []), 'formula': formula,
                          'DHf0_kJ_mol': ' '.join(span('DHf0') or []), 'DHf298_kJ_mol': ' '.join(span('DHf298') or []),
                          'uncertainty_kJ_mol': ' '.join(span('Uncert') or []), 'atct_id': ' '.join(span('ATcTID') or [])}
print(json.dumps(found, indent=1))
missing = sorted(want - set(found))
print('MISSING', missing, file=sys.stderr)
