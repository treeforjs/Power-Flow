"""Bounded primary Table I reader. Values are per H atom, in eV.

E. Sartori, JINST 19 C03026 (2024), DOI10.1088/1748-0221/19/03/C03026,
CC BY 4.0. The source is a Franck-Condon/delta calculation, not a measured
fragment spectrum. Its missing probability endpoints are never invented.
"""
from decimal import Decimal
from pathlib import Path
import hashlib
import re

SHA256='3013c05f2c0913cfa6a784d817bb7b87733af67b395e8cc273da2eaa48abff91'
NUMBER=re.compile(r'(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\Z')
def parse(raw):
    if not isinstance(raw,bytes) or len(raw)>65536:raise ValueError('bounded bytes required')
    lines=raw.decode('ascii').splitlines()
    if len(lines)<3:raise ValueError('missing headers/data')
    if lines[0]!='e,H_2-->e,H_2(b3S_u+)-->e,H,H':raise ValueError('reaction/state identity')
    if lines[1].split()!=['R[0:1]']+['v='+str(i) for i in range(15)]:raise ValueError('quantile/vibrational columns')
    rows=[]
    for line in lines[2:]:
        if not line.strip():continue
        words=line.split()
        if len(words)!=16 or any(len(w)>48 or not NUMBER.fullmatch(w) for w in words):raise ValueError('numeric columns')
        vals=tuple(Decimal(w) for w in words)
        if any(not v.is_finite() for v in vals) or not 0<vals[0]<1 or any(v<=0 for v in vals[1:]):raise ValueError('finite positive operands')
        if rows and (vals[0]<=rows[-1][0] or any(a<=b for a,b in zip(vals[1:],rows[-1][1:]))):raise ValueError('monotonic quantiles and energies')
        rows.append(vals)
        if len(rows)>4096:raise ValueError('row limit')
    if len(rows)!=100 or rows[0][0]!=Decimal('.001') or rows[-1][0]!=Decimal('.999'):raise ValueError('published support')
    return tuple((r[0],r[1]) for r in rows) # explicit v=0 only
def load(path):
    with Path(path).open('rb') as f:raw=f.read(65537)
    if hashlib.sha256(raw).hexdigest()!=SHA256:raise ValueError('primary bytes changed')
    return parse(raw)
