"""Primary OH BEB total single-ionization model; no product-branch admission."""
from dataclasses import dataclass
from decimal import Decimal as D, Context, localcontext, ROUND_HALF_EVEN, InvalidOperation, DivisionByZero, Overflow
import hashlib
import json
from pathlib import Path
import re

CTX=Context(prec=60,rounding=ROUND_HALF_EVEN,Emin=-999999,Emax=999999,flags=[],traps=[InvalidOperation,DivisionByZero,Overflow])
PI=D('3.141592653589793238462643383279502884197169399375105820974944592307816')
RYD=D('13.605') # Exact source notebook constant, not silently replaced with another convention.
A0_M=D('5.29177211e-11') # Source notebook 5.29177211e-9 cm.
SOURCE_SHA='a9791272eb0b40a104756cbec71a9a090f4660744597d9ff2fe4d09d0b759259'
NUMBER=re.compile(r'(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[Ee][+-]?[0-9]+)?\Z')
class Refusal(ValueError):pass
def need(ok,why):
    if not ok:raise Refusal(why)
def bounded_read(path,maximum):
    with Path(path).open('rb') as f:raw=f.read(maximum+1)
    need(len(raw)<=maximum,'source read bound');return raw
@dataclass(frozen=True)
class Orbital:
    spin: str
    ordinal: int
    B: D
    U: D
    N: int
    Q: D
@dataclass(frozen=True)
class Incident:
    value: D
    unit: str='eV'
    frame: str='target-rest'
    target: str='OH(X2Pi)'
@dataclass(frozen=True)
class Result:
    sigma_m2: D
    orbital_sigma_m2: tuple
    classification: str='calculated total single ionization; branches unresolved'
    event_admitted: bool=False
    network_complete: bool=False
def _parse(raw):
    need(len(raw)<=4096,'orbital size bound')
    try:lines=raw.decode('ascii').splitlines()
    except UnicodeError as exc:raise Refusal('ASCII source required') from exc
    need(len(lines)==10 and lines[0]=='#Orbital\tB\tU\tN\tQ','source header/row count')
    result=[]
    for i,line in enumerate(lines[1:]):
        fields=line.split('\t');need(len(fields)==5,'orbital fields')
        need(fields[0].isdigit() and fields[3]=='1','ordinal/occupation')
        expected=i+1 if i<5 else i-4
        need(int(fields[0])==expected,'explicit alpha/beta order')
        need(all(len(x)<=32 and NUMBER.fullmatch(x) for x in fields[1:]),'numeric source token')
        B,U,N,Q=map(D,fields[1:]);need(B>0 and U>0 and B<D(10000) and U<D(10000) and N==1 and Q==1,'orbital domain')
        result.append(Orbital('alpha' if i<5 else 'beta',expected,B,U,1,Q))
    need(sum(x.N for x in result)==9,'neutral OH electron inventory')
    return tuple(result)
CANONICAL_ORBITALS_SHA='9d991332b5bf317f83e8cfcc678e54caea8b308819f02120356629e396e1915c'
def authenticate(orbitals):
    need(type(orbitals) is tuple and len(orbitals)==9,'source orbital tuple')
    rows=[]
    for o in orbitals:
        need(type(o) is Orbital and type(o.spin) is str and type(o.ordinal) is int and type(o.N) is int,'source orbital identity types')
        need(all(type(v) is D and v.is_finite() for v in (o.B,o.U,o.Q)),'source orbital Decimal operands')
        rows.append([o.spin,o.ordinal,str(o.B),str(o.U),o.N,str(o.Q)])
    need(hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()==CANONICAL_ORBITALS_SHA,'source orbital provenance')
@dataclass(frozen=True)
class Model:
    orbitals: tuple
    def __post_init__(self):authenticate(self.orbitals)
    def diagnostic_total(self,incident,*,mode=None):
        authenticate(self.orbitals)
        need(mode=='source-BEB-total-only','explicit diagnostic source model required')
        need(type(incident) is Incident and type(incident.value) is D,'typed Decimal incident required')
        need(incident.unit=='eV' and incident.frame=='target-rest' and incident.target=='OH(X2Pi)','incident metadata')
        E=incident.value;need(E.is_finite() and D(10)<=E<=D(5000),'outside paper-supported 10..5000 eV scope')
        with localcontext(CTX):
            terms=[]
            for orbital in self.orbitals:
                if E<=orbital.B:terms.append(D(0));continue
                t=E/orbital.B;u=orbital.U/orbital.B
                S=4*PI*A0_M*A0_M*orbital.N*(RYD/orbital.B)**2
                shape=t.ln()/2*(1-1/t**2)+1-1/t-t.ln()/(t+1)
                sigma=S/(t+u+1)*shape
                need(sigma.is_finite() and sigma>=0,'model evaluation')
                terms.append(sigma)
            return Result(sum(terms),tuple(terms))
    def require_event(self,*args,**kwargs):
        raise Refusal('total ionization does not specify product branches; threshold reference and energy sharing are not event-qualified')
def load(path):
    raw=bounded_read(path,4096);need(hashlib.sha256(raw).hexdigest()==SOURCE_SHA,'orbital source identity')
    return Model(_parse(raw))
