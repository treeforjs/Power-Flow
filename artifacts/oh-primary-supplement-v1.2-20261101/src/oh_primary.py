"""Read primary OH workbook without claiming event/model qualification."""
from dataclasses import dataclass
from decimal import Decimal
import hashlib
import io
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

WORKBOOK_SHA256 = 'f3e96e924b3b6fde7583416f8ec483d708a52a982bfa80a07e61f0532c20b38e'
N = {'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
NUMBER = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\Z')


class Refusal(ValueError):
    pass


def need(condition, reason):
    if not condition:
        raise Refusal(reason)


def pinned_bytes(path):
    with Path(path).open('rb') as f:
        raw = f.read(2**20 + 1)
    need(len(raw) <= 2**20, 'workbook size bound')
    need(hashlib.sha256(raw).hexdigest() == WORKBOOK_SHA256, 'workbook source identity')
    return raw


def _xml(z, member):
    raw = z.read(member)
    need(b'<!DOCTYPE' not in raw and b'<!ENTITY' not in raw, 'DTD/entity unsupported')
    return ET.fromstring(raw)


def _cells(root, strings):
    result = {}
    for row in root.findall('m:sheetData/m:row', N):
        for cell in row:
            address = cell.get('r')
            need(address and re.fullmatch(r'[A-Z]{1,3}[1-9][0-9]*', address), 'invalid cell address')
            need(address not in result, 'duplicate cell')
            need(cell.find('m:f', N) is None, 'formulas not source values')
            value = cell.find('m:v', N)
            if value is None:
                need(cell.get('t') in (None, 'n'), 'unsupported non-value cell')
                continue
            token = value.text
            if cell.get('t') == 's':
                need(token is not None and token.isdigit() and int(token) < len(strings), 'invalid string index')
                result[address] = strings[int(token)]
            else:
                need(cell.get('t') in (None, 'n'), 'unsupported cell value type')
                need(token is not None and len(token) <= 80 and NUMBER.fullmatch(token), 'invalid numeric cell')
                number = Decimal(token)
                need(number.is_finite() and -999 <= number.adjusted() <= 999, 'numeric source bounds')
                result[address] = number
    return result


@dataclass(frozen=True)
class Channel:
    label: str
    energy_unit: str
    value_unit: str
    rows: tuple

    def summary(self):
        positive = [(e,v) for e,v in self.rows if v > 0]
        return {'label':self.label, 'row_count':len(self.rows),
                'abscissa_unit':self.energy_unit, 'value_unit_literal':self.value_unit,
                'abscissa_range':[str(self.rows[0][0]),str(self.rows[-1][0])],
                'first_positive_abscissa':str(positive[0][0]) if positive else None,
                'last_positive_abscissa':str(positive[-1][0]) if positive else None,
                'positive_rows':len(positive), 'zero_rows':len(self.rows)-len(positive),
                'below_grid':'unsupported', 'above_grid':'unsupported',
                'numeric_column_complete':True, 'event_admitted':False}


def _channel(cells, start, column, label, energy_unit, value_unit):
    rows=[]
    indices=sorted({int(re.search(r'[0-9]+$', a)[0]) for a in cells if re.fullmatch(r'A[1-9][0-9]*',a) and int(re.search(r'[0-9]+$', a)[0]) >= start})
    need(indices and indices == list(range(start,indices[-1]+1)), 'missing abscissa row')
    ordinate_indices={int(a[len(column):]) for a in cells if re.fullmatch(column+r'[1-9][0-9]*',a) and int(a[len(column):])>=start}
    need(not (ordinate_indices-set(indices)), 'orphan ordinate without abscissa')
    for index in indices:
        e,v=cells.get('A'+str(index)),cells.get(column+str(index))
        need(type(e) is Decimal and type(v) is Decimal, 'missing/non-numeric data pair at A'+str(index)+'/'+column+str(index))
        need(e.is_finite() and v.is_finite(), 'nonfinite channel operand')
        need(e > 0 and v >= 0, 'nonpositive abscissa/negative value')
        need(not rows or e > rows[-1][0], 'nonincreasing abscissa')
        rows.append((e,v))
    need(len(rows)>=2, 'insufficient rows')
    return Channel(label,energy_unit,value_unit,tuple(rows))


@dataclass(frozen=True)
class RefusedChannel:
    label: str
    reason: str

    def summary(self):
        return {'label':self.label, 'numeric_column_complete':False,
                'event_admitted':False, 'reason':self.reason,
                'missing_values_imputed':False}


def _inspect_channel(*args):
    try:
        return _channel(*args)
    except Refusal as exc:
        return RefusedChannel(args[3],str(exc))


def inspect(path):
    raw=pinned_bytes(path)
    with ZipFile(io.BytesIO(raw)) as z:
        entries=z.infolist();names=[x.filename for x in entries]
        need(len(entries)<=64 and len(set(names))==len(names), 'duplicate or excess ZIP members')
        need(sum(x.file_size for x in entries)<=8*2**20 and all(x.file_size<=2**22 for x in entries), 'inflated ZIP size bound')
        need(all(not n.startswith('/') and '..' not in Path(n).parts for n in names), 'unsafe member path')
        relations=_xml(z,'xl/_rels/workbook.xml.rels')
        need(all(r.get('TargetMode')!='External' for r in relations), 'external workbook relationships')
        targets={r.get('Id'):r.get('Target') for r in relations}
        wb=_xml(z,'xl/workbook.xml')
        sheets={s.get('name'):'xl/'+targets[s.get(REL)] for s in wb.find('m:sheets',N)}
        need(sheets=={'Elastic':'xl/worksheets/sheet1.xml', 'Excitation':'xl/worksheets/sheet2.xml',
                      'Rate':'xl/worksheets/sheet3.xml', 'Dissociation':'xl/worksheets/sheet4.xml',
                      'DEA':'xl/worksheets/sheet5.xml'}, 'workbook sheet identity')
        strings=[''.join(t.itertext()) for t in _xml(z,'xl/sharedStrings.xml')]
        data={name:_cells(_xml(z,member),strings) for name,member in sheets.items()}
    need(data['Elastic']['A2']=='Energy(eV)' and data['Elastic']['B2']=='Elastic Cross Section (bohr^2)', 'elastic units')
    need(data['Excitation']['A2']=='Energy(eV)' and data['Excitation']['B1']=='Excitation cross sections (bohr^2)', 'excitation units')
    need(data['Dissociation']['A1']=='Energy (eV)' and data['Dissociation']['B1']=='Dissociation Xsecs (bohr^2)', 'dissociation units')
    need(data['DEA']['A2']=='#E(eV)' and data['DEA']['B1']=='DEA cross srction (Angstrom^2)', 'DEA units')
    need(data['Rate']['A1']=='Temparature (K)' and data['Rate']['B1']=='Rate: X 2\\Pi --> A 2\\Sigma+(bohr^2)', 'rate literal units')
    out={}
    out['elastic']=_inspect_channel(data['Elastic'],3,'B','elastic integral','eV','bohr^2')
    for col,state in [('B','A 2\\Sigma+'),('C','A 4\\Sigma-'),('D','1 2\\Sigma-')]:
        label='X 2\\Pi --> '+state
        need(data['Excitation'][col+'2']==label, 'excitation state metadata')
        out['excitation_'+col]=_inspect_channel(data['Excitation'],3,col,label,'eV','bohr^2')
    out['dissociation']=_inspect_channel(data['Dissociation'],2,'B','dissociation; product state unspecified in workbook','eV','bohr^2')
    for col,label in [('B','Total'),('C','A 1Pi'),('D','a 3Pi'),('E','1 3Sigma+')]:
        need(data['DEA'][col+'2']==label, 'DEA resonance label')
        out['DEA_'+col]=_inspect_channel(data['DEA'],3,col,label,'eV','Angstrom^2')
    out['rate_unqualified_units']=_inspect_channel(data['Rate'],2,'B','published thermal rate (unit conflict)','K','bohr^2 as printed; rate coefficient dimension unresolved')
    return out


def load_channel(path, name):
    channels=inspect(path)
    need(name in channels, 'unknown source channel')
    channel=channels[name]
    need(type(channel) is Channel, 'source channel refused: '+channel.summary().get('reason',''))
    return channel
