from decimal import Decimal
import hashlib
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from xml.etree import ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references/PSSTab264csupp_e-OH-data.xlsx')
spec=importlib.util.spec_from_file_location('oh_primary',ROOT/'src/oh_primary.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)

class Controls(unittest.TestCase):
    def test_actual_pinned_channels(self):
        before=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        c=m.inspect(SOURCE)
        self.assertEqual(len(c),10)
        for name in ('elastic','excitation_C','excitation_D'):
            self.assertEqual(len(c[name].rows),1451)
            self.assertEqual(c[name].rows[0][0],Decimal('0.1361'))
            self.assertEqual(c[name].rows[-1][0],Decimal(10))
            self.assertEqual(c[name].value_unit,'bohr^2')
        self.assertEqual(len(c['dissociation'].rows),2186)
        self.assertEqual(c['dissociation'].rows[-1][0],Decimal(15))
        self.assertEqual(c['DEA_B'].value_unit,'Angstrom^2')
        self.assertIsInstance(c['excitation_B'],m.RefusedChannel)
        self.assertIn('B546',c['excitation_B'].reason)
        self.assertFalse(c['excitation_B'].summary()['missing_values_imputed'])
        with self.assertRaises(m.Refusal):m.load_channel(SOURCE,'excitation_B')
        self.assertEqual(m.load_channel(SOURCE,'elastic'),c['elastic'])
        with self.assertRaises(m.Refusal):m.load_channel(SOURCE,'absent')
        for channel in c.values():self.assertIs(channel.summary()['event_admitted'],False)
        self.assertEqual(before,hashlib.sha256(SOURCE.read_bytes()).hexdigest())

    def test_source_tamper_and_size_refuse(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.xlsx'
            for raw in (b'PKinvalid',SOURCE.read_bytes()+b'\n',b'x'*(2**20+1)):
                p.write_bytes(raw)
                with self.assertRaises(m.Refusal):m.inspect(p)

    def test_manufactured_channel_controls(self):
        base={'A1':Decimal(1),'B1':Decimal(0),'A2':Decimal(2),'B2':Decimal(1)}
        def read(c):return m._channel(c,1,'B','synthetic','eV','m^2')
        self.assertEqual(read(base).summary()['first_positive_abscissa'],'2')
        for update in ({'A2':Decimal(1)},{'A2':Decimal(0)},{'B2':Decimal(-1)},
                       {'B2':'1'},{'A2':'2'},{'A2':None},{'A2':Decimal('NaN')},{'B2':Decimal('Infinity')}):
            with self.subTest(update=update),self.assertRaises(m.Refusal):read(base|update)
        for base2 in ({'A1':Decimal(1),'B1':Decimal(0)},
                      {'A1':Decimal(1),'B1':Decimal(0),'A3':Decimal(3),'B3':Decimal(1)},
                      {'A1':Decimal(1),'B1':Decimal(0),'A2':Decimal(2)},
                      base|{'B3':Decimal(2)}):
            with self.subTest(base=base2),self.assertRaises(m.Refusal):read(base2)

    def test_xml_operand_guards(self):
        def cells(body):
            return m._cells(ET.fromstring('<worksheet xmlns="'+m.N['m']+'"><sheetData><row r="1">'+body+'</row></sheetData></worksheet>'),['eV'])
        self.assertEqual(cells('<c r="A1"><v>1.5</v></c><c r="B1" t="s"><v>0</v></c>'),{'A1':Decimal('1.5'),'B1':'eV'})
        bad=['<c r="A1"><v>NaN</v></c>','<c r="A1"><v>Infinity</v></c>',
             '<c r="A1"><v>1e1000</v></c>','<c r="A1"><v>1e-1000</v></c>',
             '<c r="A1" t="s"><v>1</v></c>','<c r="A1" t="e"><v>1</v></c>',
             '<c r="A1"><f>1+1</f><v>2</v></c>',
             '<c r="A1"><v>1</v></c><c r="A1"><v>2</v></c>',
             '<c r="../bad"><v>1</v></c>']
        for body in bad:
            with self.subTest(body=body),self.assertRaises(m.Refusal):cells(body)

if __name__=='__main__':unittest.main(verbosity=2)
