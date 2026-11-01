from dataclasses import replace
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sp=importlib.util.spec_from_file_location('oh_beb',ROOT/'src/oh_beb.py');m=importlib.util.module_from_spec(sp);sys.modules['oh_beb']=m;sp.loader.exec_module(m)
SOURCE=Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references/alecs-06695e72/NIST_orbitals/mp2/OH.norb')
class Tests(unittest.TestCase):
    def setUp(self):self.model=m.load(SOURCE)
    def call(self,e):return self.model.diagnostic_total(m.Incident(D(e)),mode='source-BEB-total-only')
    def test_constructor_provenance(self):
        self.assertEqual(m.Model(self.model.orbitals),self.model)
        for rows in ((),self.model.orbitals[:-1],list(self.model.orbitals)):
            with self.subTest(rows=type(rows)),self.assertRaises(m.Refusal):m.Model(rows)
        for field,value in [('B',D('13.017')),('U',D('1')),('N',2),('Q',D('0')),('spin','beta'),('ordinal',2)]:
            changed=(replace(self.model.orbitals[0],**{field:value}),)+self.model.orbitals[1:]
            with self.subTest(field=field),self.assertRaises(m.Refusal):m.Model(changed)
        # Even an object mutated outside frozen dataclass semantics cannot retain the source label.
        broken=m.Model(self.model.orbitals);object.__setattr__(broken,'orbitals',())
        with self.assertRaises(m.Refusal):broken.diagnostic_total(m.Incident(D(70)),mode='source-BEB-total-only')
        old=ROOT.with_name('oh-beb-primary-v1-20261101')/'src/oh_beb.py'
        sp=importlib.util.spec_from_file_location('original_oh_beb',old);original=importlib.util.module_from_spec(sp);sys.modules[sp.name]=original;sp.loader.exec_module(original)
        self.assertEqual(original.Model(()).diagnostic_total(original.Incident(D(70)),mode='source-BEB-total-only').sigma_m2,0)
    def test_inventory(self):
        self.assertEqual(len(self.model.orbitals),9)
        self.assertEqual([(o.spin,o.ordinal) for o in self.model.orbitals],[('alpha',i) for i in range(1,6)]+[('beta',i) for i in range(1,5)])
        self.assertEqual(sum(o.N for o in self.model.orbitals),9)
    def test_boundary_and_onset(self):
        self.assertEqual(self.call('10').sigma_m2,0)
        self.assertEqual(self.call('12.728').sigma_m2,0)
        self.assertGreater(self.call('12.72800000001').sigma_m2,0)
        self.assertGreater(self.call('5000').sigma_m2,0)
        # Source-positive but below physical ground-ion cost remains explicitly non-admitted.
        self.assertGreater(self.call('12.9').sigma_m2,0)
        self.assertFalse(self.call('12.9').event_admitted)
    def test_operands(self):
        bad=[m.Incident(D('9.999')),m.Incident(D('5000.0001')),m.Incident(D('NaN')),m.Incident(D('Infinity')),
             m.Incident(D('-Infinity')),m.Incident(D(70),unit='keV'),m.Incident(D(70),frame='lab'),m.Incident(D(70),target='OD'),m.Incident(70.0),70]
        for x in bad:
            with self.subTest(x=x),self.assertRaises(m.Refusal):self.model.diagnostic_total(x,mode='source-BEB-total-only')
        with self.assertRaises(m.Refusal):self.model.diagnostic_total(m.Incident(D(70)))
        with self.assertRaises(m.Refusal):self.model.diagnostic_total(m.Incident(D(70)),mode='production')
    def test_event_refuses(self):
        for e in ('13.1','70','5000'):
            self.assertFalse(self.call(e).event_admitted)
            with self.assertRaises(m.Refusal):self.model.require_event(m.Incident(D(e)),branch='OH+',mode='production')
    def test_context(self):
        expected=self.call('70')
        for rounding in (ROUND_FLOOR,ROUND_CEILING):
            with localcontext() as ctx:
                ctx.prec=3;ctx.rounding=rounding;ctx.Emin=-4;ctx.Emax=4
                for k in ctx.traps:ctx.traps[k]=True
                self.assertEqual(self.call('70'),expected)
    def test_tamper_and_bound(self):
        raw=SOURCE.read_bytes()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x';p.write_bytes(raw.replace(b'12.728000',b'13.017000'))
            with self.assertRaises(m.Refusal):m.load(p)
            p.write_bytes(b'0'*4097)
            with self.assertRaises(m.Refusal):m.load(p)
    def test_malformed_parser(self):
        raw=SOURCE.read_bytes()
        cases=[raw+b'\n',raw.replace(b'#Orbital',b'#Other'),raw.replace(b'561.665298',b'NaN'),raw.replace(b'561.665298',b'-1'),raw.replace(b'\t1\t1\n',b'\t2\t1\n',1),raw.replace(b'\t1\t1\n',b'\t1\t0\n',1),raw.replace(b'5\t13.786',b'4\t13.786'),raw.replace(b'12.728000',b'1'*33)]
        for case in cases:
            with self.subTest(case=case),self.assertRaises(m.Refusal):m._parse(case)
if __name__=='__main__':unittest.main(verbosity=2)
