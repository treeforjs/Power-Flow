import dataclasses
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import hashlib
import io
from unittest.mock import patch
import importlib.util
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('h2_primary', ROOT/'src/h2_primary.py')
import sys
h = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = h
spec.loader.exec_module(h)
SOURCE = Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references/MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt')
CHEM = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101')
MANIFEST = CHEM/'data/tables/chemistry-v3.2-manifest.json'
SYNTHETIC = ('# Manufactured syntax fixture; no physical data\n'+h.REACTION+'\n'+h.UNITS+
             '\n# Threshold: 2.00000E+00 eV\n2 0\n3 1\n5 2\n').encode()


class ParserControls(unittest.TestCase):
    def test_valid_manufactured_parser_control(self):
        p = h._parse_source(SYNTHETIC)
        self.assertEqual(p.threshold, Decimal(2))
        self.assertEqual(p.rows, ((Decimal(2),Decimal(0)),(Decimal(3),Decimal(1)),(Decimal(5),Decimal(2))))

    def test_read_bound_enforced_on_stream_not_prior_stat(self):
        class GrowingStream(io.BytesIO):
            requested=[]
            def read(self, size=-1):
                self.requested.append(size)
                return super().read(size)
        stream=GrowingStream(b'x'*1000)
        with patch.object(Path,'open',return_value=stream), self.assertRaises(h.Refusal):
            h._read_pinned('synthetic-growing-file',hashlib.sha256(b'x'*1000).hexdigest(),8)
        self.assertEqual(stream.requested,[9])
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'tiny';p.write_bytes(b'12345678')
            self.assertEqual(h._read_pinned(p,hashlib.sha256(b'12345678').hexdigest(),8),b'12345678')

    def test_bad_structure_refuses(self):
        cases = {
            'wrong_initial_state': SYNTHETIC.replace(b'vi=0', b'vi=1'),
            'wrong_isotope': SYNTHETIC.replace(b'H2(', b'D2('),
            'wrong_final_state': SYNTHETIC.replace(b'b3Su', b'c3Pu'),
            'wrong_units': SYNTHETIC.replace(b'a_0^2', b'cm^2'),
            'wrong_energy_units': SYNTHETIC.replace(b'Energy (eV)', b'Energy (keV)'),
            'duplicate_reaction': SYNTHETIC + (h.REACTION+'\n').encode(),
            'duplicate_units': SYNTHETIC + (h.UNITS+'\n').encode(),
            'duplicate_threshold': SYNTHETIC + b'# Threshold: 2 eV\n',
            'missing_threshold': SYNTHETIC.replace(b'# Threshold: 2.00000E+00 eV\n', b''),
            'threshold_unit': SYNTHETIC.replace(b'2.00000E+00 eV', b'2.00000E+00 keV'),
            'zero_threshold': SYNTHETIC.replace(b'2.00000E+00', b'0'),
            'threshold_anchor_mismatch': SYNTHETIC.replace(b'2.00000E+00', b'2.1'),
            'positive_anchor': SYNTHETIC.replace(b'2 0\n', b'2 1\n'),
            'extra_column': SYNTHETIC.replace(b'3 1\n', b'3 1 9\n'),
            'negative_sigma': SYNTHETIC.replace(b'3 1\n', b'3 -1\n'),
            'negative_energy': SYNTHETIC.replace(b'3 1\n', b'-3 1\n'),
            'duplicate_energy': SYNTHETIC.replace(b'3 1\n', b'2 1\n'),
            'reordered_energy': SYNTHETIC.replace(b'3 1\n5 2\n', b'5 2\n3 1\n'),
            'nan_rate': SYNTHETIC.replace(b'3 1\n', b'3 NaN\n'),
            'inf_energy': SYNTHETIC.replace(b'3 1\n', b'Infinity 1\n'),
            'overflow_energy': SYNTHETIC.replace(b'3 1\n', b'1e999 1\n'),
            'underflow_rate': SYNTHETIC.replace(b'3 1\n', b'3 1e-999\n'),
            'no_positive_rate': SYNTHETIC.replace(b'3 1\n5 2\n', b'3 0\n5 0\n'),
            'one_row': SYNTHETIC.replace(b'3 1\n5 2\n', b''),
            'non_ascii': SYNTHETIC + b'\xff',
            'oversize': b'x'*65537,
            'empty': b''}
        for label, raw in cases.items():
            with self.subTest(case=label), self.assertRaises(h.Refusal):
                h._parse_source(raw)


class ActualPinnedControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.before = {p:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      [SOURCE, MANIFEST, CHEM/'data/tables/chemistry-v3.2.blob',
                       CHEM/'data/tables/chemistry-v3.2-coverage-policy.json',
                       CHEM/'data/tables/chemistry-v3.2-spec.json']}
        cls.component = h.load_primary(SOURCE, CHEM)

    @classmethod
    def tearDownClass(cls):
        for path, digest in cls.before.items():
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise AssertionError('source mutation: '+str(path))

    def test_actual_identity_and_scope(self):
        x = self.component.identity()
        self.assertEqual(x['row_count'],54)
        self.assertEqual([Decimal(v) for v in x['source_energy_range_eV']], [Decimal('4.47713'), Decimal(150)])
        for key in ['network_complete','production_qualified','physical_model_validated','native_adopted']:
            self.assertIs(x[key],False)
        self.assertEqual(x['source_sha256'], h.SOURCE_SHA256)

    def test_exact_table_knots_and_units(self):
        # Independently parse the two columns, bypassing module parser/query.
        raw = [tuple(map(Decimal,line.split())) for line in SOURCE.read_text().splitlines()
               if line.strip() and not line.startswith('#')]
        m = json.loads(MANIFEST.read_text())
        a0 = Decimal.from_float(m['constants']['bohr_m'])
        with localcontext() as ctx:
            ctx.prec=50
            for e,s in raw:
                result=self.component.diagnostic_rate(h.SourceEnergy(e))
                self.assertEqual(result.sigma_m2, s*a0*a0)
                self.assertIs(result.event_admitted,False)
            lo,hi=raw[8],raw[9]
            mid=(lo[0]+hi[0])/2
            # Affine interpolation should recover mean at the exact midpoint.
            self.assertEqual(self.component.diagnostic_rate(h.SourceEnergy(mid)).sigma_m2,
                             ((lo[1]+hi[1])/2)*a0*a0)

    def test_no_endpoint_hold_or_unsupported_zero(self):
        for value in ['0','4.4771299999999999','150.00000000000001','1e50']:
            with self.subTest(energy=value), self.assertRaises(h.Refusal):
                self.component.diagnostic_rate(h.SourceEnergy(Decimal(value)))
        self.assertEqual(self.component.diagnostic_rate(h.SourceEnergy(Decimal('4.47713'))).sigma_m2,0)

    def test_explicit_energy_operands(self):
        values = [float('nan'), 8.0, True, Decimal('NaN'), Decimal('Infinity'), Decimal('-1')]
        for value in values:
            with self.subTest(value=str(value)), self.assertRaises(h.Refusal):
                h.SourceEnergy(value)
        for kw in [{'unit':'keV'}, {'unit':'J'}, {'coordinate':'target-rest-frame'},
                   {'target':'H2(X1Sg,vi=1)'}, {'target':'D2(X1Sg,vi=0)'}, {'target':'H2'}]:
            with self.subTest(kw=kw), self.assertRaises(h.Refusal):
                h.SourceEnergy(Decimal(8),**kw)
        with self.assertRaises(h.Refusal):
            self.component.diagnostic_rate(Decimal(8))

    def test_rechecks_corrupt_operand(self):
        energy=h.SourceEnergy(Decimal(8))
        object.__setattr__(energy,'unit','J')
        with self.assertRaises(h.Refusal):self.component.diagnostic_rate(energy)

    def test_actual_threshold_incompatibility(self):
        x=self.component.threshold_assessment()
        self.assertGreater(x['registry_reference_cost_binary64_eV'],4.47713)
        self.assertGreater(Decimal(x['target_rest_incident_threshold_eV']),Decimal.from_float(x['registry_reference_cost_binary64_eV']))
        self.assertTrue(x['gap_to_registry_exceeds_one_printed_unit'])
        self.assertTrue(x['gap_to_target_rest_threshold_exceeds_one_printed_unit'])
        for flag in ['threshold_alignment_qualified','source_energy_recoil_convention_verified',
                     'event_kinematics_qualified','production_qualified']:
            self.assertFalse(x[flag])
        with localcontext() as ctx:
            ctx.prec=60
            # Invert the invariant threshold independently using the CM energy.
            o=self.component.operands
            m,M,Q=(Decimal.from_float(v) for v in (o.electron_mass_eV,o.target_mass_eV,o.reference_cost_eV))
            T=Decimal(x['target_rest_incident_threshold_eV'])
            kcm=(2*M*T)/(((m+M)**2+2*M*T).sqrt()+m+M)
            self.assertLess(abs(kcm-Q),Decimal('1e-50'))
            conflict=(Decimal('4.47713')+T)/2
            self.assertGreater(self.component.diagnostic_rate(h.SourceEnergy(conflict)).sigma_m2,0)
            with self.assertRaises(h.Refusal):self.component.event_rate(h.SourceEnergy(conflict))

    def test_no_event_admission_even_far_above_threshold(self):
        for value in ['4.47713','4.5','8','50','150']:
            with self.subTest(energy=value),self.assertRaises(h.Refusal):
                self.component.event_rate(h.SourceEnergy(Decimal(value)))

    def test_caller_decimal_context_cannot_change_results(self):
        original=self.component.threshold_assessment()
        value=self.component.diagnostic_rate(h.SourceEnergy(Decimal('10.25')))
        for rounding in (ROUND_FLOOR, ROUND_CEILING):
            with self.subTest(rounding=rounding), localcontext() as ctx:
                ctx.prec=6
                ctx.rounding=rounding
                ctx.Emin=-9;ctx.Emax=9
                ctx.capitals=0;ctx.clamp=1
                for signal in ctx.traps:ctx.traps[signal]=True
                self.assertEqual(original,self.component.threshold_assessment())
                self.assertEqual(value,self.component.diagnostic_rate(h.SourceEnergy(Decimal('10.25'))))

    def test_registry_operand_validation(self):
        valid=self.component.operands
        for key in ['electron_mass_eV','target_mass_eV','hydrogen_mass_eV','bohr_m']:
            for value in [0.0,-1.0,float('inf'),float('nan')]:
                with self.subTest(key=key,value=str(value)),self.assertRaises(h.Refusal):
                    dataclasses.replace(valid,**{key:value})
        for key in ['target_reference_eV','hydrogen_reference_eV']:
            with self.subTest(key=key),self.assertRaises(h.Refusal):
                dataclasses.replace(valid,**{key:float('nan')})
        with self.assertRaises(h.Refusal):dataclasses.replace(valid,hydrogen_reference_eV=0.0)
        with self.assertRaises(h.Refusal):dataclasses.replace(valid,electron_mass_eV=valid.target_mass_eV)

    def test_immutable_component(self):
        with self.assertRaises(h.Refusal):self.component._data=None
        with self.assertRaises(h.Refusal):h.PrimaryRateComponent(None,None,None)

    def test_raw_source_tamper_refuses(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'local-only-source.txt'
            p.write_bytes(SOURCE.read_bytes()+b'\n')
            with self.assertRaises(h.Refusal):h.load_primary(p,CHEM)

    def test_base_blob_manifest_spec_policy_tamper_refuses(self):
        with tempfile.TemporaryDirectory() as d:
            table=Path(d)/'data/tables';table.mkdir(parents=True)
            names=['chemistry-v3.2-manifest.json','chemistry-v3.2.blob',
                   'chemistry-v3.2-spec.json','chemistry-v3.2-coverage-policy.json']
            for name in names:shutil.copyfile(CHEM/'data/tables'/name,table/name)
            for name in names:
                with self.subTest(name=name):
                    p=table/name;raw=p.read_bytes();p.write_bytes(raw+b'\n')
                    with self.assertRaises(h.Refusal):h.load_primary(SOURCE,d)
                    p.write_bytes(raw)


if __name__ == '__main__':
    unittest.main(verbosity=2)
