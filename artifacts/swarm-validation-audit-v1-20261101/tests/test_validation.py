import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


V = module(ROOT / 'source/swarm_validation.py', 'validation')
M = module(ROOT / 'source/swarm_multiseed.py', 'patched_comparator')
M.R = str(ROOT / 'baseline')
M._validation = V
OLD = module(ROOT / 'baseline/tools/swarm_multiseed.py', 'old_comparator')
OLD.R = str(ROOT / 'baseline')
C = V.Contract(ROOT / 'source/swarm-contract.json', ROOT / 'baseline')
HISTORICAL = Path('/home/treefor/projects/.workspace/outputs/pic-water-chemistry-v3.1-20261031/qual-001/SWARM-VALIDATION.json')


def fixture(label, seed):
    c = C.data
    points = []
    for en in C.energies:
        row = {q: 1.0 for q in V.POINT_NUMBERS}
        row.update({q: 0 for q in V.POINT_COUNTS})
        row['E_over_N_Td'] = float(format(en, '.6g'))
        row['process_events'] = {name: 0 for name in C.processes[label]}
        for q in V.METRICS:
            measured = C.measured[q].get(en, 1)
            row[q] = measured * (1.0 + 0.1 + 0.001*(seed-c['seeds'][0]))
        points.append(row)
    return dict(blob_sha256=c['sets'][label]['blob_sha256'], N_m3=c['N_m3'],
                T_gas_K=c['T_gas_K'], electrons=c['electrons'], seed=seed, points=points)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='swarm-audit-', dir=ROOT / 'evidence')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        for label in C.data['sets']:
            for seed in C.data['seeds']:
                self.write(label, seed, fixture(label, seed))

    def path(self, label='v32', seed=20261101):
        return self.directory / f'{label}-{seed}.json'

    def write(self, label, seed, value):
        self.path(label, seed).write_text(json.dumps(value)+'\n')

    def invalid_run(self, change):
        run = fixture('v32', 20261101)
        change(run)
        with self.assertRaises((V.Invalid, ValueError, KeyError, TypeError)):
            V.validate_run(run, 'v32', 20261101, C)

    def test_complete_all_seeds_points_metrics(self):
        for label in C.data['sets']:
            summary = V.summarize(self.directory, label, C)
            self.assertEqual(len(summary['rows']), len(C.rows))
            self.assertTrue(all(row['n'] == 6 for row in summary['rows']))
            self.assertEqual(len(summary['input_files']), 6)
        self.assertTrue(V.inspect(self.directory, C)['complete'])

    def test_baseline_loses_six_digit_energy_labels(self):
        out = self.directory / 'old-summary.json'
        OLD.summary(str(self.directory), 'v32', str(out))
        result = json.loads(out.read_text())
        self.assertLess(len(result['rows']), len(C.rows))
        (ROOT/'evidence/original-rounding-failure.json').write_text(json.dumps(
            {'original_rows':len(result['rows']), 'required_rows':len(C.rows),
             'required_judged':C.judged, 'original_judged':sum(OLD.judged(r['quantity'],r['E_over_N_Td']) for r in result['rows'])}, indent=2)+'\n')

    def test_baseline_passes_empty_comparison(self):
        a=self.directory/'empty.json';a.write_text(json.dumps(dict(rows=[],blobs=[],seeds=[])))
        prefix=self.directory/'old-empty-result'
        self.assertEqual(OLD.compare(str(a), str(a), str(HISTORICAL), str(prefix)),0)
        result=json.loads(Path(str(prefix)+'.json').read_text())
        self.assertTrue(result['passed']);self.assertEqual(result['judged'],0)
        (ROOT/'evidence/original-empty-failure.json').write_text(json.dumps(result,indent=2)+'\n')
        with self.assertRaises(V.Invalid):
            V.validate_summary(dict(rows=[],blobs=[],seeds=[]), 'v32', C)

    def test_missing_and_extra_seed_files(self):
        self.path().unlink()
        with self.assertRaises(V.Invalid): V.summarize(self.directory,'v32',C)
        self.write('v32',20261101,fixture('v32',20261101))
        self.write('v32',20261109,fixture('v32',20261101))
        with self.assertRaises(V.Invalid): V.summarize(self.directory,'v32',C)

    def test_seed_blob_and_gas_configuration(self):
        for field,value in [('seed',20261102),('blob_sha256','0'*64),('N_m3',1),('T_gas_K',301),('electrons',1999)]:
            with self.subTest(field=field): self.invalid_run(lambda r:r.__setitem__(field,value))

    def test_missing_duplicate_unexpected_energy(self):
        self.invalid_run(lambda r:r['points'].pop())
        self.invalid_run(lambda r:r['points'].__setitem__(1,copy.deepcopy(r['points'][0])))
        self.invalid_run(lambda r:r['points'][0].__setitem__('E_over_N_Td',59.65391))

    def test_missing_invalid_metrics_and_species(self):
        self.invalid_run(lambda r:r['points'][0].pop('DT_over_mu_eV'))
        self.invalid_run(lambda r:r['points'][0].__setitem__('W_bulk_m_s',True))
        self.invalid_run(lambda r:r['points'][0].__setitem__('mean_energy_eV',None))
        self.invalid_run(lambda r:r['points'][0].__setitem__('process_events',{'e + H2':0}))
        self.invalid_run(lambda r:r['points'][0]['process_events'].__setitem__(next(iter(C.processes['v32'])),1))

    def test_nonfinite_and_duplicate_json_keys(self):
        for value in ['NaN','Infinity','-Infinity','1e999']:
            with self.subTest(value=value), self.assertRaises(V.Invalid): V.decode('{"x":'+value+'}')
        with self.assertRaises(V.Invalid): V.decode('{"seed":1,"seed":2}')
        for value in [float('nan'),float('inf'),float('-inf')]:
            with self.subTest(value=value): self.invalid_run(lambda r:r['points'][0].__setitem__('W_flux_m_s',value))

    def test_partial_writer_is_incomplete_and_not_repaired(self):
        before=b'{"blob_sha256":"unfinished'
        self.path().write_bytes(before)
        result=V.inspect(self.directory,C)
        self.assertFalse(result['complete'])
        self.assertEqual(self.path().read_bytes(),before)

    def test_summary_integrity_and_empty_rows(self):
        good=V.summarize(self.directory,'v32',C)
        for mutate in [lambda r:r.__setitem__('rows',[]),lambda r:r['rows'].pop(),
                       lambda r:r['rows'].append(copy.deepcopy(r['rows'][0])),
                       lambda r:r.__setitem__('seeds',[20261101]),lambda r:r.__setitem__('blobs',['0'*64]),
                       lambda r:r.__setitem__('contract_sha256','0'*64),lambda r:r.__setitem__('target','H2'),
                       lambda r:r['rows'][0].__setitem__('mean_rel',float('nan')),
                       lambda r:r['rows'][0].__setitem__('n',5),lambda r:r['rows'][0].__setitem__('sd_rel',100)]:
            bad=copy.deepcopy(good);mutate(bad)
            with self.subTest(mutate=mutate),self.assertRaises(V.Invalid): V.validate_summary(bad,'v32',C)

    def test_pinned_sources_and_contract(self):
        c=copy.deepcopy(C.data);c['sets']['v32']['spec_sha256']='0'*64
        p=self.directory/'contract.json';p.write_text(json.dumps(c))
        with self.assertRaises(V.Invalid): V.Contract(p,ROOT/'baseline')
        c=copy.deepcopy(C.data);c['criterion']['absolute_relative_floor']=0.1;p.write_text(json.dumps(c))
        with self.assertRaises(V.Invalid): V.Contract(p,ROOT/'baseline')

    def test_overflowing_statistics_refuse(self):
        with self.assertRaises(V.Invalid): V.statistics([1e308,-1e308,1e308,-1e308,1e308,-1e308])

    def test_adoption_patch_dry_run(self):
        result=subprocess.run(['git','apply','--check',str(ROOT/'adoption/swarm-completeness.patch')],
                              cwd=ROOT/'baseline',capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_failed_summary_preserves_existing_output(self):
        out=self.directory/'summary.json';out.write_bytes(b'preserved\n');self.path().unlink()
        with self.assertRaises(V.Invalid): M.summary(str(self.directory),'v32',str(out))
        self.assertEqual(out.read_bytes(),b'preserved\n')

    def test_same_scientific_tolerance_pass_and_fail(self):
        ref=V.summarize(self.directory,'v3',C)
        rp=self.directory/'summary-v3.json';V.atomic_json(rp,ref)
        for shift,expected in [(0.00999,0),(0.01001,1)]:
            for seed in C.data['seeds']:
                run=fixture('v32',seed)
                for row in run['points']:
                    en=C.serialized[row['E_over_N_Td']]
                    for q in V.METRICS:
                        if en in C.measured[q]: row[q]+=shift*C.measured[q][en]
                self.write('v32',seed,run)
            sp=self.directory/'summary-v32.json';V.atomic_json(sp,V.summarize(self.directory,'v32',C))
            prefix=self.directory/f'comparison-{expected}'
            self.assertEqual(M.compare(str(sp),str(rp),str(HISTORICAL),str(prefix)),expected)
            result=json.loads(Path(str(prefix)+'.json').read_text())
            self.assertEqual(result['judged'],C.judged)
            self.assertEqual(len(result['rows']),len(C.rows))
            self.assertTrue(all(r['tolerance']==max(2*r['v3_sd_rel'],0.01) for r in result['rows']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
