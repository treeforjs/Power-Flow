import copy
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj)
    return obj


Q=module(ROOT/'source/qualification_integrity.py','qualification')
P=module(ROOT/'source/swarm_job_integrity.py','job_proof')
F=module(ROOT/'tests/test_validation.py','fixture_data')


class SealingTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='seal-audit-',dir=ROOT/'evidence')
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)/'component';self.out=Path(self.tmp.name)/'outputs'
        for directory in ('include','tests','device','tools','data/tables','build/bin'):
            (self.root/directory).mkdir(parents=True)
        self.out.mkdir()
        (self.root/'include/fixed.hpp').write_text('original\n')
        (self.root/'run_checks.py').write_text('# synthetic execution-identity fixture, not a simulation\n')
        for name in ('swarm_validation.py','swarm_job_integrity.py','qualification_integrity.py','make_seal.py','swarm-contract.json'):
            shutil.copy2(ROOT/'source'/name,self.root/'tools'/name)
        shutil.copytree(ROOT/'baseline/data/extracted',self.root/'data/extracted')
        shutil.copy2(ROOT/'baseline/data/tables/chemistry-v3.2-spec.json',self.root/'data/tables/chemistry-v3.2-spec.json')
        shutil.copy2(ROOT/'baseline/tests/swarm_mc.cpp',self.root/'tests/swarm_mc.cpp')
        for name in Q.BINARIES: (self.root/'build/bin'/name).write_text('synthetic binary '+name+'\n')
        self.blob=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.2-20261031/data/tables/chemistry-v3.2.blob')
        shutil.copy2(self.blob,self.root/'data/tables/chemistry-v3.2.blob')
        (self.root/'data/tables/chemistry-v3.2-coverage-policy.json').write_text('{}\n')
        (self.root/'data/tables/chemistry-v3.2-manifest.json').write_text('{"input_files_sha256":{}}\n')
        self.receipt={'schema':Q.SCHEMA,'stages':[]}
        before=Q.capture(self.root,[])
        self.binding=before['sha256']
        for directory in ('swarm','swarm-multiseed'):
            (self.out/directory).mkdir()
        jobs=[]
        for label in ('v3','v32'):
            for seed in F.C.data['seeds']:
                name=f'{label}-{seed}'
                source_blob=(Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1-20261031/data/tables/chemistry-v3.1.blob')
                             if label=='v3' else self.root/'data/tables/chemistry-v3.2.blob')
                job=self.job(name,self.out/'swarm-multiseed',source_blob,F.C.energies,seed)
                self.complete(job,F.fixture(label,seed))
                jobs.append(job)
        self.record_jobs('swarm-multiseed',jobs)
        jobs=[];self.sets={}
        for label in ('v3','v3-iso','lxcat-muroranit','lxcat-triniti','lxcat-hayashi','lxcat-phelps','lxcat-morgan'):
            legacy_blob=self.out/'swarm'/(label+'.blob')
            legacy_blob.write_text('synthetic distinct legacy blob '+label)
            self.sets[label]=Q.sha(legacy_blob)
            for k in range(3):
                name=f'{label}-{k}'; seed=20261031
                ens=[30,60,100,150,200,300,500,700,1000,1500,2000,3000][4*k:4*(k+1)]
                job=self.job(name,self.out/'swarm',legacy_blob,ens,seed)
                run=F.fixture('v32',seed);run['points']=run['points'][:4]
                run['blob_sha256']=self.sets[label]
                for row,en in zip(run['points'],ens):row['E_over_N_Td']=en
                self.complete(job,run);jobs.append(job)
        self.record_jobs('swarm',jobs)
        for label in ('v3','v32'):
            F.V.atomic_json(self.out/(label+'-summary.json'),F.V.summarize(self.out/'swarm-multiseed',label,F.C))
        with contextlib.redirect_stdout(io.StringIO()):
            F.M.compare(str(self.out/'v32-summary.json'),str(self.out/'v3-summary.json'),str(F.HISTORICAL),str(self.out/'SWARM-REVALIDATION'))
        comparison=json.loads((self.out/'SWARM-REVALIDATION.json').read_text())
        for name in Q.REQUIRED:
            before=Q.begin_stage(self.root,[],self.receipt,name)
            binding=Q.finish_stage(self.root,[],before)
            stage={'stage':name,'passed':True,'entries':[{'name':'synthetic '+name,'exit':0}],'binding':binding}
            if name=='tables': stage['blob_sha256']=Q.sha(self.root/'data/tables/chemistry-v3.2.blob')
            if name=='policy': stage['policy_sha256']=Q.sha(self.root/'data/tables/chemistry-v3.2-coverage-policy.json')
            if name=='inputs': stage['pinned_dependencies']={}
            if name=='swarm':
                stage['revalidation_sha256']=Q.sha(self.out/'SWARM-REVALIDATION.json')
                stage['scientific_criterion_status']=comparison['scientific_criterion_status']
                stage['sets']=self.sets
                stage['job_plan_sha256']={d:Q.sha(self.out/d/'jobs.json') for d in ('swarm','swarm-multiseed')}
            self.receipt['stages'].append(stage)

    def job(self,name,directory,blob,energies,seed):
        output=directory/(name+'.json')
        return {'name':name,'output':str(output),'argv':[str(self.root/'build/bin/swarm_mc'),str(blob),str(output),
                ','.join(str(x) for x in energies),'3.22e22','300','2000',str(seed)]}

    def complete(self,job,data):
        temporary=Path(job['output']+'.partial');temporary.write_text(json.dumps(data))
        P.publish(temporary,job,P.identity(job,self.binding),self.binding)

    def record_jobs(self,directory,jobs):
        p=self.out/directory
        (p/'jobs.json').write_text(json.dumps(jobs))
        (p/'jobs-receipt.json').write_text(json.dumps({'jobs':[{'name':j['name'],'exit':0,'completion_verified':True} for j in jobs]}))

    def verify(self,receipt=None):
        return Q.verify_receipt(self.root,receipt or self.receipt,self.out)

    def test_complete_fixture_is_technical_only(self):
        result=self.verify()
        self.assertFalse(result['production_qualified']);self.assertFalse(result['physics_model_validated'])
        self.assertEqual(result['required_stages'],list(Q.REQUIRED))

    def test_missing_stage_refuses(self):
        for missing in Q.REQUIRED:
            r=copy.deepcopy(self.receipt);r['stages']=[s for s in r['stages'] if s['stage']!=missing]
            with self.subTest(missing=missing),self.assertRaises(ValueError): self.verify(r)

    def test_failed_empty_duplicate_and_reordered_stages(self):
        for mode in ('failed','empty','duplicate','reorder','unbound'):
            r=copy.deepcopy(self.receipt)
            if mode=='failed':r['stages'][4]['entries'][0]['exit']=1
            if mode=='empty':r['stages'][4]['entries']=[]
            if mode=='duplicate':r['stages'].append(copy.deepcopy(r['stages'][4]))
            if mode=='reorder':r['stages'][4],r['stages'][5]=r['stages'][5],r['stages'][4]
            if mode=='unbound':r.pop('schema')
            with self.subTest(mode=mode),self.assertRaises(ValueError):self.verify(r)

    def test_current_source_and_binary_drift_refuse(self):
        path=self.root/'include/fixed.hpp';before=path.read_bytes();path.write_text('changed\n')
        with self.assertRaises(ValueError):self.verify()
        path.write_bytes(before)
        (self.root/'build/bin/test_core').write_text('different executable\n')
        with self.assertRaises(ValueError):self.verify()

    def test_stale_build_test_binding_refuses(self):
        for mode in ('source','binary','during'):
            r=copy.deepcopy(self.receipt);b=r['stages'][4]['binding']
            if mode=='source':b['before']['source_files']['include/fixed.hpp']='0'*64
            if mode=='binary':b['binaries_before']['test_core']='0'*64
            if mode=='during':b['after']['source_files']['include/fixed.hpp']='0'*64
            with self.subTest(mode=mode),self.assertRaises(ValueError):self.verify(r)

    def test_mutation_during_stage_refuses(self):
        before=Q.begin_stage(self.root,[],self.receipt,'inputs')
        (self.root/'include/fixed.hpp').write_text('modified mid-stage\n')
        with self.assertRaises(ValueError):Q.finish_stage(self.root,[],before)

    def test_parseable_file_without_completion_does_not_resume(self):
        jobs=json.loads((self.out/'swarm-multiseed/jobs.json').read_text());job=jobs[0]
        self.assertTrue(P.completed(job,self.binding))
        Path(job['output']+'.complete.json').unlink()
        self.assertFalse(P.completed(job,self.binding))
        with self.assertRaises(ValueError):self.verify()

    def test_source_output_or_command_mismatch_cannot_resume(self):
        job=json.loads((self.out/'swarm-multiseed/jobs.json').read_text())[0]
        self.assertFalse(P.completed(job,'other binding'))
        changed=copy.deepcopy(job);changed['argv'][7]='20261199'
        self.assertFalse(P.completed(changed,self.binding))
        Path(job['output']).write_text('{"valid":"json"}')
        self.assertFalse(P.completed(job,self.binding))

    def test_stale_comparison_after_new_valid_raw_output_refuses(self):
        job=json.loads((self.out/'swarm-multiseed/jobs.json').read_text())[0]
        data=json.loads(Path(job['output']).read_text());data['points'][0]['W_bulk_m_s']*=1.0001
        self.complete(job,data)
        self.assertTrue(P.completed(job,self.binding))
        with self.assertRaises(ValueError):self.verify()

    def test_legacy_wrong_blob_or_chunk_refuses(self):
        path=self.out/'swarm/jobs.json';jobs=json.loads(path.read_text());original=copy.deepcopy(jobs)
        for mode in ('blob','chunk'):
            jobs=copy.deepcopy(original)
            if mode=='blob':jobs[0]['argv'][1]=jobs[1]['argv'][1]+'-wrong'
            else:jobs[0]['argv'][3]='200,300,500,700'
            path.write_text(json.dumps(jobs))
            r=copy.deepcopy(self.receipt);r['stages'][-1]['job_plan_sha256']['swarm']=Q.sha(path)
            with self.subTest(mode=mode),self.assertRaises((ValueError,OSError)):self.verify(r)

    def test_failed_completion_keeps_authoritative_output(self):
        job=json.loads((self.out/'swarm-multiseed/jobs.json').read_text())[0]
        before=Path(job['output']).read_bytes()
        temporary=Path(job['output']+'.partial');temporary.write_text('{"points":[]}')
        with self.assertRaises(ValueError):P.publish(temporary,job,P.identity(job,self.binding),self.binding)
        self.assertEqual(Path(job['output']).read_bytes(),before)
        self.assertTrue(P.completed(job,self.binding))

    def test_sealer_refusal_preserves_existing_seal(self):
        sealer=module(self.root/'tools/make_seal.py','test_sealer')
        seal=self.root/'SEAL.json';seal.write_bytes(b'preserved historical seal\n')
        r=copy.deepcopy(self.receipt);r['stages'].pop()
        receipt=self.out/'receipt.json';receipt.write_text(json.dumps(r))
        with self.assertRaises(ValueError):sealer.main(str(receipt))
        self.assertEqual(seal.read_bytes(),b'preserved historical seal\n')


if __name__=='__main__':unittest.main(verbosity=2)
