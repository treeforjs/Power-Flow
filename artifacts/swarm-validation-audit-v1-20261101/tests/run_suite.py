from pathlib import Path
import json
import os
import unittest

suite=unittest.defaultTestLoader.discover(str(Path(__file__).parent),pattern='test_*.py')
result=unittest.TextTestRunner(verbosity=2).run(suite)
(Path(__file__).resolve().parents[1]/'evidence/test-summary.json').write_text(json.dumps({
    'scope':'synthetic validation/provenance fixtures; no physical simulation',
    'tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
    'skipped':len(result.skipped),'passed':result.wasSuccessful(),
    'logical_cpus_allowed':sorted(os.sched_getaffinity(0))},indent=2)+'\n')
raise SystemExit(not result.wasSuccessful())
