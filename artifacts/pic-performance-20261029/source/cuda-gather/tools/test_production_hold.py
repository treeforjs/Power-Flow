"""CPU-only checks of actual namespaced hold plus corrupted identity rejection."""
import copy
import json
from production_hold import paused_production, _check_identity, HOST_PROC, RESULT_ROOT, CGROUP
import os

hold = paused_production()
identity = json.loads((RESULT_ROOT / 'full60-manual-pause.json').read_text())['native_identity']
native = HOST_PROC / str(identity['pid'])
args = [identity, (native/'stat').read_text(), (native/'status').read_text(),
        hold['exe'], (native/'cgroup').read_text(),
        dict(line.split() for line in (CGROUP/'cgroup.events').read_text().splitlines())]
_check_identity(*args)
cases = []
def rejects(label, values):
    try:
        _check_identity(*values)
    except RuntimeError:
        cases.append(label)
    else:
        raise AssertionError(f'{label} bypassed the hold guard')

bad = copy.deepcopy(args); bad[0]['pid']+=1; rejects('wrong host PID',bad)
bad = copy.deepcopy(args); bad[0]['start_ticks']='0'; rejects('PID reuse/start ticks',bad)
bad = copy.deepcopy(args); bad[3]+=' (deleted)'; rejects('changed executable',bad)
bad = copy.deepcopy(args); bad[1]=bad[1].rsplit(')',1)[0]+')'+bad[1].rsplit(')',1)[1].replace(' T ',' R ',1); rejects('production not SIGSTOP',bad)
bad = copy.deepcopy(args); bad[5]['frozen']='0'; rejects('unfrozen cgroup',bad)
bad = copy.deepcopy(args); bad[4]='0::/another-job\n'; rejects('wrong cgroup membership',bad)
bad = copy.deepcopy(args); bad[2]='\n'.join(line if not line.startswith('NSpid:') else 'NSpid:\t2' for line in bad[2].splitlines()); rejects('private PID substituted for host PID',bad)
print(json.dumps({'passed':True,'gpu_initialized_or_executed':False,'production_signalled':False,
                  'positive_hold':hold,'negative_cases_rejected':cases},indent=2))
