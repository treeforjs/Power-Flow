"""Read-only production identity/hold validation through the owned host-proc bind.

The ordinary /proc belongs to run-frozen's PID namespace. Host PID identities
must only be read through /run/pic-host-proc; there is deliberately no fallback.
"""
import json
import os
from pathlib import Path
import socket

HOST_PROC = Path('/run/pic-host-proc')
RESULT_ROOT = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs')
CGROUP_NAME = '/powerflow-pic-gpu-reserved'
CGROUP = Path('/sys/fs/cgroup' + CGROUP_NAME)

def _host_guard_check(stage):
    endpoint = os.environ.get('PIC_HOST_GUARD_SOCKET', '')
    token = os.environ.get('PIC_HOST_GUARD_TOKEN', '')
    if not endpoint.startswith('@picg-') or len(token)!=64:
        raise RuntimeError('Audited host pidfd guard is absent; no CUDA load or signal allowed.')
    with socket.socket(socket.AF_UNIX,socket.SOCK_STREAM) as connection:
        connection.settimeout(3)
        connection.connect('\0'+endpoint[1:])
        connection.sendall(json.dumps({'token':token,'action':'check','stage':stage}).encode()+b'\n')
        with connection.makefile('rb') as stream:
            line=stream.readline(16385)
    if len(line)>16384 or not line.endswith(b'\n'):
        raise RuntimeError('Invalid host guard reply')
    response=json.loads(line)
    if response.get('passed') is not True or response['hold'].get('pidfd_retained') is not True:
        raise RuntimeError('Host guard has no retained production pidfd')
    return response['hold']

def _check_identity(identity, stat_text, status_text, executable, cgroup_text, events):
    before_comm, after_comm = stat_text.rsplit(')', 1)
    actual_pid = int(before_comm.split('(', 1)[0])
    process = after_comm.split()
    nspid_line = next(line for line in status_text.splitlines() if line.startswith('NSpid:'))
    namespace_pids = [int(x) for x in nspid_line.split()[1:]]
    membership = [line.split(':', 2) for line in cgroup_text.splitlines()]
    if (actual_pid != int(identity['pid']) or process[0] not in ('T', 't') or
        process[19] != str(identity['start_ticks']) or executable != identity['exe'] or
        not namespace_pids or namespace_pids[0] != int(identity['pid']) or
        ['0', '', CGROUP_NAME] not in membership or events.get('frozen') != '1' or
        events.get('populated') != '1'):
        raise RuntimeError('Production identity/manual hold is not intact; parent action required. No signal sent.')
    return {'production_host_pid': actual_pid, 'production_namespace_pids': namespace_pids,
            'start_ticks': process[19], 'state': process[0], 'exe': executable,
            'cgroup': CGROUP_NAME, 'cgroup_frozen': True}

def paused_production(stage='hold_probe'):
    # The mount must be actual read-only procfs at this exact path, rather than
    # an uploaded snapshot or the private namespace's ordinary proc mount.
    mounts = []
    for line in Path('/proc/self/mountinfo').read_text().splitlines():
        left, right = line.split(' - ', 1)
        fields = left.split()
        if fields[4] == str(HOST_PROC):
            mounts.append((set(fields[5].split(',')), right.split()[0]))
    if len(mounts) != 1 or 'ro' not in mounts[0][0] or mounts[0][1] != 'proc':
        raise RuntimeError('Owned launcher read-only host procfs bind is absent. No CUDA load or signal allowed.')
    # The caller is visible through both mounts. Compare its PID presentations
    # and NSpid chain; inspecting root-owned host PID1 namespace links would
    # require permissions unrelated to this user's production process.
    host_self = int(os.readlink(HOST_PROC / 'self'))
    private_self = int(os.readlink('/proc/self'))
    caller_status = (HOST_PROC / 'self/status').read_text()
    caller_chain = [int(x) for x in next(line for line in caller_status.splitlines() if line.startswith('NSpid:')).split()[1:]]
    if host_self == private_self or private_self != os.getpid() or len(caller_chain)<2 or caller_chain[0]!=host_self or caller_chain[-1]!=private_self:
        raise RuntimeError('Host-proc and frozen PID namespace are not distinct; refusing ambiguous PID view.')
    identity = json.loads((RESULT_ROOT / 'full60-manual-pause.json').read_text())['native_identity']
    native = HOST_PROC / str(identity['pid'])
    host_check = _host_guard_check(stage)
    events = dict(line.split() for line in (CGROUP / 'cgroup.events').read_text().splitlines())
    # The nested user namespace may read stat/status/cgroup through host procfs,
    # but ptrace permission rules prohibit reading another user namespace's exe
    # link. The retained-pidfd host guard performs that exact read itself, then
    # the child revalidates every other identity/hold field through its RO bind.
    checked = _check_identity(identity, (native / 'stat').read_text(), (native / 'status').read_text(),
                              host_check['exe'], (native / 'cgroup').read_text(), events)
    if host_check['production_host_pid']!=checked['production_host_pid'] or host_check['start_ticks']!=checked['start_ticks']:
        raise RuntimeError('Host guard and frozen host-proc identity disagree')
    # Re-read identity and frozen state to detect a change while taking the view.
    repeated_stat = (native / 'stat').read_text().rsplit(')', 1)[1].split()
    repeated_events = dict(line.split() for line in (CGROUP / 'cgroup.events').read_text().splitlines())
    if repeated_stat[0] not in ('T', 't') or repeated_stat[19] != checked['start_ticks'] or repeated_events.get('frozen') != '1':
        raise RuntimeError('Production hold changed while validating it; no signal sent.')
    checked.update(host_proc_path=str(HOST_PROC), host_proc_readonly=True,
                   distinct_pid_views=True, benchmark_namespace_pid_chain=caller_chain,
                   benchmark_host_pid=host_self, benchmark_namespace_pid=private_self,
                   host_executable_checked_with_retained_pidfd=True,host_guard=host_check)
    return checked

if __name__ == '__main__':
    print(json.dumps({'passed': True, 'gpu_initialized_or_executed': False,
                      'production_signalled': False, 'hold': paused_production()}, indent=2))
