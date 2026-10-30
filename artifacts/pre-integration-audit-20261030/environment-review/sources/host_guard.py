#!/usr/bin/env python3
"""Own the host identity/pidfd across an isolated frozen child; never signal PIC.

The host guard handles the executable-link check which nested user namespaces
cannot perform. It serves authenticated read-only checks to the child, watches
the retained production pidfd and stopped/frozen state, and can terminate only
its own new child process group if the hold fails. No GPU API is imported here.
"""
import argparse
import json
import os
from pathlib import Path
import secrets
import select
import signal
import socket
import struct
import subprocess
import time
from production_hold import _check_identity, RESULT_ROOT, CGROUP

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    launcher = ROOT / 'tools/run-frozen-with-host-proc'
    if not command or Path(command[0]).resolve() != launcher:
        parser.error('The child must use this investigation\'s owned frozen launcher')
    if os.getpriority(os.PRIO_PROCESS, 0) < 15 or os.sched_getaffinity(0) not in ({6}, {7}):
        parser.error('Host guard and child must use nice15 and one allocated CPU6 or7')
    identity = json.loads((RESULT_ROOT / 'full60-manual-pause.json').read_text())['native_identity']
    # /proc here is the caller's host view. If mistakenly nested, the exact host
    # PID/exe check below fails. The fd is retained until after child completion.
    pidfd = os.pidfd_open(int(identity['pid']), 0)
    poller = select.poll()
    poller.register(pidfd, select.POLLIN)
    native = Path('/proc') / str(identity['pid'])
    def hold():
        if poller.poll(0):
            raise RuntimeError('Production pidfd became readable: original process exited')
        result = _check_identity(identity, (native/'stat').read_text(), (native/'status').read_text(),
                                 os.readlink(native/'exe'), (native/'cgroup').read_text(),
                                 dict(line.split() for line in (CGROUP/'cgroup.events').read_text().splitlines()))
        if poller.poll(0):
            raise RuntimeError('Production exited while checking identity')
        result.update(pidfd_retained=True, host_guard_pid=os.getpid(), checked_unix_s=time.time())
        return result
    before = hold()
    # Abstract AF_UNIX avoids Linux's short pathname limit without creating
    # files outside the assigned prefix. The launcher does not unshare network.
    endpoint = '@picg-' + secrets.token_hex(16)
    token = secrets.token_hex(32)
    environment = dict(os.environ, PIC_HOST_GUARD_SOCKET=str(endpoint), PIC_HOST_GUARD_TOKEN=token)
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    child = None
    error = None
    stopped = False
    checks = 1
    requests = []
    def interrupt(signum, frame):
        nonlocal stopped
        stopped = True
    previous = {s: signal.signal(s, interrupt) for s in (signal.SIGTERM, signal.SIGINT)}
    args.out.mkdir(parents=True, exist_ok=True)
    try:
        server.bind('\0' + endpoint[1:])
        server.listen(1)
        child = subprocess.Popen(command, env=environment, start_new_session=True)
        while child.poll() is None:
            if stopped:
                raise RuntimeError('Host guard interrupted; terminating only owned benchmark child')
            hold()
            checks += 1
            ready, _, _ = select.select([server, pidfd], [], [], .1)
            if pidfd in ready:
                raise RuntimeError('Production pidfd signalled exit')
            if server in ready:
                connection, _ = server.accept()
                with connection:
                    connection.settimeout(2)
                    peer_pid,peer_uid,_=struct.unpack('3i',connection.getsockopt(socket.SOL_SOCKET,socket.SO_PEERCRED,12))
                    if peer_uid!=os.getuid() or os.getpgid(peer_pid)!=child.pid:
                        raise RuntimeError('Host guard request is outside the owned frozen child process group')
                    with connection.makefile('rb') as stream:
                        line = stream.readline(4097)
                    if len(line)>4096 or not line.endswith(b'\n'):
                        raise RuntimeError('Invalid host guard request length')
                    request = json.loads(line)
                    if not secrets.compare_digest(request.get('token', ''), token) or request.get('action') != 'check':
                        raise RuntimeError('Invalid host guard check request')
                    result = hold()
                    requests.append({'stage': request.get('stage'), 'check': result})
                    connection.sendall(json.dumps({'passed': True, 'hold': result}).encode()+b'\n')
        code = child.returncode
        after = hold()
    except BaseException as exc:
        error = str(exc)
        code = 1
        after = None
        raise
    finally:
        # The new session/group belongs exclusively to this Popen. Never use
        # production PID, cgroup members or the caller's process group here.
        if child is not None and child.poll() is None:
            assert child.pid != int(identity['pid'])
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                child.wait(timeout=2)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait(timeout=2)
        try:
            final_hold = hold()
        except Exception as exc:
            final_hold = {'passed': False, 'error': str(exc)}
        receipt = {'production_signalled': False, 'gpu_actions_in_guard': False,
                   'host_guard_pid': os.getpid(), 'logical_cpu': sorted(os.sched_getaffinity(0)),
                   'nice': os.getpriority(os.PRIO_PROCESS,0), 'production_pidfd_retained_until_child_exit': True,
                   'command': command, 'before': before, 'after': final_hold,
                   'periodic_checks': checks, 'child_requests': requests,
                   'child_host_pid': None if child is None else child.pid,
                   'child_returncode': None if child is None else child.returncode, 'guard_error': error}
        (args.out/'host-guard.json').write_text(json.dumps(receipt,indent=2)+'\n')
        server.close()
        os.close(pidfd)
        for s, handler in previous.items():
            signal.signal(s, handler)
    return code

if __name__ == '__main__':
    raise SystemExit(main())
