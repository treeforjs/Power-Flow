#!/usr/bin/python -I
"""Transient CPU6/14 shield transaction. No system access on import.

Backend supplies observed allocation, cgroup reads/writes and durable audit.
Only the canary CLI below can activate this first version. No numerical run
entry point is implemented. See CANARY.txt before privileged use.
"""
from __future__ import annotations
import json

PAIR = frozenset((6, 14))
BASE = 'user.slice/user-1000.slice/user@1000.service/simq.slice'
ANCESTORS = ('', 'user.slice', 'user.slice/user-1000.slice',
             'user.slice/user-1000.slice/user@1000.service', BASE)


class Refusal(RuntimeError):
    pass


def require(value, message):
    if not value:
        raise Refusal(message)


def cpus(value):
    result = set()
    for item in value.strip().split(','):
        if item:
            pair = item.split('-')
            require(len(pair) <= 2, 'malformed CPU list')
            a, b = int(pair[0]), int(pair[-1])
            require(0 <= a <= b <= 4095, 'invalid CPU range')
            result.update(range(a, b + 1))
    return result


def path(group, name):
    return group + '/' + name if group else name


class Shield:
    """Serial transaction; never moves, signals or changes a job's affinity.

    An undo conflict stops ancestor rollback. Continuing past a changed child
    could invalidate another controller owner's partition. Report and retain
    the remaining journal for explicit recovery instead.
    """
    def __init__(self, backend, unit):
        self.b = backend
        self.group = BASE + '/' + unit
        self.journal = []
        self.identity = None
        self.active = False

    def preflight(self):
        require('cpuset' in self.b.read('cgroup.controllers').split(), 'cpuset controller unavailable')
        identity = self.b.allocation(self.group, active=False)
        if self.identity is None:
            self.identity = identity
        require(identity == self.identity, 'live allocation identity changed')
        require(not cpus(self.b.read('cpuset.cpus.isolated')), 'another isolated partition exists')
        for group in ANCESTORS:
            if group:
                require(self.b.read(path(group, 'cgroup.type')) == 'domain', 'ancestor is not a domain')
                require(not self.b.read(path(group, 'cgroup.procs')), 'ancestor contains processes')
                for name, allowed in [('cpuset.cpus.partition', ('member',)),
                                      ('cpuset.cpus.exclusive', ('',))]:
                    p = path(group, name)
                    require(not self.b.exists(p) or self.b.read(p) in allowed,
                            'pre-existing ancestor partition/exclusive assignment')
        for name, allowed in [('cpuset.cpus.partition', ('member',)),
                              ('cpuset.cpus.exclusive', ('',)), ('cpuset.cpus', ('',))]:
            p = path(self.group, name)
            require(not self.b.exists(p) or self.b.read(p) in allowed,
                    'job already has an explicit cpuset policy')
        self.b.foreign_workers(self.group)

    def change(self, p, desired, command=None, undo_command=None):
        old, inode = self.b.read(p), self.b.inode(p)
        if old == desired:
            return
        rec = {'path': p, 'before': old, 'after': desired, 'inode': inode,
               'command': command if command is not None else desired,
               'undo': undo_command if undo_command is not None else old,
               'state': 'intent'}
        # Persist intent BEFORE writing. Even a write which changes state then
        # raises is recoverable by compare-and-restore against this intent.
        self.b.audit({'event': 'intent', 'record': rec.copy()})
        self.journal.append(rec)
        self.b.write(p, rec['command'])
        require(self.b.inode(p) == inode and self.b.read(p) == desired,
                'cgroup write did not produce its exact expected state')
        rec['state'] = 'applied'
        self.b.audit({'event': 'applied', 'record': rec.copy()})

    def setup(self):
        self.preflight()
        for group in ANCESTORS:
            p = path(group, 'cgroup.subtree_control')
            enabled = set(self.b.read(p).split())
            if 'cpuset' not in enabled:
                # Linux exposes controller names in this fixed order.
                order = self.b.read('cgroup.controllers').split()
                wanted = ' '.join(x for x in order if x in enabled | {'cpuset'})
                self.change(p, wanted, '+cpuset', '-cpuset')
        # Revalidate newly exposed cpuset policy before routing exclusives.
        for group in ANCESTORS[1:] + (self.group,):
            require(self.b.read(path(group, 'cpuset.cpus.partition')) == 'member',
                    'partition appeared while enabling controller')
            require(not self.b.read(path(group, 'cpuset.cpus.exclusive')),
                    'exclusive assignment appeared while enabling controller')
        require(not self.b.read(path(self.group, 'cpuset.cpus')), 'job cpuset appeared')
        require(self.b.allocation(self.group, active=False) == self.identity, 'allocation changed during setup')
        self.b.foreign_workers(self.group)
        for group in ANCESTORS[1:]:
            self.change(path(group, 'cpuset.cpus.exclusive'), '6,14')
        self.change(path(self.group, 'cpuset.cpus'), '6,14')
        self.change(path(self.group, 'cpuset.cpus.exclusive'), '6,14')
        self.change(path(self.group, 'cpuset.cpus.partition'), 'isolated')
        self.active = True
        self.verify()
        self.b.audit({'event': 'isolated_verified', 'allocation': self.identity})

    def verify(self):
        require(self.b.allocation(self.group, active=True) == self.identity, 'allocation changed while isolated')
        require(self.b.read(path(self.group, 'cpuset.cpus.partition')) == 'isolated', 'partition lost isolation')
        for n in ('cpuset.cpus.effective', 'cpuset.cpus.exclusive.effective'):
            require(cpus(self.b.read(path(self.group, n))) == PAIR, 'effective shield is not exactly CPU6/14')
        require(cpus(self.b.read('cpuset.cpus.isolated')) == PAIR, 'global isolated set differs')
        # Ancestors necessarily describe the descendant resource as well;
        # only non-ancestor, non-job branches must exclude the shield pair.
        self.b.outsiders(self.group)
        self.b.foreign_workers(self.group)

    def rollback(self):
        errors = []
        for rec in reversed(self.journal):
            p = rec['path']
            try:
                if not self.b.exists(p):
                    # A completed transient unit can disappear normally. Its
                    # partition is already destroyed; ancestors must persist.
                    require(p.startswith(self.group + '/'), 'ancestor disappeared during cleanup')
                    rec['state'] = 'removed_with_job'
                else:
                    require(self.b.inode(p) == rec['inode'], 'cleanup inode changed: ' + p)
                    current = self.b.read(p)
                    if current == rec['before']:
                        rec['state'] = 'already_restored'
                    else:
                        require(current == rec['after'], 'cleanup concurrent write: ' + p)
                        self.b.write(p, rec['undo'])
                        require(self.b.read(p) == rec['before'], 'cleanup verification failed: ' + p)
                        rec['state'] = 'restored'
            except Exception as error:
                errors.append(str(error))
                # Do not change an ancestor underneath unresolved child state.
                break
            # Logging failure must not strand an otherwise safely owned
            # partition. Continue verified restoration even when storage is
            # unavailable, then report every missing audit record.
            try:
                self.b.audit({'event': 'undo', 'record': rec.copy()})
            except Exception as error:
                errors.append('cleanup audit unavailable: ' + str(error))
        self.active = False
        return errors
