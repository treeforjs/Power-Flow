Generic fixed-pin component dispatcher v1

This isolated component closes the admission mechanisms identified as E1-E5.
It does not deploy a new launcher or authorize numerical/full integration work.
Existing frozen launchers, paused deployment, native/Hall sources and old test
wrappers remain unchanged. The component accepts a fixed reviewed manifest;
there is no capture-current-host, accept-current-state or refresh-pins operation.

API
---
source/test_dispatch_v1.py:
  Dispatcher(fs, context, leases, runner).dispatch(
      manifest_path, trusted_manifest_sha256, execute=False)

The four adapters are injected. The returned compact receipt is successful only
when explicit checks pass; it never marks a physics integration qualified.
Default operation verifies without calling the runner. The CLI additionally
requires a fresh receipt destination and an explicit --execute flag before a
reviewed component command can be dispatched. Any refusal returns CLI status 78.
The trusted digest is supplied separately from manifest bytes. Pin changes
produce refusal, including host changes after the child, rather than new pins.

Admission
---------
* No assert statements. Optimized interpreter execution, inherited nonzero or
  malformed PYTHONOPTIMIZE, and -O/-OO/PYTHONOPTIMIZE command overrides refuse.
* Affinity must be exactly CPU 7 and nice exactly 15 before lock acquisition.
  CPU mode acquires only the shared testing-core7.lock; explicit GPU mode
  acquires testing-gpu.lock THEN testing-core7.lock. Release is reverse order.
  Missing/symlink/nonregular lock files or ten seconds of contention refuse.
  The wrapper must not be started while already holding either of those locks.
* A fresh child environment contains fixed one-thread OMP/BLAS settings,
  one-job MAKE/CMake settings, pinned PATH/HOME/locale and unoptimized Python.
  Ambient secrets, compiler/MPI overrides, Python search paths and LD_PRELOAD
  are not copied. Only explicit RAIKOU_PIC_FIELD_BACKEND, TMPDIR and sealed
  LD_LIBRARY_PATH additions are admitted. Fixed settings cannot be overridden.
  These are execution settings, not proof that a program never creates service
  threads or ignores job flags. The reviewed command/build recipe must itself
  retain the sequential one-job contract.
* Manifest admission uses exact fields/types, rejects duplicate JSON keys,
  booleans/fractions in integer limits, unknown roles/keys, traversal and invalid
  hashes. Maximum manifest 1 MiB, 10000 identity records, 2 GiB bytes hashed, 32 sealed
  roots, 600 s child wall time and 64 MiB output. No cap may be silently increased.
* Explicit source/build/runtime/tool identities are mandatory. Regular files
  require exact bytes/hash; symlinks require the target plus a completely listed
  chain, with cycle refusal. All source files lie in sealed roots. Each sealed
  root's complete nondirectory membership must equal the listed identities;
  inserted headers/libraries therefore refuse. Every explicit library search
  directory is a sealed root with runtime pins. Hashes and memberships are
  verified before and after an execution while leases remain held.
* CPU mode requires exact host kernel release AND kernel notes. GPU mode also
  requires versions and GNU build-ID notes for nvidia/nvidia_uvm/nvidia_modeset/
  nvidia_drm, nvidia srcversion, exact GPU information and a character-device
  identity. UUID, PCI path and device minor must agree. The admitted UUID is the
  only CUDA_VISIBLE_DEVICES value passed to that child. CPU mode sets it to -1;
  this is CUDA selection, not a device-access security boundary.
* Host proc/sys pins use bounded stream length/hash rather than their commonly
  zero stat sizes. Missing, changed or inconsistent kernel/module/device pins
  refuse. No driver query executable, module reload or driver update is used.

Scope and limitations
---------------------
The engine verifies the declared expected closure; it cannot discover arbitrary
compiler includes, dlopen/plugin searches, callbacks or external configuration.
Root must derive/review that closure independently, including runtime plugin
directories and build scripts, and pin the dispatcher/interpreter/bootstrap
itself from the sealed source archive. A declaration alone is not evidence that
an omitted dependency is irrelevant. Paths are interpreted in the dispatcher's
namespace. A manifest prepared for frozen execution must name frozen-visible
paths once; do not accidentally bind host /usr libraries or nest launchers.

The generic ProcessRunner adapter has fresh bounded output, a wall timeout and
owned-session cleanup. It has NOT been qualified with real child execution.
It is not a replacement for an independently reviewed GPU RSS/VRAM/output guard,
PIDFD identity/production hold guard, or complete process-tree cleanup. Commands
that daemonize or create other sessions require that separate pinned guard.
Before/after pin checks are not continuous attestation and do not defend against
an actor changing and restoring files during a run. Use immutable isolated
copies and the sole coordinated ownership policy. No numerical guard/tolerance
or scientific source is changed by this component.

Tests
-----
tests/test_dispatch.py instantiates only MemoryFS, manufactured process context,
fake lease contexts and fake runner. It never calls RealFS/RealLeases/ProcessRunner
or reads real module/device state. The test harness itself runs inside the
unchanged frozen userspace on CPU 7/nice 15 under the shared lease. Its synthetic
GPU case performs zero GPU calls and acquires zero real GPU locks.

See outputs/pic-test-dispatch-v1-20261030 for per-version primitive receipts,
final results, source bindings and handoff. No real process adapter, numerical
integration, physical mechanism or production deployment is qualified here.
