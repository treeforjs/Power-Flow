from __future__ import annotations

import json
import os
import pwd
import stat
import subprocess
from pathlib import Path
from typing import Any

from .config import parse_cpu_list, project_root, state_dir


PRODUCTION_CONFIG_PATH = Path("/etc/simqueue/exclusive.json")
CONFIG_PATH = PRODUCTION_CONFIG_PATH
CONFIG_VERSION = 1
INSTALL_FOOTPRINT = (
    Path("/etc/systemd/system/simqueue-runner.service"),
    Path("/etc/systemd/system/simqueue.slice"),
    Path("/etc/simqueue/shield.conf"),
    Path("/usr/local/libexec/simqueue-shield"),
)


def _reject_incomplete_production_install() -> None:
    for artifact in INSTALL_FOOTPRINT:
        try:
            artifact.lstat()
        except FileNotFoundError:
            continue
        except OSError as error:
            raise RuntimeError(
                f"cannot inspect exclusive-install artifact {artifact}: {error}"
            ) from error
        raise RuntimeError(
            "exclusive simulation queue installation is incomplete: "
            f"{PRODUCTION_CONFIG_PATH} is missing while {artifact} remains"
        )


def load_config() -> dict[str, Any] | None:
    """Return the validated root-installed exclusive-mode marker, if present.

    Presence of this file is a fail-closed policy switch.  Callers must not
    fall back to an ordinary user worker when a malformed marker is present.
    """
    path = CONFIG_PATH
    if path == PRODUCTION_CONFIG_PATH:
        try:
            metadata = path.lstat()
        except FileNotFoundError:
            _reject_incomplete_production_install()
            return None
        except OSError as error:
            raise RuntimeError(
                f"cannot inspect exclusive-mode configuration {path}: {error}"
            ) from error
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid not in (0, 65534)
            or metadata.st_mode & 0o022
        ):
            raise RuntimeError(
                f"exclusive-mode configuration must be a root-owned, "
                f"non-writable regular file: {path}"
            )
    try:
        raw = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None
    except OSError as error:
        raise RuntimeError(
            f"cannot read exclusive-mode configuration {path}: {error}"
        ) from error
    try:
        config = json.loads(raw)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"invalid exclusive-mode configuration {path}: {error}"
        ) from error
    if not isinstance(config, dict):
        raise RuntimeError(
            f"invalid exclusive-mode configuration {path}: expected a JSON object"
        )

    required_strings = (
        "user",
        "project_root",
        "state_dir",
        "pcores",
        "ecores",
        "runner_unit",
        "helper",
    )
    if config.get("version") != CONFIG_VERSION or not isinstance(
        config.get("uid"), int
    ):
        raise RuntimeError(f"unsupported exclusive-mode configuration {path}")
    if any(
        not isinstance(config.get(key), str) or not config[key]
        for key in required_strings
    ):
        raise RuntimeError(f"incomplete exclusive-mode configuration {path}")

    current_user = pwd.getpwuid(os.getuid()).pw_name
    if int(config["uid"]) != os.getuid() or config["user"] != current_user:
        raise RuntimeError("exclusive simulation queue is installed for a different user")
    if Path(str(config["project_root"])).resolve() != project_root():
        raise RuntimeError(
            "exclusive simulation queue project root does not match this CLI"
        )
    if Path(str(config["state_dir"])).resolve() != state_dir():
        raise RuntimeError(
            "exclusive simulation queue state directory does not match this CLI"
        )
    if config["runner_unit"] != "simqueue-runner.service":
        raise RuntimeError("exclusive simulation queue names an unexpected runner unit")
    helper = Path(str(config["helper"]))
    if helper != Path("/usr/local/libexec/simqueue-shield"):
        raise RuntimeError("exclusive simulation queue names an unexpected shield helper")
    try:
        pcores = parse_cpu_list(str(config["pcores"]))
        ecores = parse_cpu_list(str(config["ecores"]))
    except ValueError as error:
        raise RuntimeError(
            f"exclusive simulation queue has an invalid CPU set: {error}"
        ) from error
    if not pcores:
        raise RuntimeError("exclusive simulation queue has an empty P-core set")
    if not ecores:
        raise RuntimeError("exclusive simulation queue has an empty E-core set")
    if set(pcores) & set(ecores):
        raise RuntimeError("exclusive simulation queue P-core and E-core sets overlap")
    return config


def helper_status(config: dict[str, Any]) -> dict[str, Any]:
    """Read the root helper's machine-readable status without changing state."""
    helper = str(config["helper"])
    try:
        result = subprocess.run(
            [helper, "status"],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RuntimeError(f"cannot query P-core shield helper {helper}: {error}") from error
    detail = result.stderr.strip()
    try:
        status = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        suffix = f": {detail}" if detail else ""
        raise RuntimeError(f"P-core shield returned invalid status{suffix}") from error
    if result.returncode != 0 or not isinstance(status, dict):
        message = status.get("message") if isinstance(status, dict) else None
        errors = status.get("errors") if isinstance(status, dict) else None
        if isinstance(errors, list) and errors:
            message = "; ".join(str(item) for item in errors)
        detail = str(message or detail or f"exit status {result.returncode}")
        raise RuntimeError(f"P-core shield is not in a safe state: {detail}")
    return status


def current_cgroup() -> str:
    try:
        for line in Path("/proc/self/cgroup").read_text(encoding="utf-8").splitlines():
            if line.startswith("0::"):
                return line[3:]
    except OSError as error:
        raise RuntimeError(f"cannot determine worker cgroup: {error}") from error
    raise RuntimeError("cannot determine worker's unified cgroup")


def verify_worker_context(config: dict[str, Any]) -> dict[str, Any]:
    """Fail unless this process is the shielded system runner."""
    expected = "/simqueue.slice/simqueue-runner.service"
    actual = current_cgroup()
    if actual != expected:
        raise RuntimeError(
            "exclusive mode forbids a foreground or detached worker; "
            f"expected cgroup {expected}, found {actual}"
        )
    status = helper_status(config)
    if status.get("state") != "root" or status.get("valid") is not True:
        raise RuntimeError(
            "exclusive runner started without an active, valid P-core partition"
        )
    expected_pcores = set(parse_cpu_list(str(config["pcores"])))
    reported = status.get("exclusive_effective", status.get("exclusive_cpus"))
    if not isinstance(reported, str):
        raise RuntimeError("P-core shield did not report its exclusive CPU set")
    try:
        actual_pcores = set(parse_cpu_list(reported))
    except ValueError as error:
        raise RuntimeError(
            "P-core shield reported an invalid exclusive CPU set"
        ) from error
    if actual_pcores != expected_pcores:
        raise RuntimeError(
            "exclusive runner's effective P-core partition does not match configuration"
        )
    return status
