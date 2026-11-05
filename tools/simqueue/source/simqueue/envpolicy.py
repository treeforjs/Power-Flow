"""What a submission passes to its job from the submitter's environment (0.7.0).

Up to 0.6.3 the whole environment of `simq submit` was captured, stored in
the database and given to the job: a GitHub workflow's ``GH_TOKEN`` ended up
in ``queue.sqlite3`` and ``RUNNER_TRACKING_ID`` made the runner kill the jobs
that inherited it.  Now:

* only an allowlist is taken from the environment (``ALLOWED``), plus the
  names given with ``--export NAME,...`` (``--export ALL``: everything);
* a name matching a secret pattern (``*TOKEN*``, ``*SECRET*``,
  ``*PASSWORD*``, ``*KEY*``, ``GITHUB_*``, ``ACTIONS_*``, ``RUNNER_*``) is
  never taken from the environment, whatever ``--export`` says;
* ``-e NAME=VALUE`` sets a variable explicitly; a secret given that way
  reaches the job, but the database stores ``<redacted>`` for it and the
  value waits in a private file of the state directory
  (``secrets/<job>.json``, mode 600) that the worker deletes when
  the job ends, when it also redacts the job's ``spec.json``.

``scrub_environment`` redacts secret-pattern values of stored environments
(the schema-9 migration does so for every existing job).
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Iterable, Mapping

ALLOWED = ("HOME", "USER", "LOGNAME", "LANG", "PATH", "TZ")
# SimQueue's own settings (state directory, node config, project root): a job
# that runs `simq` (an ensemble supervisor, a nested submission) stays on the
# same queue
ALLOWED_PREFIXES = ("SIMQUEUE_",)
SECRET_SUBSTRINGS = ("TOKEN", "SECRET", "PASSWORD", "KEY")
SECRET_PREFIXES = ("GITHUB_", "ACTIONS_", "RUNNER_")
REDACTED = "<redacted>"
# Provided by the user manager to every job (not from the submitter): what
# `systemctl --user` and D-Bus clients inside a job need.
SESSION = ("XDG_RUNTIME_DIR", "DBUS_SESSION_BUS_ADDRESS")
NAME_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def is_secret(name: str) -> bool:
    upper = name.upper()
    if upper.startswith("SIMQ_"):
        return False  # SimQueue's own placement variables (SIMQ_LAUNCH_TOKEN is the adoption token)
    return any(part in upper for part in SECRET_SUBSTRINGS) or upper.startswith(SECRET_PREFIXES)


def parse_export(text: str | None) -> tuple[str, ...] | str:
    """``--export``: ``None``/``NONE`` (the allowlist only), ``ALL``, or names."""
    if text is None or text.strip().upper() == "NONE":
        return ()
    if text.strip().upper() == "ALL":
        return "ALL"
    names = tuple(item.strip() for item in text.split(",") if item.strip())
    for name in names:
        if NAME_PATTERN.fullmatch(name) is None:
            raise ValueError(f"--export: invalid variable name {name!r}")
    return names


def capture(
    environ: Mapping[str, str],
    overrides: Iterable[tuple[str, str]] = (),
    export: tuple[str, ...] | str = (),
) -> tuple[dict[str, str], dict[str, str], list[str], list[str]]:
    """Returns (environment to store, secrets for the job only, override
    names, names withheld because they look secret)."""
    stored: dict[str, str] = {}
    withheld: list[str] = []
    for name, value in environ.items():
        wanted = (name in ALLOWED or name.startswith(ALLOWED_PREFIXES) or export == "ALL"
                  or (export != "ALL" and name in export))
        if not wanted:
            continue
        if is_secret(name):
            withheld.append(name)
            continue
        stored[name] = value
    secrets: dict[str, str] = {}
    keys: list[str] = []
    for name, value in overrides:
        keys.append(name)
        if is_secret(name):
            stored[name] = REDACTED
            secrets[name] = value
        else:
            stored[name] = value
            secrets.pop(name, None)
    return stored, secrets, keys, sorted(withheld)


def scrub_environment(environment: Mapping[str, str]) -> tuple[dict[str, str], list[str]]:
    """Secret-pattern values replaced by ``<redacted>``; returns the names."""
    cleaned = dict(environment)
    names = []
    for name, value in environment.items():
        if is_secret(name) and value != REDACTED:
            cleaned[name] = REDACTED
            names.append(name)
    return cleaned, names


def scrub_json(text: str | None) -> tuple[str | None, list[str]]:
    """``scrub_environment`` of a stored ``env_json`` (unchanged text when
    nothing is secret, so untouched rows keep their exact bytes)."""
    if not text:
        return text, []
    try:
        value = json.loads(text)
    except ValueError:
        return text, []
    if not isinstance(value, dict):
        return text, []
    cleaned, names = scrub_environment({str(k): str(v) for k, v in value.items()})
    if not names:
        return text, []
    return json.dumps(cleaned, ensure_ascii=False, sort_keys=True), names


# ----- the secrets file of a job ---------------------------------------------------


def secrets_path(state: Path, job_id: int, token: str | None = None) -> Path:
    """By job ID: the launch token changes at every launch."""
    return state / "secrets" / f"{int(job_id):06d}.json"


def write_secrets(state: Path, job_id: int, token: str, secrets: Mapping[str, str]) -> None:
    path = secrets_path(state, job_id, token)
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    try:
        os.write(descriptor, json.dumps(dict(secrets)).encode("utf-8"))
    finally:
        os.close(descriptor)
    os.replace(temporary, path)


def read_secrets(state: Path, job_id: int, token: str) -> dict[str, str]:
    try:
        value = json.loads(secrets_path(state, job_id, token).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {str(k): str(v) for k, v in value.items()} if isinstance(value, dict) else {}


def delete_secrets(state: Path, job_id: int, token: str) -> None:
    try:
        secrets_path(state, job_id, token).unlink()
    except FileNotFoundError:
        pass


def scrub_spec(path: Path) -> list[str]:
    """Redact secret values in a job's ``spec.json`` (after the job ended)."""
    try:
        spec = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    env = spec.get("env") if isinstance(spec, dict) else None
    if not isinstance(env, dict):
        return []
    cleaned, names = scrub_environment({str(k): str(v) for k, v in env.items()})
    if not names:
        return []
    spec["env"] = cleaned
    temporary = path.with_name(path.name + ".scrub")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    try:
        os.write(descriptor, (json.dumps(spec, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    finally:
        os.close(descriptor)
    os.replace(temporary, path)
    return names
