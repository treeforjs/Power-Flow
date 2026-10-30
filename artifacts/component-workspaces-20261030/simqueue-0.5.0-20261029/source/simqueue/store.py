from __future__ import annotations

import hashlib
import json
import re
import secrets
import sqlite3
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

from .config import (
    database_path,
    ensure_state_directories,
    format_cpu_list,
    log_dir,
    state_dir,
)


ACTIVE_STATES = ("starting", "running", "cancelling")
TERMINAL_STATES = ("succeeded", "failed", "canceled", "interrupted")
MEMBER_TERMINAL_STATES = ("succeeded", "failed", "canceled", "interrupted")
MIGRATION_SENTINEL = Path("/run/simqueue/exclusive-installing")
IDEMPOTENCY_KEY_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:/-]{0,254}$")
SCHEMA_VERSION = 6
ENSEMBLE_CONTROL_PROTOCOL = 1
ENSEMBLE_DRAIN_EXIT_CODE = 75
ENSEMBLE_DRAIN_PROTOCOL_ENV = "SIMQ_ENSEMBLE_DRAIN_PROTOCOL"


def _validate_idempotency_key(value: str | None) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or IDEMPOTENCY_KEY_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "idempotency key must start with an alphanumeric and contain at most "
            "255 letters, digits, dot, underscore, colon, slash, or dash characters"
        )
    return value


def _semantic_environment(environment: dict[str, str]) -> dict[str, str]:
    """Remove values assigned by SimQueue itself before request hashing.

    A controller replay can run under a different queue job and therefore have
    different ``SIMQ_*`` lineage values.  Those values, and ``PWD``, are
    replaced by the worker at launch and are not part of the caller's logical
    submission identity.
    """
    return {
        key: value
        for key, value in environment.items()
        if key != "PWD" and not key.startswith("SIMQ_")
    }


def _request_sha256(payload: dict[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _command_request_sha256(
    *,
    name: str,
    cwd: Path,
    argv: list[str],
    environment: dict[str, str],
    env_override_keys: Iterable[str],
    requested_cores: int,
    retry_of: int | None,
    attempt: int,
) -> str:
    semantic_environment = _semantic_environment(environment)
    return _request_sha256(
        {
            "kind": "command-v1",
            "name": name,
            "cwd": str(cwd),
            "argv": argv,
            "environment": semantic_environment,
            "env_override_keys": sorted(
                {
                    key
                    for key in env_override_keys
                    if key != "PWD" and not key.startswith("SIMQ_")
                }
            ),
            "requested_cores": requested_cores,
            "retry_of": retry_of,
            "attempt": attempt,
        }
    )


def _ensemble_request_sha256(
    *,
    name: str,
    cwd: Path,
    environment: dict[str, str],
    env_override_keys: Iterable[str],
    requested_cores: int,
    runner_prefix: list[str],
    schema_version: str,
    canonical_spec: dict[str, Any],
    members: list[dict[str, Any]],
    max_parallel: int,
    gpu_device: str,
    gpu_budget_bytes: int,
    diagnostics_interval: float,
    mps_mode: str,
    mps_pipe_directory: str | None,
    fail_fast: bool,
    allow_shared_cwd: bool,
    retry_of: int | None,
    attempt: int,
) -> str:
    semantic_environment = _semantic_environment(environment)
    return _request_sha256(
        {
            "kind": "ensemble-v1",
            "name": name,
            "cwd": str(cwd),
            "environment": semantic_environment,
            "env_override_keys": sorted(
                {
                    key
                    for key in env_override_keys
                    if key != "PWD" and not key.startswith("SIMQ_")
                }
            ),
            "requested_cores": requested_cores,
            "runner_prefix": runner_prefix,
            "schema_version": schema_version,
            "canonical_spec": canonical_spec,
            "members": members,
            "max_parallel": max_parallel,
            "gpu_device": gpu_device,
            "gpu_budget_bytes": gpu_budget_bytes,
            "diagnostics_interval": diagnostics_interval,
            "mps_mode": mps_mode,
            "mps_pipe_directory": mps_pipe_directory,
            "fail_fast": fail_fast,
            "allow_shared_cwd": allow_shared_cwd,
            "retry_of": retry_of,
            "attempt": attempt,
        }
    )


def _safe_member_filename(index: int, name: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip(".-") or "member"
    return f"{index:04d}-{slug[:80]}.log"


class _Connection(sqlite3.Connection):
    """A sqlite connection whose context manager also closes the handle."""

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> bool:
        try:
            return super().__exit__(exc_type, exc, traceback)
        finally:
            self.close()


def utc_now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")


class QueueStore:
    def __init__(self, path: Path | None = None) -> None:
        ensure_state_directories()
        self.path = path or database_path()
        self._initialize()

    def connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.path, timeout=5.0, factory=_Connection)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA busy_timeout = 5000")
        connection.execute("PRAGMA foreign_keys = ON")
        connection.execute("PRAGMA synchronous = FULL")
        return connection

    def _initialize(self) -> None:
        with self.connect() as connection:
            connection.execute("PRAGMA journal_mode = WAL")
            stored_version = int(
                connection.execute("PRAGMA user_version").fetchone()[0]
            )
            if stored_version > SCHEMA_VERSION:
                raise RuntimeError(
                    f"queue database schema {stored_version} is newer than this "
                    f"SimQueue build (maximum {SCHEMA_VERSION})"
                )
            connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
                    name                TEXT NOT NULL,
                    state               TEXT NOT NULL CHECK (state IN (
                                            'queued', 'starting', 'running',
                                            'cancelling', 'succeeded', 'failed',
                                            'canceled', 'interrupted')),
                    launch_token        TEXT NOT NULL,
                    cwd                 TEXT NOT NULL,
                    argv_json           TEXT NOT NULL,
                    env_json            TEXT NOT NULL,
                    env_override_keys   TEXT NOT NULL DEFAULT '[]',
                    requested_cores     INTEGER NOT NULL CHECK (requested_cores > 0),
                    cpu_list            TEXT,
                    log_path            TEXT NOT NULL,
                    created_at          TEXT NOT NULL,
                    started_at          TEXT,
                    finished_at         TEXT,
                    exit_code           INTEGER,
                    pid                 INTEGER,
                    pgid                INTEGER,
                    boot_id             TEXT,
                    process_start_ticks INTEGER,
                    cancel_requested_at TEXT,
                    message             TEXT,
                    retry_of            INTEGER REFERENCES jobs(id),
                    attempt             INTEGER NOT NULL DEFAULT 1,
                    idempotency_key     TEXT,
                    request_sha256      TEXT,
                    queue_order         INTEGER NOT NULL DEFAULT 0
                );

                CREATE TABLE IF NOT EXISTS queue_meta (
                    key   TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );

                INSERT OR IGNORE INTO queue_meta(key, value)
                    VALUES ('paused', '0');

                CREATE TABLE IF NOT EXISTS ensembles (
                    job_id                    INTEGER PRIMARY KEY
                                              REFERENCES jobs(id) ON DELETE CASCADE,
                    schema_version            TEXT NOT NULL,
                    spec_json                 TEXT NOT NULL,
                    max_parallel              INTEGER NOT NULL CHECK (max_parallel > 0),
                    gpu_device                TEXT NOT NULL,
                    gpu_budget_bytes          INTEGER NOT NULL CHECK (gpu_budget_bytes >= 0),
                    diagnostics_interval      REAL NOT NULL CHECK (diagnostics_interval > 0),
                    mps_mode                  TEXT NOT NULL CHECK (mps_mode IN ('off', 'required')),
                    mps_pipe_directory        TEXT,
                    fail_fast                 INTEGER NOT NULL CHECK (fail_fast IN (0, 1)),
                    allow_shared_cwd          INTEGER NOT NULL CHECK (allow_shared_cwd IN (0, 1)),
                    artifact_dir              TEXT NOT NULL,
                    launch_paused             INTEGER NOT NULL DEFAULT 0
                                              CHECK (launch_paused IN (0, 1)),
                    control_generation        INTEGER NOT NULL DEFAULT 0
                                              CHECK (control_generation >= 0),
                    control_changed_at        TEXT,
                    drained_generation        INTEGER,
                    drained_at                TEXT,
                    supervisor_protocol       INTEGER,
                    supervisor_pid            INTEGER
                );

                CREATE TABLE IF NOT EXISTS ensemble_members (
                    job_id                    INTEGER NOT NULL
                                              REFERENCES ensembles(job_id) ON DELETE CASCADE,
                    member_index              INTEGER NOT NULL CHECK (member_index >= 0),
                    name                      TEXT NOT NULL,
                    state                     TEXT NOT NULL CHECK (state IN (
                                                  'queued', 'running', 'succeeded',
                                                  'failed', 'canceled', 'interrupted')),
                    cwd                       TEXT NOT NULL,
                    argv_json                 TEXT NOT NULL,
                    env_json                  TEXT NOT NULL,
                    gpu_admission_bytes       INTEGER NOT NULL CHECK (gpu_admission_bytes >= 0),
                    uvm_arena_bytes           INTEGER NOT NULL CHECK (uvm_arena_bytes >= 0),
                    lane                      INTEGER,
                    cpu                       INTEGER,
                    log_path                  TEXT NOT NULL,
                    member_token              TEXT NOT NULL,
                    created_at                TEXT NOT NULL,
                    started_at                TEXT,
                    finished_at               TEXT,
                    exit_code                 INTEGER,
                    pid                       INTEGER,
                    pgid                      INTEGER,
                    message                   TEXT,
                    last_sample_at            TEXT,
                    last_rss_bytes            INTEGER,
                    peak_rss_bytes            INTEGER,
                    last_threads              INTEGER,
                    uvm_policy                TEXT,
                    uvm_limit_bytes           INTEGER,
                    uvm_peak_live_bytes       INTEGER,
                    uvm_peak_reserved_bytes   INTEGER,
                    uvm_pool_reclaims         INTEGER,
                    uvm_teardown_ok           INTEGER,
                    uvm_release               TEXT,
                    uvm_marker_json           TEXT,
                    PRIMARY KEY(job_id, member_index),
                    UNIQUE(job_id, name)
                );

                CREATE INDEX IF NOT EXISTS ensemble_members_state_idx
                    ON ensemble_members(job_id, state, member_index);
                """
            )
            connection.execute("BEGIN IMMEDIATE")
            # Re-read the version while holding the migration lock. Another
            # initializer may have upgraded the database while this process
            # was waiting for BEGIN IMMEDIATE.
            locked_version = int(
                connection.execute("PRAGMA user_version").fetchone()[0]
            )
            if locked_version > SCHEMA_VERSION:
                raise RuntimeError(
                    f"queue database schema {locked_version} is newer than this "
                    f"SimQueue build (maximum {SCHEMA_VERSION})"
                )
            columns = {
                str(row["name"])
                for row in connection.execute("PRAGMA table_info(jobs)").fetchall()
            }
            if "launch_token" not in columns:
                connection.execute("ALTER TABLE jobs ADD COLUMN launch_token TEXT")
            if "idempotency_key" not in columns:
                connection.execute("ALTER TABLE jobs ADD COLUMN idempotency_key TEXT")
            if "request_sha256" not in columns:
                connection.execute("ALTER TABLE jobs ADD COLUMN request_sha256 TEXT")
            if "queue_order" not in columns:
                connection.execute(
                    "ALTER TABLE jobs ADD COLUMN queue_order INTEGER NOT NULL DEFAULT 0"
                )
                # IDs are the historical FIFO order. Seed the new explicit
                # ordering from them exactly once during the v6 migration.
                connection.execute("UPDATE jobs SET queue_order = id")
            ensemble_columns = {
                str(row["name"])
                for row in connection.execute(
                    "PRAGMA table_info(ensembles)"
                ).fetchall()
            }
            ensemble_migrations = (
                (
                    "launch_paused",
                    "INTEGER NOT NULL DEFAULT 0 CHECK (launch_paused IN (0, 1))",
                ),
                (
                    "control_generation",
                    "INTEGER NOT NULL DEFAULT 0 CHECK (control_generation >= 0)",
                ),
                ("control_changed_at", "TEXT"),
                ("drained_generation", "INTEGER"),
                ("drained_at", "TEXT"),
                ("supervisor_protocol", "INTEGER"),
                ("supervisor_pid", "INTEGER"),
            )
            for name, declaration in ensemble_migrations:
                if name not in ensemble_columns:
                    connection.execute(
                        f"ALTER TABLE ensembles ADD COLUMN {name} {declaration}"
                    )
            connection.execute(
                """
                UPDATE jobs
                SET launch_token = lower(hex(randomblob(16)))
                WHERE launch_token IS NULL
                """
            )
            connection.execute(
                """
                CREATE UNIQUE INDEX IF NOT EXISTS jobs_idempotency_key_uidx
                    ON jobs(idempotency_key)
                    WHERE idempotency_key IS NOT NULL
                """
            )
            connection.execute(
                """
                CREATE INDEX IF NOT EXISTS jobs_dispatch_order_idx
                    ON jobs(state, queue_order, id)
                """
            )
            connection.execute("DROP INDEX IF EXISTS jobs_dispatch_idx")
            connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
        self._make_files_private()

    def _make_files_private(self) -> None:
        for suffix in ("", "-wal", "-shm"):
            path = Path(f"{self.path}{suffix}")
            if path.exists():
                try:
                    path.chmod(0o600)
                except OSError:
                    pass

    @staticmethod
    def _row(row: sqlite3.Row | None) -> dict[str, Any] | None:
        return dict(row) if row is not None else None

    @staticmethod
    def _ensure_not_migrating() -> None:
        if MIGRATION_SENTINEL.exists():
            raise RuntimeError(
                "exclusive queue installation is in progress; retry this command shortly"
            )

    @staticmethod
    def _existing_idempotent_job(
        connection: sqlite3.Connection,
        *,
        idempotency_key: str | None,
        request_sha256: str | None,
    ) -> int | None:
        if idempotency_key is None:
            return None
        row = connection.execute(
            """
            SELECT id, request_sha256 FROM jobs
            WHERE idempotency_key = ?
            """,
            (idempotency_key,),
        ).fetchone()
        if row is None:
            return None
        if row["request_sha256"] != request_sha256:
            raise ValueError(
                f"idempotency key {idempotency_key!r} is already bound to a "
                "different submission request"
            )
        return int(row["id"])

    def submit(
        self,
        *,
        name: str,
        cwd: Path,
        argv: list[str],
        environment: dict[str, str],
        env_override_keys: Iterable[str],
        requested_cores: int,
        retry_of: int | None = None,
        attempt: int = 1,
        idempotency_key: str | None = None,
    ) -> dict[str, Any]:
        idempotency_key = _validate_idempotency_key(idempotency_key)
        override_keys = sorted(set(env_override_keys))
        request_sha256 = (
            _command_request_sha256(
                name=name,
                cwd=cwd,
                argv=argv,
                environment=environment,
                env_override_keys=override_keys,
                requested_cores=requested_cores,
                retry_of=retry_of,
                attempt=attempt,
            )
            if idempotency_key is not None
            else None
        )
        created_at = utc_now()
        path: Path | None = None
        created_log = False
        submission_created = True
        try:
            with self.connect() as connection:
                connection.execute("BEGIN IMMEDIATE")
                self._ensure_not_migrating()
                existing_job_id = self._existing_idempotent_job(
                    connection,
                    idempotency_key=idempotency_key,
                    request_sha256=request_sha256,
                )
                if existing_job_id is not None:
                    job_id = existing_job_id
                    submission_created = False
                else:
                    cursor = connection.execute(
                        """
                        INSERT INTO jobs(
                            name, state, launch_token, cwd, argv_json, env_json,
                            env_override_keys, requested_cores, log_path, created_at,
                            retry_of, attempt, idempotency_key, request_sha256)
                        VALUES (?, 'queued', ?, ?, ?, ?, ?, ?, '', ?, ?, ?, ?, ?)
                        """,
                        (
                            name,
                            secrets.token_hex(16),
                            str(cwd),
                            json.dumps(argv, ensure_ascii=False),
                            json.dumps(environment, ensure_ascii=False, sort_keys=True),
                            json.dumps(override_keys),
                            requested_cores,
                            created_at,
                            retry_of,
                            attempt,
                            idempotency_key,
                            request_sha256,
                        ),
                    )
                    job_id = int(cursor.lastrowid)
                    path = log_dir() / f"{job_id:06d}.log"
                    with path.open("xb"):
                        pass
                    created_log = True
                    path.chmod(0o600)
                    connection.execute(
                        "UPDATE jobs SET log_path = ? WHERE id = ?",
                        (str(path), job_id),
                    )
                    connection.execute(
                        "UPDATE jobs SET queue_order = ? WHERE id = ?",
                        (job_id, job_id),
                    )
        except BaseException:
            if created_log and path is not None:
                try:
                    path.unlink()
                except OSError:
                    pass
            raise
        job = self.require(job_id)
        if idempotency_key is not None:
            job["submission_created"] = submission_created
        return job

    def submit_ensemble(
        self,
        *,
        name: str,
        cwd: Path,
        environment: dict[str, str],
        env_override_keys: Iterable[str],
        requested_cores: int,
        runner_prefix: list[str],
        schema_version: str,
        canonical_spec: dict[str, Any],
        members: list[dict[str, Any]],
        max_parallel: int,
        gpu_device: str,
        gpu_budget_bytes: int,
        diagnostics_interval: float,
        mps_mode: str,
        mps_pipe_directory: str | None,
        fail_fast: bool,
        allow_shared_cwd: bool,
        retry_of: int | None = None,
        attempt: int = 1,
        idempotency_key: str | None = None,
    ) -> dict[str, Any]:
        """Atomically persist an ensemble parent and its immutable members."""
        if not members:
            raise ValueError("an ensemble must contain at least one member")
        idempotency_key = _validate_idempotency_key(idempotency_key)
        override_keys = sorted(set(env_override_keys))
        request_sha256 = (
            _ensemble_request_sha256(
                name=name,
                cwd=cwd,
                environment=environment,
                env_override_keys=override_keys,
                requested_cores=requested_cores,
                runner_prefix=runner_prefix,
                schema_version=schema_version,
                canonical_spec=canonical_spec,
                members=members,
                max_parallel=max_parallel,
                gpu_device=gpu_device,
                gpu_budget_bytes=gpu_budget_bytes,
                diagnostics_interval=diagnostics_interval,
                mps_mode=mps_mode,
                mps_pipe_directory=mps_pipe_directory,
                fail_fast=fail_fast,
                allow_shared_cwd=allow_shared_cwd,
                retry_of=retry_of,
                attempt=attempt,
            )
            if idempotency_key is not None
            else None
        )
        created_at = utc_now()
        created_files: list[Path] = []
        created_directories: list[Path] = []
        submission_created = True
        try:
            with self.connect() as connection:
                connection.execute("BEGIN IMMEDIATE")
                self._ensure_not_migrating()
                existing_job_id = self._existing_idempotent_job(
                    connection,
                    idempotency_key=idempotency_key,
                    request_sha256=request_sha256,
                )
                if existing_job_id is not None:
                    job_id = existing_job_id
                    submission_created = False
                else:
                    cursor = connection.execute(
                        """
                        INSERT INTO jobs(
                            name, state, launch_token, cwd, argv_json, env_json,
                            env_override_keys, requested_cores, log_path, created_at,
                            retry_of, attempt, idempotency_key, request_sha256)
                        VALUES (?, 'queued', ?, ?, '[]', ?, ?, ?, '', ?, ?, ?, ?, ?)
                        """,
                        (
                            name,
                            secrets.token_hex(16),
                            str(cwd),
                            json.dumps(environment, ensure_ascii=False, sort_keys=True),
                            json.dumps(override_keys),
                            requested_cores,
                            created_at,
                            retry_of,
                            attempt,
                            idempotency_key,
                            request_sha256,
                        ),
                    )
                    job_id = int(cursor.lastrowid)
                    parent_log = log_dir() / f"{job_id:06d}.log"
                    artifact_dir = state_dir() / "ensembles" / f"{job_id:06d}"
                    members_dir = artifact_dir / "members"
                    for directory in (artifact_dir.parent, artifact_dir, members_dir):
                        if not directory.exists():
                            directory.mkdir(mode=0o700)
                            created_directories.append(directory)
                        directory.chmod(0o700)
                    for path in (
                        parent_log,
                        artifact_dir / "events.jsonl",
                        artifact_dir / "gpu_samples.jsonl",
                        artifact_dir / "process_samples.jsonl",
                    ):
                        with path.open("xb"):
                            pass
                        path.chmod(0o600)
                        created_files.append(path)

                    argv = [*runner_prefix, str(job_id)]
                    connection.execute(
                        """
                        UPDATE jobs SET argv_json = ?, log_path = ? WHERE id = ?
                        """,
                        (
                            json.dumps(argv, ensure_ascii=False),
                            str(parent_log),
                            job_id,
                        ),
                    )
                    connection.execute(
                        "UPDATE jobs SET queue_order = ? WHERE id = ?",
                        (job_id, job_id),
                    )
                    connection.execute(
                        """
                        INSERT INTO ensembles(
                            job_id, schema_version, spec_json, max_parallel,
                            gpu_device, gpu_budget_bytes, diagnostics_interval,
                            mps_mode, mps_pipe_directory, fail_fast,
                            allow_shared_cwd, artifact_dir)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            job_id,
                            schema_version,
                            json.dumps(canonical_spec, ensure_ascii=False, sort_keys=True),
                            max_parallel,
                            gpu_device,
                            gpu_budget_bytes,
                            diagnostics_interval,
                            mps_mode,
                            mps_pipe_directory,
                            int(fail_fast),
                            int(allow_shared_cwd),
                            str(artifact_dir),
                        ),
                    )
                    for index, member in enumerate(members):
                        member_log = members_dir / _safe_member_filename(
                            index, str(member["name"])
                        )
                        with member_log.open("xb"):
                            pass
                        member_log.chmod(0o600)
                        created_files.append(member_log)
                        connection.execute(
                            """
                            INSERT INTO ensemble_members(
                                job_id, member_index, name, state, cwd, argv_json,
                                env_json, gpu_admission_bytes, uvm_arena_bytes,
                                log_path, member_token, created_at)
                            VALUES (?, ?, ?, 'queued', ?, ?, ?, ?, ?, ?, ?, ?)
                            """,
                            (
                                job_id,
                                index,
                                str(member["name"]),
                                str(member["cwd"]),
                                json.dumps(member["argv"], ensure_ascii=False),
                                json.dumps(
                                    member.get("env", {}),
                                    ensure_ascii=False,
                                    sort_keys=True,
                                ),
                                int(member["gpu_admission_bytes"]),
                                int(member.get("uvm_arena_bytes", 0)),
                                str(member_log),
                                secrets.token_hex(16),
                                created_at,
                            ),
                        )
        except BaseException:
            for path in reversed(created_files):
                try:
                    path.unlink()
                except OSError:
                    pass
            for directory in reversed(created_directories):
                try:
                    directory.rmdir()
                except OSError:
                    pass
            raise
        job = self.require(job_id)
        if idempotency_key is not None:
            job["submission_created"] = submission_created
        return job

    def get(self, job_id: int) -> dict[str, Any] | None:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM jobs WHERE id = ?", (job_id,)
            ).fetchone()
        return self._row(row)

    def get_by_idempotency_key(self, idempotency_key: str) -> dict[str, Any] | None:
        key = _validate_idempotency_key(idempotency_key)
        assert key is not None
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM jobs WHERE idempotency_key = ?", (key,)
            ).fetchone()
        return self._row(row)

    def require(self, job_id: int) -> dict[str, Any]:
        job = self.get(job_id)
        if job is None:
            raise KeyError(f"no job with ID {job_id}")
        return job

    def get_ensemble(self, job_id: int) -> dict[str, Any] | None:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM ensembles WHERE job_id = ?", (job_id,)
            ).fetchone()
        return self._row(row)

    def require_ensemble(self, job_id: int) -> dict[str, Any]:
        self.require(job_id)
        ensemble = self.get_ensemble(job_id)
        if ensemble is None:
            raise ValueError(f"job {job_id} is not an ensemble")
        return ensemble

    @staticmethod
    def _ensemble_control_snapshot(
        connection: sqlite3.Connection, job_id: int
    ) -> dict[str, Any]:
        row = connection.execute(
            """
            SELECT j.state AS job_state, j.cancel_requested_at,
                   e.launch_paused, e.control_generation,
                   e.control_changed_at, e.drained_generation, e.drained_at,
                   e.supervisor_protocol, e.supervisor_pid,
                   SUM(CASE WHEN m.state = 'queued' THEN 1 ELSE 0 END) AS queued,
                   SUM(CASE WHEN m.state = 'running' THEN 1 ELSE 0 END) AS running,
                   SUM(CASE WHEN m.state = 'succeeded' THEN 1 ELSE 0 END) AS succeeded,
                   SUM(CASE WHEN m.state = 'failed' THEN 1 ELSE 0 END) AS failed,
                   SUM(CASE WHEN m.state = 'canceled' THEN 1 ELSE 0 END) AS canceled,
                   SUM(CASE WHEN m.state = 'interrupted' THEN 1 ELSE 0 END) AS interrupted
            FROM jobs AS j
            JOIN ensembles AS e ON e.job_id = j.id
            LEFT JOIN ensemble_members AS m ON m.job_id = j.id
            WHERE j.id = ?
            GROUP BY j.id
            """,
            (job_id,),
        ).fetchone()
        if row is None:
            job = connection.execute(
                "SELECT id FROM jobs WHERE id = ?", (job_id,)
            ).fetchone()
            if job is None:
                raise KeyError(f"no job with ID {job_id}")
            raise ValueError(f"job {job_id} is not an ensemble")
        snapshot = dict(row)
        snapshot["launch_paused"] = bool(snapshot["launch_paused"])
        for key in (
            "control_generation",
            "drained_generation",
            "supervisor_protocol",
            "supervisor_pid",
            "queued",
            "running",
            "succeeded",
            "failed",
            "canceled",
            "interrupted",
        ):
            if snapshot.get(key) is not None:
                snapshot[key] = int(snapshot[key])
        return snapshot

    def ensemble_control_snapshot(self, job_id: int) -> dict[str, Any]:
        with self.connect() as connection:
            return self._ensemble_control_snapshot(connection, job_id)

    def set_ensemble_paused(self, job_id: int, paused: bool) -> dict[str, Any]:
        """Persist member-launch control without signaling admitted members."""
        changed_at = utc_now()
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            self._ensure_not_migrating()
            row = connection.execute(
                """
                SELECT j.state, e.launch_paused, e.supervisor_protocol
                FROM jobs AS j
                LEFT JOIN ensembles AS e ON e.job_id = j.id
                WHERE j.id = ?
                """,
                (job_id,),
            ).fetchone()
            if row is None:
                connection.rollback()
                raise KeyError(f"no job with ID {job_id}")
            if row["launch_paused"] is None:
                connection.rollback()
                raise ValueError(f"job {job_id} is not an ensemble")
            state = str(row["state"])
            if state in TERMINAL_STATES or state == "cancelling":
                connection.rollback()
                raise ValueError(
                    f"cannot {'pause' if paused else 'resume'} ensemble {job_id} "
                    f"while its parent is {state}"
                )
            if (
                paused
                and state in ("starting", "running")
                and row["supervisor_protocol"] != ENSEMBLE_CONTROL_PROTOCOL
            ):
                connection.rollback()
                raise RuntimeError(
                    f"ensemble {job_id} is running under a supervisor that "
                    "predates member-level pause support; the existing outer "
                    "job must finish before this control can be used"
                )
            desired = int(bool(paused))
            if int(row["launch_paused"]) != desired:
                connection.execute(
                    """
                    UPDATE ensembles
                    SET launch_paused = ?,
                        control_generation = control_generation + 1,
                        control_changed_at = ?
                    WHERE job_id = ?
                    """,
                    (desired, changed_at, job_id),
                )
                if not paused and state == "queued":
                    connection.execute(
                        "UPDATE jobs SET message = NULL WHERE id = ?", (job_id,)
                    )
            # Return the state serialized with this control mutation. A second
            # connection here would let a concurrent opposite request change
            # the result between commit and the caller's response.
            snapshot = self._ensemble_control_snapshot(connection, job_id)
        return snapshot

    def register_ensemble_supervisor(
        self,
        job_id: int,
        *,
        pid: int,
        launch_token: str,
        protocol: int = ENSEMBLE_CONTROL_PROTOCOL,
    ) -> None:
        if protocol != ENSEMBLE_CONTROL_PROTOCOL:
            raise ValueError(
                f"unsupported ensemble launch-control protocol {protocol}"
            )
        if pid <= 0:
            raise ValueError("ensemble supervisor PID must be positive")
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            changed = connection.execute(
                """
                UPDATE ensembles
                SET supervisor_protocol = ?, supervisor_pid = ?
                WHERE job_id = ?
                  AND EXISTS (
                      SELECT 1 FROM jobs
                      WHERE id = ?
                        AND state IN ('starting', 'running')
                        AND cancel_requested_at IS NULL
                        AND launch_token = ?
                  )
                """,
                (protocol, pid, job_id, job_id, launch_token),
            ).rowcount
        if changed != 1:
            raise RuntimeError(
                f"ensemble {job_id} could not register its launch-control protocol"
            )

    def reserve_member_launch(
        self,
        job_id: int,
        member_index: int,
        *,
        lane: int,
        cpu: int,
    ) -> bool:
        """Linearize pause/cancel against a member launch before ``Popen``.

        ``False`` means launch admission is currently closed. A successful
        reservation is an admitted member even while its PID is briefly NULL.
        """
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            control = connection.execute(
                """
                SELECT j.state, j.cancel_requested_at, e.launch_paused,
                       e.drained_generation, e.supervisor_protocol,
                       e.supervisor_pid
                FROM jobs AS j JOIN ensembles AS e ON e.job_id = j.id
                WHERE j.id = ?
                """,
                (job_id,),
            ).fetchone()
            if control is None:
                connection.rollback()
                raise RuntimeError(f"ensemble {job_id} disappeared before launch")
            if (
                str(control["state"]) not in ("starting", "running")
                or control["cancel_requested_at"] is not None
                or int(control["launch_paused"]) == 1
                or control["drained_generation"] is not None
                or control["supervisor_protocol"] != ENSEMBLE_CONTROL_PROTOCOL
                or control["supervisor_pid"] is None
            ):
                connection.commit()
                return False
            changed = connection.execute(
                """
                UPDATE ensemble_members
                SET state = 'running', lane = ?, cpu = ?, pid = NULL, pgid = NULL,
                    started_at = ?, finished_at = NULL, exit_code = NULL,
                    message = NULL
                WHERE job_id = ? AND member_index = ? AND state = 'queued'
                """,
                (lane, cpu, utc_now(), job_id, member_index),
            ).rowcount
            if changed != 1:
                connection.rollback()
                raise RuntimeError(
                    f"ensemble {job_id} member {member_index} was not queued at reservation"
                )
        return True

    def record_member_process(
        self,
        job_id: int,
        member_index: int,
        *,
        pid: int,
        pgid: int,
    ) -> None:
        with self.connect() as connection:
            changed = connection.execute(
                """
                UPDATE ensemble_members
                SET pid = ?, pgid = ?
                WHERE job_id = ? AND member_index = ?
                  AND state = 'running' AND pid IS NULL
                """,
                (pid, pgid, job_id, member_index),
            ).rowcount
        if changed != 1:
            raise RuntimeError(
                f"ensemble {job_id} member {member_index} lost its launch reservation"
            )

    def mark_ensemble_drained(self, job_id: int, generation: int) -> bool:
        """Record a crash-recoverable drain marker for the active supervisor."""
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            control = connection.execute(
                """
                SELECT j.state, j.cancel_requested_at, e.launch_paused,
                       e.control_generation, e.supervisor_protocol,
                       e.supervisor_pid,
                       SUM(CASE WHEN m.state = 'running' THEN 1 ELSE 0 END) AS running,
                       SUM(CASE WHEN m.state = 'queued' THEN 1 ELSE 0 END) AS queued
                FROM jobs AS j
                JOIN ensembles AS e ON e.job_id = j.id
                LEFT JOIN ensemble_members AS m ON m.job_id = j.id
                WHERE j.id = ?
                GROUP BY j.id
                """,
                (job_id,),
            ).fetchone()
            if control is None:
                connection.rollback()
                raise RuntimeError(f"ensemble {job_id} disappeared while draining")
            if (
                str(control["state"]) not in ("starting", "running")
                or control["cancel_requested_at"] is not None
                or int(control["launch_paused"]) != 1
                or int(control["control_generation"]) != generation
                or control["supervisor_protocol"] != ENSEMBLE_CONTROL_PROTOCOL
                or control["supervisor_pid"] is None
                or int(control["running"] or 0) != 0
                or int(control["queued"] or 0) == 0
            ):
                connection.commit()
                return False
            connection.execute(
                """
                UPDATE ensembles
                SET drained_generation = ?, drained_at = ?
                WHERE job_id = ?
                """,
                (generation, utc_now(), job_id),
            )
        return True

    def requeue_drained_ensemble(self, job_id: int) -> bool:
        """Suspend a process-tree-quiet drained parent without touching members."""
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            control = connection.execute(
                """
                SELECT j.state, j.cancel_requested_at, e.launch_paused,
                       e.drained_generation,
                       SUM(CASE WHEN m.state = 'running' THEN 1 ELSE 0 END) AS running,
                       SUM(CASE WHEN m.state = 'queued' THEN 1 ELSE 0 END) AS queued
                FROM jobs AS j
                JOIN ensembles AS e ON e.job_id = j.id
                LEFT JOIN ensemble_members AS m ON m.job_id = j.id
                WHERE j.id = ?
                GROUP BY j.id
                """,
                (job_id,),
            ).fetchone()
            if control is None:
                connection.rollback()
                return False
            if (
                str(control["state"]) not in ACTIVE_STATES
                or control["cancel_requested_at"] is not None
                or control["drained_generation"] is None
                or int(control["running"] or 0) != 0
                or int(control["queued"] or 0) == 0
            ):
                connection.commit()
                return False
            paused = int(control["launch_paused"]) == 1
            changed = connection.execute(
                """
                UPDATE jobs
                SET state = 'queued', cpu_list = NULL, started_at = NULL,
                    finished_at = NULL, exit_code = NULL, pid = NULL, pgid = NULL,
                    boot_id = NULL, process_start_ticks = NULL,
                    queue_order = (
                        SELECT COALESCE(MIN(queued.queue_order), 1) - 1
                        FROM jobs AS queued
                        WHERE queued.state = 'queued'
                    ),
                    message = ?
                WHERE id = ? AND state IN ('starting', 'running')
                  AND cancel_requested_at IS NULL
                """,
                (
                    "ensemble paused after draining admitted members"
                    if paused
                    else "ensemble drained during a concurrent resume; queued to continue",
                    job_id,
                ),
            ).rowcount
            if changed == 1:
                connection.execute(
                    """
                    UPDATE ensembles
                    SET supervisor_protocol = NULL, supervisor_pid = NULL
                    WHERE job_id = ?
                    """,
                    (job_id,),
                )
        return changed == 1

    def list_ensemble_members(self, job_id: int) -> list[dict[str, Any]]:
        self.require_ensemble(job_id)
        with self.connect() as connection:
            rows = connection.execute(
                """
                SELECT * FROM ensemble_members
                WHERE job_id = ? ORDER BY member_index
                """,
                (job_id,),
            ).fetchall()
        return [dict(row) for row in rows]

    def require_ensemble_member(
        self, job_id: int, selector: str | int
    ) -> dict[str, Any]:
        self.require_ensemble(job_id)
        with self.connect() as connection:
            if isinstance(selector, int) or (
                isinstance(selector, str) and selector.isdecimal()
            ):
                row = connection.execute(
                    """
                    SELECT * FROM ensemble_members
                    WHERE job_id = ? AND member_index = ?
                    """,
                    (job_id, int(selector)),
                ).fetchone()
            else:
                row = connection.execute(
                    """
                    SELECT * FROM ensemble_members
                    WHERE job_id = ? AND name = ?
                    """,
                    (job_id, str(selector)),
                ).fetchone()
        if row is None:
            raise KeyError(f"ensemble {job_id} has no member {selector!r}")
        return dict(row)

    def mark_member_running(
        self,
        job_id: int,
        member_index: int,
        *,
        lane: int,
        cpu: int,
        pid: int,
        pgid: int,
    ) -> None:
        with self.connect() as connection:
            changed = connection.execute(
                """
                UPDATE ensemble_members
                SET state = 'running', lane = ?, cpu = ?, pid = ?, pgid = ?,
                    started_at = ?, finished_at = NULL, exit_code = NULL,
                    message = NULL
                WHERE job_id = ? AND member_index = ? AND state = 'queued'
                """,
                (lane, cpu, pid, pgid, utc_now(), job_id, member_index),
            ).rowcount
        if changed != 1:
            raise RuntimeError(
                f"ensemble {job_id} member {member_index} was not queued at launch"
            )

    def finish_member(
        self,
        job_id: int,
        member_index: int,
        *,
        state: str,
        exit_code: int | None,
        message: str | None = None,
        uvm_marker: dict[str, Any] | None = None,
    ) -> None:
        if state not in MEMBER_TERMINAL_STATES:
            raise ValueError(f"invalid member terminal state: {state}")
        fields: dict[str, Any] = {}
        if uvm_marker is not None:
            fields = {
                "uvm_policy": uvm_marker.get("policy"),
                "uvm_limit_bytes": uvm_marker.get("limit_bytes"),
                "uvm_peak_live_bytes": uvm_marker.get("peak_live_bytes"),
                "uvm_peak_reserved_bytes": uvm_marker.get("peak_reserved_bytes"),
                "uvm_pool_reclaims": uvm_marker.get("pool_reclaims"),
                "uvm_teardown_ok": uvm_marker.get("teardown_ok"),
                "uvm_release": uvm_marker.get("release"),
                "uvm_marker_json": json.dumps(
                    uvm_marker, ensure_ascii=False, sort_keys=True
                ),
            }
        with self.connect() as connection:
            connection.execute(
                """
                UPDATE ensemble_members
                SET state = ?, finished_at = ?, exit_code = ?, message = ?,
                    uvm_policy = COALESCE(?, uvm_policy),
                    uvm_limit_bytes = COALESCE(?, uvm_limit_bytes),
                    uvm_peak_live_bytes = COALESCE(?, uvm_peak_live_bytes),
                    uvm_peak_reserved_bytes = COALESCE(?, uvm_peak_reserved_bytes),
                    uvm_pool_reclaims = COALESCE(?, uvm_pool_reclaims),
                    uvm_teardown_ok = COALESCE(?, uvm_teardown_ok),
                    uvm_release = COALESCE(?, uvm_release),
                    uvm_marker_json = COALESCE(?, uvm_marker_json)
                WHERE job_id = ? AND member_index = ?
                  AND state IN ('queued', 'running')
                """,
                (
                    state,
                    utc_now(),
                    exit_code,
                    message,
                    fields.get("uvm_policy"),
                    fields.get("uvm_limit_bytes"),
                    fields.get("uvm_peak_live_bytes"),
                    fields.get("uvm_peak_reserved_bytes"),
                    fields.get("uvm_pool_reclaims"),
                    fields.get("uvm_teardown_ok"),
                    fields.get("uvm_release"),
                    fields.get("uvm_marker_json"),
                    job_id,
                    member_index,
                ),
            )

    def update_member_sample(
        self,
        job_id: int,
        member_index: int,
        *,
        rss_bytes: int | None,
        threads: int | None,
    ) -> None:
        with self.connect() as connection:
            connection.execute(
                """
                UPDATE ensemble_members
                SET last_sample_at = ?, last_rss_bytes = ?,
                    peak_rss_bytes = CASE
                        WHEN ? IS NULL THEN peak_rss_bytes
                        WHEN peak_rss_bytes IS NULL OR ? > peak_rss_bytes THEN ?
                        ELSE peak_rss_bytes END,
                    last_threads = ?
                WHERE job_id = ? AND member_index = ? AND state = 'running'
                """,
                (
                    utc_now(),
                    rss_bytes,
                    rss_bytes,
                    rss_bytes,
                    rss_bytes,
                    threads,
                    job_id,
                    member_index,
                ),
            )

    def finish_unstarted_members(
        self,
        job_id: int,
        *,
        state: str,
        message: str,
    ) -> None:
        if state not in MEMBER_TERMINAL_STATES:
            raise ValueError(f"invalid member terminal state: {state}")
        with self.connect() as connection:
            connection.execute(
                """
                UPDATE ensemble_members
                SET state = ?, finished_at = ?, message = ?
                WHERE job_id = ? AND state = 'queued'
                """,
                (state, utc_now(), message, job_id),
            )

    @staticmethod
    def _finish_dangling_members(
        connection: sqlite3.Connection,
        job_id: int,
        *,
        parent_state: str,
        message: str,
    ) -> None:
        member_state = (
            "canceled"
            if parent_state == "canceled"
            else "interrupted"
            if parent_state in ("interrupted", "failed")
            else "failed"
        )
        connection.execute(
            """
            UPDATE ensemble_members
            SET state = ?, finished_at = ?, message = COALESCE(message, ?)
            WHERE job_id = ? AND state IN ('queued', 'running')
            """,
            (member_state, utc_now(), message, job_id),
        )

    def list_jobs(self, *, include_all: bool = False, limit: int = 20) -> list[dict[str, Any]]:
        with self.connect() as connection:
            if include_all:
                rows = connection.execute(
                    "SELECT * FROM jobs ORDER BY id DESC"
                ).fetchall()
            else:
                rows = connection.execute(
                    """
                    SELECT * FROM jobs
                    ORDER BY
                        CASE state
                            WHEN 'running' THEN 0
                            WHEN 'cancelling' THEN 0
                            WHEN 'starting' THEN 0
                            WHEN 'queued' THEN 1
                            ELSE 2
                        END,
                        CASE WHEN state = 'queued' THEN queue_order ELSE 0 END,
                        CASE WHEN state IN ('starting', 'running', 'cancelling')
                             THEN id
                             WHEN state = 'queued' THEN id
                             ELSE -id END
                    LIMIT ?
                    """,
                    (limit,),
                ).fetchall()
        return [dict(row) for row in rows]

    def _move_queued_job(self, job_id: int, *, direction: str) -> dict[str, Any]:
        if direction not in ("front", "back"):
            raise ValueError(f"invalid queue move direction: {direction}")
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            self._ensure_not_migrating()
            target = connection.execute(
                "SELECT * FROM jobs WHERE id = ?", (job_id,)
            ).fetchone()
            if target is None:
                connection.rollback()
                raise KeyError(f"no job with ID {job_id}")
            state = str(target["state"])
            if state != "queued":
                connection.rollback()
                raise ValueError(
                    f"job {job_id} is {state}; only queued jobs can be "
                    f"{'promoted' if direction == 'front' else 'demoted'}"
                )
            rows = connection.execute(
                """
                SELECT id FROM jobs
                WHERE state = 'queued'
                ORDER BY queue_order, id
                """
            ).fetchall()
            ordered_ids = [int(row["id"]) for row in rows]
            old_position = ordered_ids.index(job_id)
            new_position = 0 if direction == "front" else len(ordered_ids) - 1
            changed = old_position != new_position
            if changed:
                ordered_ids.pop(old_position)
                ordered_ids.insert(new_position, job_id)
                connection.executemany(
                    "UPDATE jobs SET queue_order = ? WHERE id = ? AND state = 'queued'",
                    [
                        (position, queued_id)
                        for position, queued_id in enumerate(ordered_ids, start=1)
                    ],
                )
            stored = connection.execute(
                "SELECT * FROM jobs WHERE id = ?", (job_id,)
            ).fetchone()
            assert stored is not None
            result = {
                "job_id": job_id,
                "changed": changed,
                "direction": direction,
                "queue_position": new_position + 1,
                "queued_count": len(ordered_ids),
                "job": dict(stored),
            }
        return result

    def promote(self, job_id: int) -> dict[str, Any]:
        """Move one queued outer job to the dispatch front."""
        return self._move_queued_job(job_id, direction="front")

    def demote(self, job_id: int) -> dict[str, Any]:
        """Move one queued outer job to the dispatch back."""
        return self._move_queued_job(job_id, direction="back")

    def claim_next(self, pcores: tuple[int, ...]) -> dict[str, Any] | None:
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            if MIGRATION_SENTINEL.exists():
                connection.commit()
                return None
            paused = connection.execute(
                "SELECT value FROM queue_meta WHERE key = 'paused'"
            ).fetchone()
            if paused is not None and paused["value"] == "1":
                connection.commit()
                return None
            # Inspect the true dispatch head before considering launch control.
            # A paused ensemble holds its explicit queue position; dispatch
            # never skips around it implicitly.
            row = connection.execute(
                """
                SELECT j.*, e.job_id AS ensemble_job_id,
                       e.launch_paused AS ensemble_launch_paused
                FROM jobs AS j
                LEFT JOIN ensembles AS e ON e.job_id = j.id
                WHERE j.state = 'queued'
                ORDER BY j.queue_order, j.id
                LIMIT 1
                """
            ).fetchone()
            if row is None:
                connection.commit()
                return None
            if row["ensemble_launch_paused"] is not None and int(
                row["ensemble_launch_paused"]
            ) == 1:
                connection.commit()
                return None
            requested = int(row["requested_cores"])
            cpu_list = format_cpu_list(pcores[:requested])
            changed = connection.execute(
                """
                UPDATE jobs
                SET state = 'starting', started_at = ?, cpu_list = ?, message = NULL
                WHERE id = ? AND state = 'queued'
                """,
                (utc_now(), cpu_list, row["id"]),
            ).rowcount
            if changed == 1 and row["ensemble_job_id"] is not None:
                # A resumed allocation always starts with a clean drain marker.
                # Advertising the protocol at claim time closes the brief
                # starting-before-supervisor-registration window without
                # making a legacy active supervisor pause-capable during an
                # in-place schema migration.
                connection.execute(
                    """
                    UPDATE ensembles
                    SET drained_generation = NULL, drained_at = NULL,
                        supervisor_protocol = ?, supervisor_pid = NULL
                    WHERE job_id = ? AND launch_paused = 0
                    """,
                    (ENSEMBLE_CONTROL_PROTOCOL, row["id"]),
                )
            connection.commit()
        if changed != 1:
            return None
        return self.require(int(row["id"]))

    def dispatch_blocker(self) -> dict[str, Any] | None:
        """Describe an intentional reason the dispatch head cannot be claimed.

        This is a status/autostart hint only; :meth:`claim_next` repeats every
        check under its own write transaction and remains authoritative.
        """
        if MIGRATION_SENTINEL.exists():
            return {"kind": "migration"}
        with self.connect() as connection:
            paused = connection.execute(
                "SELECT value FROM queue_meta WHERE key = 'paused'"
            ).fetchone()
            if paused is not None and paused["value"] == "1":
                return {"kind": "queue-paused"}
            row = connection.execute(
                """
                SELECT j.id, j.name, e.launch_paused
                FROM jobs AS j
                LEFT JOIN ensembles AS e ON e.job_id = j.id
                WHERE j.state = 'queued'
                ORDER BY j.queue_order, j.id
                LIMIT 1
                """
            ).fetchone()
        if row is not None and row["launch_paused"] is not None and int(
            row["launch_paused"]
        ) == 1:
            return {
                "kind": "ensemble-paused",
                "job_id": int(row["id"]),
                "name": str(row["name"]),
            }
        return None

    def mark_running(
        self,
        job_id: int,
        *,
        pid: int,
        pgid: int,
        boot_id: str | None,
        process_start_ticks: int | None,
    ) -> None:
        with self.connect() as connection:
            connection.execute(
                """
                UPDATE jobs
                SET state = CASE WHEN cancel_requested_at IS NULL
                                 THEN 'running' ELSE 'cancelling' END,
                    pid = ?, pgid = ?, boot_id = ?, process_start_ticks = ?
                WHERE id = ? AND state IN ('starting', 'cancelling')
                """,
                (pid, pgid, boot_id, process_start_ticks, job_id),
            )

    def cancellation_requested(self, job_id: int) -> bool:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT cancel_requested_at FROM jobs WHERE id = ?", (job_id,)
            ).fetchone()
        return bool(row and row["cancel_requested_at"])

    def finish(
        self,
        job_id: int,
        *,
        state: str,
        exit_code: int | None,
        message: str | None = None,
    ) -> None:
        if state not in TERMINAL_STATES:
            raise ValueError(f"invalid terminal state: {state}")
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            changed = connection.execute(
                """
                UPDATE jobs
                SET state = ?, finished_at = ?, exit_code = ?, message = ?
                WHERE id = ? AND state IN ('starting', 'running', 'cancelling')
                """,
                (state, utc_now(), exit_code, message, job_id),
            ).rowcount
            if changed:
                self._finish_dangling_members(
                    connection,
                    job_id,
                    parent_state=state,
                    message=f"ensemble parent finished {state} before member completion",
                )

    def cancel(self, job_id: int) -> str:
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            self._ensure_not_migrating()
            row = connection.execute(
                "SELECT state FROM jobs WHERE id = ?", (job_id,)
            ).fetchone()
            if row is None:
                connection.rollback()
                raise KeyError(f"no job with ID {job_id}")
            state = str(row["state"])
            now = utc_now()
            if state == "queued":
                connection.execute(
                    """
                    UPDATE jobs
                    SET state = 'canceled', cancel_requested_at = ?,
                        finished_at = ?, message = 'canceled before launch'
                    WHERE id = ?
                    """,
                    (now, now, job_id),
                )
                self._finish_dangling_members(
                    connection,
                    job_id,
                    parent_state="canceled",
                    message="canceled before ensemble launch",
                )
                result = "canceled"
            elif state in ACTIVE_STATES:
                connection.execute(
                    """
                    UPDATE jobs
                    SET state = 'cancelling',
                        cancel_requested_at = COALESCE(cancel_requested_at, ?)
                    WHERE id = ?
                    """,
                    (now, job_id),
                )
                result = "cancelling"
            else:
                result = state
            connection.commit()
        return result

    def active_jobs(self) -> list[dict[str, Any]]:
        placeholders = ",".join("?" for _ in ACTIVE_STATES)
        with self.connect() as connection:
            rows = connection.execute(
                f"SELECT * FROM jobs WHERE state IN ({placeholders}) ORDER BY id",
                ACTIVE_STATES,
            ).fetchall()
        return [dict(row) for row in rows]

    def queued_jobs(self) -> list[dict[str, Any]]:
        """Return waiting jobs in explicit dispatch order."""
        with self.connect() as connection:
            rows = connection.execute(
                """
                SELECT * FROM jobs
                WHERE state = 'queued'
                ORDER BY queue_order, id
                """
            ).fetchall()
        return [dict(row) for row in rows]

    def mark_interrupted(self, job_id: int, message: str) -> None:
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            changed = connection.execute(
                """
                UPDATE jobs
                SET state = 'interrupted', finished_at = ?, message = ?
                WHERE id = ? AND state IN ('starting', 'running', 'cancelling')
                """,
                (utc_now(), message, job_id),
            ).rowcount
            if changed:
                self._finish_dangling_members(
                    connection,
                    job_id,
                    parent_state="interrupted",
                    message="ensemble interrupted with its parent job",
                )

    def retry(self, job_id: int) -> dict[str, Any]:
        original = self.require(job_id)
        if original["state"] not in TERMINAL_STATES:
            raise ValueError(
                f"job {job_id} is {original['state']}; only finished jobs can be retried"
            )
        ensemble = self.get_ensemble(job_id)
        if ensemble is not None:
            members = self.list_ensemble_members(job_id)
            original_argv = list(json.loads(original["argv_json"]))
            if not original_argv or original_argv[-1] != str(job_id):
                raise RuntimeError(
                    f"ensemble {job_id} has an invalid stored supervisor command"
                )
            return self.submit_ensemble(
                name=str(original["name"]),
                cwd=Path(str(original["cwd"])),
                environment=dict(json.loads(original["env_json"])),
                env_override_keys=list(json.loads(original["env_override_keys"])),
                requested_cores=int(original["requested_cores"]),
                runner_prefix=original_argv[:-1],
                schema_version=str(ensemble["schema_version"]),
                canonical_spec=dict(json.loads(ensemble["spec_json"])),
                members=[
                    {
                        "name": member["name"],
                        "cwd": member["cwd"],
                        "argv": list(json.loads(member["argv_json"])),
                        "env": dict(json.loads(member["env_json"])),
                        "gpu_admission_bytes": member["gpu_admission_bytes"],
                        "uvm_arena_bytes": member["uvm_arena_bytes"],
                    }
                    for member in members
                ],
                max_parallel=int(ensemble["max_parallel"]),
                gpu_device=str(ensemble["gpu_device"]),
                gpu_budget_bytes=int(ensemble["gpu_budget_bytes"]),
                diagnostics_interval=float(ensemble["diagnostics_interval"]),
                mps_mode=str(ensemble["mps_mode"]),
                mps_pipe_directory=ensemble["mps_pipe_directory"],
                fail_fast=bool(ensemble["fail_fast"]),
                allow_shared_cwd=bool(ensemble["allow_shared_cwd"]),
                retry_of=job_id,
                attempt=int(original["attempt"]) + 1,
            )
        return self.submit(
            name=str(original["name"]),
            cwd=Path(str(original["cwd"])),
            argv=list(json.loads(original["argv_json"])),
            environment=dict(json.loads(original["env_json"])),
            env_override_keys=list(json.loads(original["env_override_keys"])),
            requested_cores=int(original["requested_cores"]),
            retry_of=job_id,
            attempt=int(original["attempt"]) + 1,
        )

    def set_paused(self, paused: bool) -> None:
        with self.connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            self._ensure_not_migrating()
            connection.execute(
                """
                INSERT INTO queue_meta(key, value) VALUES ('paused', ?)
                ON CONFLICT(key) DO UPDATE SET value = excluded.value
                """,
                ("1" if paused else "0",),
            )

    def is_paused(self) -> bool:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT value FROM queue_meta WHERE key = 'paused'"
            ).fetchone()
        return bool(row and row["value"] == "1")

    def unfinished_ids(self) -> list[int]:
        with self.connect() as connection:
            rows = connection.execute(
                """
                SELECT id FROM jobs
                WHERE state IN ('queued', 'starting', 'running', 'cancelling')
                ORDER BY id
                """
            ).fetchall()
        return [int(row["id"]) for row in rows]
