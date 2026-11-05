"""A local stand-in for the GitHub REST API (SimQueue 0.7.0 tests).

Serves the endpoints the worker uses: commit statuses, issue comments
(list, create, edit).  Every request is recorded with its method, path,
body and whether the Authorization header carried the expected token.
``down`` makes it answer 503 (GitHub unreachable), ``fail_status`` a given
HTTP code for statuses only.

Run standalone: ``fake_github.py PORTFILE LOGFILE TOKEN`` (writes the port).
"""

from __future__ import annotations

import json
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any


class FakeGitHub:
    def __init__(self, token: str, log_path: str | None = None) -> None:
        self.token = token
        self.log_path = log_path
        self.requests: list[dict[str, Any]] = []
        self.statuses: list[dict[str, Any]] = []
        self.comments: dict[int, dict[str, Any]] = {}
        self.next_comment = 1000
        self.down = False
        self.fail_status: int | None = None
        self.lock = threading.Lock()
        fake = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args: Any) -> None:
                pass

            def _reply(self, code: int, body: Any) -> None:
                data = json.dumps(body).encode()
                self.send_response(code)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)

            def _handle(self, method: str) -> None:
                length = int(self.headers.get("Content-Length") or 0)
                raw = self.rfile.read(length) if length else b""
                body = json.loads(raw) if raw else None
                auth = self.headers.get("Authorization") or ""
                record = {"method": method, "path": self.path, "body": body,
                          "auth_ok": auth == f"Bearer {fake.token}", "user_agent": self.headers.get("User-Agent")}
                with fake.lock:
                    fake.requests.append(record)
                    if fake.log_path:
                        with open(fake.log_path, "a") as handle:
                            handle.write(json.dumps(record) + "\n")
                    if fake.down:
                        return self._reply(503, {"message": "Service Unavailable"})
                    if not record["auth_ok"]:
                        return self._reply(401, {"message": "Bad credentials"})
                    parts = self.path.split("?")[0].strip("/").split("/")
                    if method == "POST" and len(parts) == 5 and parts[0] == "repos" and parts[3] == "statuses":
                        if fake.fail_status:
                            return self._reply(fake.fail_status, {"message": "Validation Failed"})
                        fake.statuses.append({"repo": f"{parts[1]}/{parts[2]}", "sha": parts[4], **(body or {})})
                        return self._reply(201, {"id": len(fake.statuses), "state": (body or {}).get("state")})
                    if len(parts) == 6 and parts[0] == "repos" and parts[3] == "issues" and parts[5] == "comments":
                        pr = int(parts[4])
                        if method == "GET":
                            items = [c for c in fake.comments.values() if c["pr"] == pr]
                            return self._reply(200, [{"id": c["id"], "body": c["body"]} for c in items])
                        if method == "POST":
                            fake.next_comment += 1
                            fake.comments[fake.next_comment] = {"id": fake.next_comment, "pr": pr,
                                                                "body": body["body"], "edits": 0}
                            return self._reply(201, {"id": fake.next_comment})
                    if (method == "PATCH" and len(parts) == 6 and parts[3] == "issues" and parts[4] == "comments"):
                        comment = fake.comments.get(int(parts[5]))
                        if comment is None:
                            return self._reply(404, {"message": "Not Found"})
                        comment["body"] = body["body"]
                        comment["edits"] += 1
                        return self._reply(200, {"id": comment["id"]})
                    return self._reply(404, {"message": "Not Found"})

            def do_GET(self) -> None:  # noqa: N802
                self._handle("GET")

            def do_POST(self) -> None:  # noqa: N802
                self._handle("POST")

            def do_PATCH(self) -> None:  # noqa: N802
                self._handle("PATCH")

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.port = self.server.server_address[1]
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)

    @property
    def api_base(self) -> str:
        return f"http://127.0.0.1:{self.port}"

    def start(self) -> "FakeGitHub":
        self.thread.start()
        return self

    def stop(self) -> None:
        self.server.shutdown()
        self.server.server_close()

    def statuses_for(self, context: str) -> list[dict[str, Any]]:
        with self.lock:
            return [s for s in self.statuses if s.get("context") == context]


if __name__ == "__main__":
    port_file, log_file, token = sys.argv[1:4]
    fake = FakeGitHub(token, log_file)
    with open(port_file, "w") as handle:
        handle.write(str(fake.port))
    fake.server.serve_forever()
