from __future__ import annotations

import json
import os
from http.server import BaseHTTPRequestHandler, HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path


class Handler(SimpleHTTPRequestHandler):
    """Serve the interactive map and keep a machine-readable health endpoint."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(Path(__file__).parents[2] / "docs"), **kwargs)

    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/healthz":
            body = {"service": "cannabis-carbon", "status": "ok", "version": "0.1.0"}
            payload = (json.dumps(body) + "\n").encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        super().do_GET()

    def log_message(self, fmt: str, *args: object) -> None:
        print(json.dumps({"message": fmt % args}))


def main() -> None:
    port = int(os.environ.get("PORT", "8080"))
    HTTPServer(("0.0.0.0", port), Handler).serve_forever()


if __name__ == "__main__":
    main()
