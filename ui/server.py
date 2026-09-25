"""Serve the saved orbitForge baseline. Run with the project's Python environment."""
import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
AU_KM = 149597870.7


def load_data(directory):
    """Validate aligned samples, then convert positions to Sun-relative AU."""
    report = json.loads((directory / "summary.json").read_text())
    with np.load(directory / "trajectory_samples.npz", allow_pickle=False) as archive:
        times = archive["times_s"]
        if times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all() or not np.all(np.diff(times) > 0):
            raise ValueError("times_s must contain at least two finite, increasing timestamps")
        arrays = {}
        for name, width in {"no_burn": 6, "correction": 6, "sun_positions_km": 3,
                            "earth_positions_km": 3, "jupiter_positions_km": 3,
                            "saturn_positions_km": 3}.items():
            value = archive[name]
            if value.shape != (len(times), width) or not np.isfinite(value).all():
                raise ValueError(f"{name} must be finite and have shape ({len(times)}, {width})")
            arrays[name] = value
        sun = arrays["sun_positions_km"]
        positions = {key: ((value[:, :3] - sun) / AU_KM).tolist() for key, value in arrays.items()}
        # Speeds retain the saved barycentric frame; no inferred Sun velocity.
        speeds = {key: np.linalg.norm(arrays[key][:, 3:], axis=1).tolist() for key in ("no_burn", "correction")}
    return {"times_s": times.tolist(), "positions": positions, "speeds_km_s": speeds, "report": report}


def make_handler(payload):
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            routes = {"/": ("index.html", "text/html"), "/app.js": ("app.js", "text/javascript"),
                      "/style.css": ("style.css", "text/css")}
            if self.path == "/api/trajectory":
                body, mime = payload, "application/json"
            elif self.path in routes:
                filename, mime = routes[self.path]
                body = (ROOT / filename).read_bytes()
            else:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Type", mime + "; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(body)
    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT.parent / "outputs/interplanetary_baseline")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    try:
        data = load_data(args.data_dir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Cannot load saved baseline: {error}\nRun python -m experiments.interplanetary_baseline --no-show first.\n")
    payload = json.dumps(data, allow_nan=False, separators=(",", ":")).encode()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(payload))
    print(f"orbitForge viewer → http://127.0.0.1:{server.server_port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
