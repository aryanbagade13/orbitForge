"""Launch orbitForge in a native desktop window, without a browser or web server."""
import argparse
import json
from pathlib import Path

from server import ROOT, load_data


def make_document(directory):
    """Bundle local styles, scripts and saved results into an offline document."""
    payload = json.dumps(load_data(directory), allow_nan=False).replace("</", "<\\/")
    document = (ROOT / "index.html").read_text()
    document = document.replace('<link rel="stylesheet" href="/style.css">',
                                "<style>" + (ROOT / "style.css").read_text() + "</style>")
    document = document.replace('<script src="/app.js" defer></script>', "")
    document = document.replace('href="/" aria-label="orbitForge home"', 'href="#" aria-label="orbitForge home"')
    return document.replace("</body>", "<script>window.orbitForgeData=" + payload + ";</script><script>" +
                            (ROOT / "app.js").read_text() + "</script></body>")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT.parent / "outputs/interplanetary_baseline")
    args = parser.parse_args()
    try:
        import webview
    except ImportError:
        parser.exit(1, "Desktop dependency missing. Install with: python -m pip install -r ui/requirements.txt\n")
    try:
        document = make_document(args.data_dir)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Cannot load saved baseline: {error}\nRun the baseline experiment first.\n")
    webview.create_window("orbitForge — Trajectory viewer", html=document,
                          width=1440, height=960, min_size=(900, 650), background_color="#0c1012")
    webview.start()


if __name__ == "__main__":
    main()
