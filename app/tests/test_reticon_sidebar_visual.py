"""Visual regression of the navigation sidebar (Standard level, collapsed categories).

Runs only when RETICON_VISUAL=1 (its own CI job; it needs Playwright and an installed Microsoft Edge).
It starts the app on a free port, captures the left 300 px with Playwright (msedge channel), and compares
the picture with tests/golden/sidebar_visual/sidebar_standard.png.

Baseline rule: the baseline changes only with RETICON_VISUAL_UPDATE=1, and the commit that changes it must
say why. A changed baseline is a visible UI change and needs the coordinator's and the user's approval.
Tolerance: at most MAX_DIFF_SHARE of the pixels may differ by more than CHANNEL_TOL per channel.
These two values are a first setting; calibrate them on the CI machine before relying on the test.
"""
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parents[1]
BASELINE = Path(__file__).resolve().parent / "golden" / "sidebar_visual" / "sidebar_standard.png"
CHANNEL_TOL = 8
MAX_DIFF_SHARE = 0.005
CLIP = {"x": 0, "y": 0, "width": 300, "height": 900}

pytestmark = pytest.mark.skipif(os.environ.get("RETICON_VISUAL") != "1",
                                reason="visual sidebar test runs only in its own job (set RETICON_VISUAL=1)")


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _capture(out: Path) -> None:
    from playwright.sync_api import sync_playwright
    port = _free_port()
    env = dict(os.environ, MFS_NO_UPDATE_CHECK="1", PYTHONUTF8="1", MFS_DEFAULT_LEVEL="Standard")
    server = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", str(APP / "app.py"), "--server.headless", "true",
         "--server.port", str(port), "--browser.gatherUsageStats", "false"],
        cwd=str(APP), env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        import urllib.request
        deadline = time.time() + 120
        while time.time() < deadline:
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=5)
                break
            except Exception:
                time.sleep(1)
        with sync_playwright() as p:
            browser = p.chromium.launch(channel="msedge", headless=True)
            page = browser.new_page(viewport={"width": 1280, "height": 900})
            page.goto(f"http://127.0.0.1:{port}/", wait_until="domcontentloaded", timeout=180000)
            page.wait_for_function(
                "() => document.querySelectorAll('section[data-testid=\"stSidebar\"] button').length > 5",
                timeout=240000)
            page.wait_for_timeout(3000)
            page.screenshot(path=str(out), clip=CLIP)
            browser.close()
    finally:
        server.terminate()
        try:
            server.wait(timeout=20)
        except subprocess.TimeoutExpired:
            server.kill()


def _diff_share(a: Path, b: Path) -> float:
    from PIL import Image
    ia, ib = Image.open(a).convert("RGB"), Image.open(b).convert("RGB")
    assert ia.size == ib.size, f"size differs: {ia.size} vs {ib.size}"
    da, db = ia.tobytes(), ib.tobytes()
    px = ia.size[0] * ia.size[1]
    differing = 0
    for i in range(0, len(da), 3):
        if max(abs(da[i + k] - db[i + k]) for k in range(3)) > CHANNEL_TOL:
            differing += 1
    return differing / px


def test_sidebar_matches_baseline(tmp_path):
    out = tmp_path / "sidebar_now.png"
    _capture(out)
    if os.environ.get("RETICON_VISUAL_UPDATE") == "1":
        BASELINE.parent.mkdir(parents=True, exist_ok=True)
        BASELINE.write_bytes(out.read_bytes())
        pytest.skip("baseline written; commit it with the reason")
    assert BASELINE.is_file(), "no baseline: run once with RETICON_VISUAL_UPDATE=1 and review it"
    share = _diff_share(out, BASELINE)
    assert share <= MAX_DIFF_SHARE, f"sidebar differs from baseline: {share:.4%} of pixels (max {MAX_DIFF_SHARE:.2%})"
