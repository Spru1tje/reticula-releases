"""v2.11 update check: asks once, background fetch, banner Download / Later / Skip, silent offline, right asset per OS."""
import os, sys, json, time, threading, http.server, socketserver
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import pytest
from mfs import update as U

REL = {"tag_name": "v9.1.0", "html_url": "https://github.com/Spru1tje/reticula-releases/releases/tag/v9.1.0", "body": "* faster",
       "assets": [{"name": "Reticula_Setup-1.bin", "browser_download_url": "https://x/1.bin", "size": 1},
                  {"name": "Reticula_Setup_9.1.0.exe", "browser_download_url": "https://x/Reticula_Setup_9.1.0.exe", "size": 180_000_000},
                  {"name": "Reticula_9.1.0_macOS.dmg", "browser_download_url": "https://x/R.dmg", "size": 2}]}


def test_versions_and_assets():
    assert U.is_newer("v2.11.1", "2.11.0") and not U.is_newer("v2.11.0", "2.11.0") and U.is_newer("2.12", "2.11.9")
    assert U.pick_asset(REL, "Windows")["name"] == "Reticula_Setup_9.1.0.exe"
    assert U.pick_asset(REL, "Darwin")["name"].endswith(".dmg")
    assert U.pick_asset(REL, "Linux") is None
    s = U.summarize(REL, "2.11.0", "Windows")
    assert s["newer"] and s["version"] == "9.1.0" and s["url"].endswith(".exe")
    assert U.summarize(dict(REL, prerelease=True), "2.11.0") is None
    assert U.DEFAULT_URL == "https://api.github.com/repos/Spru1tje/reticula-releases/releases/latest"


@pytest.fixture
def feed():
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            b = json.dumps(REL).encode(); self.send_response(200); self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(b))); self.end_headers(); self.wfile.write(b)
        def log_message(self, *a):
            pass
    srv = socketserver.TCPServer(("127.0.0.1", 0), H); t = threading.Thread(target=srv.serve_forever, daemon=True); t.start()
    yield f"http://127.0.0.1:{srv.server_address[1]}/releases/latest"
    srv.shutdown()


def _app(tmp_path, monkeypatch, url):
    from streamlit.testing.v1 import AppTest
    monkeypatch.delenv("MFS_NO_UPDATE_CHECK", raising=False)
    monkeypatch.setenv("MFS_UPDATE_URL", url)
    os.environ.update(MFS_STORE=str(tmp_path / "s.json"), MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path / "ex"))
    U._STATE.update(result=None, error=None, running=False, t=0.0, url=None)
    return AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300)


def _wait_done():
    for _ in range(100):
        if not U.status()["running"]:
            return
        time.sleep(0.1)


def test_asks_once_then_banner_later_skip(tmp_path, monkeypatch, feed):
    at = _app(tmp_path, monkeypatch, feed); at.run()
    assert any("Check for Reticula updates" in m.value for m in at.markdown)
    at.button(key="upd_yes").click().run()
    assert json.load(open(tmp_path / "update.json", encoding="utf-8"))["opt_in"] is True
    _wait_done(); at.run()
    assert not at.exception, [e.value for e in at.exception][:1]
    assert any("Reticula 9.1.0 is available" in m.value for m in at.markdown)
    at.button(key="upd_later").click().run()
    assert not any("is available" in m.value for m in at.markdown)
    b = _app(tmp_path, monkeypatch, feed); U._STATE.update(t=0.0); b.run(); _wait_done(); b.run()
    assert any("is available" in m.value for m in b.markdown)          # 'Later' = this session only
    b.button(key="upd_skip").click().run()
    assert json.load(open(tmp_path / "update.json", encoding="utf-8"))["skip"] == "9.1.0"
    assert not any("is available" in m.value for m in b.markdown)
    assert not any("Check for Reticula updates" in m.value for m in b.markdown)   # never asked again


def test_offline_is_silent(tmp_path, monkeypatch):
    (tmp_path).mkdir(exist_ok=True); (tmp_path / "update.json").write_text(json.dumps({"opt_in": True}))
    at = _app(tmp_path, monkeypatch, "http://127.0.0.1:9/nothing"); at.run(); _wait_done(); at.run()
    assert not at.exception and not any("is available" in m.value for m in at.markdown)
    assert U.status()["error"]


def test_declined_never_fetches(tmp_path, monkeypatch, feed):
    (tmp_path / "update.json").write_text(json.dumps({"opt_in": False}))
    at = _app(tmp_path, monkeypatch, feed); at.run()
    assert not U.status()["running"] and U.status()["url"] is None
    assert not any("Check for Reticula updates" in m.value for m in at.markdown)


def test_pick_asset_macos_arch_v212():
    from mfs import update as U
    rel = {"assets": [{"name": n} for n in ("Reticula_Setup_2.12.0.exe", "Reticula_2.12.0_macOS_x86_64.dmg", "Reticula_2.12.0_macOS_arm64.dmg")]}
    assert U.pick_asset(rel, "Darwin", "arm64")["name"].endswith("macOS_arm64.dmg")
    assert U.pick_asset(rel, "Darwin", "x86_64")["name"].endswith("macOS_x86_64.dmg")
    only_arm = {"assets": [{"name": "Reticula_2.12.0_macOS_arm64.dmg"}]}
    assert U.pick_asset(only_arm, "Darwin", "x86_64") is None          # Intel Mac -> release page, not an arm64 image
    only_intel = {"assets": [{"name": "Reticula_2.12.0_macOS_x86_64.dmg"}]}
    assert U.pick_asset(only_intel, "Darwin", "arm64")["name"].endswith("x86_64.dmg")   # Rosetta
    untagged = {"assets": [{"name": "Reticula_2.12.0.dmg"}]}
    assert U.pick_asset(untagged, "Darwin", "arm64")["name"] == "Reticula_2.12.0.dmg"
    assert U.pick_asset(rel, "Windows")["name"] == "Reticula_Setup_2.12.0.exe"
