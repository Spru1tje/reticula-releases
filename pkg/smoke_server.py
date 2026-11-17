#!/usr/bin/env python3
"""Start Reticula the way a user starts it and check that the app is really served.

    python pkg/smoke_server.py [--timeout 900] [--client-python PY] [--cwd DIR] -- <command to start Reticula ...>

1. runs the command (browser and control window switched off: MFS_NO_BROWSER=1, MFS_NO_GUI=1);
2. waits until http://127.0.0.1:<8501-8600>/_stcore/health answers "ok";
3. fetches the index page and its JavaScript bundle (a bundle served with a wrong MIME type = white page);
4. opens a session over the websocket exactly like a browser does, lets the app script run once and fails if the
   script raised an exception or showed nothing.

Used by CI for: the run_* scripts (from a folder with a space in its name), `python mfs_launcher.py`, and the built
Windows / macOS apps. The orchestration needs only the standard library; step 4 runs in --client-python (default: this
interpreter), which needs streamlit (for its protobuf messages and the websockets client)."""
from __future__ import annotations

import argparse, json, os, re, signal, subprocess, sys, tempfile, time, urllib.request

HOST = "127.0.0.1"
PORTS = range(8501, 8601)


def http_get(url, timeout=5.0):
    req = urllib.request.Request(url, headers={"Cache-Control": "no-cache"})
    with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req, timeout=timeout) as r:
        return r.status, r.headers.get("Content-Type", ""), r.read()


def healthy_port():
    for p in PORTS:
        try:
            st, _, body = http_get(f"http://{HOST}:{p}/_stcore/health", 0.5)
            if st == 200 and body.strip() == b"ok":
                return p
        except Exception:
            pass
    return None


def frontend(port):
    st, ct, body = http_get(f"http://{HOST}:{port}/", 10)
    m = re.search(rb'src="(\.?/?static/js/[^"]+\.js)"', body)
    if st != 200 or not m:
        return False, f"index.html incomplete (HTTP {st}, {len(body)} bytes, no JS bundle)"
    src = m.group(1).decode().lstrip(".").lstrip("/")
    st2, ct2, js = http_get(f"http://{HOST}:{port}/{src}", 20)
    if st2 != 200 or len(js) < 1000:
        return False, f"JS bundle missing ({src}: HTTP {st2}, {len(js)} bytes)"
    if "javascript" not in ct2:
        return False, f"JS bundle served as '{ct2}' (expected application/javascript)"
    return True, f"frontend ok ({src.rsplit('/', 1)[-1]}: {ct2})"


# ----------------------------------------------------------------------------- step 4: one browser-like session
def session(port, timeout):
    """Connect like the Streamlit frontend, request one script run, collect what the script sent."""
    # Streamlit >= 1.5x serves through uvicorn/starlette and no longer installs tornado; the `websockets` package is a
    # Streamlit dependency, so the client uses that.
    import asyncio
    from websockets.asyncio.client import connect
    from websockets.exceptions import ConnectionClosed
    from streamlit.proto.BackMsg_pb2 import BackMsg
    from streamlit.proto.ForwardMsg_pb2 import ForwardMsg

    async def go():
        conn = await connect(f"ws://{HOST}:{port}/_stcore/stream", origin=f"http://{HOST}:{port}", subprotocols=["streamlit"],
                             max_size=1 << 28, open_timeout=30)
        back = BackMsg()
        back.rerun_script.SetInParent()
        await conn.send(back.SerializeToString())
        out = {"elements": 0, "kinds": {}, "exceptions": [], "texts": [], "finished": None, "runs": 0}
        t_end = time.time() + timeout
        while time.time() < t_end:
            try:
                raw = await asyncio.wait_for(conn.recv(), max(1.0, t_end - time.time()))
            except ConnectionClosed:
                raw = None
            if raw is None:
                out["finished"] = out["finished"] or "connection closed by the server"
                break
            if isinstance(raw, str):
                continue
            fm = ForwardMsg()
            fm.ParseFromString(raw)
            kind = fm.WhichOneof("type")
            if kind == "delta" and fm.delta.WhichOneof("type") == "new_element":
                el = fm.delta.new_element
                ek = el.WhichOneof("type")
                out["elements"] += 1
                out["kinds"][ek] = out["kinds"].get(ek, 0) + 1
                if ek == "exception":
                    out["exceptions"].append({"type": el.exception.type, "message": el.exception.message,
                                              "trace": list(el.exception.stack_trace)[-12:]})
                elif ek in ("markdown", "heading") and len(out["texts"]) < 40:
                    out["texts"].append(getattr(el, ek).body[:120])
            elif kind == "script_finished":
                name = ForwardMsg.ScriptFinishedStatus.Name(fm.script_finished)
                out["runs"] += 1
                if name != "FINISHED_EARLY_FOR_RERUN":
                    out["finished"] = name
                    break
        await conn.close()
        return out

    return asyncio.run(go())


# ----------------------------------------------------------------------------- orchestration
def kill_tree(proc):
    try:
        if os.name == "nt":
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True, timeout=30)
        else:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            time.sleep(1.5)
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except Exception:
                pass
    except Exception as e:
        print(f"(could not stop the process tree: {e})")


def launcher_logs():
    home = os.path.expanduser("~")
    return [os.path.join(os.environ.get("LOCALAPPDATA") or home, "Reticula", "server_log.txt"),
            os.path.join(home, "Library", "Application Support", "Reticula", "server_log.txt"),
            os.path.join(os.environ.get("XDG_STATE_HOME") or os.path.join(home, ".local", "state"), "Reticula", "server_log.txt")]


def show(path, title, tail=120):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return
    print(f"----- {title} ({path}; last {min(tail, len(lines))} of {len(lines)} lines)")
    for l in lines[-tail:]:
        print("   " + l)


def main():
    if len(sys.argv) >= 4 and sys.argv[1] == "--session":
        print(json.dumps(session(int(sys.argv[2]), float(sys.argv[3]))))
        return 0
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout", type=float, default=900, help="seconds to wait for the server (includes a first-run pip install)")
    ap.add_argument("--client-python", default=sys.executable)
    ap.add_argument("--cwd", default=None)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd[:1] == ["--"] else a.cmd
    if not cmd:
        ap.error("no command given")
    if healthy_port():
        print("A Streamlit server is already running on this machine; stop it first."); return 2
    for lg in launcher_logs():
        try:
            os.remove(lg)
        except OSError:
            pass
    env = dict(os.environ, MFS_NO_BROWSER="1", MFS_NO_GUI="1", MFS_NO_UPDATE_CHECK="1", PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")
    out_path = os.path.join(tempfile.mkdtemp(), "start_output.txt")
    print(f"$ {' '.join(cmd)}   (cwd: {a.cwd or os.getcwd()})", flush=True)
    kw = {"creationflags": 0x00000200} if os.name == "nt" else {"start_new_session": True}
    with open(out_path, "w", encoding="utf-8") as out:
        proc = subprocess.Popen(cmd, cwd=a.cwd, env=env, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT, **kw)
    ok, t0, port = False, time.time(), None
    try:
        while time.time() - t0 < a.timeout:
            port = healthy_port()
            if port:
                break
            if proc.poll() is not None:
                time.sleep(3)                          # a launcher may hand over to a child process
                port = healthy_port()
                if port:
                    break
                print(f"FAIL: the start command ended (exit code {proc.returncode}) and no server is running.")
                break
            time.sleep(1)
        else:
            print(f"FAIL: no server answered within {a.timeout:.0f} s.")
        if port:
            print(f"ok   server healthy on port {port} after {time.time() - t0:.0f} s", flush=True)
            f_ok, f_msg = frontend(port)
            print(("ok   " if f_ok else "FAIL: ") + f_msg, flush=True)
            cp = subprocess.run([a.client_python, os.path.abspath(__file__), "--session", str(port), "300"],
                                capture_output=True, text=True, timeout=420, env=env)
            try:
                s = json.loads(cp.stdout.strip().splitlines()[-1])
            except Exception:
                s = None
                print("FAIL: the browser-like session could not be run:\n" + cp.stdout[-2000:] + cp.stderr[-4000:])
            if s is not None:
                print(f"     script run: {s['finished']}; {s['elements']} elements {json.dumps(s['kinds'], sort_keys=True)}")
                for t in s["texts"][:6]:
                    print("       | " + t.replace("\n", " ")[:110])
                for e in s["exceptions"]:
                    print(f"FAIL: the app raised {e['type']}: {e['message']}")
                    for l in e["trace"]:
                        print("       " + l)
                s_ok = s["finished"] == "FINISHED_SUCCESSFULLY" and not s["exceptions"] and s["elements"] >= 10
                if not s_ok and not s["exceptions"]:
                    print(f"FAIL: unexpected script result ({s['finished']}, {s['elements']} elements)")
                ok = f_ok and s_ok
    finally:
        kill_tree(proc)
        time.sleep(1)
        show(out_path, "output of the start command")
        for lg in launcher_logs():
            show(lg, "launcher server log")
    print("SMOKE TEST " + ("PASSED" if ok else "FAILED"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
