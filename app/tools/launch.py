"""Cross-platform launcher: start Streamlit and open the browser once the server answers."""
import os, socket, subprocess, sys, time, webbrowser

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def free_port(start=8501):
    for p in range(start, start + 50):
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", p)) != 0:
                return p
    return start


def main():
    port = free_port()
    cmd = [sys.executable, "-m", "streamlit", "run", os.path.join(ROOT, "app.py"), "--server.port", str(port),
           "--server.headless", "true", "--browser.gatherUsageStats", "false", *sys.argv[1:]]
    proc = subprocess.Popen(cmd, cwd=ROOT)
    url = f"http://localhost:{port}"
    for _ in range(120):
        time.sleep(0.5)
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                break
    print(f"\n  Reticula is running at {url}  (close this window or press Ctrl+C to stop)\n", flush=True)
    if os.environ.get("MFS_NO_BROWSER") != "1":
        webbrowser.open(url)
    try:
        proc.wait()
    except KeyboardInterrupt:
        proc.terminate()


if __name__ == "__main__":
    main()
