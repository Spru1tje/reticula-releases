"""Reproducible check of the PL window fit (mfs/plwindows.py) for reviewer and author.

Prints, in one place: the git head, the sha256 and the CR count of both files, and the literal pytest summary.
Run from the app folder:  python tools/check_plwindows.py [logfile]
Exit code: 0 when pytest passes, otherwise the pytest exit code. Nothing is written into the repository except the optional log file.
"""
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
FILES = ["mfs/plwindows.py", "tests/test_v2220_plwindows.py"]


def _git_head():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=APP, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def main(argv):
    lines = [f"git HEAD: {_git_head()}"]
    for rel in FILES:
        path = os.path.join(APP, rel)
        data = open(path, "rb").read()
        cr = data.count(b"\r")
        nl = data.count(b"\n")
        lines.append(f"{rel}: sha256 {hashlib.sha256(data).hexdigest()}  CR {cr}  lines {nl}")
    proc = subprocess.run([sys.executable, "-m", "pytest", "tests/test_v2220_plwindows.py", "-q", "-p", "no:cacheprovider"],
                          cwd=APP, capture_output=True, text=True)
    summary = (proc.stdout.strip().splitlines() or ["(no pytest output)"])[-1]
    lines.append(f"pytest: {summary}")
    lines.append(f"pytest exit code: {proc.returncode}")
    report = "\n".join(lines)
    print(report)
    if len(argv) > 1:
        with open(argv[1], "w", encoding="utf-8") as fh:
            fh.write(report + "\n")
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv))
