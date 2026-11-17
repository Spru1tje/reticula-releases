#!/bin/bash
# Build the unsigned Reticula macOS app (Reticula.app) and the disk image (.dmg) on THIS Mac.
# Normally started by double-clicking build_mac.command in the top folder. Or:   bash pkg/build_mac.sh
# Needs: a Mac (Apple Silicon M1-M4 is the target), internet the first time, ~6 GB free disk, 20-40 minutes.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
VER="$(sed -n 's/^__version__ *= *"\([^"]*\)".*/\1/p' "$ROOT/app/mfs/__init__.py")"
ARCH="$(uname -m)"; [[ "$ARCH" == "i386" ]] && ARCH="x86_64"
LOG="$ROOT/build_mac.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo " Reticula ${VER} - macOS build"
echo " $(date)"
echo " folder: $ROOT"
echo " chip:   $ARCH      macOS: $(sw_vers -productVersion 2>/dev/null || echo unknown)"
echo " log:    $LOG"
echo "============================================================"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This has to run on a Mac (it uses iconutil, codesign and hdiutil)." >&2; exit 1
fi
if [[ -z "$VER" ]]; then echo "Could not read the version from app/mfs/__init__.py" >&2; exit 1; fi
if [[ "$ARCH" != "arm64" ]]; then
  echo "NOTE: this Mac is $ARCH (Intel). The build is meant for Apple Silicon; it should still work but is untested on Intel."
  echo "      The disk image will be named ...macOS_${ARCH}.dmg."
fi
FREE_GB=$(df -g "$ROOT" | awk 'NR==2{print $4}')
if [[ "${FREE_GB:-0}" -lt 6 ]]; then echo "WARNING: only ${FREE_GB} GB free on this disk; the build needs about 6 GB."; fi

have() { command -v "$1" >/dev/null 2>&1; }

# ---------------- Python 3.11 - 3.13 with Tcl/Tk (3.12.0 is refused: it breaks the frozen app) ----------------
PY=""
good_python() {
  "$1" - <<'PYCHECK' 2>/dev/null
import sys
v = sys.version_info[:3]
assert (3, 11) <= v[:2] <= (3, 13) and v != (3, 12, 0), v
import tkinter  # the control window needs Tk
PYCHECK
}
pick_python() {
  local c
  for c in /Library/Frameworks/Python.framework/Versions/3.13/bin/python3.13 \
           /Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12 \
           /Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 \
           /opt/homebrew/bin/python3.13 /opt/homebrew/bin/python3.12 /opt/homebrew/bin/python3.11 \
           /usr/local/bin/python3.13 /usr/local/bin/python3.12 /usr/local/bin/python3.11 \
           "$(command -v python3.13 || true)" "$(command -v python3.12 || true)" "$(command -v python3.11 || true)" "$(command -v python3 || true)"; do
    [[ -n "$c" && -x "$c" ]] || continue
    if good_python "$c"; then PY="$c"; return 0; fi
  done
  return 1
}
install_python_org() {
  echo
  echo "No suitable Python (3.11-3.13 with Tcl/Tk) found."
  echo "I can download the official installer from python.org and run it (macOS asks for your password)."
  read -r -p "Install Python 3.12 from python.org now? [y/N] " ans
  [[ "$ans" == "y" || "$ans" == "Y" ]] || return 1
  local v url pkg
  for v in 3.12.12 3.12.11 3.12.10 3.12.9 3.12.8 3.12.7; do
    url="https://www.python.org/ftp/python/${v}/python-${v}-macos11.pkg"
    if curl -fsI "$url" >/dev/null; then
      pkg="/tmp/reticula-python-${v}.pkg"
      echo "Downloading Python $v ..."; curl -L --retry 5 -o "$pkg" "$url"
      sudo installer -pkg "$pkg" -target /
      return 0
    fi
  done
  return 1
}
if ! pick_python; then
  if install_python_org && pick_python; then :; else
    echo
    echo "Install Python 3.12 from https://www.python.org/downloads/macos/ (the official installer includes Tcl/Tk)"
    echo "or:  brew install python@3.12 python-tk@3.12      and run this again."
    exit 1
  fi
fi
echo "Python: $("$PY" -c 'import sys; print(sys.executable, sys.version.split()[0])')"

# ---------------- Tesseract (reads the scale bars of SEM images) - optional ----------------
if have brew; then
  have tesseract || { echo "Installing tesseract (Homebrew) ..."; brew install tesseract || true; }
  have dylibbundler || brew install dylibbundler || true
else
  echo "Homebrew (https://brew.sh) is not installed: continuing without Tesseract (SEM scale-bar OCR will not work in this build)."
fi
if have tesseract; then tesseract --version | head -1; fi

# ---------------- virtual environment + packages ----------------
VENV="$ROOT/.venv-build"
[[ -x "$VENV/bin/python" ]] || "$PY" -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install --upgrade pip wheel setuptools
echo "Installing the Python packages (several minutes the first time) ..."
python -m pip install -r "$ROOT/requirements-build.txt"
python -c "import streamlit, numpy, scipy, matplotlib, pymatgen, rdkit, tkinter, h5py; print('package imports ok')"

# ---------------- icon, Tesseract bundle, PyInstaller ----------------
bash "$ROOT/pkg/make_icns.sh"
if have tesseract; then
  bash "$ROOT/pkg/bundle_tesseract.sh" || echo "WARNING: could not bundle Tesseract; continuing without it"
fi
rm -rf "$ROOT/build" "$ROOT/dist"
echo "Running PyInstaller (10-20 minutes is normal) ..."
python -m PyInstaller --noconfirm --clean "$ROOT/pkg/mfs_mac.spec"
[[ -d "$ROOT/dist/Reticula.app" ]] || { echo "PyInstaller did not make dist/Reticula.app" >&2; exit 1; }
codesign --force --deep --sign - "$ROOT/dist/Reticula.app" || true     # ad-hoc signature (not notarised)

# ---------------- self-test of the built app ----------------
BIN="$ROOT/dist/Reticula.app/Contents/MacOS/Reticula"
OUT="$ROOT/selftest.json"; rm -f "$OUT"
echo "Self-test of the built app (up to 15 minutes) ..."
python - "$BIN" "$OUT" <<'PYRUN' || true
import subprocess, sys
try:
    subprocess.run([sys.argv[1], "--selftest", sys.argv[2]], timeout=900)
except subprocess.TimeoutExpired:
    print("self-test timed out")
PYRUN
SELFTEST_OK="not run"
if [[ -s "$OUT" ]]; then
  if python "$ROOT/pkg/check_selftest.py" "$OUT"; then SELFTEST_OK="all checks ok"; else SELFTEST_OK="SOME CHECKS FAILED (see selftest.json / build_mac.log)"; fi
fi

# ---------------- disk image ----------------
RETICULA_VERSION="$VER" bash "$ROOT/pkg/make_dmg.sh"
DMG="$ROOT/Reticula_${VER}_macOS_${ARCH}.dmg"
echo
echo "============================================================"
if [[ -f "$DMG" ]]; then
  ls -lh "$DMG"; cat "$DMG.sha256"
  echo
  echo "DONE.  Disk image:   $DMG"
  echo "Self-test:           $SELFTEST_OK"
  echo
  echo "Install:  open the .dmg, drag Reticula to Applications. First start: right-click -> Open (the app is unsigned)."
  echo "Share:    send the .dmg (and the .sha256 file) - or attach it to the GitHub release."
else
  echo "The app was built but the .dmg is missing. You can still use:  $ROOT/dist/Reticula.app"
  exit 1
fi
echo "============================================================"
