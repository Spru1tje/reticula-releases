#!/bin/bash
# Copy the Homebrew tesseract binary + dylibs + eng tessdata into ./tesseract
# so PyInstaller can ship it inside the .app (same layout as the Windows build).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/tesseract"
BIN="$(command -v tesseract || true)"
if [[ -z "$BIN" ]]; then
  echo "tesseract not on PATH" >&2
  exit 1
fi
rm -rf "$OUT"
mkdir -p "$OUT/lib" "$OUT/tessdata"
cp "$BIN" "$OUT/tesseract"
chmod +x "$OUT/tesseract"

# Prefer dylibbundler when available
if command -v dylibbundler >/dev/null 2>&1; then
  dylibbundler -od -b -x "$OUT/tesseract" -d "$OUT/lib" -p "@executable_path/lib/"
else
  python3 - "$OUT" << 'PY'
import os, shutil, subprocess, sys
from pathlib import Path
out = Path(sys.argv[1])
bin_path = out / "tesseract"
libdir = out / "lib"

def deps(p: Path):
    text = subprocess.check_output(["otool", "-L", str(p)], text=True)
    libs = []
    for line in text.splitlines()[1:]:
        path = line.strip().split(" (")[0].strip()
        if not path or path.startswith(("/usr/lib/", "/System/", "@")):
            continue
        if Path(path).resolve() == Path(p).resolve():
            continue
        libs.append(path)
    return libs

seen = {}
queue = [bin_path]
while queue:
    cur = queue.pop()
    for lib in deps(cur):
        name = os.path.basename(lib)
        if name in seen:
            continue
        if not os.path.exists(lib):
            print("missing dylib", lib)
            continue
        dest = libdir / name
        shutil.copy2(lib, dest)
        os.chmod(dest, 0o755)
        seen[name] = lib
        queue.append(dest)

# Rewrite install names
def rewrite(target: Path, is_bin: bool):
    for lib in deps(target) if False else []:
        pass
    # re-read original deps via otool before we change anything — use seen map by basename
    text = subprocess.check_output(["otool", "-L", str(target)], text=True)
    for line in text.splitlines()[1:]:
        path = line.strip().split(" (")[0].strip()
        name = os.path.basename(path)
        if name not in seen and not is_bin:
            continue
        if name not in seen:
            continue
        new = f"@executable_path/lib/{name}" if is_bin else f"@loader_path/{name}"
        subprocess.check_call(["install_name_tool", "-change", path, new, str(target)])
    if not is_bin:
        subprocess.check_call(["install_name_tool", "-id", f"@loader_path/{target.name}", str(target)])

rewrite(bin_path, True)
for dest in libdir.iterdir():
    rewrite(dest, False)
print("bundled", len(seen), "dylibs")
PY
fi

# tessdata: Homebrew share, then the checked-in eng fallback
TESSDATA_SRC=""
for c in \
  "$(brew --prefix tesseract 2>/dev/null)/share/tessdata" \
  "$(brew --prefix)/share/tessdata" \
  "/opt/homebrew/share/tessdata" \
  "/usr/local/share/tessdata" \
  "$ROOT/pkg/tessdata"
do
  if [[ -n "$c" && -d "$c" ]]; then
    TESSDATA_SRC="$c"
    break
  fi
done
if [[ -z "$TESSDATA_SRC" ]]; then
  echo "no tessdata directory found" >&2
  exit 1
fi
# Copy eng (and osd if present); keep the payload small
cp "$TESSDATA_SRC/eng.traineddata" "$OUT/tessdata/" 2>/dev/null || cp "$ROOT/pkg/tessdata/eng.traineddata" "$OUT/tessdata/"
if [[ -f "$TESSDATA_SRC/osd.traineddata" ]]; then
  cp "$TESSDATA_SRC/osd.traineddata" "$OUT/tessdata/" || true
fi
# If brew tessdata is a configs/tessconfigs tree, ignore it.

# Smoke-test the bundled binary (TESSDATA_PREFIX must be set)
export TESSDATA_PREFIX="$OUT/tessdata"
export DYLD_LIBRARY_PATH="$OUT/lib:${DYLD_LIBRARY_PATH:-}"
"$OUT/tesseract" --version
ls -la "$OUT" "$OUT/lib" "$OUT/tessdata"
