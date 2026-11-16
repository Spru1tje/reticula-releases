#!/bin/bash
# Package dist/Reticula.app + README_macOS.txt into a UDZO dmg.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
APP="$ROOT/dist/Reticula.app"
README="$ROOT/pkg/README_macOS.txt"
VER="${RETICULA_VERSION:?set RETICULA_VERSION, e.g. 2.13.0}"
ARCH="$(uname -m)"
# Normalise: arm64 stays arm64, i386/x86_64 -> x86_64
if [[ "$ARCH" == "x86_64" || "$ARCH" == "i386" ]]; then
  ARCH="x86_64"
fi
OUT="$ROOT/Reticula_${VER}_macOS_${ARCH}.dmg"
STAGE="$ROOT/dmg_root"

if [[ ! -d "$APP" ]]; then
  echo "missing $APP" >&2
  exit 1
fi

rm -rf "$STAGE" "$OUT"
mkdir -p "$STAGE"
# Copy the bundle (ditto preserves resource forks / Finder info)
ditto "$APP" "$STAGE/Reticula.app"
cp "$README" "$STAGE/README_macOS.txt"
ln -s /Applications "$STAGE/Applications"

# Clear leftover extended attributes from the build tree
xattr -cr "$STAGE/Reticula.app" || true

hdiutil create -volname "Reticula ${VER}" -srcfolder "$STAGE" -ov -format UDZO -imagekey zlib-level=9 "$OUT"
# detach any accidental mount
ls -lh "$OUT"
shasum -a 256 "$OUT" | tee "$OUT.sha256"
rm -rf "$STAGE"
echo "DMG=$OUT"
