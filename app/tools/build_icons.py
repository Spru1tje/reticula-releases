r"""Rebuild the Reticon icon library that Reticula uses (mfs/icons_runtime and mfs/components/scheme_canvas).

The icon toolkit (reticon-v4.2, build_library.py) is kept outside this repository. Point RETICON_TOOLKIT at its folder, then run:
    set RETICON_TOOLKIT=C:\path\to\reticon          (Windows)   /  export RETICON_TOOLKIT=/path/to/reticon   (macOS / Linux)
    python tools/build_icons.py
The result is copied into mfs/icons_runtime and mfs/components/scheme_canvas, and the vendored code in mfs/reticon is refreshed."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent          # this file lives in app/tools/


def main() -> int:
    tk = os.environ.get("RETICON_TOOLKIT")
    if not tk or not (Path(tk) / "build_library.py").exists():
        print("Set RETICON_TOOLKIT to the reticon toolkit folder (the one that contains build_library.py).")
        return 2
    tk = Path(tk)
    r = subprocess.run([sys.executable, "build_library.py"], cwd=tk)
    if r.returncode:
        return r.returncode
    dist = tk / "dist" / "runtime"
    dst = APP / "mfs" / "icons_runtime"
    comp = APP / "mfs" / "components" / "scheme_canvas"
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(dist, dst, ignore=shutil.ignore_patterns("component_scheme_canvas", "__pycache__"))
    shutil.rmtree(comp, ignore_errors=True)
    shutil.copytree(dist / "component_scheme_canvas", comp, ignore=shutil.ignore_patterns("__pycache__"))
    for f in ("__init__.py", "runtime.py", "scene.py", "raster.py", "core.py", "generators.py", "anchors.py", "variants.py", "external.py",
              "taxonomy.py", "library.py", "overrides.py"):
        if (tk / "reticon" / f).exists():
            shutil.copy(tk / "reticon" / f, APP / "mfs" / "reticon" / f)
    print("icon library copied into mfs/icons_runtime and mfs/components/scheme_canvas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
