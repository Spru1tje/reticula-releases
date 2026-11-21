"""v2.18.11: start-up check of the system tools that some readers need (Tesseract OCR with the English data file).
Returns the missing items with the install hint for the platform, so the app can say what to install instead of failing later."""
from __future__ import annotations

import os
import shutil
import sys


def missing_tools(which=shutil.which, tessdata_prefix: str | None = None, platform: str | None = None) -> list[dict]:
    """A list of {tool, why, hint} for every missing item. An empty list means all checks pass."""
    plat = platform or sys.platform
    out = []
    if which("tesseract") is None:
        hint = {"win32": "install Tesseract for Windows (or use the Reticula installer, which bundles it)",
                "darwin": "brew install tesseract tesseract-lang",
                }.get(plat, "sudo apt install tesseract-ocr tesseract-ocr-eng")
        out.append({"tool": "tesseract", "why": "OCR of scale bars and gel images", "hint": hint})
    prefix = tessdata_prefix if tessdata_prefix is not None else os.environ.get("TESSDATA_PREFIX", "")
    if prefix and not os.path.isfile(os.path.join(prefix, "eng.traineddata")):
        out.append({"tool": "eng.traineddata", "why": "English OCR language data missing in TESSDATA_PREFIX", "hint": f"copy eng.traineddata into {prefix}"})
    return out
