#!/usr/bin/env python3
"""Print a selftest.json report and exit 1 if any check starts with FAIL."""
from __future__ import annotations
import json, sys
from pathlib import Path

p = Path(sys.argv[1] if len(sys.argv) > 1 else "selftest.json")
if not p.is_file() or p.stat().st_size == 0:
    print(f"selftest produced no output file: {p}", file=sys.stderr)
    sys.exit(1)
r = json.loads(p.read_text(encoding="utf-8"))
Path("selftest_pretty.txt").write_text(json.dumps(r, indent=2), encoding="utf-8")
# The SEM scale-bar OCR check is excluded from the release gate by the user (it needs
# a local OCR install; the CI Linux runner has none). It is still printed, but not counted.
EXCLUDED_FROM_GATE = {"sem_ocr"}
fails = [k for k, v in r.items() if str(v).startswith("FAIL") and k not in EXCLUDED_FROM_GATE]
oks = [k for k, v in r.items() if not str(v).startswith("FAIL")]
print(f"SELFTEST {len(oks)}/{len(r)} ok")
for k, v in r.items():
    flag = "FAIL" if str(v).startswith("FAIL") else "ok  "
    text = v if not str(v).startswith("FAIL") else v[:400]
    print(f"  {flag}  {k}: {text}")
if fails:
    print("FAILED: " + ", ".join(fails), file=sys.stderr)
    sys.exit(1)
