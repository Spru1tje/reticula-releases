"""Build docs/icons/INDEX.md from the proof SVGs in assets/icons_proof.

Domain, source and licence come from METADATA below. A new SVG without metadata is an error,
so the index can never list an icon without a declared source and licence.
LEGACY_MAP lists which 24-grid navigation icon (assets/nav_icons, unchanged) each 64-grid icon replaces.
The index text is generated from LEGACY_MAP, so the mapping cannot drift from the index.
Usage (from app/):  python tools/build_icon_index.py [--check]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
PROOF = APP / "assets" / "icons_proof"
NAV_LEGACY = APP / "assets" / "nav_icons"
INDEX = APP.parent / "docs" / "icons" / "INDEX.md"
ORDER = [
    "pxrd_pattern", "tem_instrument", "tem_lattice", "sem_particle", "mof_pore_frame", "sorption_isotherm",
    "nmr_spectrum", "ftir_spectrum", "zeta_potential", "structure_cell",
    "tga_curve", "dsc_curve", "icp_bars", "uvvis_spectrum", "xps_peaks",
]
OWN = {"source": "own drawing", "license": "own work (Reticula)"}
METADATA = {
    "pxrd_pattern": {"domain": "PXRD", **OWN},
    "tem_instrument": {"domain": "TEM", **OWN},
    "tem_lattice": {"domain": "TEM", **OWN},
    "sem_particle": {"domain": "SEM", **OWN},
    "mof_pore_frame": {"domain": "MOF / porous materials", **OWN},
    "sorption_isotherm": {"domain": "Sorption", **OWN},
    "nmr_spectrum": {"domain": "NMR", **OWN},
    "ftir_spectrum": {"domain": "FTIR", **OWN},
    "zeta_potential": {"domain": "Zetasizer", **OWN},
    "structure_cell": {"domain": "Structure (unit cell)", **OWN},
    "tga_curve": {"domain": "TGA", **OWN},
    "dsc_curve": {"domain": "DSC", **OWN},
    "icp_bars": {"domain": "ICP", **OWN},
    "uvvis_spectrum": {"domain": "UV/Vis", **OWN},
    "xps_peaks": {"domain": "XPS", **OWN},
}
# legacy 24-grid nav icon -> 64-grid icon in this set
LEGACY_MAP = {
    "nmr": "nmr_spectrum",
    "ftir": "ftir_spectrum",
    "zeta": "zeta_potential",
    "structure_cif": "structure_cell",
    "sorption": "sorption_isotherm",
    "sem": "sem_particle",
    "pxrd": "pxrd_pattern",
    "tga": "tga_curve",
    "dsc": "dsc_curve",
    "icp": "icp_bars",
    "uvvis": "uvvis_spectrum",
    "xps": "xps_peaks",
}
HEADER = """# Reticon proof set: source index

Branch `reticon/biorender-style`. All icons below are drawn for Reticula from scratch.
No external artwork, no BioRender or other commercial sets, no CC0 sources used.
"""
FOOTER_TAIL = """
Exports: `app/assets/icons_proof/png/<name>_<size>.png` at 24, 48, 96 and 192 px,
made with `app/tools/export_icons.py` (headless Edge).
Contact sheet: `app/assets/icons_proof/contact_sheet.png`.
"""


def mapping_text() -> str:
    lines = [
        "",
        "Migration of the 24-grid navigation icons (`app/assets/nav_icons`, still used by the app and left unchanged):",
    ]
    for legacy, new in LEGACY_MAP.items():
        lines.append(f"- {legacy} -> {new}")
    return "\n".join(lines) + "\n"


def grid_of(svg: Path) -> str:
    m = re.search(r'viewBox="0 0 (\d+) \1"', svg.read_text(encoding="utf-8"))
    if not m:
        raise SystemExit(f"{svg.name}: viewBox is not a square starting at 0 0")
    return m.group(1)


def build() -> str:
    names = sorted(p.stem for p in PROOF.glob("*.svg"))
    missing = [n for n in names if n not in METADATA]
    if missing:
        raise SystemExit(f"no metadata for: {', '.join(missing)}")
    ordered = [n for n in ORDER if n in names] + [n for n in names if n not in ORDER]
    rows = ["| Icon | Domain | Source | License | Grid |", "|---|---|---|---|---|"]
    for n in ordered:
        m = METADATA[n]
        rows.append(f"| {n} | {m['domain']} | {m['source']} | {m['license']} | {grid_of(PROOF / (n + '.svg'))} |")
    return HEADER + "\n" + "\n".join(rows) + "\n" + FOOTER_TAIL + mapping_text()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="fail if INDEX.md is out of date")
    args = ap.parse_args()
    text = build()
    if args.check:
        if INDEX.read_bytes().decode("utf-8") != text:
            print("INDEX.md is out of date; run tools/build_icon_index.py")
            return 1
        print("INDEX.md is up to date")
        return 0
    INDEX.write_bytes(text.encode("utf-8"))  # LF bytes, not the platform newline
    print(f"wrote {INDEX}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
