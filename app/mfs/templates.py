"""Journal figure templates (ACS, RSC, Wiley) and a pre-submission check of a finished figure."""
from __future__ import annotations
import matplotlib.pyplot as plt

CM = 1 / 2.54
TEMPLATES = {
    "ACS": {"single": 3.25, "double": 7.0, "font_size": 7.0, "min_font": 4.5, "label_style": "(a)", "dpi": 600,
            "fonts": "Arial / Helvetica", "note": "ACS: 3.25 in single, 7 in double column; Arial/Helvetica; text ≥ 4.5 pt at print size; "
                                                   "line art ≥ 600 dpi (1200 dpi recommended), halftones ≥ 300 dpi."},
    "RSC": {"single": 8.3 * CM, "double": 17.1 * CM, "font_size": 7.0, "min_font": 6.0, "label_style": "(a)", "dpi": 600,
            "fonts": "Arial / Helvetica", "note": "RSC: 8.3 cm single, 17.1 cm double column; sans-serif text ≈ 7 pt (≥ 6 pt); "
                                                   "≥ 600 dpi line art, ≥ 300 dpi images."},
    "Wiley": {"single": 8.5 * CM, "double": 17.5 * CM, "font_size": 7.0, "min_font": 6.0, "label_style": "a)", "dpi": 600,
              "fonts": "Arial / Helvetica", "note": "Wiley-VCH: 8.5 cm single, 17.5 cm double column; sans-serif text ≈ 7 pt (≥ 6 pt); "
                                                     "panel labels a), b) …; ≥ 600 dpi line art, ≥ 300 dpi images."},
}


def width(template: str, columns: str = "single") -> float:
    t = TEMPLATES[template]
    return t["double"] if columns.startswith("double") else (t["single"] + t["double"]) / 2 if columns.startswith("1.5") else t["single"]


def check(fig, template: str, columns: str = "single", tol_in: float = 0.05):
    """Returns a list of (level, message): width against the column width, every visible text against the minimum size."""
    t = TEMPLATES[template]; out = []
    w_target = width(template, columns)
    w = fig.get_size_inches()[0]
    if abs(w - w_target) > tol_in:
        out.append(("warn", f"figure width {w:.2f} in ≠ {template} {columns} column {w_target:.2f} in"))
    small = []
    for txt in fig.findobj(plt.Text):
        if not txt.get_visible() or not txt.get_text().strip():
            continue
        if txt.get_fontsize() < t["min_font"] - 1e-6:
            small.append(f"'{txt.get_text()[:18]}' {txt.get_fontsize():.1f} pt")
    if small:
        out.append(("warn", f"{len(small)} text item(s) below {t['min_font']:g} pt: " + ", ".join(small[:6])))
    if not out:
        out.append(("ok", f"meets the {template} {columns}-column width and minimum text size ({t['min_font']:g} pt)"))
    return out
