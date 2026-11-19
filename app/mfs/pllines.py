"""v2.18.4: approximate emission ranges of common lanthanide (and Cr3+) transitions, for labelling PL spectra, and the thermally coupled pairs that are used
for ratiometric thermometry with their usual energy gap. The ranges are typical values from the luminescence literature; the exact positions depend on the host
(crystal field) and are only a guide - the program never assigns a band on this basis alone."""
from __future__ import annotations

# ion -> [(transition, lo nm, hi nm)]
LINES = {
    "Er³⁺": [("²H₉/₂→⁴I₁₅/₂", 403, 415), ("²H₁₁/₂→⁴I₁₅/₂", 515, 535), ("⁴S₃/₂→⁴I₁₅/₂", 535, 560), ("⁴F₉/₂→⁴I₁₅/₂", 640, 680), ("⁴I₁₁/₂→⁴I₁₅/₂", 960, 1010),
             ("⁴I₁₃/₂→⁴I₁₅/₂", 1480, 1620)],
    "Ho³⁺": [("⁵F₄,⁵S₂→⁵I₈", 535, 555), ("⁵F₅→⁵I₈", 640, 670), ("⁵F₄,⁵S₂→⁵I₇", 740, 760), ("⁵I₆→⁵I₈", 1150, 1200), ("⁵I₇→⁵I₈", 1900, 2100)],
    "Yb³⁺": [("²F₅/₂→²F₇/₂", 960, 1060)],
    "Tm³⁺": [("¹D₂→³F₄", 445, 455), ("¹G₄→³H₆", 470, 490), ("¹G₄→³F₄", 640, 660), ("³H₄→³H₆", 780, 820), ("³H₄→³F₄", 1450, 1500), ("³F₄→³H₆", 1700, 1900)],
    "Nd³⁺": [("⁴F₅/₂,²H₉/₂→⁴I₉/₂", 795, 815), ("⁴F₃/₂→⁴I₉/₂", 870, 910), ("⁴F₃/₂→⁴I₁₁/₂", 1050, 1070), ("⁴F₃/₂→⁴I₁₃/₂", 1320, 1350)],
    "Eu³⁺": [("⁵D₀→⁷F₀", 577, 581), ("⁵D₀→⁷F₁", 585, 600), ("⁵D₀→⁷F₂", 610, 625), ("⁵D₀→⁷F₃", 645, 660), ("⁵D₀→⁷F₄", 680, 710)],
    "Tb³⁺": [("⁵D₄→⁷F₆", 485, 495), ("⁵D₄→⁷F₅", 540, 550), ("⁵D₄→⁷F₄", 580, 590), ("⁵D₄→⁷F₃", 615, 625)],
    "Dy³⁺": [("⁴I₁₅/₂→⁶H₁₅/₂", 450, 460), ("⁴F₉/₂→⁶H₁₅/₂", 475, 485), ("⁴F₉/₂→⁶H₁₃/₂", 570, 580)],
    "Sm³⁺": [("⁴G₅/₂→⁶H₅/₂", 560, 566), ("⁴G₅/₂→⁶H₇/₂", 595, 605), ("⁴G₅/₂→⁶H₉/₂", 640, 650)],
    "Cr³⁺": [("²E→⁴A₂ (R lines)", 690, 700)],
}

# thermally coupled pairs: (ion, upper-level band, lower-level band, typical gap range in cm-1)
COUPLED = [
    ("Er³⁺", "²H₁₁/₂→⁴I₁₅/₂", "⁴S₃/₂→⁴I₁₅/₂", (650, 900)),
    ("Nd³⁺", "⁴F₅/₂,²H₉/₂→⁴I₉/₂", "⁴F₃/₂→⁴I₉/₂", (900, 1100)),
    ("Dy³⁺", "⁴I₁₅/₂→⁶H₁₅/₂", "⁴F₉/₂→⁶H₁₅/₂", (850, 1150)),
]


def ranges(ion):
    return LINES.get(ion, [])


def find(nm):
    """All (ion, transition) whose typical range contains the wavelength."""
    return [(ion, tr) for ion, ls in LINES.items() for tr, lo, hi in ls if lo <= nm <= hi]


def coupled_pair(c1_nm, c2_nm):
    """A known thermally coupled pair that matches the two band centres (in either order): (ion, upper, lower, gap range) or None."""
    for ion, up, low, gap in COUPLED:
        u = next(((lo, hi) for tr, lo, hi in LINES[ion] if tr == up), None)
        l_ = next(((lo, hi) for tr, lo, hi in LINES[ion] if tr == low), None)
        if not u or not l_:
            continue
        for a, b in ((c1_nm, c2_nm), (c2_nm, c1_nm)):
            if u[0] - 3 <= a <= u[1] + 3 and l_[0] - 3 <= b <= l_[1] + 3:
                return ion, up, low, gap
    return None


def draw_labels(ax, ions, lo, hi, ymax, colour="#555555"):
    """Light bands and short labels of the selected ions' transitions inside lo..hi (nm) on an existing axis."""
    import itertools
    cols = itertools.cycle(["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"])
    for ion in ions:
        c = next(cols)
        for tr, a, b in LINES.get(ion, []):
            if b < lo or a > hi:
                continue
            ax.axvspan(max(a, lo), min(b, hi), color=c, alpha=0.08, lw=0)
            ax.text((max(a, lo) + min(b, hi)) / 2, ymax, f"{ion} {tr}", rotation=90, fontsize=4.5, ha="center", va="top", color=c)
