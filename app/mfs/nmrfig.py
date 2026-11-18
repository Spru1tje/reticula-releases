"""Labelled publication ¹H NMR figure + numbered structure depiction (v2.13)."""
from __future__ import annotations

import io

import numpy as np

OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]


def mol_drawing(mol, groups=None, size=(360, 280), fmt="png", notes=True, highlight=None):
    """Structure with atom numbers on the H-bearing atoms (same labels as the assignment). Returns (bytes|str, coords)
    where coords = {atom index: (x, y) px}."""
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit import Chem
    m = Chem.Mol(mol)
    if not m.GetNumConformers():
        rdDepictor.Compute2DCoords(m)
    if notes:
        lab = {a: ("H" + g["tag"] if g.get("tag") else l) for g in (groups or []) for a, l in zip(g["atoms"], g["labels"])}
        for a in m.GetAtoms():
            if a.GetIdx() in lab:
                a.SetProp("atomNote", str(lab[a.GetIdx()]))
    W, H = size
    d = rdMolDraw2D.MolDraw2DCairo(W, H) if fmt == "png" else rdMolDraw2D.MolDraw2DSVG(W, H)
    o = d.drawOptions()
    o.annotationFontScale = 0.75
    o.addStereoAnnotation = False
    o.padding = 0.08
    o.clearBackground = fmt == "png"
    hl = list(highlight or [])
    if hl:
        d.DrawMolecule(m, highlightAtoms=hl, highlightAtomColors={i: (0.98, 0.6, 0.25) for i in hl}, highlightBonds=[])
    else:
        d.DrawMolecule(m)
    d.FinishDrawing()
    coords = {}
    for a in m.GetAtoms():
        p = d.GetDrawCoords(a.GetIdx())
        coords[a.GetIdx()] = (float(p.x), float(p.y))
    return d.GetDrawingText(), coords


def auto_xlim(sp, mults, pad=0.35):
    use = [m for m in mults if not m.note or m.label]
    if not use:
        return float(sp.ppm.max()), float(sp.ppm.min())
    hi = max(m.hi for m in use) + pad; lo = min(m.lo for m in use) - pad
    if lo > 0 and lo < 0.8:
        lo = -0.2
    return min(hi, float(sp.ppm.max())), max(lo, float(sp.ppm.min()))


def _spread(xs, minsep):
    """Label x positions (descending ppm order) pushed apart by at least minsep, staying close to the targets."""
    xs = list(xs)
    if not xs:
        return []
    order = np.argsort(xs)[::-1]
    pos = [xs[i] for i in order]
    for _ in range(200):
        moved = False
        for k in range(1, len(pos)):
            if pos[k - 1] - pos[k] < minsep:
                mid = (pos[k - 1] + pos[k]) / 2
                pos[k - 1] = mid + minsep / 2; pos[k] = mid - minsep / 2; moved = True
        if not moved:
            break
    out = [0.0] * len(xs)
    for k, i in enumerate(order):
        out[i] = pos[k]
    return out


def _free_box(x, y, labelled, hi, lo, w=0.3):
    """Inset box (axes fraction) over the widest part of the window without signal or labels."""
    span = hi - lo
    occ = np.zeros(200, bool)
    f = (hi - x) / span
    occ[np.clip((f[y > 0.04] * 199).astype(int), 0, 199)] = True
    for m in labelled:
        a = int(np.clip((hi - m.hi) / span * 199, 0, 199)); b = int(np.clip((hi - m.lo) / span * 199, 0, 199))
        occ[max(0, a - 6):b + 7] = True
    best, bs, cur = (0, 0), None, None
    for i, o in enumerate(np.r_[occ, True]):
        if not o and cur is None:
            cur = i
        if o and cur is not None:
            if i - cur > best[1] - best[0]:
                best = (cur, i)
            cur = None
    a, b = best[0] / 199, best[1] / 199
    ww = min(w, max(b - a - 0.02, 0.18))
    x0 = a + (b - a - ww) / 2
    return (max(0.0, x0), 0.22, ww, 0.62)


def _h_sub(label):
    """'Ha, Hb' becomes H with a subscript letter in the figure (mathtext)."""
    import re
    return re.sub(r"H([a-z]+)", lambda m: "H$_{\\mathrm{" + m.group(1) + "}}$", label)


def _plain_len(txt):
    """Length of a label without the mathtext markup (for the head-room estimate)."""
    import re
    return len(re.sub(r"\$_\{\\mathrm\{|\}\}\$", "", txt))


def _ax_size_pt(ax):
    """Size of the axes in points (width, height); a sensible guess when the figure has no renderer yet."""
    try:
        b = ax.get_window_extent()
        d = ax.figure.dpi
        w, h = b.width / d * 72, b.height / d * 72
        if w > 20 and h > 20:
            return w, h
    except Exception:
        pass
    fw, fh = ax.figure.get_size_inches()
    return fw * 72 * 0.8, fh * 72 * 0.6


def _stagger(centres, widths, gap):
    """Row index (0, 1, 2 ...) for each text so that texts on one row do not overlap. Positions in descending ppm order."""
    order = np.argsort(centres)[::-1]
    rows, last = [0] * len(centres), []
    for i in order:
        for r in range(len(last)):
            if centres[i] + widths[i] / 2 + gap <= last[r]:      # fits left of the previous text on row r
                rows[i] = r; last[r] = centres[i] - widths[i] / 2
                break
        else:
            rows[i] = len(last); last.append(centres[i] - widths[i] / 2)
    return rows


def draw(ax, sp, mults, mol=None, groups=None, impurities=None, xlim=None, show_integrals=True, show_labels=True,
         show_structure=True, show_mult=True, inset_box="auto", color="black", lw=0.6,
         label_fmt="{label}", fontsize=None, mark_unknown=True, show_int_numbers=True, simple=False):
    """Spectrum with multiplet/assignment labels (leader lines), integrals under the axis, solvent (*) and impurity (#)
    marks and an optional numbered structure inset. impurities: {mult index: name}.
    simple=True is the 'simple mode' for the final figure: every ink (spectrum, labels, leaders, integrals, marks) is black."""
    def ink(c):
        return "black" if simple else c
    if simple:
        color = "black"
    hi, lo = xlim or auto_xlim(sp, mults)
    m_ = (sp.ppm <= hi) & (sp.ppm >= lo)
    heights = []
    for m in mults:
        mm = (sp.ppm <= m.hi) & (sp.ppm >= m.lo) & m_
        if mm.any() and not m.note:
            heights.append(float(np.max(sp.y[mm])))
    heights = sorted(heights, reverse=True)
    norm = heights[0] if heights else (np.max(sp.y[m_]) if m_.any() else 1)
    clipped = False
    if len(heights) > 1 and heights[0] > 3.0 * heights[1]:          # one dominant line (CH3, solvent): clip it
        norm = 1.6 * heights[1]; clipped = True
    y = sp.y / max(norm, 1e-12)
    yc = np.minimum(y, 1.0)
    ax.plot(sp.ppm[m_], yc[m_], color=color, lw=lw)
    if clipped:
        for x in sp.ppm[m_][np.diff(np.r_[0, (y[m_] > 1.0).astype(int)]) == 1]:
            ax.text(x, 1.0, "≈", ha="center", va="center", fontsize=8, color=color, rotation=90,
                    bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
    y = yc
    ax.set_xlim(hi, lo)
    ymax = 1.0
    ax.set_xlabel("δ (ppm)" if sp.nucleus == "1H" else f"δ$_{{{sp.nucleus}}}$ (ppm)")
    ax.set_yticks([]); [ax.spines[s].set_visible(False) for s in ("left", "right", "top")]
    ax.tick_params(axis="x", which="both", top=False, labeltop=False)
    import matplotlib as mpl
    fs = fontsize or max(5.0, 0.8 * float(mpl.rcParams.get("font.size", 8)))
    impurities = impurities or {}
    vis = [(i, m) for i, m in enumerate(mults) if lo <= m.centre <= hi]
    lab = [(i, m) for i, m in vis if (m.label or (not m.note and i not in impurities and (m.nH or 0) >= 0.5))]
    top = 1.08
    w_pt, h_pt = _ax_size_pt(ax)
    span = hi - lo
    lab_fs = max(4.5, 0.8 * fs)                                   # v2.14: smaller, regular weight (was heavy and clipped)
    t_frac = 0.03                                                 # axes fraction reserved above the labels
    lab_items = []
    if show_labels and lab:
        xs = [m.centre for _, m in lab]
        tx = _spread(xs, minsep=max(span * 0.03, lab_fs * 1.35 / w_pt * span))
        for (i, m), x0 in zip(lab, tx):
            mm = (sp.ppm <= m.hi) & (sp.ppm >= m.lo)
            ypk = float(np.max(y[mm])) if mm.any() else 0.5
            txt = label_fmt.format(label=_h_sub(m.label) if m.label else f"{m.centre:.2f}", mult=m.mult, nH=m.nH, d=m.centre)
            lab_items.append((m, x0, ypk, txt))
        longest = max(_plain_len(t[3]) for t in lab_items)
        t_frac = min(0.45, longest * lab_fs * 0.62 / h_pt + 0.05)
    # integrals: bars + numbers in their own band under the baseline (horizontal, staggered when crowded)
    ints = []
    if show_integrals:
        for i, m in vis:
            if (m.note and not m.label) or m.nH is None or (m.nH < 0.5 and not m.label):
                continue
            ints.append((m, f"{m.nH:.2f}"))
    int_fs = max(4.5, 0.75 * fs)
    row_frac = int_fs * 1.35 / h_pt
    n_rows = 0
    rows = []
    if ints and show_int_numbers:
        wid = [len(t) * 0.58 * int_fs / w_pt * span for _, t in ints]
        rows = _stagger([m.centre for m, _ in ints], wid, gap=0.15 * int_fs / w_pt * span)
        n_rows = max(rows) + 1
    b_frac = (0.03 + 0.035 + n_rows * row_frac + 0.01) if (ints and show_int_numbers) else (0.075 if ints else 0.03)
    t_data = top + 0.04
    span_y = t_data / max(0.2, 1.0 - t_frac - b_frac)
    y_lo, y_hi = -b_frac * span_y, -b_frac * span_y + span_y
    ax.set_ylim(y_lo, y_hi)
    xt = ax.get_xaxis_transform()                                 # x in data, y in axes fraction (stable against the ylim)
    for m, x0, ypk, txt in lab_items:
        ax.plot([m.centre, m.centre, x0], [ypk + 0.03, top - 0.06, top - 0.02], color=ink("#6B7F99"), lw=0.45, solid_capstyle="round")
        ax.text(x0, top, txt, rotation=90, ha="center", va="bottom", fontsize=lab_fs, color=ink("#33475B"), fontweight="normal")
    ymax = y_hi
    for i, m in vis:
        mm = (sp.ppm <= m.hi) & (sp.ppm >= m.lo)
        ypk = float(np.max(y[mm])) if mm.any() else 0.5
        if m.note and not m.label:
            ax.text(m.centre, min(ypk, 1.0) + 0.03, "*", ha="center", va="bottom", fontsize=fs + 1, color=ink("#7A8CA5"))
        elif i in impurities:
            ax.text(m.centre, min(ypk, 1.0) + 0.03, "#", ha="center", va="bottom", fontsize=fs, color=ink("#D55E00"))
        elif not m.label and (m.nH or 0) >= 0.05 and mark_unknown and mol is not None:
            ax.text(m.centre, min(ypk, 1.0) + 0.03, "?", ha="center", va="bottom", fontsize=fs, color=ink("#7A8CA5"))
    for k_, (m, t) in enumerate(ints):
        a_, b_ = (max(m.lines) + 0.012, min(m.lines) - 0.012) if m.lines else (m.hi, m.lo)
        ax.plot([min(a_, m.hi), max(b_, m.lo)], [b_frac - 0.03] * 2, color=ink("#0072B2"), lw=0.9, transform=xt, clip_on=False)
        if show_int_numbers:
            ax.text(m.centre, b_frac - 0.055 - rows[k_] * row_frac, t, transform=xt, ha="center", va="top", fontsize=int_fs, color=ink("#0072B2"))
    if show_structure and mol is not None:
        png, _ = mol_drawing(mol, groups, size=(520, 400))
        import matplotlib.image as mpimg
        img = mpimg.imread(io.BytesIO(png), format="png")
        if inset_box is None or inset_box == "auto":
            inset_box = _free_box(sp.ppm[m_], y[m_], [m for _, m in lab] if show_labels else [], hi, lo)
        ia = ax.inset_axes(list(inset_box))
        ia.imshow(img); ia.axis("off")
    return (hi, lo)
