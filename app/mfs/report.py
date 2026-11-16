"""Defect-calculator report: JSON round-trip, CSV and a multi-page PDF with rendered equations."""
from __future__ import annotations
import re, textwrap, datetime
from dataclasses import asdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .defects import DefectInput, Species, results_table, SUMMARY_PITFALLS
from . import plots, style


def input_to_dict(inp: DefectInput) -> dict:
    d = asdict(inp)
    return d


def input_from_dict(d: dict) -> DefectInput:
    d = dict(d)
    d["species"] = [Species(**s) for s in d.get("species", [])]
    return DefectInput(**d)


def tex_to_mathtext(s: str) -> str:
    s = re.sub(r"\\text\{([^}]*)\}", lambda m: "{" + m.group(1) + "}" if m.group(1) in ("-",) else r"\mathrm{" + m.group(1) + "}", s)
    s = re.sub(r"\\xrightarrow\{[^}]*\}", r"\\rightarrow", s)
    s = s.replace(r"\tfrac", r"\frac").replace(r"\left(", "(").replace(r"\right)", ")").replace(r"\left[", "[").replace(r"\right]", "]")
    s = s.replace(r"\arcsin", r"\mathrm{arcsin}")
    return s


def _render_ok(expr):
    from matplotlib.mathtext import MathTextParser
    try:
        MathTextParser("path").parse(f"${expr}$")
        return True
    except Exception:
        return False


class _Page:
    def __init__(self, pdf, title=None):
        self.pdf = pdf
        self.fig = plt.figure(figsize=(8.27, 11.69))
        self.y = 0.95
        if title:
            self.text(title, size=12, bold=True, gap=0.03)

    def need(self, h):
        if self.y - h < 0.05:
            self.flush()
            self.fig = plt.figure(figsize=(8.27, 11.69)); self.y = 0.95

    def text(self, s, size=8.5, bold=False, gap=None, indent=0.06, width=118, color="black"):
        lines = []
        for para in str(s).split("\n"):
            lines += textwrap.wrap(para, width) or [""]
        h = (size / 72 / 11.69) * 1.35
        for l in lines:
            self.need(h)
            self.fig.text(indent, self.y, l, fontsize=size, fontweight="bold" if bold else "normal", va="top", color=color)
            self.y -= h
        self.y -= gap if gap is not None else h * 0.3

    def math(self, expr, size=10):
        e = tex_to_mathtext(expr)
        h = (size / 72 / 11.69) * 2.6
        self.need(h)
        if _render_ok(e):
            self.fig.text(0.09, self.y, f"${e}$", fontsize=size, va="top")
        else:
            self.fig.text(0.09, self.y, expr[:160], fontsize=7, va="top", family="monospace")
        self.y -= h

    def figure(self, fig, h=0.3):
        self.need(h)
        import io
        from PIL import Image
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=220); buf.seek(0)
        im = Image.open(buf)
        ax = self.fig.add_axes([0.05, self.y - h, 0.9, h]); ax.imshow(im); ax.axis("off")
        self.y -= h + 0.01

    def flush(self):
        self.pdf.savefig(self.fig); plt.close(self.fig)


def pdf_report(res, path):
    inp = res["input"]; cons = res["consensus"]
    with PdfPages(path) as pdf:
        pg = _Page(pdf, f"Defect analysis report - {inp.sample}")
        pg.text(f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} by Reticula. Framework {inp.framework} "
                f"(CN = {inp.cn()}, topicity t = {inp.t()}, ideal L0 = {inp.L0():g}); metal {inp.metal} (f_Hf = {inp.fhf():.3f})"
                + ("; extra-framework Zr treated as ZrO2" if inp.extraframework_Zr else "") + ".")
        if cons:
            pg.text("Consensus", size=10, bold=True)
            pg.text(f"L = {cons['L']:.2f} linkers/node (range {cons['L_min']:.2f}–{cons['L_max']:.2f}, SD {cons['L_sd']:.2f}, n = {cons['n_methods']}); "
                    f"missing linkers x = {cons['x']:.2f}; caps m = {cons['m']:.2f}; Cl = {cons['Cl']:.2f}; open sites OS = {cons['OS']:.2f} per node.")
            pg.text(f"Minimal formula (coordinated H2O on terminal OH/Cl): {cons['formula']}   M = {cons['M']:.1f} g/mol; "
                    f"metal {cons['metal_wt_pct']:.1f} wt%.")
            pg.text(f"Anhydrous: {cons['formula_anhydrous']}   M = {cons['M_anhydrous']:.1f} g/mol.")
        fig = plots.fig_defect_comparison(res, 7.0)
        pg.figure(fig, 0.27); plt.close(fig)
        tab = results_table(res)
        pg.text("Method comparison (per node)", size=10, bold=True)
        for _, r in tab.iterrows():
            pg.text(f"{r['method']}  [{r['status']}{', semi-quant.' if r['semi-quantitative'] else ''}]  L={r['linkers/node L']}  x={r['missing linkers/node x']}  "
                    f"m={r['caps/node m']}  Cl={r['Cl/node']}  OS={r['open sites/node']}  M={r['M (g/mol per node)']}  {r['note'] or ''}", size=7.5)
        pg.text("When methods disagree", size=10, bold=True)
        for a, b in SUMMARY_PITFALLS:
            pg.text(f"{a}: {b}", size=7.5)
        pg.text("Inputs", size=10, bold=True)
        for k, v in input_to_dict(inp).items():
            if k == "species":
                for s in v:
                    pg.text(f"  species: {s['name']}  {s['formula']}  z={s['z']}  role={s['role']}  nH={s['nH']}  integral={s['integral']}  at plateau={s['at_plateau']}  main={s['main']}", size=7)
            elif v not in (None, [], ""):
                pg.text(f"  {k} = {v}", size=7)
        # open-site rule
        h = res["h"]
        pg.text(h["title"], size=10, bold=True)
        for e in h["equations"]:
            pg.math(e)
        for k, v in h["symbols"].items():
            pg.text(f"{k}: {v}", size=7.5)
        for a in h["assumptions"]:
            pg.text("• " + a, size=7.5)
        pg.flush()
        for m in res["methods"]:
            pg = _Page(pdf, m.title)
            pg.text(f"Status: {m.status}{' (semi-quantitative)' if m.semi_quant else ''}. {m.message}", size=8.5)
            pg.text("References", bold=True, size=9.5)
            for c, doi in m.refs:
                pg.text(f"{c} https://doi.org/{doi}", size=7.5)
            pg.text("Equations", bold=True, size=9.5)
            for e in m.equations:
                pg.math(e)
            pg.text("Symbols", bold=True, size=9.5)
            for k, v in m.symbols.items():
                pg.text(f"{k}: {v}", size=7.5)
            pg.text("Assumptions", bold=True, size=9.5)
            for a in m.assumptions:
                pg.text("• " + a, size=7.5)
            if m.steps:
                pg.text("Worked substitution with your numbers", bold=True, size=9.5)
                for s in m.steps:
                    pg.math(s, size=9)
            if m.L is not None:
                pg.text(f"Result: L = {m.L:.3f}, x = {m.x:.3f}" + (f", m = {m.m:.3f}" if m.m is not None else "") +
                        (f", Cl = {m.Cl:.3f}" if m.Cl is not None else "") + (f", OS = {m.OS:.3f}" if m.OS is not None else ""), bold=True)
            if m.formula:
                pg.text(f"Formula: {m.formula}  (M = {m.M_formula:.1f} g/mol)")
            pg.text("Pitfalls", bold=True, size=9.5)
            for p in m.pitfalls:
                pg.text("• " + p, size=7.5)
            pg.flush()
    return path


def csv_report(res) -> str:
    tab = results_table(res)
    lines = [tab.to_csv(index=False)]
    lines.append("\n# step-by-step\nmethod,step\n")
    for m in res["methods"]:
        for s in m.steps:
            lines.append(f"\"{m.key}\",\"{s}\"\n")
    return "".join(lines)
