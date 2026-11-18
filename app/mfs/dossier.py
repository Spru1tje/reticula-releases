"""v2.15.1: a dossier per material - everything about ONE material in one ZIP: its figures (PDF / PNG / SVG), the plotted data, the settings and
provenance, the results of that material from the sample store, the methods text of the techniques used and one Excel with the figure data.
Figures are assigned to a material in the SI export page (column *material*, guessed from the sample name of a graph that shows one sample)."""
from __future__ import annotations

import io
import json
import re
import zipfile

from mfs import samplestore as SS
from mfs import si_bundle as SB
from mfs import summary as SUM


def guess_material(series) -> str:
    """The material of a figure: the sample name when the graph shows exactly one sample, else '' (the user fills it in)."""
    s = [str(x) for x in (series or []) if str(x).strip()]
    return s[0] if len(set(s)) == 1 else ""


def group_by_material(figs: dict, assignment: dict) -> dict:
    """{material: {stem: rec}} for the non-empty materials; `assignment` = {stem: material}."""
    out: dict = {}
    for stem, rec in figs.items():
        m = str(assignment.get(stem, "") or "").strip()
        if m:
            out.setdefault(m, {})[stem] = rec
    return out


def store_for(material: str, store: dict) -> dict:
    """The part of the sample store that belongs to the material (same normalised id, or one contains the other)."""
    sid = SS.find({"sample_store": store}, material)
    return {sid: store[sid]} if sid else {}


def build_dossier(material: str, rendered: dict, store: dict, version: str = "", include=("pdf", "png", "svg", "csv"), provenance: dict | None = None) -> bytes:
    """`rendered` = {stem: {page, time, params, cfg, title, n_rows, pdf?, png?, svg?, csv?, data?(DataFrame)}} for this material."""
    from mfs import journal as J
    base = SB.build({k: {kk: v for kk, v in r.items() if kk != "data"} for k, r in rendered.items()}, store_for(material, store), version, include)
    buf = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(base)) as zin, zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for it in zin.infolist():
            z.writestr(it, zin.read(it.filename))
        z.writestr("figure_data.xlsx", J.figure_data_xlsx({k: {"data": r.get("data"), "page": r.get("page", ""), "title": r.get("title")} for k, r in rendered.items()}, version))
        pages = [r.get("page", "") for r in rendered.values()]
        z.writestr("methods.md", f"# {material} - methods\n\n" + SUM.methods_text(pages, version) + "\n")
        sm = SUM.sample_summary(store_for(material, store))
        if len(sm):
            z.writestr("sample_summary.csv", sm.to_csv(index=False))
            z.writestr("sample_summary.tex", SUM.summary_latex(sm, f"Characterisation of {material}."))
            z.writestr("sample_summary_sources.csv", SUM.summary_long(store_for(material, store)).to_csv(index=False))
            z.writestr("sample_summary_sourced.tex", SUM.summary_latex_sourced(store_for(material, store), f"Characterisation of {material}."))
        if provenance:
            z.writestr("provenance.json", json.dumps(provenance, indent=1, ensure_ascii=False, default=str))
        z.writestr("README.txt", f"Reticula {version} - dossier of '{material}'\n{len(rendered)} figure(s) from: {', '.join(sorted({p for p in pages if p}))}.\n"
                                 "figures/ and data/ as in the SI bundle; figure_data.xlsx = all plotted data (a sheet per figure); methods.md = starting text "
                                 "for the Methods (fill the [brackets]); sample_summary.* = results of this material from the sample store.\n")
    return buf.getvalue()


def build_all(materials: dict, store: dict, version: str = "", include=("pdf", "png", "svg", "csv")) -> bytes:
    """One ZIP with a folder per material. `materials` = {material: {stem: rendered}}."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for m, rend in materials.items():
            folder = re.sub(r"[^A-Za-z0-9._@+-]+", "_", m)
            inner = build_dossier(m, rend, store, version, include)
            with zipfile.ZipFile(io.BytesIO(inner)) as zi:
                for it in zi.infolist():
                    z.writestr(f"{folder}/{it.filename}", zi.read(it.filename))
    return buf.getvalue()
