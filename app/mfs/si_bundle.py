"""SI export bundle: every figure collected in this session (vector PDF + PNG + data CSV), a parameters table (figure
settings and the analysis widgets on the page at capture time), the sample-store results and the defect SI table."""
from __future__ import annotations
import datetime, io, json, re, zipfile
import pandas as pd


def _flat(d, pre=""):
    out = {}
    for k, v in (d or {}).items():
        if isinstance(v, dict):
            out.update(_flat(v, f"{pre}{k}.")); continue
        out[f"{pre}{k}"] = v if isinstance(v, (int, float, str, bool)) or v is None else json.dumps(v, ensure_ascii=False, default=str)[:300]
    return out


def tables(figs: dict, store: dict):
    rows, par = [], []
    for i, (stem, f) in enumerate(figs.items(), start=1):
        rows.append({"Figure": f"S{i}", "file stem": stem, "page": f.get("page", ""), "title": f.get("title") or "",
                     "captured": f.get("time", ""), "data points": f.get("n_rows")})
        for k, v in sorted((f.get("params") or {}).items()):
            par.append({"Figure": f"S{i}", "file stem": stem, "parameter": k, "value": str(v)[:300]})
        for k, v in sorted(_flat(f.get("cfg") or {}).items()):
            par.append({"Figure": f"S{i}", "file stem": stem, "parameter": f"figure.{k}", "value": str(v)[:300]})
    res = []
    for sid, rec in (store or {}).items():
        for sec, d in rec.get("sections", {}).items():
            for k, v in _flat(d.get("values") or {}).items():
                res.append({"sample": rec.get("name", sid), "technique": sec, "quantity": k, "value": v, "source": d.get("source", "")})
    return pd.DataFrame(rows), pd.DataFrame(par, columns=["Figure", "file stem", "parameter", "value"]), \
        pd.DataFrame(res, columns=["sample", "technique", "quantity", "value", "source"])


def build(figs: dict, store: dict, version: str = "", include=("pdf", "png", "svg", "csv")) -> bytes:
    from . import si_export as SI
    fig_t, par_t, res_t = tables(figs, store)
    buf = io.BytesIO(); z = zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED)
    for i, (stem, f) in enumerate(figs.items(), start=1):
        base = f"Figure_S{i}_{re.sub(r'[^A-Za-z0-9._@+-]+', '_', stem)}"
        for ext in ("pdf", "png", "svg"):
            if ext in include and f.get(ext):
                z.writestr(f"figures/{base}.{ext}", f[ext])
        if "csv" in include and f.get("csv"):
            z.writestr(f"data/{base}.csv", f["csv"])
    xb = io.BytesIO()
    with pd.ExcelWriter(xb, engine="openpyxl") as w:
        fig_t.to_excel(w, sheet_name="Figures", index=False); par_t.to_excel(w, sheet_name="Parameters", index=False)
        res_t.to_excel(w, sheet_name="Results per sample", index=False)
    z.writestr("SI_parameters.xlsx", xb.getvalue())
    z.writestr("SI_parameters.csv", par_t.to_csv(index=False))
    try:
        df = SI.defect_rows(store or {})
        if len(df):
            z.writestr("SI_Table_defects.docx", SI.to_docx(df)); z.writestr("SI_Table_defects.xlsx", SI.to_xlsx(df))
    except Exception as e:
        z.writestr("SI_Table_defects_ERROR.txt", str(e))
    z.writestr("README.txt", f"Reticula {version} SI bundle – {datetime.datetime.now():%Y-%m-%d %H:%M}\n"
               f"{len(figs)} figure(s): figures/ (PDF vector, PNG, SVG), data/ (plotted data as CSV).\n"
               "SI_parameters.xlsx: Figures (list), Parameters (analysis settings on the page when the figure was captured + "
               "figure customisation), Results per sample (sample store).\nSI_Table_defects.docx/.xlsx: defect table for all samples "
               "with defect results (if any).\n")
    z.close(); return buf.getvalue()
