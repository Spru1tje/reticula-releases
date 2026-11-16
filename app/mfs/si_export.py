"""Supporting-information tables (one consistent defect report per sample) as XLSX and DOCX."""
from __future__ import annotations
import io
import pandas as pd

COLS = [("sample", "Sample"), ("framework", "Framework"), ("metal", "Metal"), ("L", "Linkers / node"), ("x", "Missing linkers / node"),
        ("m", "Caps / node"), ("OS", "Open sites / node"), ("M", "M per node (g mol⁻¹)"), ("formula", "Minimal formula"),
        ("methods", "Methods (consensus)"), ("L_range", "L range (methods)"),
        ("W200", "TGA W(200 °C) (%)"), ("W_plateau", "TGA W(plateau) (%)"), ("W_end", "TGA residue (%)"),
        ("icp_Zr", "ICP Zr (wt%)"), ("icp_Hf", "ICP Hf (wt%)"), ("nmr_ratio", "¹H NMR cap : linker"),
        ("S_BET", "S_BET (m² g⁻¹)"), ("V_total", "V_total (cm³ g⁻¹)"), ("V_micro", "V_micro t-plot (cm³ g⁻¹)")]


def defect_rows(store: dict):
    """One row per sample that has a 'defects' section; other sections fill the input columns."""
    rows = []
    for sid, rec in store.items():
        sec = rec.get("sections", {})
        d = (sec.get("defects") or {}).get("values")
        if not d:
            continue
        t = (sec.get("tga") or {}).get("values", {}); i = (sec.get("icp") or {}).get("values", {})
        n = (sec.get("nmr") or {}).get("values", {}); s = (sec.get("n2") or {}).get("values", {})
        ratio = next((v for k, v in n.items() if k.startswith("ratio ") and v is not None), None)
        rows.append({"sample": rec.get("name", sid), "framework": d.get("framework"), "metal": d.get("metal"),
                     "L": d.get("L"), "x": d.get("x"), "m": d.get("m"), "OS": d.get("OS"), "M": d.get("M"),
                     "formula": d.get("formula"), "methods": d.get("methods"), "L_range": d.get("L_range"),
                     "W200": t.get("W200"), "W_plateau": t.get("W_plateau"), "W_end": t.get("W_end"),
                     "icp_Zr": i.get("Zr_wt"), "icp_Hf": i.get("Hf_wt"), "nmr_ratio": ratio,
                     "S_BET": s.get("S_BET"), "V_total": s.get("V_total"), "V_micro": s.get("V_micro")})
    df = pd.DataFrame(rows, columns=[c for c, _ in COLS])
    return df.rename(columns=dict(COLS))


def _round(df):
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].round(3 if "cm³" in c else 2)
    return out


def to_xlsx(df: pd.DataFrame, methods: pd.DataFrame | None = None, notes: str = "") -> bytes:
    b = io.BytesIO()
    with pd.ExcelWriter(b, engine="openpyxl") as xw:
        _round(df).to_excel(xw, sheet_name="Table S-defects", index=False)
        if methods is not None and len(methods):
            methods.to_excel(xw, sheet_name="Methods per sample", index=False)
        if notes:
            pd.DataFrame({"notes": notes.split("\n")}).to_excel(xw, sheet_name="Notes", index=False)
        ws = xw.sheets["Table S-defects"]
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = max(10, min(45, max(len(str(c.value or "")) for c in col) + 2))
    return b.getvalue()


def _x(t):
    from xml.sax.saxutils import escape
    return escape(str(t))


def _run(text, bold=False, sz=14):
    return (f'<w:r><w:rPr>{"<w:b/>" if bold else ""}<w:sz w:val="{sz}"/></w:rPr>'
            f'<w:t xml:space="preserve">{_x(text)}</w:t></w:r>')


def to_docx(df: pd.DataFrame, title="Table S1. Defect analysis per sample.", notes: str = "") -> bytes:
    """Word table (landscape, 7 pt, grid). Written directly as WordprocessingML (zip + XML) so it needs no python-docx/lxml
    (the Windows build ships neither); opens in Word, LibreOffice and python-docx."""
    import zipfile
    d = _round(df).dropna(axis=1, how="all")
    fmt = lambda v: "" if v is None or (isinstance(v, float) and pd.isna(v)) else (f"{v:g}" if isinstance(v, float) else str(v))
    cell = lambda t, b=False: f'<w:tc><w:tcPr><w:tcW w:w="0" w:type="auto"/></w:tcPr><w:p>{_run(t, b)}</w:p></w:tc>'
    rows = ['<w:tr>' + "".join(cell(c, True) for c in d.columns) + '</w:tr>']
    rows += ['<w:tr>' + "".join(cell(fmt(v)) for v in r) + '</w:tr>' for r in d.itertuples(index=False)]
    bd = "".join(f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="000000"/>' for e in ("top", "left", "bottom", "right", "insideH", "insideV"))
    body = (f'<w:p>{_run(title, True, 20)}</w:p>'
            f'<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblBorders>{bd}</w:tblBorders></w:tblPr>'
            f'<w:tblGrid>{"".join("<w:gridCol/>" for _ in d.columns)}</w:tblGrid>{"".join(rows)}</w:tbl>'
            + "".join(f'<w:p>{_run(line, False, 16)}</w:p>' for line in (notes.split("\n") if notes else []))
            + '<w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>'
              '<w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134" w:header="708" w:footer="708" w:gutter="0"/></w:sectPr>')
    W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    doc = f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="{W}"><w:body>{body}</w:body></w:document>'
    ct = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/>'
          '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>')
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", ct); z.writestr("_rels/.rels", rels); z.writestr("word/document.xml", doc)
    return b.getvalue()
