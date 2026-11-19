"""v2.18.3 page: Help and info - one info document (what it does, terms, tutorial, formulas, what is checked, limits) per module, readable here
and downloadable as a PDF in the Reticula colours (one module, or all as a ZIP)."""
from __future__ import annotations

import io
import re
import zipfile

import streamlit as st

from mfs import helpdocs as HD
from mfs import ui as UI


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


@st.cache_data(show_spinner=False, max_entries=64)
def _pdf(key: str) -> bytes:
    from mfs import helppdf
    return helppdf.render_pdf(HD.DOCS[key])


def info_button(page: str):
    """A download button for the info PDF of this module, at the top of the module (the PDF is made only when it is clicked)."""
    if page in ("Home", "Help and info") or page not in HD.DOCS:
        return
    label = UI.DISPLAY.get(page, page)
    c1, c2 = st.columns([5, 2])
    c2.download_button("📖 Info about this module (PDF)", lambda: _pdf(page), file_name=f"Reticula_info_{_safe(label)}.pdf",
                       mime="application/pdf", key=f"info_pdf_{_safe(page)}", on_click="ignore", width="stretch",
                       help="What the module does, the terms, a tutorial, the formulas, what is checked and the limits. More in 'Help and info'.")


def _show(doc: dict):
    step = 0
    for kind, text in doc["sections"]:
        if kind == "h":
            step = 0
            st.subheader(text)
        elif kind == "p":
            st.markdown(text)
        elif kind == "b":
            st.markdown("- " + text)
        elif kind == "n":
            step += 1
            st.markdown(f"**{step}.** {text}")
        elif kind == "f":
            st.latex(text)
        elif kind == "k":
            st.info(text)


def page_help():
    st.header("Help and info")
    st.caption("An info document for every module: what it does, the terms (written out before they are abbreviated), a short tutorial, the formulas, "
               "what the program checks, and the known limits. Read it here or download it as a PDF.")
    keys = list(HD.DOCS)
    labels = {k: (UI.DISPLAY.get(k, k) if k in UI.ICON else k) for k in keys}
    q = st.text_input("Search the documents", "", key="help_q", placeholder="e.g. Scherrer, Rouquerol, zeta, polynomial …").strip().lower()
    if q:
        hits = [k for k in keys if q in " ".join(t for _, t in HD.DOCS[k]["sections"]).lower() or q in labels[k].lower()]
        if not hits:
            st.info("No document mentions this."); return
        st.caption(f"{len(hits)} document(s) mention “{q}”: " + ", ".join(labels[k] for k in hits))
        keys = hits
    sel = st.selectbox("Module", keys, format_func=lambda k: f"{UI.ICON.get(k, '📖')}  {labels[k]}", key="help_sel")
    doc = HD.DOCS[sel]
    c1, c2 = st.columns([1, 1])
    c1.download_button("⬇️ Download this document (PDF)", _pdf(sel), file_name=f"Reticula_info_{_safe(labels[sel])}.pdf", mime="application/pdf",
                       key="help_dl", on_click="ignore", width="stretch")

    def _zip():
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for i, k in enumerate(HD.DOCS, 1):
                z.writestr(f"{i:02d}_{_safe(labels.get(k, k))}.pdf", _pdf(k))
        return buf.getvalue()
    c2.download_button("⬇️ Download all documents (ZIP of PDFs)", _zip, file_name="Reticula_info_all_modules.zip", mime="application/zip",
                       key="help_dl_all", on_click="ignore", width="stretch")
    st.markdown(f"### {doc['title']}")
    if doc.get("subtitle"):
        st.caption(doc["subtitle"])
    _show(doc)
