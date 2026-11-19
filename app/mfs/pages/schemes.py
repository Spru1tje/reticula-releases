"""v2.17 page: Reticon – a BioRender-style scheme editor with the Reticon icon library (offline, CC0), publication export,
and the content that comes from Reticula itself: a methods schematic, your own data as an icon, a molecule, a crystal structure, and a
scheme as a panel of the combined figure."""
from __future__ import annotations

import streamlit as st

from mfs import levels as LV
from mfs import schemes as SC


def _sample_records() -> dict[str, dict]:
    """Samples of the project: name -> record (with its sections of measurements)."""
    out = {}
    for sid, rec in (st.session_state.get("sample_store") or {}).items():
        out[rec.get("name", sid)] = rec
    return out



def _sources(store: dict, pick: str):
    """The 'add to this scheme' tabs: methods schematic, own data, molecule, crystal structure, combined figure."""
    t_meth, t_data, t_mol, t_cif, t_fig = st.tabs(["🧪 Methods schematic", "📈 My data as icon", "⚗️ Molecule", "💎 Crystal structure", "🧩 Combined figure"])
    uid = st.session_state.get("sch_uid", 0)

    with t_meth:
        recs = _sample_records()
        if not recs:
            st.caption("No samples yet. Analyse a sample on a technique page first; its measurements then show up here.")
        else:
            from mfs.summary import TECH_NAMES
            name = st.selectbox("Sample", list(recs), key="sch_meth_sample")
            secs = list((recs[name].get("sections") or {}).keys())
            chosen = st.multiselect("Techniques to show", secs, default=secs, key="sch_meth_secs",
                                    format_func=lambda k: TECH_NAMES.get(k, k))
            st.button("Make the methods schematic", key="sch_meth_go", disabled=not chosen, on_click=_methods_scheme, args=(name, chosen))
            st.caption("The sample and every technique it went through, as a flow with icons. Open it, then edit.")

    with t_data:
        up = st.file_uploader("Your data (two columns: x, y; CSV or text)", type=["csv", "txt", "dat", "xy"], key=f"sch_xy_{uid}")
        kind = st.selectbox("Drawing", ["line", "sticks", "dots"], key="sch_xy_kind")
        label = st.text_input("Label under the icon", "My data", key="sch_xy_label")
        if up is not None and st.button("Add the plot to the scheme", key="sch_xy_go"):
            try:
                from mfs.reticon import generators as G
                x, y = SC.parse_xy(up.getvalue())
                ic = G.icon_from_xy(x, y, kind=kind, name="my_data", label=label)
                SC.add_svgart(store[pick], SC.icon_svg(ic, prefix=f"xy{uid}"), f"xy{uid}")
                st.session_state["sch_uid"] = uid + 1
                st.rerun()
            except ValueError as e:
                st.error(str(e))

    with t_mol:
        smi = st.text_input("SMILES of the molecule", "OC(=O)c1ccc(cc1)C(=O)O", key="sch_smi",
                            help="For example benzene-1,4-dicarboxylic acid (BDC). Drawn from the SMILES by RDKit; nothing is sent anywhere.")
        if st.button("Add the molecule to the scheme", key="sch_smi_go", disabled=not smi.strip()):
            try:
                from mfs.reticon import generators as G
                ic = G.icon_from_smiles(smi.strip(), name="molecule")
                SC.add_svgart(store[pick], SC.icon_svg(ic, prefix=f"mol{uid}"), f"mol{uid}")
                st.session_state["sch_uid"] = uid + 1
                st.rerun()
            except Exception as e:      # an invalid SMILES must not break the page
                st.error(f"Could not draw this molecule: {e}")

    with t_cif:
        cif = st.file_uploader("Crystal structure (CIF)", type=["cif"], key=f"sch_cif_{uid}")
        mode = st.radio("Picture", ["net", "atoms"], horizontal=True, key="sch_cif_mode",
                        help="net: the framework as a topology net; atoms: the atoms and bonds of the cell.")
        if cif is not None and st.button("Add the structure to the scheme", key="sch_cif_go"):
            try:
                from mfs.reticon import generators as G
                ic = G.icon_from_cif(cif.getvalue().decode("utf-8", "replace"), name="cif_structure", mode=mode)
                SC.add_svgart(store[pick], SC.icon_svg(ic, prefix=f"cif{uid}"), f"cif{uid}")
                st.session_state["sch_uid"] = uid + 1
                st.rerun()
            except Exception as e:
                st.error(f"Could not read this CIF: {e}")

    with t_fig:
        st.caption("Puts this scheme as one panel into the combined figure (sidebar → 🧩 Combined figure), so it sits next to the graphs "
                   "of the other techniques. The scheme is included as an image at 600 dpi.")
        if st.button("🧩 Add this scheme to the combined figure", key="sch_fig_go"):
            png, _info = SC.export(store[pick], "png", width_mm=178, dpi=600)
            key = "scheme|" + pick
            ss = st.session_state
            ss.setdefault("fig_registry", {})[key] = {"fn": SC.figure_panel(png), "cfg": {}, "fams": [], "page": "Reticon", "stem": pick}
            comp = ss.setdefault("composer", {})
            comp[key] = ss["fig_registry"][key]["fn"]
            st.toast(f"Added to the combined figure ({len(comp)} panel{'s' if len(comp) != 1 else ''}).", icon="🧩")


def _new_scheme(tpl: str):
    """Button callback: create a scheme and select it before the selectbox is drawn (setting a widget key after it is drawn raises)."""
    base = SC.empty_scene() if tpl == "empty" else SC.templates()[tpl]
    st.session_state["sch_pick"] = SC.add_scene(tpl.replace("_", " ") if tpl != "empty" else "scheme", base)


def _methods_scheme(name: str, chosen: list):
    from mfs.reticon import scene as RS
    st.session_state["sch_pick"] = SC.add_scene(f"methods {name}", RS.methods_schematic(name, chosen))


def _import_scheme(up):
    try:
        st.session_state["sch_pick"] = SC.add_scene(up.name.rsplit(".", 1)[0], SC.load_json_scene(up.getvalue()))
        st.session_state["sch_import_err"] = ""
    except (ValueError, UnicodeDecodeError) as e:
        st.session_state["sch_import_err"] = f"Could not read this file: {e}"


@st.cache_data(show_spinner=False)
def _tpl_preview(name: str) -> str:
    import base64
    data, _info = SC.export(SC.templates()[name], "svg")
    return base64.b64encode(data).decode("ascii")


def esc_alt(text: str) -> str:
    return text.replace('"', "&quot;").replace("<", "&lt;")


def page_schemes():
    st.header("Reticon – schemes with lab and materials icons")
    st.caption("Drag icons from the library onto the canvas, add arrows and text, then export to SVG, PDF, PNG or TIFF at the width and DPI of your journal. "
               "The icons are CC0 and bundled: nothing is downloaded.")
    # v2.18.2: the editor and its icon library load when the page opens (no extra click)
    store = SC.scenes()
    left, right = st.columns([3, 2])
    with left:
        opts = list(store) or ["(no scheme yet)"]
        pick = st.selectbox("Scheme", opts, key="sch_pick")
    with right:
        tpl_names = list(SC.templates())
        tpl = st.selectbox("Start from", ["empty"] + tpl_names, key="sch_tpl",
                           format_func=lambda t: {"empty": "empty canvas"}.get(t, t.replace("_", " ")))
        st.button("➕ New scheme", key="sch_new", width="stretch", on_click=_new_scheme, args=(tpl,))
        if tpl != "empty":                     # preview of the chosen starting scheme (our own bundled template, shown as SVG)
            try:
                st.markdown(f'<img alt="{esc_alt(tpl)} preview" src="data:image/svg+xml;base64,{_tpl_preview(tpl)}" '
                            f'style="max-width:100%;border:1px solid #CBD5E1;border-radius:6px;background:#fff">', unsafe_allow_html=True)
            except Exception:                    # a broken preview must not break the page
                pass

    if pick not in store:
        st.caption("Create a scheme to start. Schemes are saved with the project file.")
        return

    with st.expander("Import a scheme (JSON)", expanded=False):
        up = st.file_uploader("Scheme JSON", type=["json"], key="sch_up")
        if up is not None and st.button("Import as a new scheme", key="sch_imp", on_click=_import_scheme, args=(up,)):
            pass

    edited = SC.canvas(store[pick], key=f"canvas_{pick}")
    if edited is not None and edited != store[pick]:
        store[pick] = edited

    with st.expander("➕ Add content to this scheme (methods, your data, molecule, structure, combined figure)", expanded=False):
        _sources(store, pick)

    with st.expander("⬇️ Export this scheme (publication)", expanded=True):
        from mfs.reticon.scene import JOURNAL_WIDTHS_MM
        c1, c2, c3 = st.columns(3)
        fmt = c1.selectbox("Format", ["PDF", "SVG", "PNG", "TIFF"], key="sch_fmt")
        width_name = c2.selectbox("Width", list(JOURNAL_WIDTHS_MM), index=0, key="sch_w",
                                  format_func=lambda k: f"{k} ({JOURNAL_WIDTHS_MM[k]} mm)")
        dpi = c3.selectbox("DPI (PNG / TIFF)", [300, 600, 1200], index=1, key="sch_dpi", disabled=fmt in ("PDF", "SVG"))
        try:
            data, info = SC.export(store[pick], fmt, width_mm=JOURNAL_WIDTHS_MM[width_name], dpi=dpi)
        except Exception as e:      # a failed export must not break the page
            st.error(f"Export failed: {e}")
            return
        for w in info.get("warnings", []):
            st.warning(w)
        if info.get("vector_only"):
            st.caption("All icons are vector: the PDF / SVG have no bitmaps.")
        base = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in pick) or "scheme"
        ext = fmt.lower()
        mime = {"pdf": "application/pdf", "svg": "image/svg+xml", "png": "image/png", "tiff": "image/tiff"}[ext]
        st.download_button(f"⬇️ Download {fmt}", data=data, file_name=f"{base}.{ext}", mime=mime, key="sch_dl")
        if fmt in ("PNG", "TIFF"):
            st.caption(f"{info.get('width_px', '')} × {info.get('height_px', '')} px at {dpi} dpi.")
