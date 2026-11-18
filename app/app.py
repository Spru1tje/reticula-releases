"""Reticula (formerly MOF Figure Studio) - publication figures for framework materials: one app for the common characterisation techniques of framework materials + a defect calculator.
Run:  streamlit run app.py"""
from __future__ import annotations
import time as _time
_T0 = _time.perf_counter()     # v2.11: per-run timing (MFS_TIMING=<file> appends one JSON line per script run)
import contextlib, io, os, re, json, glob, zipfile, warnings, datetime, hashlib
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from mfs import readers as R, analysis as A, plots as P, style as S, defects as D, composer as C, report as REP
from mfs import figctl as F, figctl_ui as FU, ui as UI, __version__
from mfs.chem import SPECIES, FRAMEWORKS, molar_mass

HERE = os.path.dirname(os.path.abspath(__file__))


# v2.11: parsed instrument files are cached on their content hash (st.cache_data), so reruns / page revisits / the same
# file on two pages do not parse again. R.read_* and DSC.read_dsc go through this; results are copies (safe to modify).
class _NamedBytes(io.BytesIO):
    def __init__(self, b, name):
        super().__init__(b); self.name = name


@st.cache_data(max_entries=400, show_spinner=False)
def _cached_read_impl(mod, fn, data, name, args, kw):
    import importlib
    return getattr(importlib.import_module(mod), fn)(_NamedBytes(data, name), *args, **dict(kw))


def _plain(v):
    return v is None or isinstance(v, (str, int, float, bool)) or (isinstance(v, tuple) and all(_plain(x) for x in v))


def _cached_reader(mod, fn, f, *args, **kw):
    func = getattr(__import__(mod, fromlist=[fn]), fn)
    if not (all(_plain(a) for a in args) and all(_plain(v) for v in kw.values())):
        return func(f, *args, **kw)
    if isinstance(f, (bytes, bytearray)):
        data, name = bytes(f), ""
    elif isinstance(f, str):
        with open(f, "rb") as fh:
            data = fh.read()
        name = os.path.basename(f)
    elif hasattr(f, "getvalue"):
        data, name = f.getvalue(), getattr(f, "name", "") or ""
    else:
        return func(f, *args, **kw)
    return _cached_read_impl(mod, fn, data, name, tuple(args), tuple(sorted(kw.items())))


class _CachedReaders:
    """Module proxy: read_* functions cached (except those in `skip`), everything else passed through."""
    def __init__(self, mod, skip=()):
        self._m, self._skip = mod, set(skip)

    def __getattr__(self, n):
        a = getattr(self._m, n)
        if n.startswith("read_") and callable(a) and n not in self._skip:
            import functools
            return functools.partial(_cached_reader, self._m.__name__, n)
        return a


R = _CachedReaders(R, skip=("read_sem",))
EX = os.path.join(HERE, "example_data")
REF = os.path.join(HERE, "references")
_ICON = os.path.join(HERE, "assets", "reticula_icon.png")
st.set_page_config(page_title="Reticula", layout="wide", page_icon=_ICON if os.path.exists(_ICON) else "📈")
from mfs import pagestate as PS   # v2.9: keep every page's work across navigation + unsaved-work dialog
PS.configure(os.path.join(os.path.dirname(os.environ.get("MFS_STORE") or os.path.join(os.path.dirname(
    os.environ.get("MFS_EXPORT_DIR") or os.path.join(HERE, "exports")), "sample_store.json")), "workspace"))
PS.install_uploader_hook()
from mfs import perf as _PF0
_PF0.install_fast_widget_state()      # v2.11: O(1) widget state lookups
PS.install_autokey_hook()
PS.begin()
from mfs import project as PRJ     # v2.11 project files (.reticula)
PRJ.open_from_argv()
PRJ.apply_pending(os.environ.get("MFS_STORE") or os.path.join(os.path.dirname(os.environ.get("MFS_EXPORT_DIR") or os.path.join(HERE, "exports")), "sample_store.json"))
try:  # presets write widget values through session_state; don't show Streamlit's one-off warning about it
    from streamlit.elements.lib import policies as _pol
    _pol._shown_default_value_warning = True
except Exception:
    pass
S.set_series_colours(st.session_state.get("series_colours", {}))

# ============================================================== sidebar: global style
UI.inject(HERE)
UI.sidebar_brand(HERE, __version__)
PAGES = UI.PAGES
if "composer" not in st.session_state:
    st.session_state.composer = {}
from mfs import levels as LV            # v2.15: Standard / Advanced / Developer
LV.sidebar()                            # v2.15.1: coloured box at the TOP of the sidebar (reads the current page from the nav widget)
page = st.sidebar.radio("Technique", PAGES, key="nav", format_func=UI.nav_label, label_visibility="collapsed", on_change=PS.on_nav)
with st.sidebar.container(border=True, key="composer_tray"):          # v2.14: always visible: where the combined figure lives
    _nc = len(st.session_state.composer)
    st.markdown("**🧩 Combined figure**" + (f" · {_nc} panel{'s' if _nc != 1 else ''}" if _nc else ""))
    if _nc:
        st.caption(" · ".join(str(k) for k in list(st.session_state.composer)[:4]) + (" …" if _nc > 4 else ""))
        st.button("Open and arrange", key="tray_open", on_click=UI._go, args=("Figure composer",), width="stretch", type="primary")
    else:
        st.caption("Put graphs of several techniques in ONE figure (A, B, C …): click **🧩 Add to combined figure** under any graph.")
        st.button("How it works", key="tray_open", on_click=UI._go, args=("Figure composer",), width="stretch")
st.sidebar.markdown("---")
from mfs import globalstyle as GS   # v2.11: one figure style for every page (journal presets, fonts, lines, palette)
WIDTH, ASPECT, FONT, DPI = GS.sidebar(os.path.join(os.path.dirname(os.environ.get("MFS_STORE") or os.path.join(os.path.dirname(
    os.environ.get("MFS_EXPORT_DIR") or os.path.join(HERE, "exports")), "sample_store.json")), "figure_style.json"))
width_name = st.session_state.get("gs_width")
with (st.sidebar.expander("Sample → colour rules") if LV.at_least("Advanced") else contextlib.nullcontext()):
  if LV.at_least("Advanced"):
    st.caption("Samples take the palette colours in list order. Add a rule to give every sample whose name matches a pattern "
               "(regex, case-insensitive) a fixed colour in all figures, e.g. `uio-66` → `#0072B2`.")
    if "_cmap_base" not in st.session_state:
        st.session_state["_cmap_base"] = [list(r_) for r_ in S.SAMPLE_COLOURS]
    _cb1, _cb2 = st.columns(2)
    if _cb1.button("Example-set rules", key="_cmap_ex", help="Colour rules for the bundled MOF-801 example data."):
        st.session_state["_cmap_base"] = [list(r_) for r_ in S.EXAMPLE_SAMPLE_COLOURS]; st.session_state.pop("cmap", None)
    if _cb2.button("Clear rules", key="_cmap_clear"):
        st.session_state["_cmap_base"] = []; st.session_state.pop("cmap", None)
    cmap_df = pd.DataFrame(st.session_state["_cmap_base"], columns=["pattern", "colour"]).astype(str)
    cmap_df = st.data_editor(cmap_df, num_rows="dynamic", key="cmap")
    S.SAMPLE_COLOURS[:] = [(str(r.pattern), str(r.colour)) for r in cmap_df.itertuples()
                           if isinstance(r.pattern, str) and r.pattern.strip() and isinstance(r.colour, str) and r.colour.strip().startswith("#")]

if LV.at_least("Advanced"):
    FU.sidebar_house_style()
    st.sidebar.checkbox("📦 Collect figures for the SI bundle", True, key="si_collect",
                        help="Every figure shown is kept (with its data and settings) for Tools → SI export bundle.")
PS.mark_global()
PS.maybe_dialog()

if "composer" not in st.session_state:
    st.session_state.composer = {}


# ============================================================== helpers
import threading
_RENDER_LOCK = threading.Lock()   # matplotlib rendering is not thread-safe; deferred downloads render in a worker thread
EXPORT_DIR = os.environ.get("MFS_EXPORT_DIR") or os.path.join(HERE, "exports")  # the Windows .exe sets Documents\Reticula\exports
STORE_PATH = os.environ.get("MFS_STORE") or os.path.join(os.path.dirname(EXPORT_DIR), "sample_store.json")
MIME = {"xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "docx": "application/vnd.openxmlformats-officedocument.wordprocessingml.document", "png": "image/png", "tiff": "image/tiff", "tif": "image/tiff", "svg": "image/svg+xml", "pdf": "application/pdf",
        "csv": "text/csv", "json": "application/json", "zip": "application/zip", "txt": "text/plain"}


def open_with_default_app(path):
    """Open a local file with the OS default program (Windows / macOS / Linux)."""
    import sys as _sys, subprocess as _sp
    try:
        if _sys.platform.startswith("win"):
            os.startfile(path)              # noqa  (Windows only)
        elif _sys.platform == "darwin":
            _sp.Popen(["open", str(path)])
        else:
            _sp.Popen(["xdg-open", str(path)])
    except Exception as e:
        st.toast(f"Could not open {os.path.basename(str(path))}: {e}")


def safe_name(s: str) -> str:
    """File names that work on Windows/macOS/Linux (no $ \\ / : * ? " < > | or mathtext)."""
    s = re.sub(r"\$|\\mathrm|[{}]", "", str(s))
    s = re.sub(r"[^A-Za-z0-9._@+-]+", "_", s).strip("._")
    return (s or "figure")[:120]


def store_put(sample, section, values, source=""):
    """Put one page's per-sample results into the shared sample store (session + Documents/Reticula/sample_store.json)."""
    from mfs import samplestore as SS
    try:
        SS.put(st.session_state, sample, section, values, source, STORE_PATH if st.session_state.get("store_persist", True) else None,
               files=list(dict.fromkeys(SOURCES)), version=__version__)
    except Exception as e:      # the store must never break a page
        st.session_state["_store_err"] = str(e)


def _store_boot():
    from mfs import samplestore as SS
    if not st.session_state.get("_store_loaded"):
        st.session_state["_store_loaded"] = True
        try:
            SS.load(st.session_state, STORE_PATH)
        except Exception as e:
            st.session_state["_store_err"] = str(e)


def newfig(aspect=None, width=None):
    S.apply_style(FONT)
    return plt.subplots(figsize=(width or WIDTH, (width or WIDTH) * (aspect or ASPECT)))


def export_bytes(fig, fmt, prov=None):
    with _RENDER_LOCK:
        return S.export(fig, fmt, DPI, prov)


SOURCES: list[str] = []          # v2.15: names of the data files of this page run (for the provenance record)


def provenance(stem, cfg=None):
    """v2.15: what is needed to reproduce a figure: program version, page, source files, figure settings, time.
    Built in the main thread (deferred downloads render in a worker thread without session access)."""
    try:
        st_cfg = json.loads(json.dumps(cfg, default=str)) if cfg is not None else None
    except Exception:
        st_cfg = None
    return {"program": "Reticula", "version": __version__, "page": page, "figure": safe_name(stem), "sources": list(dict.fromkeys(SOURCES)),
            "export": {"width_in": WIDTH, "aspect": ASPECT, "font": FONT, "dpi": DPI}, "customisation": st_cfg,
            "created": datetime.datetime.now().isoformat(timespec="seconds")}


def dl(container, label, data, file_name, key, mime=None, help=None, width="content"):
    """Download button that works reliably:
    - data may be bytes/str or a zero-argument callable -> generated only when clicked (deferred), so reruns do not
      have to re-render every format and old files can't expire while the page is busy;
    - on_click='ignore' -> clicking does NOT rerun the app (the v1 rerun invalidated the other buttons' files);
    - explicit MIME type and a file name that is valid on Windows."""
    ext = file_name.rsplit(".", 1)[-1].lower()
    return container.download_button(label, data, file_name=file_name, mime=mime or MIME.get(ext, "application/octet-stream"),
                                     key=key, on_click="ignore", help=help, width=width)


def save_to_folder(fig, stem, data=None, formats=("png", "tiff", "svg", "pdf"), prov=None):
    os.makedirs(EXPORT_DIR, exist_ok=True)
    paths = []
    for fmt in formats:
        p = os.path.join(EXPORT_DIR, f"{stem}.{'tif' if fmt == 'tiff' else fmt}")
        with open(p, "wb") as fh:
            fh.write(export_bytes(fig, fmt, prov))
        paths.append(p)
    if prov:
        p = os.path.join(EXPORT_DIR, f"{stem}_provenance.json")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(prov, indent=1, ensure_ascii=False, default=str))
        paths.append(p)
    if data is not None:
        p = os.path.join(EXPORT_DIR, f"{stem}_processed.csv"); data.to_csv(p, index=False); paths.append(p)
    return paths


def downloads(fig, stem, data: pd.DataFrame | None = None, key=None, cfg=None):
    stem = safe_name(stem); key = key or stem
    prov = provenance(stem, cfg)
    bar = st.container(border=True, key=f"exportbar_{safe_name(key)}")
    bar.markdown('<div class="rt-sec">⬇️ Export</div>', unsafe_allow_html=True)
    cols = bar.columns([1, 1, 1, 1, 1, 1.15, 1.7])
    for c, fmt in zip(cols, ["png", "tiff", "svg", "pdf"]):
        dl(c, fmt.upper(), (lambda f=fmt: export_bytes(fig, f, prov)), f"{stem}.{'tif' if fmt == 'tiff' else fmt}", f"{key}_{fmt}", width="stretch",
           help=f"{'%d dpi ' % DPI if fmt in ('png', 'tiff') else 'vector, '}{fmt.upper()}")
    if data is not None:
        dl(cols[4], "CSV", (lambda: data.to_csv(index=False).encode("utf-8-sig")), f"{stem}_processed.csv", f"{key}_csv", width="stretch",
           help="processed data (UTF-8, opens in Excel)")

    def _zip():
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for fmt in ("png", "tiff", "svg", "pdf"):
                z.writestr(f"{stem}.{'tif' if fmt == 'tiff' else fmt}", export_bytes(fig, fmt, prov))
            z.writestr(f"{stem}_provenance.json", json.dumps(prov, indent=1, ensure_ascii=False, default=str))
            if data is not None:
                z.writestr(f"{stem}_processed.csv", data.to_csv(index=False).encode("utf-8-sig"))
        return buf.getvalue()
    dl(cols[5], "ZIP (all)", _zip, f"{stem}_all_formats.zip", f"{key}_zip", width="stretch")
    if cols[6].button("💾 Save to folder", key=f"{key}_save", help=f"Writes all formats to {EXPORT_DIR}", width="stretch"):
        try:
            paths = save_to_folder(fig, stem, data, prov=prov)
            bar.success("Saved:\n\n" + "\n\n".join(f"`{p}`" for p in paths))
            st.session_state["last_export_dir"] = EXPORT_DIR
        except Exception as e:
            bar.error(f"Could not write to {EXPORT_DIR}: {e}")


CURRENT_SERIES: list[str] = []


def show(fig, stem, data=None, composer_fn=None, title=None, units=None, series=None, controls=True):
    """Every figure of every page goes through here: customisation panel -> preview -> downloads -> composer.
    v2.11: runs as an st.fragment on a pristine copy of the figure, so changing a customisation setting re-renders only
    this figure instead of rerunning the whole page (fallback: plain call if the figure can't be pickled)."""
    series = series if series is not None else (list(CURRENT_SERIES) or None)
    if os.environ.get("MFS_NO_FRAGMENT") != "1":
        import pickle
        try:
            blob = pickle.dumps(fig)
        except Exception:
            blob = None
        if blob is not None:
            plt.close(fig)
            _show_fragment(blob, stem, data, composer_fn, title, units, series, controls)
            return
    _show_impl(fig, stem, data, composer_fn, title, units, series, controls)


def _label_editor(fig, stem, controls=True):
    """v2.12: apply the stored label positions / texts of this figure; with the toggle on, the preview becomes the
    drag-and-drop label editor (returns True then). Overrides live in st.session_state['lbl|<stem>'] (page state,
    project files and undo)."""
    from mfs import labelmove as LM
    sn = safe_name(stem); lk = LM.PREFIX + sn; ss = st.session_state
    PS.track(lk); PS.track_undo(lk)
    ov = ss.get(lk) or {}
    try:
        LM.apply(fig, ov)
    except Exception as e:
        st.warning(f"Label positions not applied: {e}")
    if not controls or not LM.items(fig):
        return False
    on = st.toggle("🏷️ Move / edit labels on the figure", key=f"lmon_{sn}",
                   help="Drag sample labels (and the legend) directly on the figure, double-click to change a label's text. "
                        "Positions are kept relative to each trace's offset; exports use them.")
    if not on:
        if ov:
            st.caption(f"🏷️ {len(ov)} label(s) moved or renamed – turn the switch on to change or reset them.")
        return False
    ev = LM.editor(fig, sn, component_key=f"lm_{sn}")
    if isinstance(ev, dict) and ev.get("key") == sn and ev.get("t") and ev.get("t") != ss.get(lk + "|t"):
        ss[lk + "|t"] = ev["t"]
        ss[lk] = LM.handle(fig, ov, ev)
        try:
            PS.history_record(ss.get("_ps_page"))
        except Exception:
            pass
        try:
            st.rerun(scope="fragment")
        except Exception:
            st.rerun()
    return True


@st.fragment
def _show_fragment(blob, stem, data, composer_fn, title, units, series, controls):
    import pickle
    _show_impl(pickle.loads(blob), stem, data, composer_fn, title, units, series, controls)


def _show_impl(fig, stem, data, composer_fn, title, units, series, controls):
    from mfs import labelmove as LM
    LM.tag_originals(fig)              # v2.12: stable keys for movable labels (before renaming)
    if controls:
        try:
            cfg, fams = FU.panel(fig, safe_name(stem), units=units, series=series)
            F.apply(fig, cfg, fams)
        except Exception as e:  # a bad setting must never cost the user the figure
            st.warning(f"Customisation not applied: {e}")
    from mfs import perf as _PF       # v2.11: unchanged figures are not re-rendered; cheap preview interpolation for big images
    if not _label_editor(fig, stem, controls):
        st.image(_PF.preview_png(fig, 200, _RENDER_LOCK))
    if st.session_state.get("si_collect", True):
        _si_capture(fig, stem, data, title, locals().get("cfg"))
    downloads(fig, stem, data, cfg=locals().get("cfg"))
    if composer_fn is not None:
        _composer_controls(stem, title or stem, composer_fn, locals().get("cfg"), locals().get("fams"))
    plt.close(fig)


def _composer_controls(stem, key, fn, cfg, fams):
    """v2.14: the button under every graph that puts it into the combined figure (and takes it out again). The latest drawing and
    settings of every graph are also remembered in fig_registry, so the composer page can offer ALL graphs seen in this session."""
    ss = st.session_state
    comp = ss.composer
    ss.setdefault("fig_registry", {})[key] = {"fn": fn, "cfg": cfg, "fams": fams, "page": page, "stem": stem}
    if key in comp:
        comp[key] = fn                                     # keep an added panel in step with its data
    inside = key in comp
    c1, c2 = st.columns([2, 5], vertical_alignment="center")
    if c1.button("✔ In the combined figure – remove" if inside else "🧩 Add to combined figure", key=f"add_{safe_name(stem)}",
                 type="secondary" if inside else "primary", width="stretch",
                 help="Collect graphs of different techniques into ONE multi-panel figure (A, B, C …). Open it in the sidebar: 🧩 Combined figure."):
        if inside:
            comp.pop(key, None)
        else:
            comp[key] = fn
        st.toast(("Removed from" if inside else "Added to") + f" the combined figure ({len(comp)} panel{'s' if len(comp) != 1 else ''}) – open it in the sidebar.", icon="🧩")
        st.rerun()
    n = len(comp)
    c2.caption(f"🧩 {n} panel{'s' if n != 1 else ''} in the combined figure – sidebar → *Open and arrange*." if n
               else "🧩 Combine graphs of several techniques into one figure (A, B, C …): add this one, then go to another technique.")


_SI_SKIP = re.compile(r"^(fc_|sc_|data_|exportbar_|add_|si_|db_|bp_|FormSubmitter|nav$|cmap$|house_|store_|_)")


def _si_capture(fig, stem, data, title, cfg):
    """Keeps the finished figure + data + the scalar widget values of this run for the SI bundle (overwrites by stem)."""
    try:
        _pg = "ak|" + PS.slug(page) + "|"
        params = {(k.split("|")[3] if k.startswith("ak|") and k.count("|") >= 4 else k): v for k, v in st.session_state.items()
                  if isinstance(k, str) and not _SI_SKIP.match(k) and (not k.startswith("ak|") or k.startswith(_pg))
                  and isinstance(v, (int, float, str, bool)) and not (isinstance(v, str) and len(v) > 200)}
        st.session_state.setdefault("si_figs", {})[safe_name(stem)] = {
            "fig": fig, "data": data if isinstance(data, pd.DataFrame) else None, "page": page, "title": title, "params": params,
            "cfg": cfg, "series": list(CURRENT_SERIES), "time": datetime.datetime.now().strftime("%H:%M:%S"), "n_rows": int(len(data)) if isinstance(data, pd.DataFrame) else None}
    except Exception as e:      # never break a page
        st.session_state["_si_err"] = str(e)


GUESS = [(r"mof-?801hf_3flex_full", "MOF-801(Hf) batch 3 (3Flex report)"), (r"zro2.?@|composite", "ZrO$_2$@Sulfo-MOF-801(Hf)"), (r"^pm|physical", "PM"), (r"sulfo-?mof|sulfo_mof", "Sulfo-MOF-801(Hf)"),
         (r"^zro2", "ZrO$_2$"), (r"mof-?801", "MOF-801(Hf)"), (r"sulfosucc", "2-Sulfosuccinic acid")]


_EXAMPLE_FILES: set | None = None


def _example_files():
    global _EXAMPLE_FILES
    if _EXAMPLE_FILES is None:
        _EXAMPLE_FILES = {f_ for _r, _d, fs_ in os.walk(EX) for f_ in fs_}
    return _EXAMPLE_FILES


def guess_name(fn):
    """Sample name from a file name: the file stem. Only the bundled example files get their display names
    (e.g. 'MOF801_STOE.xy' -> 'MOF-801(Hf)'); a user's own files are never renamed by guesswork."""
    base = os.path.basename(fn)
    stem = os.path.splitext(base)[0]
    if base in _example_files():
        for p, n in GUESS:
            if re.search(p, stem, re.I):
                return n
    return stem


class _F(io.BytesIO):
    def __init__(self, path):
        super().__init__(open(path, "rb").read()); self.name = os.path.basename(path)


def files_input(label, types, example_glob, key, multiple=True):
    box = st.container(border=True, key=f"data_{key}")
    box.markdown('<div class="rt-sec">📁 Data</div>', unsafe_allow_html=True)
    up = box.file_uploader(label, type=types, accept_multiple_files=multiple, key=key)
    use_ex = box.checkbox("Use the bundled EXAMPLE data (synthetic)", value=not up, key=key + "_ex")
    files = (up if multiple else [up]) if up else []
    SOURCES.extend(getattr(f, "name", "") for f in files if f is not None)
    if use_ex and not files:
        files = [_F(p) for p in sorted(glob.glob(os.path.join(EX, example_glob)))]
        SOURCES.extend(f.name + " (bundled example)" for f in files)
    return [f for f in files if f is not None]


def sample_table(files, key):
    """v2.11: the unified sample list (show/hide, colour, rename, drag to reorder; order = stacking + legend order).
    Returns [(sample name, file)] of the shown samples in list order."""
    from mfs import samplelist as SL
    items = SL.sample_list(files, key, [guess_name(f.name) for f in files])
    lut = {f.name: f for f in files}
    out = [(n, lut[i]) for i, n in SL.visible(items) if i in lut]
    CURRENT_SERIES.extend(n for n, _ in out if n not in CURRENT_SERIES)
    return out


def sample_names_list(names, key, label="🧪 Samples – show/hide, colour, rename, drag to reorder", max_shown=None):
    """The same list for pages whose samples are rows/columns of one table (kinetics, zeta, ICP …).
    Returns [(original name, display name)] of the shown samples in list order."""
    from mfs import samplelist as SL
    names = [str(n) for n in names]
    items = SL.sample_list(names, key, names, label=label, max_shown=max_shown)
    out = SL.visible(items)
    CURRENT_SERIES.extend(n for _, n in out if n not in CURRENT_SERIES)
    return out


def _apply_sample_list(df, col, key, sel=None):
    """Filter / order / rename the rows of a table by the unified sample list over the values of `col`
    (pass `sel` = sample_names_list(...) to apply one list to several tables)."""
    if sel is None:
        sel = sample_names_list(list(pd.unique(df[col].astype(str))), key)
    ren = dict(sel); order = {k: i for i, (k, _) in enumerate(sel)}
    out = df[df[col].astype(str).isin(ren)].copy()
    out["_sl_ord"] = out[col].astype(str).map(order)
    out = out.sort_values("_sl_ord", kind="stable").drop(columns="_sl_ord")
    out[col] = out[col].astype(str).map(ren)
    return out


def offset_controls(col, key, default, lo, hi, step, label="Offset between traces", even=None, help=None):
    """Offset slider + '⇅ Stack evenly' (sets the offset to the value `even` computed by the page on the last run)."""
    k = key
    ev = st.session_state.get("_even|" + k)

    def _even():
        v = st.session_state.get("_even|" + k)
        if v is not None:
            st.session_state[k] = float(min(max(v, lo), hi))
    v = col.slider(label, lo, hi, default, step, key=k, help=help)
    col.button("⇅ Stack evenly", key=k + "_even", on_click=_even, disabled=ev is None,
               help="Offset = the largest trace height (+5 %), so the traces never overlap and are equally spaced")
    return v


def set_even_offset(key, traces, pad=1.05, xlim=None):
    """Remember the even-stacking offset for the next '⇅ Stack evenly' click.
    traces: [(name, x, y)] / [(x, y)] / [y]; only the part inside xlim counts."""
    try:
        h = []
        for t in traces:
            if isinstance(t, (tuple, list)) and len(t) >= 3:
                x, y = np.asarray(t[1], float), np.asarray(t[2], float)
            elif isinstance(t, (tuple, list)) and len(t) == 2:
                x, y = np.asarray(t[0], float), np.asarray(t[1], float)
            else:
                y = np.asarray(t, float); x = np.arange(len(y))
            if xlim is not None and len(x) == len(y):
                m = (x >= min(xlim)) & (x <= max(xlim))
                if m.any():
                    y = y[m]
            if y.size and np.isfinite(y).any():
                h.append(float(np.nanmax(y) - np.nanmin(y)))
        h = [v for v in h if np.isfinite(v) and v > 0]
        if h:
            st.session_state["_even|" + key] = round(max(h) * pad, 3)
    except Exception:
        pass


def proc_controls(key, default_baseline="none", default_norm="max", norm_opts=("max", "min-max", "area", "none")):
    c1, c2, c3, c4 = st.expander("⚙️ Processing – baseline, smoothing, normalisation, offset", expanded=True).columns(4)
    bl = c1.selectbox("Baseline", ["none", "ALS", "rolling minimum", "polynomial", "linear (ends)"],
                      index=["none", "ALS", "rolling minimum", "polynomial", "linear (ends)"].index(default_baseline), key=key + "bl")
    lam = c1.number_input("ALS λ", 1e2, 1e9, 1e6, format="%.0e", key=key + "lam") if bl == "ALS" else 1e6
    sm = c2.selectbox("Smoothing", ["none", "Savitzky-Golay", "moving average"], key=key + "sm")
    win = c2.slider("Window (points)", 3, 51, 11, 2, key=key + "win") if sm != "none" else 11
    nm = c3.selectbox("Normalisation", list(norm_opts), index=list(norm_opts).index(default_norm), key=key + "nm")
    off = offset_controls(c4, key + "off", 1.05, 0.0, 2.0, 0.05)
    return dict(baseline=bl, lam=lam, smooth_method=sm, window=win, norm=nm), off


# ============================================================== pages
def page_home():
    from mfs import home as HM           # v2.16: guided, coloured home page
    HM.render(HERE, __version__)


def page_nmr():
    from mfs import nmr_page as NMRP          # v2.13: characterisation + digestion modes
    NMRP.render(globals())


# -------------------------------------------------------------- defect calculator


# ============================================================== Batch D: batch processing, SI bundle, sample database


def page_uvvis():
    from mfs import uvvis_page as UVP
    UVP.render(globals())


def page_pl():
    from mfs import pl_page as PLP          # v2.14: photoluminescence (steady-state, temperature-dependent, thermometry)
    PLP.render(globals())


UI.header(HERE, page)
try:                                   # v2.11: opt-in update check (asks once; background thread; silent offline)
    from mfs import update as _UPD
    _UPD.banner(os.path.dirname(STORE_PATH), __version__)
except Exception:
    pass
_store_boot()
PS.start_page(page)
from mfs import pages as PG        # v2.13: the page code lives in mfs/pages/*.py (moved verbatim out of this file)
from mfs.pages.spectra import page_pxrd, page_ftir
from mfs.pages.nmr import page_nmr_digest
from mfs.pages.sorption import page_sorption
from mfs.pages.thermal import page_tga, page_dsc
from mfs.pages.defects import page_defects
from mfs.pages.kinetics import page_kinetics
from mfs.pages.elemental import page_icp
from mfs.pages.microscopy import page_sem
from mfs.pages.bioassays import page_sds, page_zeta, page_xps
from mfs.pages.structure import page_structure
from mfs.pages.workflow import page_composer, page_batch, page_si, page_db, page_guide
PG.bind(globals())                 # hand the helpers and this run's sidebar values (WIDTH, FONT …) to the page modules
_PAGEFN = {"Home": page_home, "PXRD": page_pxrd, "FTIR / Raman": page_ftir, "1H NMR (digest)": page_nmr, "N2 sorption (BET/PSD)": page_sorption,
 "TGA": page_tga, "DSC": page_dsc, "Defect calculator": page_defects, "Kinetics & adsorption": page_kinetics, "ICP-OES / XRF": page_icp, "SEM / TEM / EDX": page_sem,
 "SDS-PAGE (Image Lab)": page_sds, "Zeta / DLS": page_zeta, "XPS": page_xps, "Figure composer": page_composer, "Batch processing": page_batch,
 "SI export bundle": page_si, "Sample database": page_db, "Structure (CIF)": page_structure, "UV/Vis": page_uvvis, "Photoluminescence": page_pl, "Export guide": page_guide}
try:
    if os.environ.get("MFS_PROFILE"):        # dev: cProfile of the page function -> <MFS_PROFILE>_<page>.prof
        import cProfile as _cp
        _pr = _cp.Profile(); _pr.enable()
        try:
            _PAGEFN[page]()
        finally:
            _pr.disable(); _pr.dump_stats(f"{os.environ['MFS_PROFILE']}_{PS.slug(page)}.prof")
    else:
        try:
            _PAGEFN[page]()
        except Exception as _exc:                # v2.16: remember the error for the problem report, then show it as before
            from mfs import problem_report as _PR
            _PR.record_error(st.session_state, page, _exc)
            raise
finally:
    if os.environ.get("MFS_TIMING"):
        try:
            with open(os.environ["MFS_TIMING"], "a", encoding="utf-8") as _fh:
                _fh.write(json.dumps({"page": page, "s": round(_time.perf_counter() - _T0, 3), "at": round(_time.time(), 2)}) + "\n")
        except Exception:
            pass
    PS.end(page)
    PS.history_record(page)       # v2.11 undo/redo
    PS.history_box(page)
    PS.status_badge(page)
    UI.about_box(__version__, os.path.dirname(STORE_PATH))
    PRJ.sidebar(page, __version__)
    with st.sidebar.expander("🐞 Report a problem", expanded=bool(st.session_state.get("_last_error"))):
        from mfs import problem_report as _PR
        _le = st.session_state.get("_last_error")
        if _le:
            st.error(f"Last error on **{_le['page']}**: {_le['type']}: {_le['message'][:160]}")
        st.caption("Creates a small ZIP with the version, your settings and the last error - **no data**. Send it to the developer.")
        _note = st.text_area("What were you doing? (optional)", key="_pr_note", height=70)
        dl(st, "⬇️ Create the problem report", lambda: _PR.build(st.session_state, page, __version__, SOURCES, _note),
           f"Reticula_problem_report_{datetime.datetime.now():%Y%m%d_%H%M}.zip", "_pr_dl", width="stretch")
    with st.sidebar.expander("💾 Saved page work", expanded=False):
        sp = PS.saved_pages()
        st.caption("Saved with the dialog or 💾 Save; restored on the next start. " + (", ".join(f"{k} ({v[5:16].replace('T', ' ')})"
                   for k, v in sp.items() if v) if sp else "Nothing saved yet."))
        if sp and st.button("Forget all saved page work", key="_ps_forget_all"):
            PS.forget_all(); st.rerun()

