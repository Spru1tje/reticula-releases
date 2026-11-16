"""Reticula (formerly MOF Figure Studio) - publication figures for framework materials: one app for every Paper-1 characterisation technique + the defect calculator.
Run:  streamlit run app.py"""
from __future__ import annotations
import time as _time
_T0 = _time.perf_counter()     # v2.11: per-run timing (MFS_TIMING=<file> appends one JSON line per script run)
import io, os, re, json, glob, zipfile, warnings, datetime, hashlib
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
page = st.sidebar.radio("Technique", PAGES, key="nav", format_func=UI.nav_label, label_visibility="collapsed", on_change=PS.on_nav)
st.sidebar.markdown("---")
from mfs import globalstyle as GS   # v2.11: one figure style for every page (journal presets, fonts, lines, palette)
WIDTH, ASPECT, FONT, DPI = GS.sidebar(os.path.join(os.path.dirname(os.environ.get("MFS_STORE") or os.path.join(os.path.dirname(
    os.environ.get("MFS_EXPORT_DIR") or os.path.join(HERE, "exports")), "sample_store.json")), "figure_style.json"))
width_name = st.session_state.get("gs_width")
with st.sidebar.expander("Sample → colour mapping"):
    st.caption("Fixed Okabe–Ito colours, matched on the sample name (regex). Edit to override.")
    cmap_df = pd.DataFrame(S.SAMPLE_COLOURS, columns=["pattern", "colour"])
    cmap_df = st.data_editor(cmap_df, num_rows="dynamic", key="cmap")
    S.SAMPLE_COLOURS[:] = [(r.pattern, r.colour) for r in cmap_df.itertuples() if r.pattern]

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
        SS.put(st.session_state, sample, section, values, source, STORE_PATH if st.session_state.get("store_persist", True) else None)
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


def export_bytes(fig, fmt):
    with _RENDER_LOCK:
        return S.export(fig, fmt, DPI)


def dl(container, label, data, file_name, key, mime=None, help=None, width="content"):
    """Download button that works reliably:
    - data may be bytes/str or a zero-argument callable -> generated only when clicked (deferred), so reruns do not
      have to re-render every format and old files can't expire while the page is busy;
    - on_click='ignore' -> clicking does NOT rerun the app (the v1 rerun invalidated the other buttons' files);
    - explicit MIME type and a file name that is valid on Windows."""
    ext = file_name.rsplit(".", 1)[-1].lower()
    return container.download_button(label, data, file_name=file_name, mime=mime or MIME.get(ext, "application/octet-stream"),
                                     key=key, on_click="ignore", help=help, width=width)


def save_to_folder(fig, stem, data=None, formats=("png", "tiff", "svg", "pdf")):
    os.makedirs(EXPORT_DIR, exist_ok=True)
    paths = []
    for fmt in formats:
        p = os.path.join(EXPORT_DIR, f"{stem}.{'tif' if fmt == 'tiff' else fmt}")
        with open(p, "wb") as fh:
            fh.write(export_bytes(fig, fmt))
        paths.append(p)
    if data is not None:
        p = os.path.join(EXPORT_DIR, f"{stem}_processed.csv"); data.to_csv(p, index=False); paths.append(p)
    return paths


def downloads(fig, stem, data: pd.DataFrame | None = None, key=None):
    stem = safe_name(stem); key = key or stem
    bar = st.container(border=True, key=f"exportbar_{safe_name(key)}")
    bar.markdown('<div class="rt-sec">⬇️ Export</div>', unsafe_allow_html=True)
    cols = bar.columns([1, 1, 1, 1, 1, 1.15, 1.7])
    for c, fmt in zip(cols, ["png", "tiff", "svg", "pdf"]):
        dl(c, fmt.upper(), (lambda f=fmt: export_bytes(fig, f)), f"{stem}.{'tif' if fmt == 'tiff' else fmt}", f"{key}_{fmt}", width="stretch",
           help=f"{'%d dpi ' % DPI if fmt in ('png', 'tiff') else 'vector, '}{fmt.upper()}")
    if data is not None:
        dl(cols[4], "CSV", (lambda: data.to_csv(index=False).encode("utf-8-sig")), f"{stem}_processed.csv", f"{key}_csv", width="stretch",
           help="processed data (UTF-8, opens in Excel)")

    def _zip():
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for fmt in ("png", "tiff", "svg", "pdf"):
                z.writestr(f"{stem}.{'tif' if fmt == 'tiff' else fmt}", export_bytes(fig, fmt))
            if data is not None:
                z.writestr(f"{stem}_processed.csv", data.to_csv(index=False).encode("utf-8-sig"))
        return buf.getvalue()
    dl(cols[5], "ZIP (all)", _zip, f"{stem}_all_formats.zip", f"{key}_zip", width="stretch")
    if cols[6].button("💾 Save to folder", key=f"{key}_save", help=f"Writes all formats to {EXPORT_DIR}", width="stretch"):
        try:
            paths = save_to_folder(fig, stem, data)
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


@st.fragment
def _show_fragment(blob, stem, data, composer_fn, title, units, series, controls):
    import pickle
    _show_impl(pickle.loads(blob), stem, data, composer_fn, title, units, series, controls)


def _show_impl(fig, stem, data, composer_fn, title, units, series, controls):
    if controls:
        try:
            cfg, fams = FU.panel(fig, safe_name(stem), units=units, series=series)
            F.apply(fig, cfg, fams)
        except Exception as e:  # a bad setting must never cost the user the figure
            st.warning(f"Customisation not applied: {e}")
    from mfs import perf as _PF       # v2.11: unchanged figures are not re-rendered; cheap preview interpolation for big images
    st.image(_PF.preview_png(fig, 200, _RENDER_LOCK))
    if st.session_state.get("si_collect", True):
        _si_capture(fig, stem, data, title, locals().get("cfg"))
    downloads(fig, stem, data)
    if composer_fn is not None:
        if st.button("➕ Add this panel to the figure composer", key=f"add_{safe_name(stem)}"):
            st.session_state.composer[title or stem] = composer_fn
            st.success(f"Added '{title or stem}' ({len(st.session_state.composer)} panels in composer)")
    plt.close(fig)


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
            "cfg": cfg, "time": datetime.datetime.now().strftime("%H:%M:%S"), "n_rows": int(len(data)) if isinstance(data, pd.DataFrame) else None}
    except Exception as e:      # never break a page
        st.session_state["_si_err"] = str(e)


GUESS = [(r"mof-?801hf_3flex_full", "MOF-801(Hf) batch 3 (3Flex report)"), (r"zro2.?@|composite", "ZrO$_2$@Sulfo-MOF-801(Hf)"), (r"^pm|physical", "PM"), (r"sulfo-?mof|sulfo_mof", "Sulfo-MOF-801(Hf)"),
         (r"^zro2", "ZrO$_2$"), (r"mof-?801", "MOF-801(Hf)"), (r"sulfosucc", "2-Sulfosuccinic acid")]


def guess_name(fn):
    stem = os.path.splitext(os.path.basename(fn))[0]
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
    if use_ex and not files:
        files = [_F(p) for p in sorted(glob.glob(os.path.join(EX, example_glob)))]
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
    UI.home(HERE)
    return
    st.markdown("""
Turn raw characterisation exports into **journal-ready figures** (ACS / Wiley-VCH style: Arial/Helvetica 7–9 pt,
3.25 in or 7 in wide, inward ticks, no grid, Okabe–Ito colour-blind-safe palette, units in parentheses) and quantify
**defects** in any MOF (Zr/Hf/Ce/Th₆, paddlewheels, Zn₄O, Al-OH chains, M₃O trimers, MOF-74, Ti₈, ZIFs or a custom node) with every established method side by side.

**Workflow**: pick a technique on the left → upload your files (or tick *Use the bundled EXAMPLE data*) → adjust
processing → download PNG/TIFF (600 dpi), SVG (editable text), PDF and the processed CSV → optionally *Add to the
figure composer* to build a multi-panel figure with (a), (b), (c) labels.

The bundled example datasets are **synthetic EXAMPLE DATA** written in each instrument's export format for testing only.
See **Export guide** for exactly which file to export from each instrument's software.
""")
    st.info("Fixed colours: MOF-801 parent = black, Sulfo = blue, ZrO$_2$ = orange, physical mixture (PM) = purple, "
            "ZrO$_2$@Sulfo composite = vermillion, simulated/reference = grey.")


def page_pxrd():
    st.header("Powder X-ray diffraction")
    st.caption("STOE WinXPOW ASCII export (.xy / .xye / .dat), any 2θ–intensity text. Binary STOE .raw only if the optional `xylib` is installed.")
    files = files_input("PXRD files", ["xy", "xye", "dat", "txt", "csv", "raw", "asc", "rasx", "ras"], "PXRD/*.xy", "pxrd")
    if not files:
        return
    samples = sample_table(files, "pxrd_tab")
    opts, off = proc_controls("pxrd", "ALS", "max")
    c1, c2, c3 = st.columns(3)
    lo, hi = c1.slider("2θ range (°)", 2.0, 90.0, (4.0, 50.0), 0.5)
    ref_mode = c2.selectbox("Reference pattern", ["Simulate from bundled MOF-801 CIF (COD 4121459)", "Simulate from bundled UiO-66 CIF (COD 4512074)",
                                                  "Upload CIF", "Upload .xy reference", "none"])
    subst = c3.selectbox("Metal in simulation", ["Hf (replace Zr)", "Zr (as in CIF)"])
    fwhm = c3.number_input("Simulated FWHM (°)", 0.03, 1.0, 0.12, 0.01)
    sticks_on = c2.checkbox("Show reflection sticks", True)
    ref, sticks = None, None
    try:
        if ref_mode.startswith("Simulate") or ref_mode == "Upload CIF":
            if ref_mode == "Upload CIF":
                cf = st.file_uploader("CIF", type=["cif"], key="cif")
                cif = cf.read() if cf else None
            else:
                cif = os.path.join(REF, "MOF-801_Zr_COD4121459.cif" if "801" in ref_mode else "UiO-66_Zr_COD4512074.cif")
            if cif:
                sim, sticks, formula = A.simulate_pxrd(cif, (lo, hi), fwhm=fwhm, substitute={"Zr": "Hf"} if subst.startswith("Hf") else None)
                name = st.text_input("Reference label", "Simulated " + ("MOF-801" if "801" in ref_mode else "UiO-66" if "UiO" in ref_mode else formula) + ("(Hf)" if subst.startswith("Hf") else "(Zr)"))
                ref = (name, sim.x.values, sim.y.values)
        elif ref_mode == "Upload .xy reference":
            rf = st.file_uploader("Reference .xy", key="refxy")
            if rf:
                d = R.read_xy(rf); ref = (st.text_input("Reference label", "Reference"), d.x.values, d.y.values)
    except Exception as e:
        st.warning(f"Reference simulation failed: {e}")
    traces, rows = [], []
    for nm, f in samples:
        try:
            d = R.read_pxrd(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        rg = d.attrs.get("rigaku")
        extra = []
        if rg:                                   # Rigaku SmartLab .rasx / .ras: metadata → store, flag non-2θ/θ scans
            from mfs import rigaku as RG
            allsc = d.attrs.get("rigaku_all") or [rg]
            m_ = rg["meta"]
            st.caption(f"{f.name}: Rigaku SmartLab – {len(allsc)} scan(s); target {m_.get('target', '?')}, Kα1 "
                       f"{m_.get('wavelength_ka1', '?')} Å, axis {rg['axis']}, {m_.get('start', '?')}–{m_.get('end', '?')}°, step "
                       f"{m_.get('step', '?')}°, speed {m_.get('speed', '?')} {m_.get('speed_unit', '')}; intensity = counts × attenuator "
                       f"(max factor {float(np.max(rg['att'])):.4g}).")
            for s_ in allsc:
                if not s_["is_2theta_theta"]:
                    st.warning(f"{f.name} · {s_['name']}: scan axis is **{s_['axis']}**, not 2θ/θ – "
                               + ("not plotted as a powder pattern." if s_ is not rg else "plotted anyway (no 2θ/θ scan in the file), check before use."))
            extra = [s_ for s_ in allsc if s_ is not rg and s_["is_2theta_theta"]]
            store_put(nm, "pxrd", RG.store_values(rg), "PXRD (Rigaku)")
        y = A.process(d.x.values, d.y.values, **opts)
        traces.append((nm, d.x.values, y))
        rows.append(pd.DataFrame({"sample": nm, "two_theta": d.x, "I_raw": d.y, "I_processed": y}))
        for k_, s_ in enumerate(extra, 2):
            ye = A.process(s_["x"], s_["y"], **opts)
            traces.append((f"{nm} · scan {k_}", s_["x"], ye))
            rows.append(pd.DataFrame({"sample": f"{nm} · scan {k_}", "two_theta": s_["x"], "I_raw": s_["y"], "I_processed": ye}))
    st.session_state["_even|pxrdoff"] = 1.05      # traces are min-max scaled to 1 in the window -> even = 1.05
    def draw(ax, traces=traces, ref=ref, sticks=sticks):
        P.draw_stacked(ax, traces, "pxrd", offset=off, xlim=(lo, hi), reference=ref, sticks=sticks if sticks_on else None)
    fig, ax = newfig(); draw(ax)
    show(fig, "PXRD", pd.concat(rows) if rows else None, draw, "PXRD", units={"x": "2theta"})
    from mfs import xrd as X
    with st.expander("📏 Peak width / Scherrer crystallite size (with instrumental correction)"):
        if traces:
            c1, c2, c3, c4 = st.columns(4)
            pks = c1.text_input("Peaks (°2θ, comma separated)", "8.6, 9.9", key="sch_pk")
            Kf = c2.number_input("Shape factor K", 0.5, 1.5, 0.9, 0.01, key="sch_K")
            shape = c3.selectbox("Width correction", ["gaussian", "lorentzian", "voigt"], key="sch_shape",
                                 help="gaussian: β² = B² − b² · lorentzian: β = B − b · voigt: β = B − b²/B")
            st.caption("Defaults: MOF-801 (111) 8.6° and (200) 9.9° (Cu Kα1); UiO-66: 7.4° and 8.5°.")
            inst_src = c4.radio("Instrumental FWHM", ["constant", "from a standard (LaB₆ / Si) pattern"], key="sch_src")
            pk_list = [float(v) for v in pks.replace(";", ",").split(",") if v.strip()]
            inst = lambda tt: 0.05
            if inst_src == "constant":
                b0 = st.number_input("Instrumental FWHM (°)", 0.0, 1.0, 0.05, 0.005, key="sch_b")
                inst = lambda tt, b0=b0: b0
            else:
                sf = st.file_uploader("Standard pattern (.xy)", key="sch_std")
                spk = st.text_input("Standard peaks (°2θ)", "21.36, 30.38, 37.44, 43.51" , key="sch_spk",
                                    help="LaB₆ SRM 660 (Cu Kα1): 21.36 (100), 30.38 (110), 37.44 (111), 43.51 (200) …")
                if sf:
                    try:
                        dstd = R.read_xy(sf)
                        inst, pos, fw = X.instrument_fwhm_from_standard(dstd.x.values, dstd.y.values, [float(v) for v in spk.split(",") if v.strip()])
                        st.caption("Standard fitted: " + ", ".join(f"{p_:.2f}° → {f_:.3f}°" for p_, f_ in zip(pos, fw)) + " (linear FWHM(2θ))")
                    except Exception as e:
                        st.warning(f"Standard fit failed: {e}")
                else:
                    st.caption("No standard loaded – 0.05° used.")
            srows = []
            for nm, x, y in traces:
                for g in pk_list:
                    try:
                        r = X.fit_peak(x, y, g)
                        sz = X.scherrer_size(r["center"], r["fwhm"], inst(r["center"]), K=Kf, shape=shape)
                        srows.append({"sample": nm, "2θ (°)": r["center"], "FWHM obs (°)": r["fwhm"], "± (°)": r["fwhm_err"],
                                      "FWHM inst (°)": inst(r["center"]), "β sample (°)": sz["beta_sample"], "D (nm)": sz["size_nm"], "note": sz["note"]})
                    except Exception as e:
                        srows.append({"sample": nm, "2θ (°)": g, "note": f"fit failed ({e})"})
            st.dataframe(pd.DataFrame(srows).style.format(precision=3), width="stretch", hide_index=True)
            st.caption("pseudo-Voigt fit with linear background; D = Kλ/(β cos θ), λ = 1.540598 Å (Cu Kα1). Sizes > ~100 nm are not reliable.")
    with st.expander("💎 Crystallinity index vs a reference sample"):
        if len(traces) >= 2:
            names = [t[0] for t in traces]
            c1, c2 = st.columns(2)
            rname = c1.selectbox("Reference (100 % crystalline)", names, key="ci_ref")
            win = c2.text_input("Reflection windows (°2θ, lo-hi; …)", "8.2-9.0; 9.5-10.3; 14.6-15.4", key="ci_win")
            ws = [tuple(float(v) for v in w.split("-")) for w in win.split(";") if "-" in w]
            rt = traces[names.index(rname)]
            raw = {r_["sample"].iloc[0]: r_ for r_ in rows}
            use_raw = st.checkbox("use the raw intensities (recommended: same mass / holder / time); untick to use the processed curves", True, key="ci_raw")
            ci = []
            for nm, x, y in traces:
                yy = raw[nm].I_raw.values if use_raw else y; yr = raw[rname].I_raw.values if use_raw else rt[2]
                r_ = X.crystallinity_index(x, yy, rt[1], yr, ws)
                ci.append({"sample": nm, "CI (%)": r_["CI_pct"], "Σ area sample": r_["A_sample"], "Σ area reference": r_["A_ref"]})
            st.dataframe(pd.DataFrame(ci).style.format(precision=1), width="stretch", hide_index=True)
        else:
            st.caption("Load at least two patterns (sample(s) + a reference).")
    with st.expander("↔️ Peak-shift tracking (Δ2θ, d-spacing, cubic lattice parameter)"):
        if traces:
            c1, c2, c3 = st.columns(3)
            tp = c1.text_input("Peaks (°2θ)", "8.6, 9.9", key="sh_pk")
            hk = c2.text_input("hkl (optional, same order)", "111, 200", key="sh_hkl")
            rsel = c3.selectbox("Reference trace", [t[0] for t in traces], key="sh_ref")
            pl = [float(v) for v in tp.split(",") if v.strip()]
            hl = [tuple(int(ch) for ch in h.strip()) if h.strip().isdigit() and len(h.strip()) == 3 else None for h in hk.split(",")]
            tr = X.track_shifts(traces, pl, ref_index=[t[0] for t in traces].index(rsel), hkl=hl)
            tdf = pd.DataFrame(tr)
            st.dataframe(tdf.style.format(precision=4), width="stretch", hide_index=True)

            def draw_sh(ax, tdf=tdf):
                for i_, (g, gg) in enumerate(tdf.groupby("peak (°2θ, guess)", sort=False)):
                    ax.errorbar(range(len(gg)), 1000 * gg["Δ2θ vs ref (°)"], yerr=1000 * gg["± (°)"], fmt="o-", ms=4, lw=0.8, capsize=2, label=f"{g:g}°")
                ax.axhline(0, color="0.6", lw=0.5)
                ax.set_xticks(range(len(traces))); ax.set_xticklabels([t[0] for t in traces], rotation=30, ha="right")
                ax.set_ylabel(r"$\Delta 2\theta$ vs reference (m°)"); ax.legend(frameon=False)
            fig, ax = newfig(); draw_sh(ax)
            show(fig, "PXRD_peak_shifts", tdf, draw_sh, "Peak shifts")


def page_ftir():
    st.header("FTIR (ATR) / FT-Raman")
    st.caption("Bruker OPUS binary (.0, .1, …) via brukeropusreader, or OPUS *Save As → Data Point Table* (.dpt) / CSV.")
    mode = st.radio("Spectrum type", ["ATR-FTIR (transmittance)", "ATR-FTIR (absorbance)", "FT-Raman (1064 nm)"], horizontal=True)
    files = files_input("Spectra", None, "FTIR/*.dpt" if "FTIR" in mode else "Raman/*.dpt", "ftir")
    if not files:
        return
    samples = sample_table(files, "ftir_tab")
    opts, off = proc_controls("ftir", "none" if "FTIR" in mode else "ALS", "min-max")
    c1, c2 = st.columns(2)
    lo, hi = c1.slider("Range (cm⁻¹)", 100.0, 4000.0, (400.0, 2000.0) if "FTIR" in mode else (150.0, 1800.0), 10.0)
    conv = c2.checkbox("Convert transmittance → absorbance", False) if "transmittance" in mode else False
    hl = c2.text_input("Highlight bands (lo-hi:label; …)", "1700-1740:C=O; 1190-1230:S=O; 1030-1055:" if "FTIR" in mode else "")
    highlight = []
    for part in [p for p in hl.split(";") if p.strip()]:
        try:
            rng_, lab = part.split(":") if ":" in part else (part, "")
            a, b = [float(v) for v in rng_.split("-")]; highlight.append((a, b, lab.strip()))
        except Exception:
            pass
    traces, rows = [], []
    for nm, f in samples:
        try:
            d = R.read_ftir(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        y = A.transmittance_to_absorbance(d.y.values) if conv else d.y.values
        y = A.process(d.x.values, y, **opts)
        traces.append((nm, d.x.values, y)); rows.append(pd.DataFrame({"sample": nm, "wavenumber": d.x, "y_processed": y}))
    kind = "raman" if "Raman" in mode else ("ftir_A" if ("absorbance" in mode or conv) else "ftir_T")
    from mfs import spectra as SP
    is_T = kind == "ftir_T"
    with st.expander("🧪 Reference band, peak labels, zoomed inset, difference spectra", expanded=False):
        c1, c2, c3, c4 = st.columns(4)
        nb = c1.checkbox("Normalise to a reference band", False, key="ft_nb")
        dlo, dhi, dlab = (1540.0, 1620.0, "νas(COO⁻)") if kind != "raman" else (1630.0, 1690.0, "ν(C=C)")
        blo = c2.number_input("band from (cm⁻¹)", 50.0, 4000.0, dlo, key="ft_blo")
        bhi = c3.number_input("band to (cm⁻¹)", 50.0, 4000.0, dhi, key="ft_bhi")
        bmode = c4.radio("by", ["height", "area"], horizontal=True, key="ft_bm")
        c1, c2, c3, c4 = st.columns(4)
        pl = c1.checkbox("Label peaks", False, key="ft_pl")
        prom = c2.slider("min. prominence (fraction of range)", 0.01, 0.5, 0.08, 0.01, key="ft_prom")
        npk = c3.number_input("max labels per trace", 1, 40, 8, key="ft_npk")
        npk0 = {t[0]: len(SP.pick_peaks(t[1], t[2], prom, 40, (lo, hi), minima=is_T)) for t in traces}
        pl_which = c4.multiselect("traces", [t[0] for t in traces], default=[max(npk0, key=npk0.get)] if npk0 else [], key="ft_plw")
        c1, c2, c3 = st.columns(3)
        ins = c1.checkbox("Zoomed inset", False, key="ft_ins")
        zr = c2.text_input("inset range (cm⁻¹, lo-hi)", "1000-1250" if kind != "raman" else "200-700", key="ft_zr")
        zloc = c3.selectbox("inset position", ["upper left", "upper center", "upper right", "center left", "lower left"], key="ft_zloc")
        st.caption(f"Reference band suggestion: {dlab}. Normalisation keeps the relative intensities between traces (one common scale).")
    if nb:
        nt, weak = [], []
        for nm, x, y in traces:
            m_ = (x >= lo) & (x <= hi)
            try:
                v_ = SP.band_value(x, -y if is_T else y, blo, bhi, bmode if bmode == "height" else "height")
                if v_ < 0.05 * np.nanmax(np.abs(np.ptp(y[m_]))):
                    raise ValueError("band absent (< 5 % of the trace range)")
                nt.append((nm, x, SP.normalise_to_band(x, y, blo, bhi, bmode, minimum=is_T)))
            except Exception as e:
                weak.append(f"{nm} ({e})")
                r0 = np.ptp(y[m_]) or 1.0          # fallback: own range = 1 so the trace stays readable
                nt.append((nm, x, (y - np.nanmax(y[m_])) / r0 + 1 if is_T else (y - np.nanmin(y[m_])) / r0))
        if weak:
            st.caption("Not band-normalised (scaled to their own range instead): " + "; ".join(weak))
        traces = nt
    peaks_ = {nm: SP.pick_peaks(x, y, prom, int(npk), (lo, hi), minima=is_T) for nm, x, y in traces if pl and nm in pl_which}
    try:
        zlo, zhi = sorted(float(v) for v in zr.split("-")[:2])
    except Exception:
        zlo, zhi = lo, hi

    def draw(ax, traces=traces):
        st.session_state["_even|ftiroff"] = 1.05
        P.draw_stacked(ax, traces, kind, offset=off, xlim=(lo, hi), highlight=highlight, invert=(kind != "raman"),
                       scale="common" if nb else "each")
        info = {e[0]: e for e in getattr(ax, "_mfs_stack", [])}
        for nm, x, y in traces:
            if nm in peaks_ and nm in info:
                _, y0, r_, base, _k = info[nm]
                yl = (np.asarray(y) - y0) / r_ + base
                if is_T:     # transmittance: write the labels just above the trace's top, at the dip positions
                    m_ = (np.asarray(x) >= lo) & (np.asarray(x) <= hi); yl = np.full_like(yl, np.nanmax(yl[m_]) + 0.01)
                SP.label_peaks(ax, x, yl, peaks_[nm])
        if ins and zhi > zlo:
            SP.add_inset(ax, lambda ia: P.draw_stacked(ia, traces, kind, offset=off, xlim=(zlo, zhi), invert=(kind != "raman"),
                                                      scale="common" if nb else "each"), (zhi, zlo) if kind != "raman" else (zlo, zhi), loc=zloc)
    fig, ax = newfig(1.1); draw(ax)
    show(fig, "FTIR" if kind != "raman" else "Raman", pd.concat(rows) if rows else None, draw, "FTIR" if kind != "raman" else "Raman",
         units={"x": "wavenumber"} if kind != "raman" else None)
    if peaks_:
        st.dataframe(pd.DataFrame([{"sample": k, "peaks (cm⁻¹)": ", ".join(f"{v:.0f}" for v in vv)} for k, vv in peaks_.items()]),
                     hide_index=True, width="stretch")
    if len(traces) >= 2:
        with st.expander("➖ Difference spectra (sample − f × reference)", expanded=False):
            names = [t[0] for t in traces]
            c1, c2, c3 = st.columns(3)
            refn = c1.selectbox("reference", names, index=0, key="ft_dref")
            sams = c2.multiselect("samples", [n for n in names if n != refn], default=[n for n in names if n != refn][:2], key="ft_dsam")
            fmode = c3.radio("factor", ["1", "least squares in a window", "manual"], key="ft_dfm")
            fw = c3.text_input("window (cm⁻¹)", "1350-1450", key="ft_dfw") if fmode.startswith("least") else None
            fman = c3.number_input("f", 0.0, 10.0, 1.0, 0.05, key="ft_dfman") if fmode == "manual" else 1.0
            rt = traces[names.index(refn)]
            diffs, drows = [], []
            for sn in sams:
                stt = traces[names.index(sn)]
                f_ = fman
                if fw:
                    try:
                        a_, b_ = sorted(float(v) for v in fw.split("-")[:2]); f_ = SP.auto_factor((rt[1], rt[2]), (stt[1], stt[2]), a_, b_)
                    except Exception as e:
                        st.caption(f"factor: {e}")
                g_, d_ = SP.difference((rt[1], rt[2]), (stt[1], stt[2]), f_)
                diffs.append((f"{sn} − {f_:.2f}×{refn}", g_, d_)); drows.append(pd.DataFrame({"difference": diffs[-1][0], "x": g_, "y": d_}))

            def draw_d(ax, diffs=diffs):
                for k, (nm, g_, d_) in enumerate(diffs):
                    sc = np.nanmax(np.abs(d_[(g_ >= lo) & (g_ <= hi)])) or 1
                    ax.plot(g_, d_ / sc + 2.2 * k, lw=0.9, color=P.colours_for([n for n, _, _ in diffs])[k])
                    ax.axhline(2.2 * k, color="0.7", lw=0.5, ls=":")
                    ax.text(0.02, (2.2 * k + 0.6), nm, transform=ax.get_yaxis_transform(), fontsize=6)
                ax.set_xlim((hi, lo) if kind != "raman" else (lo, hi)); ax.set_yticks([])
                ax.set_xlabel(P.LABELS[kind][0]); ax.set_ylabel("Difference (norm.)")
            if diffs:
                fig, ax = newfig(0.8); draw_d(ax)
                show(fig, ("FTIR" if kind != "raman" else "Raman") + "_difference", pd.concat(drows), draw_d, "Difference spectra")


def page_nmr():
    st.header("¹H NMR of digested MOFs (composition) – Bruker Avance 400")
    st.caption("TopSpin experiment folder as .zip (pdata/1/1r + procs), TopSpin *totxt* export, or ppm/intensity text.")
    files = files_input("Spectra", ["zip", "txt", "csv", "dat"], "NMR/*_totxt.txt", "nmr")
    if not files:
        return
    samples = sample_table(files, "nmr_tab")
    c1, c2, c3 = st.columns(3)
    lo, hi = c1.slider("δ range (ppm)", -1.0, 14.0, (1.5, 9.0), 0.1)
    cut = c2.text_input("Cut regions (e.g. solvent HDO 4.55-5.05)", "4.55-5.05")
    off = offset_controls(c3, "nmr_off", 1.1, 0.2, 2.0, 0.05, "Offset")
    st.session_state["_even|nmr_off"] = 1.05
    ann = c3.text_input("Annotations (ppm:label; …)", "8.44:HCOO$^-$; 6.52:fumarate; 3.62:SS CH; 2.68:SS CH$_2$")
    cuts = []
    for p in cut.split(";"):
        try:
            a, b = [float(v) for v in p.split("-")]; cuts.append((a, b))
        except Exception:
            pass
    traces, specs = [], {}
    for nm, f in samples:
        try:
            d = R.read_nmr(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        specs[nm] = d
        m = np.ones(len(d), bool)
        for a, b in cuts:
            m &= ~((d.x > a) & (d.x < b))
        y = np.where(m, d.y, np.nan)
        traces.append((nm, d.x.values, y))
    annotations = []
    for p in ann.split(";"):
        if ":" in p:
            try:
                a, b = p.split(":", 1); annotations.append((float(a), b.strip()))
            except Exception:
                pass
    def draw(ax, traces=traces):
        P.draw_stacked(ax, traces, "nmr", offset=off, xlim=(hi, lo), annotations=annotations)
    fig, ax = newfig(0.9); draw(ax)
    show(fig, "NMR_digest", None, draw, "1H NMR")
    st.subheader("Integration → linker : modulator ratios (feeds the defect calculator)")
    st.caption("Integrals are baseline-corrected (linear between region ends). Internal standard TMSP-d4 = 9H at 0.00 ppm.")
    regions = pd.DataFrame({"species": ["TMSP-d4 (standard)", "Fumarate", "Formate", "2-Sulfosuccinate (dicarboxylate, SO3H kept)", "Acetate"],
                            "from_ppm": [-0.15, 6.40, 8.35, 2.40, 1.80], "to_ppm": [0.15, 6.65, 8.55, 2.90, 2.00], "nH": [9, 2, 1, 2, 3]})
    regions = st.data_editor(regions, num_rows="dynamic", key="nmr_regions")
    out = []
    for nm, d in specs.items():
        row = {"sample": nm}
        for r in regions.itertuples():
            row[f"I {r.species}"] = A.integrate(d.x.values, d.y.values, r.from_ppm, r.to_ppm)
        std = row.get("I TMSP-d4 (standard)") or np.nan
        for r in regions.itertuples():
            if r.species.startswith("TMSP"):
                continue
            row[f"I {r.species}"] = row[f"I {r.species}"] / std * 9 if std == std and std else row[f"I {r.species}"]
            row[f"n {r.species} (rel. to 9H std)"] = row[f"I {r.species}"] / r.nH
        main = "n Fumarate (rel. to 9H std)"
        if main in row and row[main]:
            for r in regions.itertuples():
                k = f"n {r.species} (rel. to 9H std)"
                if k in row and k != main:
                    row[f"ratio {r.species}:linker"] = row[k] / row[main]
            parts = [("linker", 1.0)] + [(r.species.split(" (")[0], row.get(f"ratio {r.species}:linker")) for r in regions.itertuples()
                                         if not r.species.startswith("TMSP") and f"ratio {r.species}:linker" in row]
            row["linker : modulators"] = " : ".join(f"{v:.2f}" for _, v in parts if v is not None) + "  (" + " : ".join(n_ for n_, v in parts if v is not None) + ")"
        out.append(row)
    if out:
        df = pd.DataFrame(out)
        st.dataframe(df, width="stretch")
        st.caption("Ratios per linker (main linker = 1). The integrals go to the sample store and are pre-filled in the defect "
                   "calculator ('Pull data for sample') for methods (b), (c), (d) and (e).")
        st.session_state["nmr_integrals"] = df
        for row in out:
            store_put(row["sample"], "nmr", {k: v for k, v in row.items() if k != "sample"}, "1H NMR digest page")
        dl(st, "CSV integrals", df.to_csv(index=False).encode("utf-8-sig"), "NMR_integrals.csv", "nmr_int_csv")


MARKERS = {"circle ●": "o", "square ■": "s", "triangle ▲": "^", "down-triangle ▼": "v", "diamond ◆": "D", "pentagon ⬟": "p",
           "hexagon ⬢": "h", "star ★": "*", "plus +": "P", "cross ×": "X", "none": "none"}
LINESTYLES = {"solid —": "-", "dashed – –": "--", "dotted ···": ":", "dash-dot –·": "-.", "none": "none"}


def page_sorption():
    st.header("N₂ physisorption at 77 K – Micromeritics 3Flex / TriStar / ASAP")
    st.caption("MicroActive 'Export' (.xls via xlrd / .xlsx) of the report, or CSV/TXT with p/p₀ and quantity adsorbed. "
               "v2 splits the export into its individual reports and uses only the Isotherm Tabular Report "
               "(the BET, Langmuir, t-plot and plot-overlay tables are no longer merged into the isotherm). TriStar/ASAP "
               "'Full Report Set' exports with all reports side by side on one sheet are split into their report columns.")
    files = files_input("Isotherm files", ["xls", "xlsx", "csv", "txt"], "Sorption/*_3Flex_*.xls*", "sorp")
    if not files:
        return
    samples = sample_table(files, "sorp_tab")
    data = {}
    for nm, f in samples:
        try:
            data[nm] = R.read_sorption(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}")
    if not data:
        return
    with st.expander("🔎 What was read from each file (tables found, which one is used, units, sample name)"):
        for nm, d in data.items():
            m = d["meta"]
            st.markdown(f"**{nm}** – sample in file: `{m.get('sample', '–')}` · mass {m.get('mass_g', '–')} g · "
                        f"{int((d['iso'].branch == 'ads').sum())} adsorption + {int((d['iso'].branch == 'des').sum())} desorption points · "
                        f"quantity unit in file: {m.get('q_unit_in_file')} (plotted as cm³ g⁻¹ STP) · isotherm from: {m.get('isotherm_table', '–')}")
            if d.get("tables"):
                st.dataframe(pd.DataFrame(d["tables"]), width="stretch", hide_index=True)
    # ---------------------------------------------------------------- per-sample style
    with st.expander("🎨 Series style per sample (adsorption / desorption colour, marker, line)", expanded=False):
        styles = {}
        defaults = S.colours_for(list(data))
        for nm, dc in zip(data, defaults):
            h = hashlib.md5(nm.encode()).hexdigest()[:10]
            c = st.columns([2, 1, 1, 2, 2, 1, 2, 1])
            c[0].markdown(f"**{re.sub(r'[$_{{}}]|mathrm', '', nm)}**")
            ca = c[1].color_picker("ads", dc, key=f"iso_ca_{h}")
            cd = c[2].color_picker("des", ca, key=f"iso_cd_{h}")
            mk = c[3].selectbox("marker", list(MARKERS), key=f"iso_mk_{h}")
            mkd = c[4].selectbox("des marker", ["same"] + list(MARKERS), key=f"iso_mkd_{h}")
            ms = c[5].number_input("size", 0.5, 12.0, 3.0, 0.5, key=f"iso_ms_{h}")
            ls = c[6].selectbox("line", list(LINESTYLES), key=f"iso_ls_{h}")
            lw = c[7].number_input("width", 0.0, 4.0, 0.8, 0.1, key=f"iso_lw_{h}")
            styles[nm] = {"ads_color": ca, "des_color": cd, "marker": MARKERS[mk], "des_marker": None if mkd == "same" else MARKERS[mkd],
                          "ms": ms, "ls": LINESTYLES[ls], "lw": lw}
        st.caption("Adsorption = filled symbols, desorption = open symbols (MOF convention). Lines join the points in measurement order.")
    c1, c2, c3, c4 = st.columns(4)
    logx = c1.checkbox("log p/p₀ axis (micropore filling)", False)
    both = c1.checkbox("also show linear + log side by side", False)
    show_des = c2.checkbox("show desorption", True)
    psd_kind = c3.selectbox("PSD representation", ["dvdw", "dvdlogw", "cum"], format_func={"dvdw": "dV/dw", "dvdlogw": "dV/dlog(w)", "cum": "cumulative"}.get)
    wmax = c4.number_input("PSD max width (nm)", 0.5, 50.0, 2.0)
    # ---------------------------------------------------------------- BET
    rows = []
    st.subheader("BET (Rouquerol consistency criteria, ISO 9277)")
    for nm, d in data.items():
        a = d["iso"][d["iso"].branch == "ads"]
        bet = A.bet_rouquerol(a.p_rel.values, a.q.values)
        tp = A.t_plot(a.p_rel.values, a.q.values)
        vt = A.total_pore_volume(a.p_rel.values, a.q.values, 0.95)
        m = d["meta"]
        rr = m.get("bet_reported_range")
        fixed = A.bet_fixed_range(a.p_rel.values, a.q.values, *(rr or (0.05, 0.30)))
        for w in A.isotherm_quality(a.p_rel.values, a.q.values):
            st.warning(f"{nm}: {w}")
        if m.get("bet_C_reported"):
            try:
                if float(m["bet_C_reported"]) <= 0:
                    st.info(f"{nm}: the file's BET fit has C = {float(m['bet_C_reported']):.1f} (negative). A negative C has no physical "
                            "meaning (ISO 9277 requires C > 0), so the reported area is not a valid BET area – the app's own value is used.")
            except ValueError:
                pass
        if bet is None:
            st.error(f"{nm}: fewer than 3 adsorption points below p/p₀ 0.35 – no BET area can be computed from this file.")
        if bet:
            rows.append({"sample": nm, "S_BET (m2/g) – this app": round(bet["area"], 0), "C": round(bet["C"], 0),
                         "p/p0 range": f"{bet['p_lo']:.2e}–{bet['p_hi']:.2e}", "points": bet["n"], "R2": round(bet["R2"], 5),
                         "all 4 Rouquerol criteria": bet["all_criteria"],
                         "reported in file (m2/g)": float(m["bet_reported"]) if m.get("bet_reported") else None,
                         "C in file": float(m["bet_C_reported"]) if m.get("bet_C_reported") else None,
                         "file's p/p0 range": f"{rr[0]:.3g}–{rr[1]:.3g}" if rr else None,
                         "same data, file's range (m2/g)": round(fixed["area"], 0) if fixed else None,
                         "V_total @0.95 (cm3/g)": round(vt, 3), "V_micro t-plot (cm3/g)": round(tp["V_micro"], 3) if tp else None,
                         "S_ext t-plot (m2/g)": round(tp["S_ext"], 0) if tp else None})
            d["bet"] = bet
            store_put(nm, "n2", {"S_BET": bet["area"], "C": bet["C"], "p_lo": bet["p_lo"], "p_hi": bet["p_hi"], "points": bet["n"],
                                 "all_criteria": bet["all_criteria"], "V_total": vt, "V_micro": tp["V_micro"] if tp else None,
                                 "S_ext": tp["S_ext"] if tp else None}, "N2 sorption page")
            if bet.get("few_points"):
                st.warning(f"{nm}: only {bet['n']} adsorption points lie in the Rouquerol range (p/p₀ {bet['p_lo']:.3g}–{bet['p_hi']:.3g}); "
                           "the BET area is computed from these 3 points (the app normally requires ≥ 4). Treat it as an estimate and "
                           "measure more points below p/p₀ 0.05.")
            if not bet["all_criteria"]:
                miss = [n for n, k in (("1 (n(1-p) increasing)", "c1"), ("3 (n_m inside the fitted range)", "c3"), ("4 (1/(√C+1) ≈ p at n_m, ±20 %)", "c4")) if not bet[k]]
                st.warning(f"{nm}: no p/p0 range satisfies every Rouquerol criterion – failed: {', '.join(miss)}. "
                           "Common for ultramicroporous MOFs: report a BET-equivalent area with the p/p0 range and show the Rouquerol plot in the SI.")
    if rows:
        tab = pd.DataFrame(rows); st.dataframe(tab, width="stretch", hide_index=True)
        st.caption("The value **reported in the file** is shown for reference only. Instrument software fits a fixed window "
                   "(typically p/p₀ 0.05–0.30); for microporous MOFs that window lies after micropore filling and usually gives a lower, "
                   "Rouquerol-inconsistent area. This app's value uses the Rouquerol-consistent range (ISO 9277:2022). The column "
                   "'same data, file's range' shows that the difference is the pressure window, not the data.")
        dl(st, "CSV BET/t-plot table", tab.to_csv(index=False).encode("utf-8-sig"), "BET_table.csv", "bet_tab_csv")
        crit = []
        for nm, d in data.items():
            b = d.get("bet")
            if not b:
                continue
            ok = lambda v: "✅" if v else "❌"
            crit.append({"sample": nm, "p/p₀ range (auto)": f"{b['p_lo']:.2e}–{b['p_hi']:.2e}", "points": b["n"],
                         "1 · n(1−p/p₀) increases": ok(b["c1"]), "2 · C > 0": ok(b["C"] > 0) + f" C = {b['C']:.0f}",
                         "3 · n_m inside the range": ok(b["c3"]) + (f" p(n_m) = {b['p_nm']:.2e}" if b["c3"] else ""),
                         "4 · 1/(√C+1) ≈ p(n_m) ±20 %": ok(b["c4"]) + f" {b['p_th']:.2e}",
                         "S_BET (m²/g)": round(b["area"], 0)})
        if crit:
            st.markdown("**Rouquerol criteria panel** – the BET range is chosen automatically as the range with the highest R² among "
                        "all ranges that satisfy every criterion (ISO 9277:2022, Rouquerol et al. 2007)")
            st.dataframe(pd.DataFrame(crit), width="stretch", hide_index=True)
    sbet_leg = c2.checkbox("S_BET in the legend", False, key="sorp_sbet_leg")
    lab_ = {nm: (f"{nm} ({d['bet']['area']:.0f} m$^2$ g$^{{-1}}$)" if sbet_leg and d.get("bet") else nm) for nm, d in data.items()}
    styles = {lab_[k]: v for k, v in styles.items()}
    isos = [(lab_[nm], d["iso"] if show_des else d["iso"][d["iso"].branch == "ads"]) for nm, d in data.items()]
    psds = [(nm, d["psd"]) for nm, d in data.items() if d["psd"] is not None]

    def draw_iso(ax, isos=isos, logx=logx, styles=styles):
        P.draw_isotherms(ax, isos, logx=logx, styles=styles)

    def draw_psd(ax, psds=psds):
        P.draw_psd(ax, psds, kind=psd_kind, xlim=(0.3, wmax), legend=len(psds) > 1)
    csv = pd.concat([df.assign(sample=n) for n, df in isos])
    if both:
        fig, axs = C.compose([lambda ax: draw_iso(ax, logx=False), lambda ax: P.draw_isotherms(ax, isos, logx=True, styles=styles, legend=False)],
                             ncols=2, width=max(WIDTH, 6.0), panel_aspect=0.85, wspace=0.4)
        show(fig, "N2_isotherms_lin_log", csv, None, units={"y": "uptake"})
    else:
        fig, ax = newfig(); draw_iso(ax)
        show(fig, "N2_isotherms", csv, draw_iso, "N2 isotherm", units={"y": "uptake"})
    no_psd = [nm for nm, d in data.items() if d["psd"] is None]
    if no_psd:
        st.info("No pore-size-distribution table in: " + ", ".join(no_psd) + ". These files contain at most BJH/D-H summary values; "
                "for a PSD panel export the DFT, HK or BJH *tabular* report as well.")
    if psds:
        st.subheader("Pore-size distribution (separate panel)")
        fig, ax = newfig(); draw_psd(ax)
        show(fig, "DFT_PSD", pd.concat([df.assign(sample=n) for n, df in psds]), draw_psd, "PSD", units={"x": "pore_nm"})
    sel = st.selectbox("BET plots for", [n for n in data if "bet" in data[n]])
    if sel:
        bet = data[sel]["bet"]
        fig, axs = C.compose([lambda ax: P.draw_rouquerol(ax, bet), lambda ax: P.draw_bet(ax, bet, sel)], ncols=2, width=max(WIDTH, 5.0), panel_aspect=0.8, wspace=0.5)
        show(fig, f"BET_{re.sub(r'[^A-Za-z0-9]+', '_', sel)}")
        st.subheader("t-plot (Harkins–Jura) – micropore volume")
        trng = st.slider("t range for the linear fit (Å)", 3.0, 10.0, (3.5, 5.0), 0.1, key="sorp_trng")
        a_ = data[sel]["iso"][data[sel]["iso"].branch == "ads"]
        tp = A.t_plot(a_.p_rel.values, a_.q.values, trng)
        if tp:
            st.write(f"**{sel}**: V_micro = **{tp['V_micro']:.3f} cm³ g⁻¹** (intercept {tp['intercept']:.1f} cm³ STP g⁻¹), "
                     f"S_ext = {tp['S_ext']:.0f} m² g⁻¹ (t {trng[0]}–{trng[1]} Å). A negative intercept means no micropores in this t window.")

            def draw_t(ax, tp=tp, trng=trng):
                ax.plot(tp["t"], tp["q"], "o", ms=3, color="k", mfc="white")
                tt = np.linspace(0, max(tp["t"].max(), trng[1]), 50)
                ax.plot(tt, tp["slope"] * tt + tp["intercept"], "-", color="#D55E00", lw=0.9)
                ax.axvspan(*trng, color="#D55E00", alpha=0.08, lw=0)
                ax.set_xlabel("t (Å)"); ax.set_ylabel("Quantity adsorbed (cm$^3$ g$^{-1}$ STP)"); ax.set_xlim(0, None); ax.set_ylim(0, None)
            fig, ax = newfig(); draw_t(ax)
            show(fig, f"tplot_{safe_name(sel)}", pd.DataFrame({"t (A)": tp["t"], "q (cm3/g STP)": tp["q"]}), draw_t, "t-plot")
        else:
            st.info("Fewer than 3 points in this t range.")
    st.subheader("Pore-size distribution computed from the isotherm")
    from mfs import psd as PSDm
    METH = list(PSDm.KERNELS) + ["Horváth–Kawazoe slit (carbon parameters)", "Saito–Foley cylinder (oxide parameters)",
                                 "BJH, desorption branch (mesopores)", "BJH, adsorption branch (mesopores)"]
    c1, c2, c3 = st.columns([3, 1, 1])
    meth = c1.selectbox("Method", METH, key="sorp_psdm",
                        help="The DFT kernel is the openly licensed pyGAPS N2/77 K carbon-slit kernel (MIT) – no free NLDFT/QSDFT "
                             "kernel exists for oxides or cylindrical pores, so for those Saito–Foley (micropores) and BJH (mesopores) are "
                             "offered. Commercial NLDFT/QSDFT PSDs: import the instrument's DFT tabular report (panel above).")
    pkind = c2.selectbox("Show", ["dV/dw", "cumulative"], key="sorp_psdk")
    wmax2 = c3.number_input("max width (nm)", 0.5, 100.0, 5.0 if "BJH" not in meth else 50.0, key="sorp_psdw_" + ("b" if "BJH" in meth else "m"))
    comp = []
    for nm, d in data.items():
        iso = d["iso"]; a_ = iso[iso.branch == "ads"]; de_ = iso[iso.branch == "des"]
        try:
            if meth in PSDm.KERNELS:
                r_ = PSDm.dft_fit(a_.p_rel.values, a_.q.values, meth)
            elif meth.startswith("Horv"):
                r_ = PSDm.hk_psd(a_.p_rel.values, a_.q.values, "slit", "Carbon (HK)")
            elif meth.startswith("Saito"):
                r_ = PSDm.hk_psd(a_.p_rel.values, a_.q.values, "cylinder", "AlSi oxide ion (SF)")
            else:
                br = de_ if ("desorption" in meth and len(de_) > 5) else a_
                r_ = PSDm.bjh_psd(br.p_rel.values, br.q.values, "desorption" if br is de_ else "adsorption")
            comp.append((lab_[nm], r_))
        except Exception as e:
            st.warning(f"{nm}: PSD failed ({e})")
    if comp:
        st.caption("Method: **" + comp[0][1]["method"] + "**" + (f" · fit residual " + ", ".join(f"{n}: {100*r['residual_rel']:.1f} %" for n, r in comp)
                   if "residual_rel" in comp[0][1] else "") + ". Shown on the figure legend title.")
        cols_ = S.colours_for([n for n, _ in comp])

        def draw_cpsd(ax, comp=comp, pkind=pkind, wmax2=wmax2, cols_=cols_):
            for (n, r_), col in zip(comp, cols_):
                yv = r_["dVdw"] if pkind == "dV/dw" else r_["V_cum"]
                ax.plot(r_["w"], yv, "-", color=col, lw=1.0, label=n)
            ax.set_xlim(0, wmax2); ax.set_ylim(0, None)
            ax.set_xlabel("Pore width (nm)")
            ax.set_ylabel("d$V$/d$w$ (cm$^3$ g$^{-1}$ nm$^{-1}$)" if pkind == "dV/dw" else "Cumulative pore volume (cm$^3$ g$^{-1}$)")
            ax.legend(frameon=False, title=comp[0][1]["method"].split(" (")[0], title_fontsize=FONT - 2, fontsize=FONT - 2)
        fig, ax = newfig(); draw_cpsd(ax)
        show(fig, "PSD_computed", pd.concat([pd.DataFrame({"sample": n, "w (nm)": r_["w"], "dV/dw": r_["dVdw"], "V_cum": r_["V_cum"],
                                                            "method": r_["method"]}) for n, r_ in comp]), draw_cpsd, "PSD (computed)", units={"x": "pore_nm"})
        if meth in PSDm.KERNELS:
            with st.expander("Kernel fit of the isotherm (check the fit quality)"):
                fig, ax = newfig()
                for (n, r_), col in zip(comp, cols_):
                    ax.semilogx(r_["p"], r_["q"], "o", ms=2.5, mfc="white", color=col, label=n)
                    ax.semilogx(r_["p"], r_["fit"], "-", color=col, lw=0.9)
                ax.set_xlabel("$p/p_0$"); ax.set_ylabel("Quantity adsorbed (cm$^3$ g$^{-1}$ STP)"); ax.legend(frameon=False, fontsize=FONT - 2)
                show(fig, "PSD_kernel_fit")


def page_tga():
    st.header("Thermogravimetric analysis")
    st.caption("Netzsch Proteus ASCII, TA Universal Analysis/TRIOS text, Mettler STARe export, or any CSV with temperature + mass/weight columns.")
    files = files_input("TGA files", ["txt", "csv", "dat", "xlsx"], "TGA/*.txt", "tga")
    if not files:
        return
    samples = sample_table(files, "tga_tab")
    c1, c2, c3 = st.columns(3)
    dtg = c1.selectbox("DTG", ["panel", "twin", "none"])
    norm_end = c2.checkbox("Normalise to residue = 100 % (Shearer plot)", False)
    win = c3.slider("DTG smoothing window", 5, 101, 21, 2)
    res = []
    for nm, f in samples:
        try:
            t = R.read_tga(f, f.name); a = A.tga_analysis(t["T"], t["w"], win); res.append((nm, a))
        except Exception as e:
            st.error(f"{f.name}: {e}")
    tab = pd.DataFrame([{"sample": nm, "W(200 °C) %": round(a["W200"], 2) if a["W200"] else None, "plateau T (°C)": round(a["T_plateau"], 0) if a["T_plateau"] else None,
                         "W(plateau) %": round(a["W_plateau"], 2) if a["W_plateau"] else None, "W(end) %": round(a["W_end"], 2),
                         "plateau/end": round(a["W_plateau"] / a["W_end"], 4) if a["W_plateau"] else None,
                         "DTG peaks (°C)": ", ".join(f"{p:.0f}" for p in a["peaks"])} for nm, a in res])
    st.dataframe(tab, width="stretch")
    st.session_state["tga_table"] = tab
    for nm, a in res:
        store_put(nm, "tga", {"W_initial": a["W_initial"], "W200": a["W200"], "W_plateau": a["W_plateau"], "T_plateau": a["T_plateau"],
                              "W_end": a["W_end"], "DTG_peaks": list(a["peaks"])}, "TGA page")
    from mfs import tga as TG
    c1, c2, c3 = st.columns(3)
    lab_on = c1.checkbox("Automatic step labels (mass loss per DTG step)", True, key="tga_steps")
    min_loss = c2.number_input("min. step (%-points)", 0.1, 20.0, 1.0, 0.1, key="tga_minloss")
    lab_which = c3.multiselect("label which samples", [nm for nm, _ in res], default=[nm for nm, _ in res][:1], key="tga_labw")
    stp = {nm: TG.steps(a, min_loss) for nm, a in res}
    with st.expander("📉 Steps per sample", expanded=False):
        st.dataframe(pd.DataFrame([dict(sample=nm, step=i + 1, **{"T start (°C)": s_["T_start"], "DTG peak (°C)": s_["T_peak"],
                                                                   "T end (°C)": s_["T_end"], "loss (%-points)": s_["loss"],
                                                                   "loss (% of residue)": s_["loss"] * 100 / a["W_end"]})
                                   for nm, a in res for i, s_ in enumerate(stp[nm])]).style.format(precision=1), hide_index=True, width="stretch")
    with st.expander("⚖️ Plateaus normalised to the oxide residue (= 100 %) → linkers per node (any node / metal)", expanded=False):
        gfw, _gm = guess_framework(res[0][0] if res else "")
        c1, c2 = st.columns([1, 2])
        fw_ = c1.selectbox("Framework", list(FRAMEWORKS), index=list(FRAMEWORKS).index(gfw or "UiO-66"), key="tga_fw")
        nsel_ = node_picker(c2, fw_, "tga", metal_default=_gm)
        met = nsel_["metal"]
        st.caption(f"Ideal plateau: node with all ligand charge compensated by O²⁻ + L linkers → TGA residue {_res_lab(met, nsel_['def'])} "
                   "(M₆O₆(L)₆ → 6 MO₂ for Zr/Hf: Valenzano 2011; Shearer 2014). The defect calculator uses the same equations.")
        pr = []
        for nm, a in res:
            if not a["W_plateau"]:
                continue
            lk = FRAMEWORKS[fw_]["linker"]
            inp = D.DefectInput(sample=nm, framework=fw_, metal=met, species=[D.Species(lk, SPECIES[lk][0], SPECIES[lk][1], "linker", SPECIES[lk][3], None, False, True)],
                                tga_W_plateau=a["W_plateau"], tga_W_end=a["W_end"], node=nsel_["node"], node_custom=nsel_["node_custom"])
            try:
                r_ = D.method_a(inp)
                pr.append({"sample": nm, "plateau T (°C)": a["T_plateau"], "W(plateau) (% of residue)": 100 * a["W_plateau"] / a["W_end"],
                           "ideal (%)": r_.extra.get("W_plat_ideal_%") if r_.extra else None, "linkers / node L": r_.L, "missing x": r_.x})
            except Exception as e:
                pr.append({"sample": nm, "note": str(e)[:80]})
        st.dataframe(pd.DataFrame(pr).style.format(precision=2), hide_index=True, width="stretch")

    def _labels(ax):
        if lab_on:
            for nm, a in res:
                if nm in lab_which:
                    TG.label_steps(ax, a, stp[nm], normalise_end=norm_end)
    if dtg == "panel":
        fig = P.fig_tga(res, WIDTH, "panel", normalise_end=norm_end); _labels(fig.axes[0])
    else:
        fig, ax = newfig(0.8)
        tw = ax.twinx() if dtg == "twin" else None
        if tw is not None:
            ax.tick_params(right=False, which="both"); tw.tick_params(direction="in")
        P.draw_tga(ax, res, tw, normalise_end=norm_end); _labels(ax)
    def draw(ax, res=res):
        P.draw_tga(ax, res, normalise_end=norm_end); _labels(ax)
    show(fig, "TGA", pd.concat([pd.DataFrame({"sample": nm, "T": a["T"], "w": a["w"], "DTG": a["dtg"]}) for nm, a in res]), draw, "TGA",
         units={0: {"x": "temperature", "y": "mass_pct"}, 1: {"x": "temperature"}})
    st.subheader("TGA–MS overlay")
    c1, c2 = st.columns(2)
    msf = c1.file_uploader("Ion-current table (CSV/TXT: temperature or time + m/z columns)", type=["csv", "txt", "dat"], key="tga_ms")
    ms_ex = c2.checkbox("Use the TGA–MS EXAMPLE (MOF-801, synthetic m/z 18, 28, 44)", value=msf is None, key="tga_ms_ex")
    src = msf if msf is not None else (_F(os.path.join(EX, "TGA", "MOF801_QMS_ion_currents_EXAMPLE.csv")) if ms_ex else None)
    if src is not None:
        try:
            ms = TG.read_ms(src)
            chans = [c_ for c_ in ms.columns if c_.startswith("m/z")]
            c1, c2, c3, c4 = st.columns(4)
            tsel = c1.selectbox("TGA trace", [nm for nm, _ in res], key="tga_ms_s")
            chs = c2.multiselect("m/z channels", chans, default=chans[:4], key="tga_ms_ch")
            logy = c3.checkbox("log ion current", False, key="tga_ms_log")
            if "T" not in ms:
                T0 = c4.number_input("start T (°C)", -100.0, 1000.0, 30.0, key="tga_ms_T0"); beta = c4.number_input("rate (K/min)", 0.1, 100.0, 10.0, key="tga_ms_b")
                ms["T"] = T0 + beta * ms["t"]
            a = dict(res)[tsel]

            def draw_ms(ax, a=a, ms=ms, chs=chs):
                ax.plot(a["T"], a["w"] / a["W_end"] * 100 if norm_end else a["w"], color="k", lw=1.0)
                ax.set_ylabel("Weight (%)"); ax.set_xlabel("Temperature (°C)")
                tw = ax.twinx(); TG.draw_ms(tw, ms, chs, log=logy, norm=not logy)
                ax.set_xlim(max(a["T"].min(), ms["T"].min()), min(a["T"].max(), ms["T"].max()))
            fig, ax = newfig(0.8); draw_ms(ax)
            show(fig, "TGA_MS", ms, draw_ms, "TGA–MS")
        except Exception as e:
            st.warning(f"MS file: {e}")


# -------------------------------------------------------------- DSC
def page_dsc():
    from mfs import dsc as DSC
    DSC = _CachedReaders(DSC)
    st.header("Differential scanning calorimetry")
    st.caption("Mettler-Toledo STARe ‘Curve Values’ export, TA Instruments TRIOS / Universal Analysis text, Netzsch Proteus ASCII, "
               "PerkinElmer Pyris text, or any CSV/XLSX with temperature + heat-flow columns (decimal commas ok). Heating/cooling "
               "segments are split automatically. Peaks: onset / peak / endset and ΔH by integration over time; Tg by ISO 11357-2.")
    files = files_input("DSC files", ["txt", "csv", "dat", "xlsx", "xls"], "DSC/MOF-PEG_TA-TRIOS_EXAMPLE.txt", "dsc")
    if st.checkbox("➕ Use the full-cycle example instead: REAL heating (SG_DSC_1.1, STARe) + SIMULATED cooling (separate export)",
                   False, key="dsc_cycle_ex", help="Two files of one sample that are combined into Heat 1 → Cool 1. The cooling "
                   "curve and its crystallisation exotherm are synthetic (labelled EXAMPLE)."):
        files = [_F(p_) for p_ in sorted(glob.glob(os.path.join(EX, "DSC", "full_cycle", "SG_DSC_*ing_Mettler_*.txt")))]
    if st.checkbox("➕ Use the heat–cool–heat EXAMPLE (REAL heating 1; SIMULATED cooling and heating 2) – cycle-to-cycle comparison",
                   False, key="dsc_hch_ex"):
        files = [_F(os.path.join(EX, "DSC", "full_cycle", "SG_DSC_heat-cool-heat_EXAMPLE.csv"))]
    if not files:
        return
    samples = sample_table(files, "dsc_tab")
    data = []
    for nm, f in samples:
        try:
            data.append((nm, f, DSC.read_dsc(f, f.name)))
        except Exception as e:
            st.error(f"{f.name}: could not read ({type(e).__name__}: {e})")
    if not data:
        return
    def _base(nm):   # 'X_heating' / 'X cooling' / 'X_H1' -> 'X' : files of one sample share a cycle group
        return DSC.cycle_group_name(nm)
    with st.expander("⚖️ Files – instrument, sample mass, sign, cycle", expanded=True):
        info = pd.DataFrame([{"sample": nm, "instrument / format": d["meta"]["vendor"],
                              "mass (mg)": d["meta"].get("mass_mg"), "curves": len(d["curves"]),
                              "file time": d["meta"].get("timestamp").strftime("%Y-%m-%d %H:%M") if d["meta"].get("timestamp") else "",
                              "cycle group": _base(nm), "flip sign": False} for nm, f, d in data])
        info = st.data_editor(info, key="dsc_info", hide_index=True, disabled=["sample", "instrument / format", "curves", "file time"],
                              column_config={"flip sign": st.column_config.CheckboxColumn(help="Tick if exothermic events point the wrong way for this file"),
                                             "cycle group": st.column_config.TextColumn(help="Files with the same group are one sample; with 'Combine' they are merged into one heating/cooling cycle")})
        groups = list(dict.fromkeys(info["cycle group"].astype(str)))
        multi = any((info["cycle group"].astype(str) == g).sum() > 1 for g in groups)
        k1, k2 = st.columns(2)
        combine = k1.checkbox("Combine files of one sample into one cycle (e.g. separate heating and cooling exports)", multi, key="dsc_combine_" + hashlib.md5("|".join(info["sample"].astype(str)).encode()).hexdigest()[:6],
                              help="Merges the files of each cycle group in order, makes the time axis continuous and numbers the segments Heat 1, Cool 1, Heat 2 …")
        order = k2.radio("Order", ["by file time stamp", "as in the sample table"], horizontal=True, key="dsc_order",
                         help="Time stamps come from the file header; if a file has none, the table order is used.") if combine else None
    curves = []
    packs = []
    for (nm, f, d), (_, row) in zip(data, info.iterrows()):
        mass = float(row["mass (mg)"]) if pd.notna(row["mass (mg)"]) and row["mass (mg)"] else None
        cs = []
        for c in d["curves"]:
            c = dict(c)
            if row["flip sign"]:
                c["q"] = -c["q"]
            c["mass"] = mass; c["file"] = nm
            cs.append(c)
        packs.append((str(row["cycle group"]) if combine else nm, {"curves": cs, "meta": d["meta"]}))
    for g in dict.fromkeys(p_[0] for p_ in packs):
        ds = [p_[1] for p_ in packs if p_[0] == g]
        segs = DSC.combine_cycle(ds, "time" if order == "by file time stamp" else "user") if combine else DSC.label_segments(
            [c for d_ in ds for c in d_["curves"]])
        for c in segs:
            c["sample"] = g; c["uid"] = f"{g}|{c['seg']}"
            curves.append(c)
    one = len({c["sample"] for c in curves}) == 1
    for c in curves:
        c["label"] = c["seg"] if one else f"{c['sample']} · {c['seg']}"
    has_cycle = any(c["kind"] == "cooling" for c in curves) and any(c["kind"] == "heating" for c in curves)
    with st.expander("🔁 Segments (detected from the sign of dT/dt)", expanded=False):
        st.dataframe(pd.DataFrame([{"sample": c["sample"], "segment": c["seg"], "file": c.get("file", ""), "kind": c["kind"],
                                    "T start (°C)": round(float(c["T"][0]), 2), "T end (°C)": round(float(c["T"][-1]), 2),
                                    "rate (K/min)": round(DSC.heating_rate(c), 3) if DSC.heating_rate(c) is not None else None,
                                    "t start (min)": round(float(c["t"][0]), 2) if c.get("t") is not None else None,
                                    "points": len(c["T"])} for c in curves]), hide_index=True, width="stretch")
        if not has_cycle:
            st.caption("Only one direction found. For the full cycle upload the cooling export too (same cycle group) – or "
                       "export all segments from the instrument software (STARe: all curves of the method).")
    with st.expander("⚙️ Plot and evaluation settings", expanded=True):
        c1, c2, c3, c4 = st.columns(4)
        labels = [c["label"] for c in curves]
        default = [l for l, c in zip(labels, curves) if c["kind"] != "isothermal"]
        import hashlib as _hl   # key tied to the loaded curves: new files -> fresh selection (no stale/empty choice)
        sel = c1.multiselect("Segments / cycles", labels, default=default,
                             key="dsc_sel_" + _hl.md5("|".join(labels).encode()).hexdigest()[:8])
        exo_up = c1.radio("Convention", ["exo up", "exo down"], horizontal=True, key="dsc_exo") == "exo up"
        VIEWS = ["Curves", "Full cycle vs temperature", "Full cycle vs time", "Stacked: heating above cooling"]
        view = c1.selectbox("View", VIEWS, index=1 if has_cycle else 0, key="dsc_view_" + ("c" if has_cycle else "s"),
                            help="Full cycle vs T: heating and cooling as one closed loop with direction arrows. Vs time: with the "
                                 "temperature programme on a twin axis and shaded segments. Stacked: heating panel above the cooling panel.")
        yunit = c2.selectbox("Heat flow", ["W/g (= mW/mg)", "mW/mg", "mW"], key="dsc_unit",
                             help="W/g needs normalised data or a sample mass; mW multiplies normalised data by the mass.")
        xax = c2.radio("x axis", ["temperature", "time"], horizontal=True, key="dsc_x")
        off = offset_controls(c3, "dsc_off", 0.0, 0.0, 5.0, 0.05, "Vertical offset", help="Stacks the curves (y ticks are hidden when > 0)")
        bl = c3.selectbox("Baseline", ["linear", "sigmoidal"], key="dsc_bl")
        do_pk = c4.checkbox("Detect peaks (onset / peak / endset / ΔH)", True, key="dsc_pk")
        prom = c4.number_input("Min. peak prominence (0 = auto)", 0.0, 100.0, 0.0, format="%g", key="dsc_prom")
        ann = c4.checkbox("Annotate peaks on the figure", True, key="dsc_ann")
        showbl = c4.checkbox("Show integration baselines", True, key="dsc_showbl")
        tg_on = c1.checkbox("Glass transition Tg", False, key="dsc_tg")
        if tg_on:
            Tall = np.concatenate([c["T"] for c in curves])
            lo, hi = float(np.nanmin(Tall)), float(np.nanmax(Tall))
            heat = [c for c in curves if c["kind"] == "heating"] or curves
            try:
                h0 = heat[0]
                sug = DSC.auto_tg_window(h0, DSC.evaluate_all(dict(h0, id=0), "linear", h0["mass"]))
            except Exception:
                sug = None
            sug = sug or (lo + 0.15 * (hi - lo), lo + 0.45 * (hi - lo))
            tgw = c2.slider("Tg window (°C)", lo, hi, (round(sug[0], 1), round(sug[1], 1)), key="dsc_tgw",
                            help="Auto-suggested around the steepest step outside the peaks; tangents use the outer quarters.")
            tg_cool = c2.checkbox("Tg also on cooling segments", False, key="dsc_tgcool")
    styles = DSC.default_styles(list(dict.fromkeys(c["sample"] for c in curves)))
    if view != "Curves":
        with st.expander("🎨 Segment styles – colour and line style for heating / cooling, per sample", expanded=False):
            inv = {v: k for k, v in DSC.LINESTYLES.items()}
            sty = pd.DataFrame([{"sample": smp, "segment kind": kd, "colour": col, "line style": inv.get(ls, "solid")}
                                for (smp, kd), (col, ls) in styles.items() if any(c["sample"] == smp and c["kind"] == kd for c in curves)])
            sty = st.data_editor(sty, key="dsc_sty_" + str(abs(hash(tuple(sty["sample"]))) % 10**6), hide_index=True,
                                 disabled=["sample", "segment kind"],
                                 column_config={"colour": st.column_config.TextColumn(help="Hex colour, e.g. #C0392B"),
                                                "line style": st.column_config.SelectboxColumn(options=list(DSC.LINESTYLES))})
            for _, r in sty.iterrows():
                col = str(r["colour"]).strip()
                if re.fullmatch(r"#?[0-9a-fA-F]{6}", col):
                    col = col if col.startswith("#") else "#" + col
                    styles[(r["sample"], r["segment kind"])] = (col, DSC.LINESTYLES.get(r["line style"], "-"))
            q1, q2, q3 = st.columns(3)
            arrows = q1.checkbox("Direction arrows on the curves", True, key="dsc_arrows")
            connect = q2.checkbox("Close the loop (join consecutive segments)", True, key="dsc_connect")
            prog = q3.radio("Temperature trace (vs time)", ["Ts (sample)", "Tr (programme)"], horizontal=True, key="dsc_prog")
    else:
        arrows, connect, prog = True, True, "Ts (sample)"
    use = [c for c in curves if c["label"] in sel]
    if not use:
        st.info("Select at least one segment."); return
    items, peaks, tgs, rows, warn = [], {}, {}, [], []
    for c in use:
        spec = DSC.to_specific(c, c["mass"])
        if yunit == "mW":
            y = c["q"] if c["unit"] in ("mW", "?") else (spec * c["mass"] if (spec is not None and c["mass"]) else None)
            if y is None:
                warn.append(f"{c['label']}: no sample mass – shown in W/g"); y = spec
        else:
            y = spec if spec is not None else c["q"]
            if spec is None:
                warn.append(f"{c['label']}: heat flow in {c['unit']} and no sample mass – enter it in the table to normalise")
        cc = dict(c); cc["id"] = c["uid"]
        if do_pk:
            try:
                prs = DSC.evaluate_all(cc, bl, c["mass"], prominence=prom or None)
            except Exception as e:
                prs = []; warn.append(f"{c['label']}: peak evaluation failed ({e})")
            if yunit == "mW" and c["mass"] and spec is not None:
                for r in prs:
                    r = r  # baseline drawn in the plotted unit
                    r["baseline"] = r["baseline"] * (c["mass"] if c["unit"] == "W/g" else 1)
            peaks[cc["id"]] = prs
            for r in prs:
                rows.append({"curve": c["label"], "segment": c["seg"], "file": c.get("file", ""), "event": r["kind"], "onset (°C)": r["T_onset"], "peak (°C)": r["T_peak"],
                             "endset (°C)": r["T_endset"], "limits (°C)": f"{r['T_start']:.1f}–{r['T_end']:.1f}",
                             "ΔH (J/g)": r["dH"], "area (mJ)": r["area_mJ"], "baseline": bl})
        if tg_on and (c["kind"] == "heating" or tg_cool):
            g = DSC.glass_transition(cc, *tgw, mass_mg=c["mass"])
            if g:
                tgs[cc["id"]] = [g]
                rows.append({"curve": c["label"], "segment": c["seg"], "file": c.get("file", ""), "event": "Tg", "onset (°C)": g["Tg_onset"], "peak (°C)": g["Tg_mid"],
                             "endset (°C)": g["Tg_endset"], "limits (°C)": f"{g['window'][0]:.1f}–{g['window'][1]:.1f}",
                             "ΔH (J/g)": None, "area (mJ)": None, "baseline": f"Tg midpoint; inflection {g['Tg_inflection']:.1f} °C"
                             + (f"; Δcp {g['dcp_J_gK']:.3f} J g⁻¹ K⁻¹" if g["dcp_J_gK"] else "")})
        items.append((c["label"], cc, y))
    for w in warn:
        st.warning(w)
    if rows:
        tab = pd.DataFrame(rows)
        st.dataframe(tab.round(2), width="stretch", hide_index=True)
        st.session_state["dsc_table"] = tab
    for nm, f, d in data:
        rep = d["meta"].get("reported")
        if rep:
            mine = [r for r in rows if r.get("file") == nm and r["event"] != "Tg"]
            best = min(mine, key=lambda r: abs((r["peak (°C)"] or 0) - (rep.get("peak_C") or 0))) if mine else None
            st.info(f"**{nm}** – values reported by {d['meta']['vendor']}: "
                    + ", ".join(f"{k.replace('_C', ' (°C)').replace('_Jg', ' (J/g)').replace('_mJ', ' (mJ)').replace('_', ' ')} {v:g}" for k, v in rep.items())
                    + (f".  Reticula: onset {best['onset (°C)']:.2f} °C, peak {best['peak (°C)']:.2f} °C, ΔH {best['ΔH (J/g)']:.2f} J/g." if best and best["ΔH (J/g)"] is not None else ""))
    ylab = {"W/g (= mW/mg)": "Heat flow (W g$^{-1}$)", "mW/mg": "Heat flow (mW mg$^{-1}$)", "mW": "Heat flow (mW)"}[yunit]
    pk_ = peaks if do_pk else None
    if any(c.get("seg") for _, c, _ in items) and do_pk:
        summ = DSC.cycle_summary([c for _, c, _ in items], peaks)
        if len(summ) and summ[["Tm peak (°C)", "Tc peak (°C)"]].notna().any().any():
            st.markdown("**Cycle summary – melting (heating) vs crystallisation (cooling), supercooling ΔT = Tm − Tc**")
            st.dataframe(summ.round(2), width="stretch", hide_index=True)
            st.session_state["dsc_cycle_summary"] = summ
        cmp_ = DSC.cycle_compare([c for _, c, _ in items], peaks)
        kinds_ = pd.Series([r_["kind"] for r_ in cmp_]).value_counts() if cmp_ else pd.Series(dtype=int)
        if len(kinds_) and kinds_.max() >= 2:
            with st.expander("🔄 Cycle-to-cycle comparison (Heat 1 vs Heat 2 …, Cool 1 vs Cool 2 …)", expanded=True):
                st.dataframe(pd.DataFrame(cmp_).round(3), width="stretch", hide_index=True)
                ck = st.radio("overlay", [k_ for k_ in ("heating", "cooling") if kinds_.get(k_, 0) >= 2], horizontal=True, key="dsc_cmp_k")
                its_ = [it for it in items if it[1].get("kind") == ck]
                pks_ = peaks

                def draw_cmp(ax, its_=its_, pks_=pks_):
                    DSC.draw_dsc(ax, its_, "T", exo_up, 0.0, peaks=pks_ if do_pk else None, annotate=ann, show_baseline=False, ylabel=ylab,
                                 linestyles=["-", "--", ":", "-."] * 3)
                    ax.legend(frameon=False)
                figc, axc = newfig(0.75); draw_cmp(axc)
                show(figc, f"DSC_cycle_compare_{ck}", None, draw_cmp, f"DSC {ck} cycles")
    if view == "Full cycle vs time" and any(c.get("t") is None for _, c, _ in items):
        st.info("Some segments have no time column – showing the full cycle vs temperature instead."); view = "Full cycle vs temperature"
    progk = "Tr" if prog.startswith("Tr") else "Ts"
    if view == "Stacked: heating above cooling":
        S.apply_style(FONT)
        fig, axs = plt.subplots(2, 1, figsize=(WIDTH, WIDTH * 1.25), sharex=True)
        hi_ = [it for it in items if it[1]["kind"] == "heating"]; co_ = [it for it in items if it[1]["kind"] == "cooling"]
        DSC.draw_cycle_stacked(axs, hi_, co_, exo_up, styles, pk_, ann, showbl, ylab, arrows)
        for a_, its in zip(axs, (hi_, co_)):
            if len(its) > 1:
                a_.legend(frameon=False, fontsize=FONT - 1)
        fig.tight_layout()
        draw = None
    else:
        set_even_offset("dsc_off", [(it[1], it[2]) if not isinstance(it[1], dict) else (it[1].get("T"), it[2]) for it in items])
        fig, ax = newfig(0.75)
        if view == "Full cycle vs temperature":
            DSC.draw_cycle_T(ax, items, exo_up, styles, pk_, ann, showbl, ylab, arrows, connect)
        elif view == "Full cycle vs time":
            tw = DSC.draw_cycle_time(ax, items, exo_up, styles, pk_, ann, showbl, ylab, program=progk)
            tw.tick_params(which="both", direction="in")
        else:
            DSC.draw_dsc(ax, items, "T" if xax == "temperature" else "t", exo_up, off, peaks=pk_,
                         tgs=tgs, annotate=ann, show_baseline=showbl, ylabel=ylab)
        if len(items) > 1:
            ax.legend(frameon=False, loc="best")
        def draw(ax, items=items, peaks=peaks, tgs=tgs, view=view):
            if view == "Full cycle vs temperature":
                DSC.draw_cycle_T(ax, items, exo_up, styles, peaks if do_pk else None, ann, showbl, ylab, arrows, connect)
            elif view == "Full cycle vs time":
                DSC.draw_cycle_time(ax, items, exo_up, styles, peaks if do_pk else None, ann, showbl, ylab, program=progk)
            else:
                DSC.draw_dsc(ax, items, "T" if xax == "temperature" else "t", exo_up, off, peaks=peaks if do_pk else None, tgs=tgs,
                             annotate=ann, show_baseline=showbl, ylabel=ylab)
    exp = pd.concat([pd.DataFrame({"curve": lab, "sample": c.get("sample", ""), "segment": c.get("seg", ""), "kind": c.get("kind", ""),
                                   "T (°C)": c["T"], "Tr (°C)": c["Tr"] if c.get("Tr") is not None else np.nan,
                                   "t (min)": c["t"] if c.get("t") is not None else np.nan,
                                   "heat flow (exo +)": y}) for lab, c, y in items])
    stem = {"Curves": "DSC", "Full cycle vs temperature": "DSC_cycle_T", "Full cycle vs time": "DSC_cycle_time"}.get(view, "DSC_cycle_stacked")
    show(fig, stem, exp, draw, "DSC", units={0: {"x": "temperature" if view != "Full cycle vs time" else "time"}},
         series=[lab for lab, _, _ in items])
    sta = [c for c in use if c.get("w") is not None]
    if sta:
        st.subheader("TGA–DSC overlay (STA / SDT)")
        c = sta[0] if len(sta) == 1 else next(cc for cc in sta if cc["label"] == st.selectbox("Curve", [s_["label"] for s_ in sta], key="dsc_sta_sel"))
        y = DSC.to_specific(c, c["mass"]); y = y if y is not None else c["q"]
        fig2, ax2 = newfig(0.75)
        tw = DSC.draw_sta(ax2, c, y, exo_up)
        ax2.tick_params(right=False, which="both")
        show(fig2, "TGA_DSC_overlay", pd.DataFrame({"T (°C)": c["T"], "heat flow (exo +)": y, "mass (%)": c["w"]}), None, "TGA–DSC",
             units={0: {"x": "temperature"}})
    _dsc_cp()


def _dsc_cp():
    from mfs import dsc as DSC
    DSC = _CachedReaders(DSC)
    st.subheader("Specific heat capacity – sapphire method (ASTM E1269)")
    with st.expander("🌡️ cp from blank + sapphire + sample runs (same pans, rate and purge)", expanded=False):
        st.caption("cp(T) = cp,sapphire(T) · (q_sample − q_blank)/(q_sapphire − q_blank) · m_sapphire/m_sample. Sapphire values: "
                   "ASTM E1269 Table 1 (150–600 K), Ginnings & Furukawa above. Use the heating segment of each run (isothermal holds removed).")
        ex_ = st.checkbox("Use the Cp EXAMPLE (synthetic: blank, 25.30 mg sapphire, 9.80 mg MOF-801; truth cp = 0.62 + 0.0016 (T − 25))",
                          True, key="cp_ex")
        c1, c2, c3 = st.columns(3)
        fb = c1.file_uploader("blank (empty pans)", key="cp_b"); fs_ = c2.file_uploader("sapphire", key="cp_s")
        fx = c3.file_uploader("sample(s)", accept_multiple_files=True, key="cp_x")
        D_ = os.path.join(EX, "DSC", "Cp_sapphire")
        fb = fb or (_F(os.path.join(D_, "Cp_blank_EXAMPLE.csv")) if ex_ else None)
        fs_ = fs_ or (_F(os.path.join(D_, "Cp_sapphire_EXAMPLE.csv")) if ex_ else None)
        fx = fx or ([_F(os.path.join(D_, "Cp_sample_MOF-801_EXAMPLE.csv"))] if ex_ else [])
        if not (fb and fs_ and fx):
            st.caption("Load the three runs."); return
        try:
            db, ds = DSC.read_dsc(fb, fb.name), DSC.read_dsc(fs_, fs_.name)
            heat = lambda d: max([c_ for c_ in d["curves"] if (DSC.heating_rate(c_) or 1) > 0], key=lambda c_: len(c_["T"]))
            msap = c2.number_input("sapphire mass (mg)", 0.1, 200.0, float(ds["meta"].get("mass_mg") or 25.0), key="cp_ms_" + fs_.name)
            res_ = []
            for f_ in fx:
                dx = DSC.read_dsc(f_, f_.name)
                mx = c3.number_input(f"mass {f_.name} (mg)", 0.01, 200.0, float(dx["meta"].get("mass_mg") or 10.0), key="cp_mx_" + f_.name)
                r_ = DSC.cp_from_sapphire(heat(dx), heat(ds), heat(db), mx, msap)
                res_.append((guess_name(f_.name), r_))
            show_ref = st.checkbox("show the sapphire reference curve", False, key="cp_ref")

            def draw_cp(ax, res_=res_):
                DSC.draw_cp(ax, res_, show_ref=show_ref)
            fig, ax = newfig(0.75); draw_cp(ax)
            exp_ = pd.concat([pd.DataFrame({"sample": nm, "T (°C)": r_["T"], "cp (J/g/K)": r_["Cp"], "cp direct (J/g/K)": r_.get("Cp_direct"),
                                            "cp sapphire (J/g/K)": r_["Cp_sapphire"]}) for nm, r_ in res_])
            show(fig, "DSC_cp", exp_, draw_cp, "cp")
            Tt = [25, 50, 100, 150, 200, 250]
            st.dataframe(pd.DataFrame([{"sample": nm, **{f"cp {T_} °C": float(np.interp(T_, r_["T"], r_["Cp"])) for T_ in Tt if r_["T"].min() <= T_ <= r_["T"].max()},
                                        "sapphire check (direct/table)": float(np.median(r_["Cp_sapphire_check"] / r_["Cp_sapphire"])) if r_.get("Cp_sapphire_check") is not None else None}
                                       for nm, r_ in res_]).round(4), hide_index=True, width="stretch")
            for nm, r_ in res_:
                store_put(nm, "dsc", {"cp_25C": float(np.interp(25, r_["T"], r_["Cp"])) if r_["T"].min() <= 25 else None,
                                      "cp_100C": float(np.interp(100, r_["T"], r_["Cp"])) if r_["T"].min() <= 100 <= r_["T"].max() else None}, "DSC cp")
        except Exception as e:
            st.warning(f"cp evaluation failed: {e}")


# -------------------------------------------------------------- defect calculator
def ND_sites(nd):
    from mfs import nodes as N
    return N.sites(nd)


def _res_lab(metal, nd):
    from mfs import nodes as N
    if metal in ("Zr", "Hf", "mixed"):
        return "MO₂"
    r = nd.get("residue")
    return r[0] if r else N.residue_label(metal)


def guess_framework(name):
    """v2.11: framework + metal from a sample name ('MOF-801(Hf)', 'UiO-66-NH2 Zr', 'HKUST-1') -> (framework, metal|None)."""
    import re as _re
    n = (name or "").lower().replace("_", "-").replace(" ", "")
    fw = next((k for k in sorted(FRAMEWORKS, key=len, reverse=True) if k != "custom framework"
               and k.lower().replace(" ", "").replace("(", "").replace(")", "") in n.replace("(", "").replace(")", "")), None)
    met = None
    for m in ("Hf", "Zr", "Ce", "Th", "Cu", "Zn", "Al", "Fe", "Cr", "Mg", "Co", "Ni", "Ti"):
        if _re.search(r"(^|[^a-z])" + m.lower() + r"($|[^a-z])", (name or "").lower()):
            met = m; break
    return fw, met


def default_defect_input(sample="my sample", fw=None):
    """Neutral defaults (UiO-66 / Zr); a sample name that names a framework/metal pre-selects them."""
    gfw, gmet = guess_framework(sample)
    fw = fw or gfw or "UiO-66"
    d = FRAMEWORKS[fw]; lk = d["linker"]; f_, z_, _r, nH_, _c = SPECIES[lk]
    from mfs import nodes as N
    nd = N.get(d.get("node"))
    metal = gmet if gmet in nd["metals"] else d.get("metal", nd["metals"][0])
    sp = [D.Species(lk, f_, z_, "linker", nH_, None, False, True)]
    if z_ >= 2 and nd["core"] != "{M}":
        sp.append(D.Species("Formate", "CHO2", 1, "cap", 1, None, False, False))
    return REP.input_to_dict(D.DefectInput(sample=sample, framework=fw, metal=metal, species=sp))


def node_picker(where, fw, key, node_default="", metal_default=None, custom_default=None):
    """v2.11 node library picker shared by the defect calculator and the TGA page. Returns
    {"node": key, "metal": str, "node_custom": dict|None, "def": node definition}."""
    from mfs import nodes as N
    from mfs.chem import AW as _AW
    fwd = FRAMEWORKS.get(fw, {})
    keys = list(N.NODES) + [N.CUSTOM]
    nd0 = node_default or fwd.get("node") or N.DEFAULT_NODE
    c1, c2 = where.columns([2, 1])
    node = c1.selectbox("Node (SBU)", keys, index=keys.index(nd0) if nd0 in keys else 0, key=f"{key}_node|{fw}",
                        help="Metal-oxo/metal node per formula unit. All counts (linkers, caps, Cl, open sites) are per node.")
    custom = None
    if node == N.CUSTOM:
        cd = dict(custom_default or {})
        e = where.container(border=True)
        e.caption("Custom node: metal atoms per node, oxidation state, core O²⁻ / OH⁻ and the TGA residue. "
                  "Charge for ligands = n·val − 2·O − OH.")
        k1, k2, k3, k4, k5, k6 = e.columns(6)
        cm = k1.text_input("metal", cd.get("metal", "Zn"), key=f"{key}_cm|{fw}").strip() or "Zn"
        if cm not in _AW:
            e.error(f"Unknown element '{cm}'"); cm = "Zn"
        cn_ = k2.number_input("metals/node", 1, 24, int(cd.get("nM", 2)), key=f"{key}_cn|{fw}")
        cv = k3.number_input("oxidation state", 1, 6, int(cd.get("val", 2)), key=f"{key}_cv|{fw}")
        co = k4.number_input("core O²⁻", 0.0, 48.0, float(cd.get("core_O", 0)), 1.0, key=f"{key}_co|{fw}")
        coh = k5.number_input("core OH⁻", 0.0, 48.0, float(cd.get("core_OH", 0)), 1.0, key=f"{key}_coh|{fw}")
        rs = (cd.get("residue") or [N.residue_label(cm), N.RESIDUE.get(cm, (None, 1))[1]])
        rf = k6.text_input("TGA residue", rs[0], key=f"{key}_crf|{fw}|{cm}", help="oxide formula, e.g. Co3O4")
        try:
            from mfs.chem import molar_mass as _mm
            _mm(rf); rn = max(1, int(round(_mm_count(rf, cm))))
        except Exception:
            e.error(f"Cannot parse residue '{rf}'"); rf, rn = N.residue_label(cm), N.RESIDUE.get(cm, (None, 1))[1]
        custom = dict(nM=int(cn_), val=int(cv), core_O=co, core_OH=coh, metal=cm, residue=[rf, rn])
        nd = N.custom_node(**custom)
        e.caption(f"→ core {nd['core'].format(M=cm)}, {N.sites(nd):g} charges for ligands, residue {rf} ({_mm(rf)/rn:.2f} g mol⁻¹ per {cm})")
        metal = cm
    else:
        nd = N.NODES[node]
        ms = nd["metals"]
        m0 = metal_default if metal_default in ms else fwd.get("metal") if fwd.get("metal") in ms else ms[0]
        metal = c2.selectbox("Metal", ms, index=ms.index(m0), key=f"{key}_met|{fw}|{node}")
    return {"node": node, "metal": metal, "node_custom": custom, "def": nd}


def _mm_count(formula, el):
    """Number of atoms of element `el` in a simple formula (residue oxides)."""
    import re as _re
    tot = 0.0
    for sym, n in _re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula):
        if sym == el:
            tot += float(n) if n else 1.0
    return tot or 1.0


def page_defects():
    st.header("Defect calculator – shotgun every method, compare, get the minimal formula")
    st.caption("Per node (M₆O₄(OH)₄, paddlewheel, Zn₄O, M₃O, Al-OH chain, MOF-74, Ti₈, ZIF or your own – see the node library). Enter data once; every applicable method runs in parallel. Equations, symbols, assumptions, your own substitution and "
               "pitfalls are shown for each method. Unit tests reproduce the published worked examples (see README).")
    with st.expander("Load / save inputs (JSON)", expanded=False):
        exs = sorted(glob.glob(os.path.join(EX, "Defects", "defect_input_*.json")))
        exn = st.selectbox("Example input (synthetic EXAMPLE data)", ["(none)"] + [os.path.basename(p) for p in exs])
        upj = st.file_uploader("…or load your own JSON", type=["json"], key="dj")
        if st.button("Load"):
            src = json.load(upj) if upj else (json.load(open(os.path.join(EX, "Defects", exn))) if exn != "(none)" else None)
            if src:
                st.session_state["def_in"] = src; st.rerun()
    base = st.session_state.get("def_in") or default_defect_input("my sample")
    from mfs import samplestore as SS, si_export as SI
    with st.expander("🗂️ Sample store – pull TGA / ¹H NMR / ICP / XRF / N₂ results of one sample from the other pages", expanded=True):
        ids_ = SS.ids(st.session_state)
        if not ids_:
            st.caption("Empty. Open the TGA, NMR digestion, ICP/XRF and N₂ pages with your files (same sample names): every page "
                       f"stores its per-sample results here and in `{STORE_PATH}`.")
        if ids_:
            st.dataframe(pd.DataFrame(SS.summary_rows(st.session_state)), hide_index=True, width="stretch")
        pick = st.selectbox("Pull data for sample", ["(none – type the inputs)"] + list(ids_), format_func=lambda k: ids_.get(k, k), key="def_pull")
        k1, k2, k3 = st.columns(3)
        k1.download_button("⬇️ Store (JSON)", json.dumps(st.session_state.get("sample_store", {}), indent=1, ensure_ascii=False).encode(),
                           "sample_store.json", key="store_dl")
        upst = k2.file_uploader("Load a store (JSON)", type=["json"], key="store_up")
        if upst is not None and k2.button("Merge", key="store_merge"):
            st.session_state.setdefault("sample_store", {}).update(json.load(upst)); st.rerun()
        if k3.button("🗑️ Clear the store", key="store_clear"):
            SS.clear(st.session_state); SS.save(st.session_state, STORE_PATH); st.rerun()
    pulled, pull_note = {}, []
    if pick in ids_:
        rec = SS.get(st.session_state, ids_[pick]) or {"sections": {}}
        secs = {k: v["values"] for k, v in rec["sections"].items()}
        t_ = secs.get("tga") or {}
        if t_:
            pulled.update(tga_W_initial=t_.get("W_initial"), tga_W_200=t_.get("W200"), tga_W_plateau=t_.get("W_plateau"), tga_W_end=t_.get("W_end"))
            pull_note.append("TGA")
        i_ = secs.get("icp") or {}
        if i_:
            pulled.update(icp_Zr_wt=i_.get("Zr_wt"), icp_Hf_wt=i_.get("Hf_wt")); pull_note.append("ICP")
        x_ = secs.get("xrf") or {}
        if x_.get("xrf_S_per_M") is None and i_.get("xrf_S_per_M") is not None:
            x_ = {**x_, "xrf_S_per_M": i_["xrf_S_per_M"], "xrf_S_per_M_sd": i_.get("xrf_S_per_M_sd")}
        if x_.get("xrf_S_per_M") is not None:
            pulled["xrf_S_per_M"] = x_["xrf_S_per_M"]
            pull_note.append("S/M" + (f" = {x_['xrf_S_per_M']:.3f} ± {x_['xrf_S_per_M_sd']:.3f}" if x_.get("xrf_S_per_M_sd") is not None else ""))
        if (x_.get("Hf_frac") or i_.get("Hf_frac")) is not None:
            hf_ = x_.get("Hf_frac") if x_.get("Hf_frac") is not None else i_.get("Hf_frac")
            sd_ = x_.get("Hf_frac_sd") if x_.get("Hf_frac") is not None else i_.get("Hf_frac_sd")
            pull_note.append(f"Hf/(Zr+Hf) = {hf_:.3f}" + (f" ± {sd_:.3f}" if sd_ is not None else ""))
        n2_ = secs.get("n2") or {}
        if n2_.get("V_total") is not None:
            pulled["V_meas"] = n2_["V_total"]; pull_note.append("N₂ V_total")
        nm_ = secs.get("nmr") or {}
        if nm_:
            pulled["_nmr"] = nm_; pull_note.append("¹H NMR integrals")
        pulled["sample"] = rec.get("name", pick)
        if not st.session_state.get("def_in"):          # v2.11: framework/metal/linker defaults follow the sample name
            base = default_defect_input(pulled["sample"])
        st.success(f"Pulled for **{pulled['sample']}**: " + (", ".join(pull_note) or "nothing (no results stored yet)")
                   + ". Fields below are pre-filled and can still be edited.")
    base = dict(base); base.update({k: v for k, v in pulled.items() if not k.startswith("_") and v is not None})
    if pulled.get("_nmr"):      # integrals into the species table, matched by species name
        sp_ = [dict(x) for x in base.get("species", [])]
        for x in sp_:
            key = next((k for k in pulled["_nmr"] if k.startswith("I ") and x["name"].lower() in k.lower()), None)
            if key and pulled["_nmr"][key] is not None:
                x["integral"] = pulled["_nmr"][key]
        if not any(x.get("name", "").lower().startswith("tmsp") for x in sp_):
            pass
        base["species"] = sp_
        std = next((v for k, v in pulled["_nmr"].items() if k.startswith("I TMSP")), None)
        if std:
            base["nmr_std_integral"] = std
    sfx = "" if pick not in ids_ else "_" + pick[:20]
    g = lambda k, d=None: base.get(k, d)
    st.subheader("1 · Framework, node and metal")
    c1, c2 = st.columns([1, 1])
    sample = c1.text_input("Sample", g("sample"), key="def_sample" + sfx)
    fw0 = g("framework", "UiO-66"); fw0 = fw0 if fw0 in FRAMEWORKS else "custom framework"
    fw = c2.selectbox("Framework", list(FRAMEWORKS), index=list(FRAMEWORKS).index(fw0), key="def_fw" + sfx,
                      help="Picking a framework sets its node, metal and main linker (all editable below).")
    nsel = node_picker(st, fw, "def", g("node", "") if fw == fw0 else "", g("metal") if fw == fw0 else None,
                       g("node_custom") if fw == fw0 else None)
    metal, ND = nsel["metal"], nsel["def"]
    c1, c2, c3 = st.columns(3)
    cn_def = int(g("CN") or FRAMEWORKS[fw]["CN"]) if fw == fw0 else int(FRAMEWORKS[fw]["CN"] if not nsel["node_custom"] else ND["CN"])
    CN = c1.number_input("Ideal connectivity CN (charge units)", 1, 24, cn_def, key=f"def_cn|{fw}|{nsel['node']}" + sfx,
                         help="Linker donor groups per node in the defect-free framework, counted in charge units "
                              f"(12 for UiO-66, 8 NU-1000, 6 MOF-808 / Zn4O / M3O, 4 paddlewheel / MOF-74, 2 Al-OH chain / ZIF). "
                              f"This node offers {ND_sites(ND):g} charges for ligands.")
    f_hf = c2.number_input("Hf/(Hf+Zr) (mixed; blank → from ICP)", 0.0, 1.0, float(g("f_hf") or 0.0)) if metal == "mixed" else None
    extra = c3.checkbox("All Zr is extra-framework ZrO₂ (ZrO₂@Hf-MOF composites) – correct TGA residue, ICP, porosity",
                        bool(g("extraframework_Zr"))) if metal in ("Hf", "mixed") else False
    if fw != fw0 and not st.session_state.get("def_in"):
        base = dict(base); base["species"] = default_defect_input(sample, fw)["species"]
    with st.expander("🧱 Node library – metals, charges, TGA residues (pick above or build a custom node)", expanded=False):
        from mfs import nodes as _N
        st.dataframe(pd.DataFrame(_N.table_rows()), hide_index=True, width="stretch")
        st.caption("Charge for ligands = metals × oxidation state − 2·core O²⁻ − core OH⁻. TGA residue in air per metal "
                   "(Co → Co₃O₄, Fe/Cr/Al → M₂O₃, Mn → Mn₂O₃). Al-OH chains keep their µ-OH at the plateau; M₆O₄(OH)₄ and Ti₈O₈(OH)₄ "
                   "dehydroxylate first. Linker/cap molar masses are computed from the formula you type.")
    st.subheader("2 · Species (main linker, caps, guests) and ¹H NMR integrals")
    st.caption("Bound (deprotonated) formula; z = carboxylates bound (= node sites occupied, charge); role cap = monocarboxylate modulator. "
               "'at plateau' = still present at the ~400 °C TGA plateau (e.g. benzoate, sulfonated linkers). Integrals on any common scale "
               "(set the TMSP-d4 standard to 9 for the qNMR method).")
    sp = pd.DataFrame(g("species") if fw == fw0 else base["species"])
    sp = st.data_editor(sp, num_rows="dynamic", key="def_sp" + sfx + ("" if fw == fw0 else "|" + fw), column_config={
        "role": st.column_config.SelectboxColumn(options=["linker", "cap", "neutral"]),
        "name": st.column_config.TextColumn(help="Pick from the library below or type your own")})
    _mm = []
    for r_ in sp.to_dict("records"):
        try:
            _mm.append(f"{r_['name']}: {molar_mass(str(r_['formula'])):.2f}" if r_.get("name") and r_.get("formula") else None)
        except Exception as e_:
            _mm.append(f"⚠️ {r_.get('name')}: {e_}")
    st.caption("Molar masses (bound form, g mol⁻¹, computed from the formula – type any custom linker/cap): " + " · ".join(x for x in _mm if x))
    with st.expander("Species library"):
        lib = pd.DataFrame([(k, *v) for k, v in SPECIES.items()], columns=["name", "formula", "z", "role", "nH", "¹H signal"])
        lib["M (g/mol)"] = [round(molar_mass(f), 2) for f in lib.formula]
        st.dataframe(lib, width="stretch", hide_index=True)
    nmr = st.session_state.get("nmr_integrals")
    if nmr is not None:
        st.caption("Integrals from the NMR tab are available:"); st.dataframe(nmr.filter(regex="sample|^I "), width="stretch")
    st.subheader("3 · TGA (in air, % of initial mass)")
    up = st.file_uploader("Optional: TGA file → auto-read W(200 °C), plateau, residue", key="def_tga")
    auto = {}
    if up:
        t = R.read_tga(up, up.name); a = A.tga_analysis(t["T"], t["w"])
        auto = {"tga_W_initial": a["W_initial"], "tga_W_200": a["W200"], "tga_W_plateau": a["W_plateau"], "tga_W_end": a["W_end"]}
        st.info(f"Auto: W_i = {a['W_initial']:.2f}, W(200) = {a['W200']:.2f}, plateau {a['W_plateau']:.2f} % at {a['T_plateau']:.0f} °C (DTG minimum 250–480 °C), residue {a['W_end']:.2f} %")
    gv = lambda k: auto.get(k, g(k))
    c1, c2, c3, c4 = st.columns(4)
    num = lambda col, lab, key, v, **kw: col.number_input(lab, value=None if v is None else float(v), **kw, key="dn_" + key + sfx)
    W_i = num(c1, "W initial (RT)", "wi", gv("tga_W_initial") or 100.0)
    W200 = num(c2, "W at 200 °C (M-200)", "w200", gv("tga_W_200"))
    Wpl = num(c3, "W at plateau (M-400)", "wpl", gv("tga_W_plateau"))
    Wend = num(c4, f"W residue ({_res_lab(metal, ND)})", "wend", gv("tga_W_end"))
    st.subheader("4 · qNMR with internal standard (De Roo protocol)")
    c1, c2, c3, c4 = st.columns(4)
    std_c = num(c1, "[TMSP-d4] (mM)", "stdc", g("nmr_std_conc_mM"))
    std_I = num(c2, "TMSP integral", "stdi", g("nmr_std_integral") or 9.0)
    mdig = num(c3, "MOF mass digested (mg)", "mdig", g("nmr_sample_mg"), format="%.4f")
    vol = num(c4, "NMR volume (mL)", "vol", g("nmr_volume_mL"))
    st.subheader("5 · Chloride, ICP-OES, XRF")
    c1, c2, c3, c4 = st.columns(4)
    cl_node = num(c1, f"Cl per node (XRF/EDX: {ND['nM']}·Cl/M)", "cln", g("cl_per_node"))
    cl_mM = num(c2, "[Cl⁻] in digest (mM, UV-vis)", "clm", g("cl_conc_mM"))
    cl_wt = num(c3, "Cl wt%", "clw", g("cl_wt_pct"))
    s_m = num(c4, "XRF S/M atomic ratio (check)", "sm", g("xrf_S_per_M"))
    c1, c2, c3, c4 = st.columns(4)
    if metal in ("Zr", "Hf", "mixed"):
        zr = num(c1, "ICP Zr (wt%)", "zr", g("icp_Zr_wt"))
        hf = num(c2, "ICP Hf (wt%)", "hf", g("icp_Hf_wt")); icpM = None
    else:
        zr = hf = None
        icpM = num(c1, f"ICP {metal} (wt%)", "icpm_" + metal, g("icp_M_wt"))
        c2.caption(f"framework metal {metal}: {ND['nM']} per node")
    state = c3.selectbox("ICP sample weighed", ["as-measured (RT)", "dried 200 °C"], index=0 if "RT" in g("icp_state", "RT") else 1)
    lum = num(c4, "qNMR linker content (µmol mg⁻¹)", "lum", g("linker_umol_per_mg"), format="%.5f")
    st.subheader("6 · Porosity (semi-quantitative) and titration")
    c1, c2, c3, c4 = st.columns(4)
    Vm = num(c1, "V measured (cm³ g⁻¹)", "vm", g("V_meas"), format="%.4f")
    Vi = num(c2, "V ideal (cm³ g⁻¹)", "vi", g("V_ideal"), format="%.4f")
    vL = num(c3, "Volume freed per missing linker (Å³)", "vl", g("v_linker_A3"))
    if metal in ("Zr", "Hf") and fw in ("UiO-66", "MOF-801") and c4.button("Compute ideal V from bundled CIF"):
        from mfs.porosity import probe_pore_volume
        cif = os.path.join(REF, "UiO-66_Zr_COD4512074.cif" if fw.startswith("UiO") else "MOF-801_Zr_COD4121459.cif")
        r = probe_pore_volume(cif, substitute={"Zr": "Hf"} if metal == "Hf" else None)
        st.success(f"Probe-occupiable (N₂, 1.82 Å) pore volume: {r['V_cm3_g']:.3f} cm³ g⁻¹ (void fraction {r['void_fraction']:.2f})")
    c1, c2, c3, c4 = st.columns(4)
    tm = num(c1, "Titrated mass (mg)", "tm", g("tit_mass_mg"))
    tc = num(c2, "NaOH (M)", "tc", g("tit_c_NaOH") or 0.1)
    veq = c3.text_input("Equivalence volumes (mL, comma)", ", ".join(f"{v:.4f}" for v in (g("tit_V_eq") or [])))
    hcl = num(c4, "HCl added (mmol)", "hcl", g("tit_HCl_mmol"), format="%.4f")
    tup = st.file_uploader("Optional: titration curve (V, pH) → equivalence points from dpH/dV", key="tit")
    if tup:
        df = R.read_table(tup); eq = A.equivalence_points(df.iloc[:, 0].values, df.iloc[:, 1].values, 3)
        st.info("Equivalence volumes: " + ", ".join(f"{v:.3f}" for v in eq) + " mL – copy into the field above if correct.")
        fig, ax = newfig(0.6)
        ax.plot(df.iloc[:, 0], df.iloc[:, 1], "k-", lw=0.9); ax.set_xlabel("V(NaOH) (mL)"); ax.set_ylabel("pH")
        tw = ax.twinx(); tw.plot(df.iloc[:, 0], np.gradient(A.smooth(df.iloc[:, 1].values, "Savitzky-Golay", 15, 2), df.iloc[:, 0].values), color="#0072B2", lw=0.8)
        tw.set_ylabel("dpH/dV (mL$^{-1}$)", color="#0072B2"); ax.tick_params(right=False, which="both")
        for v in eq: ax.axvline(v, color="#D55E00", lw=0.6, ls="--")
        show(fig, "titration")
    reo = st.selectbox("Low-angle reo reflections (missing-cluster check) observed in PXRD?", ["not checked", "yes", "no"])
    try:
        species = [D.Species(**{k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in r.items()}) for r in sp.to_dict("records") if r.get("name")]
        for s in species:
            s.z = int(s.z); s.nH = int(s.nH or 0); s.at_plateau = bool(s.at_plateau); s.main = bool(s.main)
        inp = D.DefectInput(sample=sample, framework=fw, CN=int(CN), metal=metal, f_hf=f_hf, extraframework_Zr=extra, species=species,
                            node=nsel["node"] if nsel["node"] != (FRAMEWORKS[fw].get("node") or "") else "", node_custom=nsel["node_custom"], icp_M_wt=icpM,
                            tga_W_initial=W_i, tga_W_200=W200, tga_W_plateau=Wpl, tga_W_end=Wend, nmr_std_conc_mM=std_c, nmr_std_integral=std_I,
                            nmr_sample_mg=mdig, nmr_volume_mL=vol, cl_per_node=cl_node, cl_conc_mM=cl_mM, cl_wt_pct=cl_wt, icp_Zr_wt=zr, icp_Hf_wt=hf,
                            icp_state=state, linker_umol_per_mg=lum, xrf_S_per_M=s_m, V_meas=Vm, V_ideal=Vi, v_linker_A3=vL, tit_mass_mg=tm,
                            tit_c_NaOH=tc, tit_V_eq=[float(v) for v in veq.split(",") if v.strip()] or None, tit_HCl_mmol=hcl,
                            reo_peaks_seen=None if reo == "not checked" else reo == "yes")
        res = D.run_all(inp)
    except Exception as e:
        st.error(f"Calculation error: {e}"); return
    st.markdown("---")
    st.subheader("Results")
    cons = res["consensus"]
    if cons:
        c = st.columns(5)
        c[0].metric("Linkers / node (median)", f"{cons['L']:.2f}", f"range {cons['L_min']:.2f}–{cons['L_max']:.2f}", delta_color="off")
        c[1].metric("Missing linkers / node", f"{cons['x']:.2f}")
        c[2].metric("Caps / node", f"{cons['m']:.2f}")
        c[3].metric("Open sites / node", f"{cons['OS']:.2f}")
        c[4].metric("M per node (g mol⁻¹)", f"{cons['M']:.1f}")
        if not inp.is_zrhf() or ND["core"] != "{M}6O4(OH)4":
            st.caption(f"Node {ND['core'].format(M=inp.M_symbol())}: {inp.nM()} {inp.M_symbol()} per node, {inp.Q():g} charges for ligands, "
                       f"TGA residue {inp.residue_label()} ({inp.MMO2():.2f} g mol⁻¹ per metal). The equations below are written for "
                       "M₆O₄(OH)₄; the substitutions use these node values.")
        st.markdown(f"**Minimal formula** (one H₂O per terminal OH⁻/Cl⁻): `{cons['formula']}`  \n**Anhydrous**: `{cons['formula_anhydrous']}` "
                    f"({cons['M_anhydrous']:.1f} g mol⁻¹) · metal {cons['metal_wt_pct']:.1f} wt%")
        if s_m:
            st.caption(f"XRF cross-check: S per node = 6·S/M = {6*s_m:.2f} (compare with the sulfonated species per node in the formula).")
    tab = D.results_table(res)
    st.dataframe(tab, width="stretch", hide_index=True)
    if cons:
        used = [m.title.split(" ·")[0].split(" –")[0] for m in res["methods"] if m.status == "ok" and m.L is not None and not m.semi_quant]
        store_put(sample, "defects", {"framework": fw, "metal": metal, "node": ND["core"].format(M=inp.M_symbol()), "L": cons["L"], "x": cons["x"], "m": cons["m"], "OS": cons["OS"],
                                      "M": cons["M"], "formula": cons["formula"], "L_range": f"{cons['L_min']:.2f}–{cons['L_max']:.2f}",
                                      "methods": "; ".join(used), "metal_wt_pct": cons.get("metal_wt_pct")}, "Defect calculator")
    with st.expander("📑 SI table – one consistent defect report per sample (all samples in the store)", expanded=False):
        sit = SI.defect_rows(st.session_state.get("sample_store", {}))
        if len(sit):
            st.dataframe(SI._round(sit), width="stretch", hide_index=True)
            notes = ("Defect analysis per node (Reticula " + __version__ + "). L = linkers per node (median of the methods listed), "
                     "x = missing linkers, caps = monocarboxylate modulators, open sites = 2x - caps. TGA in air; S_BET from the "
                     "Rouquerol-consistent range (ISO 9277); V_total at p/p0 = 0.95.")
            q1, q2 = st.columns(2)
            dl(q1, "⬇️ SI table (XLSX)", SI.to_xlsx(sit, tab.assign(sample=sample), notes), "SI_defects.xlsx", "si_xlsx")
            try:
                dl(q2, "⬇️ SI table (DOCX)", SI.to_docx(sit, notes=notes), "SI_defects.docx", "si_docx")
            except Exception as e:
                q2.caption(f"DOCX export unavailable ({e})")
        else:
            st.caption("No defect results stored yet – every computed sample is added automatically.")
    fig = P.fig_defect_comparison(res, max(WIDTH, 7.0))
    show(fig, f"defects_{re.sub(r'[^A-Za-z0-9]+', '_', sample)}")
    c1, c2, c3 = st.columns(3)
    sn = safe_name(sample)

    def _pdf(res=res):
        import tempfile
        with _RENDER_LOCK, tempfile.TemporaryDirectory() as td:
            pth = os.path.join(td, "defect_report.pdf")
            REP.pdf_report(res, pth)
            return open(pth, "rb").read()
    dl(c1, "📄 PDF report (full calculation)", _pdf, f"defect_report_{sn}.pdf", "def_pdf", help="generated when you click (takes a few s)")
    dl(c2, "CSV (table + steps)", REP.csv_report(res).encode("utf-8-sig"), f"defect_report_{sn}.csv", "def_csv")
    dl(c3, "JSON (inputs, reload later)", json.dumps(REP.input_to_dict(inp), indent=1).encode(), f"defect_input_{sn}.json", "def_json")
    if st.button("💾 Save report to folder", key="def_save"):
        os.makedirs(EXPORT_DIR, exist_ok=True)
        pth = os.path.join(EXPORT_DIR, f"defect_report_{sn}.pdf")
        open(pth, "wb").write(_pdf())
        open(os.path.join(EXPORT_DIR, f"defect_report_{sn}.csv"), "wb").write(REP.csv_report(res).encode("utf-8-sig"))
        st.success(f"Saved to `{EXPORT_DIR}`")
    st.subheader("Method by method")
    h = res["h"]
    for m in res["methods"]:
        icon = {"ok": "✅", "warn": "⚠️", "n/a": "▫️"}[m.status]
        with st.expander(f"{icon} {m.title}" + (f" — L = {m.L:.2f}, x = {m.x:.2f}" if m.L is not None else f" — {m.message}")):
            if m.semi_quant:
                st.warning("Semi-quantitative – excluded from the consensus.")
            if m.message:
                st.caption(m.message)
            st.markdown("**Equations**")
            for e in m.equations:
                st.latex(e)
            st.markdown("**Symbols**")
            st.table(pd.DataFrame(list(m.symbols.items()), columns=["symbol", "meaning"]))
            st.markdown("**Assumptions**\n" + "\n".join(f"- {a}" for a in m.assumptions))
            if m.steps:
                st.markdown("**Your numbers, step by step**")
                for s_ in m.steps:
                    st.latex(s_)
            if m.formula:
                st.markdown(f"**Formula**: `{m.formula}` (M = {m.M_formula:.1f} g mol⁻¹)")
            st.markdown("**Pitfalls**\n" + "\n".join(f"- {p}" for p in m.pitfalls))
            st.markdown("**References**\n" + "\n".join(f"- {c} [doi:{d}](https://doi.org/{d})" for c, d in m.refs))
    with st.expander(h["title"], expanded=True):
        for e in h["equations"]:
            st.latex(e)
        st.markdown("\n".join(f"- {a}" for a in h["assumptions"]))
        st.markdown("References: " + "; ".join(f"[{c.split('.')[0]}](https://doi.org/{d})" for c, d in h["refs"]))
    with st.expander("When and why the methods disagree"):
        st.table(pd.DataFrame(D.SUMMARY_PITFALLS, columns=["pattern", "likely cause / what to do"]))


def page_kinetics():
    st.header("Hydrolysis kinetics (¹H NMR) and adsorption isotherms")
    st.caption("Table with columns sample, time, concentration (or integrals vs TMSP-d4). Replicates are averaged (mean ± SD) and all points are fitted.")
    files = files_input("Kinetics table (CSV/XLSX)", ["csv", "txt", "xlsx"], "Kinetics/GlyGly_hydrolysis_60C*.csv", "kin", multiple=False)
    if st.checkbox("➕ Use the temperature-series EXAMPLE instead (40–70 °C, for Arrhenius / Eyring)", False, key="kin_tser"):
        files = [_F(os.path.join(EX, "Kinetics", "GlyGly_hydrolysis_Tseries_EXAMPLE.csv"))]
    if st.session_state.get("_uv_kin_csv") and st.checkbox("➕ Use the table sent from UV/Vis (A at one λ vs time)", True, key="kin_from_uv"):
        _b = io.BytesIO(st.session_state["_uv_kin_csv"]); _b.name = "uvvis_kinetics.csv"; files = [_b]
        st.caption(f"Time unit from UV/Vis: {st.session_state.get('_uv_kin_unit', 'min')} – set the same unit below.")
    from mfs import kinetics as K
    if files:
        df = R.read_csv_table(files[0], (r"sample", r"time"))
        st.dataframe(df.head(), width="stretch")
        cols = list(df.columns)
        c1, c2, c3, c4 = st.columns(4)
        sc = c1.selectbox("sample column", cols, index=cols.index(R.find_col(df, "sample") or cols[0]))
        tc = c2.selectbox("time column", cols, index=cols.index(R.find_col(df, "time") or cols[1]))
        df = _apply_sample_list(df, sc, "kin_samples")
        mode = c3.selectbox("y is", ["concentration", "integral vs internal standard"])
        yc = c4.selectbox("y column", cols, index=len(cols) - 1)
        if mode != "concentration":
            c1, c2, c3 = st.columns(3)
            stdc = c1.selectbox("standard integral column", cols, index=len(cols) - 1)
            nh = c2.number_input("protons in peptide signal", 1, 20, 4)
            cstd = c3.number_input("[TMSP-d4] (mM)", 0.01, 100.0, 1.0)
            df["_c"] = A.nmr_concentration(df[yc], nh, df[stdc], 9, cstd); yc = "_c"
        c1, c2, c3, c4 = st.columns(4)
        tunit = c1.selectbox("time unit", ["h", "min", "s", "d"])
        offset = c2.checkbox("fit with plateau C∞", False)
        norm = c3.checkbox("plot C/C₀", True)
        errbar = c4.radio("error bars on k", ["95 % CI", "1 SE"], horizontal=True, key="kin_eb")
        Tcol = R.find_col(df, r"^t\s*\(°?c\)|temp|^t_c$"); pDcol = R.find_col(df, r"^pd\b|^ph\b|pd$")
        gcols = [sc] + [c_ for c_ in (Tcol, pDcol) if c_ is not None and df[c_].nunique() > 1]
        sets, rows, fits = [], [], []
        for key_, g_ in df.groupby(gcols, sort=False):
            key_ = key_ if isinstance(key_, tuple) else (key_,)
            lab = str(key_[0]) + "".join(f", {v:g} °C" if c_ == Tcol else f", pD {v:g}" for c_, v in zip(gcols[1:], key_[1:]))
            try:
                fit = K.first_order_ci(g_[tc].values, g_[yc].values, offset=offset)
            except Exception as e:
                st.warning(f"{lab}: fit failed ({e})"); fit = None
            agg = g_.groupby(tc)[yc].agg(["mean", "std"]).reset_index()
            sets.append((lab, agg[tc].values, agg["mean"].values, agg["std"].fillna(0).values, fit))
            if fit:
                Tv = float(g_[Tcol].iloc[0]) if Tcol is not None else None; pDv = float(g_[pDcol].iloc[0]) if pDcol is not None else None
                fits.append({"label": lab, "sample": str(key_[0]), "T": Tv, "pD": pDv, "fit": fit})
                rows.append({"sample": lab, f"k_obs ({tunit}-1)": fit["k"], "± 95 % CI": fit["k_ci"], "SE": fit["k_se"],
                             f"t1/2 ({tunit})": fit["t_half"], "± t1/2 (95 %)": fit["t_half_ci"], "C0": fit["C0"], "R2": fit["R2"], "n": fit["n"]})
                store_put(str(key_[0]), "kinetics", {"k": fit["k"], "k_ci95": fit["k_ci"], "unit": f"{tunit}-1", "T_C": Tv, "pD": pDv,
                                                      "R2": fit["R2"], "label": lab}, "Kinetics page")
        tab = pd.DataFrame(rows); st.dataframe(tab.style.format(precision=4), width="stretch")
        st.caption("95 % CI = t(n − p, 0.975) × standard error from the fit covariance (replicates fitted point by point).")
        def draw(ax, sets=sets):
            P.draw_kinetics(ax, sets, ylabel=r"$C/C_0$" if norm else "[peptide] (mM)", time_label=f"Time ({tunit})", normalise=norm,
                            annotate_k=False, k_unit=f"{tunit}$^{{-1}}$")
        fig, ax = newfig(); draw(ax)
        show(fig, "kinetics", tab, draw, "Kinetics", units={"x": "time_h"})
        def drawk(ax, sets=sets):
            ok = [s for s in sets if s[4]]
            P.draw_k_bars(ax, [s[0] for s in ok], [s[4]["k"] for s in ok], [s[4]["k_ci" if errbar.startswith("95") else "k_se"] for s in ok],
                          ylabel=rf"$k_\mathrm{{obs}}$ (10$^{{-3}}$ {tunit}$^{{-1}}$)")
        fig, ax = newfig(); drawk(ax)
        show(fig, "k_obs_bars", None, drawk, "k_obs bars")
        with st.expander("⏱️ Initial rates (linear fit up to a chosen conversion)"):
            frac = st.slider("use points up to conversion", 0.02, 0.5, 0.10, 0.01, key="kin_ir_frac")
            ir = []
            for key_, g_ in df.groupby(gcols, sort=False):
                try:
                    r_ = K.initial_rate(g_[tc].values, g_[yc].values, frac)
                    ir.append({"sample": str(key_ if not isinstance(key_, tuple) else ", ".join(map(str, key_))), f"r0 (mM {tunit}-1)": r_["r0"],
                               "± 95 % CI": r_["r0_ci"], "points": r_["n"], f"up to t ({tunit})": r_["t_max"], "R2": r_["R2"]})
                except Exception as e:
                    st.caption(f"{key_}: {e}")
            st.dataframe(pd.DataFrame(ir).style.format(precision=5), width="stretch", hide_index=True)
        with st.expander("⚗️ TOF per open site (open sites per node from the defect calculator)"):
            from mfs import samplestore as SS
            dstore = {k_: v_ for k_, v_ in SS.ids(st.session_state).items() if SS.get(st.session_state, v_, "defects")}
            q1, q2, q3, q4 = st.columns(4)
            src = q1.selectbox("Catalyst (defect result)", ["(type values)"] + list(dstore), format_func=lambda k_: dstore.get(k_, k_), key="tof_src")
            dv = SS.get(st.session_state, dstore[src], "defects") if src in dstore else {}
            osn = q2.number_input("open sites per node", 0.0, 12.0, float(dv.get("OS") or 1.0), key="tof_os_" + str(src))
            Mn = q3.number_input("M per node (g/mol)", 100.0, 5000.0, float(dv.get("M") or 1650.0), key="tof_M_" + str(src))
            mcat = q4.number_input("catalyst (mg)", 0.01, 1000.0, 2.0, key="tof_m")
            q1, q2 = st.columns(2)
            c0 = q1.number_input("substrate C0 (mM)", 0.001, 1000.0, 2.0, key="tof_c0")
            vol = q2.number_input("volume (mL)", 0.01, 100.0, 0.6, key="tof_v")
            per_s = {"h": 3600, "min": 60, "s": 1, "d": 86400}[tunit]
            tr = []
            for f_ in fits:
                r_ = K.tof(f_["fit"]["k"] / per_s, c0, vol, mcat, osn, Mn, f_["fit"]["k_ci"] / per_s)
                tr.append({"sample": f_["label"], "open sites (µmol)": r_["n_open_sites_umol"], "TOF (h-1)": r_["TOF_h"],
                           "± 95 % (h-1)": r_["TOF_err_s"] * 3600 if r_["TOF_err_s"] else None})
            st.dataframe(pd.DataFrame(tr).style.format(precision=4), width="stretch", hide_index=True)
            st.caption("TOF = k_obs·C0·V / (m_cat / M_node × open sites per node) – initial turnover per open site. "
                       "Pull OS and M from the defect calculator (sample store) or type them.")
        with st.expander("ln(C/C₀) linearisation"):
            fig, ax = newfig()
            cols_ = S.colours_for([s[0] for s in sets])
            for c_, (nm, t, cc, e, fit) in zip(cols_, sets):
                ax.plot(t, np.log(np.asarray(cc) / (fit["C0"] if fit else cc[0])), "o", color=c_, ms=3, label=nm)
                if fit: ax.plot(t, -fit["k"] * np.asarray(t), "-", color=c_, lw=0.8)
            ax.set_xlabel(f"Time ({tunit})"); ax.set_ylabel(r"ln($C/C_0$)"); ax.legend(fontsize=FONT - 2)
            show(fig, "ln_kinetics")
    st.subheader("Adsorption isotherm (Langmuir / Freundlich)")
    files = files_input("Adsorption table (sample, Ce, qe)", ["csv", "txt", "xlsx"], "Kinetics/GlyGly_adsorption_*.csv", "ads", multiple=False)
    if files:
        df = R.read_csv_table(files[0], (r"sample", r"ce"))
        cols = list(df.columns)
        c1, c2, c3 = st.columns(3)
        sc = c1.selectbox("sample", cols, key="as"); cc = c2.selectbox("Ce", cols, index=1, key="ac"); qc = c3.selectbox("qe", cols, index=2, key="aq")
        sets, rows = [], []
        for s_, g_ in df.groupby(sc, sort=False):
            fit = A.langmuir_fit(g_[cc], g_[qc]); fr = A.freundlich_fit(g_[cc].values, g_[qc].values)
            sets.append((str(s_), g_[cc].values, g_[qc].values, fit))
            rows.append({"sample": s_, "q_max": fit["q_max"], "± q_max": fit["q_max_err"], "K_L": fit["K"], "± K_L": fit["K_err"], "R2 Langmuir": fit["R2"],
                         "K_F": fr["K_F"], "n_F": fr["n"]})
        st.dataframe(pd.DataFrame(rows).style.format(precision=4), width="stretch")
        def drawl(ax, sets=sets):
            P.draw_langmuir(ax, sets)
        fig, ax = newfig(); drawl(ax)
        show(fig, "langmuir", pd.DataFrame(rows), drawl, "Langmuir", units={"x": "conc_mM"})
    st.subheader("Arrhenius and Eyring analysis (several temperatures)")
    fT = [f_ for f_ in (locals().get("fits") or []) if f_["T"] is not None]
    if len({f_["T"] for f_ in fT}) >= 2:
        st.caption("Rates from the fitted kinetics table above (temperature column found). Edit if needed.")
        ey0 = pd.DataFrame({"T_C": [f_["T"] for f_ in fT], "k": [f_["fit"]["k"] for f_ in fT], "k_err (95 %)": [f_["fit"]["k_ci"] for f_ in fT]})
    else:
        ey0 = pd.read_csv(os.path.join(EX, "Kinetics", "Eyring_composite.csv")).rename(columns={"k_h-1": "k"})
        ey0["k_err (95 %)"] = np.nan
    ey = st.data_editor(ey0, num_rows="dynamic", key="ey2_" + str(len(ey0)))
    if len(ey) >= 2:
        unit = st.selectbox("k unit", ["h-1", "min-1", "s-1"], key="ar_unit")
        per = {"h-1": 3600, "min-1": 60, "s-1": 1}[unit]
        kerr = ey["k_err (95 %)"].values if ey["k_err (95 %)"].notna().all() else None
        ar = K.arrhenius(ey["T_C"].values, ey["k"].values / per, None if kerr is None else kerr / per)
        eyr = K.eyring(ey["T_C"].values, ey["k"].values / per, None if kerr is None else kerr / per)
        st.markdown(f"**Arrhenius**: Ea = **{ar['Ea_kJ']:.1f} ± {ar['Ea_ci']:.1f} kJ mol⁻¹** (95 % CI), ln(A/s⁻¹) = {ar['lnA']:.1f} ± {ar['lnA_ci']:.1f}, R² {ar['R2']:.4f}  \n"
                    f"**Eyring**: ΔH‡ = **{eyr['dH_kJ']:.1f} ± {eyr['dH_ci']:.1f} kJ mol⁻¹**, ΔS‡ = **{eyr['dS_J']:.0f} ± {eyr['dS_ci']:.0f} J mol⁻¹ K⁻¹**, "
                    f"ΔG‡(298 K) = {eyr['dG298_kJ']:.1f} kJ mol⁻¹")

        def draw_ar(ax, r=ar, lab=r"ln($k$/s$^{-1}$)"):
            ax.plot(1000 * r["X"], r["Y"], "o", color="#0072B2", ms=4)
            xx = np.linspace(r["X"].min(), r["X"].max(), 10); ax.plot(1000 * xx, r["slope"] * xx + r["intercept"], "-", color="#0072B2", lw=0.9)
            ax.set_xlabel(r"1000/$T$ (K$^{-1}$)"); ax.set_ylabel(lab)

        def draw_ey(ax, r=eyr):
            ax.plot(1000 * r["X"], r["Y"], "o", color="#D55E00", ms=4)
            xx = np.linspace(r["X"].min(), r["X"].max(), 10); ax.plot(1000 * xx, r["slope"] * xx + r["intercept"], "-", color="#D55E00", lw=0.9)
            ax.set_xlabel(r"1000/$T$ (K$^{-1}$)"); ax.set_ylabel(r"ln($k$/$T$)")
        fig, ax = newfig(); draw_ar(ax)
        show(fig, "arrhenius", ey, draw_ar, "Arrhenius")
    if len(ey) >= 2:
        fig, ax = newfig(); draw_ey(ax)
        show(fig, "eyring", ey, draw_ey, "Eyring")
    st.subheader("pD (pH)–rate profile")
    fpd = [f_ for f_ in (locals().get("fits") or []) if f_["pD"] is not None]
    if len({f_["pD"] for f_ in fpd}) >= 3:
        pdt0 = pd.DataFrame({"pD": [f_["pD"] for f_ in fpd], "k": [f_["fit"]["k"] for f_ in fpd], "± k": [f_["fit"]["k_ci"] for f_ in fpd]})
    else:
        pdt0 = pd.read_csv(os.path.join(EX, "Kinetics", "pD_rate_profile_EXAMPLE.csv"), comment="#")
        pdt0.columns = ["pD", "k", "± k"]
        st.caption("EXAMPLE values (synthetic) – replace with your own k_obs vs pD (or add a pD column to the kinetics table).")
    pdt = st.data_editor(pdt0, num_rows="dynamic", key="pdt_" + str(len(pdt0)))
    c1, c2 = st.columns(2)
    pmodel = c1.radio("Model", ["sigmoidal (one pKa)", "bell-shaped (two pKa)"], horizontal=True, key="pd_model")
    if len(pdt.dropna()) >= 3:
        try:
            pdt_ = pdt.dropna(subset=["pD", "k"])
            fr = K.fit_pd_profile(pdt_["pD"].values, pdt_["k"].values, "sigmoidal" if pmodel.startswith("sig") else "bell",
                                  pdt_["± k"].values if pdt_["± k"].notna().all() and (pdt_["± k"] > 0).all() else None)
            st.markdown(" · ".join(f"**{k_} = {v_:.3g} ± {fr['ci'][k_]:.2g}**" for k_, v_ in fr["params"].items()) + f" (95 % CI), R² {fr['R2']:.4f}. "
                        "pKa in D₂O; pKa(H₂O) ≈ pKa(D₂O) − 0.4–0.5 for carboxylic acids / Zr–OH₂.")

            def draw_pd(ax, pdt_=pdt_, fr=fr):
                ax.errorbar(pdt_["pD"], pdt_["k"], yerr=pdt_["± k"] if pdt_["± k"].notna().all() else None, fmt="o", color="k", ms=4, capsize=2, lw=0.8)
                xx = np.linspace(pdt_["pD"].min() - 0.3, pdt_["pD"].max() + 0.3, 200); ax.plot(xx, fr["curve"](xx), "-", color="#D55E00", lw=1.0)
                for k_, v_ in fr["params"].items():
                    if k_.startswith("pKa"):
                        ax.axvline(v_, color="#D55E00", lw=0.6, ls="--")
                ax.set_xlabel("pD"); ax.set_ylabel(r"$k_\mathrm{obs}$")
            fig, ax = newfig(); draw_pd(ax)
            show(fig, "pD_rate_profile", pdt_, draw_pd, "pD–rate profile")
        except Exception as e:
            st.warning(f"Fit failed: {e}")


def page_icp():
    st.header("ICP-OES (PerkinElmer Optima 8300 / Syngistix) and XRF")
    kind = st.radio("Data", ["ICP-OES (Syngistix CSV)", "XRF / any table (sample, element, value, sd)"], horizontal=True, key="icp_kind")
    if kind.startswith("ICP"):
        files = files_input("Syngistix export", ["csv", "txt", "xlsx"], "ICP/*.csv", "icp", multiple=False)
        if not files:
            return
        df = R.read_icp(files[0], files[0].name)
        c1, c2, c3 = st.columns(3)
        mass = c1.number_input("Digested sample mass (mg)", 0.01, 1000.0, 10.0)
        vol = c2.number_input("Final volume (mL)", 0.1, 1000.0, 50.0)
        dil = c3.number_input("Extra dilution factor", 1.0, 1e4, 1.0)
        df["wt%"] = df.conc * dil * vol / 1000 / mass * 100
        df["sd wt%"] = df.sd * dil * vol / 1000 / mass * 100
        val, err, ylabel = "wt%", "sd wt%", "wt%"
    else:
        files = files_input("Table", ["csv", "txt", "xlsx"], "XRF/*.csv", "xrf", multiple=False)
        if not files:
            return
        df = R.read_csv_table(files[0], (r"sample", r"element"))
        cols = list(df.columns)
        c1, c2, c3, c4 = st.columns(4)
        sc = c1.selectbox("sample", cols); ec = c2.selectbox("element", cols, index=1); val = c3.selectbox("value", cols, index=2); err = c4.selectbox("error", cols, index=min(3, len(cols) - 1))
        df = df.rename(columns={sc: "sample", ec: "element"}); ylabel = val
    st.dataframe(df, width="stretch")
    els = sorted(df.element.dropna().unique())
    chosen = st.multiselect("Elements", els, default=els)
    layout = st.radio("Layout", ["one panel per element", "grouped bars"], horizontal=True)
    # molar ratios
    piv = df.pivot_table(index="sample", columns="element", values=val, aggfunc="mean")
    from mfs.chem import AW
    mol = piv.apply(lambda c: c / AW.get(c.name, np.nan))
    if {"Hf", "Zr"} <= set(mol.columns):
        piv["Hf/(Hf+Zr) mol"] = mol.Hf / (mol.Hf + mol.Zr)
    if "S" in mol.columns and ({"Hf"} <= set(mol.columns)):
        piv["S per M6 node"] = 6 * mol.S / (mol.Hf + (mol.Zr.fillna(0) if "Zr" in mol.columns else 0))
    st.dataframe(piv, width="stretch")
    ratio_vals = _icp_ratios(df, val, err, "icp" if kind.startswith("ICP") else "xrf")
    for smp, r_ in piv.iterrows():
        vals = {f"{el}_wt" if kind.startswith("ICP") else f"{el}": (None if pd.isna(v) else float(v)) for el, v in r_.items()}
        if not kind.startswith("ICP") and {"S", "Hf"} <= set(mol.columns) and pd.notna(mol.loc[smp, "S"]):
            M_ = (mol.loc[smp, "Hf"] if pd.notna(mol.loc[smp, "Hf"]) else 0) + (mol.loc[smp, "Zr"] if "Zr" in mol.columns and pd.notna(mol.loc[smp, "Zr"]) else 0)
            vals["xrf_S_per_M"] = float(mol.loc[smp, "S"] / M_) if M_ else None
        vals.update(ratio_vals.get(str(smp), {}))
        store_put(str(smp), "icp" if kind.startswith("ICP") else "xrf", vals, "ICP page" if kind.startswith("ICP") else "XRF page")
    sub = df[df.element.isin(chosen)].copy()
    sub["sample"] = sub["sample"].str.replace("ZrO2", "ZrO$_2$")
    if layout.startswith("one"):
        fns = [(lambda ax, el=el: P.draw_grouped_bars(ax, sub[sub.element == el], val, err, ylabel=f"{el} ({ylabel})", legend=False)) for el in chosen]
        fig, axs = C.compose(fns, ncols=min(3, len(fns)), width=max(WIDTH, 2.4 * min(3, len(fns))), panel_aspect=0.95)
    else:
        fig, ax = newfig(); P.draw_grouped_bars(ax, sub, val, err, ylabel=ylabel, colors_by_group=False)
    show(fig, "elemental", piv.reset_index())


def _icp_ratios(df, val, err, sec):
    """Zr:Hf:S molar ratios with propagated uncertainties; Hf fraction and S per M6 node go to the defect calculator."""
    from mfs import ratios as RT
    st.subheader("Zr : Hf : S molar ratios ± σ")
    st.caption("σ propagated (Gaussian, uncorrelated) from the per-element SD; with replicate rows the larger of the replicate SD "
               "and the pooled instrument SD is used. Hf/(Zr+Hf) and S per M₆ node are stored for the defect calculator.")
    d = df.copy(); d[val] = pd.to_numeric(d[val], errors="coerce"); d[err] = pd.to_numeric(d[err], errors="coerce").fillna(0)
    rows, out = [], {}
    for smp, g in d.dropna(subset=[val]).groupby("sample", sort=False):
        wt, sd = {}, {}
        for el, ge in g.groupby("element"):
            el = str(el).split()[0]
            if el not in RT.M:
                continue
            wt[el] = float(ge[val].mean())
            sd[el] = float(max(np.sqrt((ge[err] ** 2).mean()), ge[val].std(ddof=1) if len(ge) > 1 else 0))
        r = RT.zr_hf_s(wt, sd)
        if not r:
            continue
        row = {"sample": smp}
        for k in ("Zr:Hf", "Hf/(Zr+Hf)", "S/(Zr+Hf)", "S per M6 node"):
            if k in r:
                row[k], row[f"σ {k}"] = r[k]
        if "Zr:Hf:S (per M6)" in r:
            row["Zr:Hf:S (per M6)"] = " : ".join(f"{v:.2f}" for v in r["Zr:Hf:S (per M6)"])
        rows.append(row)
        sv = {"Hf_frac": row.get("Hf/(Zr+Hf)"), "Hf_frac_sd": row.get("σ Hf/(Zr+Hf)"), "S_per_M6": row.get("S per M6 node"),
              "S_per_M6_sd": row.get("σ S per M6 node")}
        if "S/(Zr+Hf)" in row:
            sv.update(xrf_S_per_M=row["S/(Zr+Hf)"], xrf_S_per_M_sd=row["σ S/(Zr+Hf)"])
        out[str(smp)] = {k: (None if v is None or (isinstance(v, float) and not np.isfinite(v)) else float(v)) for k, v in sv.items()}
    if not rows:
        st.info("Needs at least two of Zr, Hf, S."); return out
    tab = pd.DataFrame(rows)
    st.dataframe(tab.style.format(precision=3), hide_index=True, width="stretch")
    dl(st, "CSV ratios ± σ", tab.to_csv(index=False).encode("utf-8-sig"), f"{sec.upper()}_ratios_sigma.csv", f"{sec}_ratio_csv")
    if "S per M6 node" in tab:
        t2 = tab.dropna(subset=["S per M6 node"])

        def draw_r(ax, t2=t2):
            xs = np.arange(len(t2)); cols = S.colours_for(list(t2["sample"]))
            ax.bar(xs, t2["S per M6 node"], 0.65, yerr=t2["σ S per M6 node"], color=cols, edgecolor="k", lw=0.6, error_kw=dict(lw=0.7, capsize=2.5))
            ax.set_xticks(xs, [str(x).replace("ZrO2", "ZrO$_2$") for x in t2["sample"]], rotation=30, ha="right")
            ax.tick_params(axis="x", which="minor", bottom=False, top=False); ax.set_ylabel(r"S per M$_6$ node")
        if len(t2):
            fig, ax = newfig(); draw_r(ax)
            show(fig, f"{sec.upper()}_S_per_node", t2, draw_r, "S per M6 node")
    return out


BAR_PRESETS = ["auto", "50 nm", "100 nm", "200 nm", "250 nm", "500 nm", "1 µm", "2 µm", "5 µm", "10 µm", "20 µm", "50 µm", "100 µm", "custom"]


def _sync_len(src, dst, w, to_frac):
    v = st.session_state.get(src)
    if v is None:
        return
    st.session_state[dst] = float(v / w) if to_frac else float(v * w)


def scalebar_controls(k, pix, width_px, shared=None):
    """Scale-bar widgets for one image; returns the bar settings dict for plots.draw_scalebar."""
    b = {}
    c1, c2, c3 = st.columns([2, 2, 3])
    mode = c1.radio("Scale bar set by", ["number (value + unit)", "length (bar size)", "auto"], key=f"{k}_mode",
                    help="number: the bar length is recomputed from the pixel size so it stays physically correct. "
                         "length: you set the bar size, the number is computed from the pixel size.")
    if mode.startswith("number"):
        pre = c2.selectbox("Scale value", BAR_PRESETS, index=0, key=f"{k}_pre")
        if pre == "auto":
            b["mode"] = "auto"
        elif pre == "custom":
            cc1, cc2 = c2.columns(2)
            v = cc1.number_input("value", 0.001, 1e6, 300.0, format="%g", key=f"{k}_val")
            u = cc2.selectbox("unit", ["nm", "µm", "mm", "Å", "pm"], key=f"{k}_unit")
            b.update(mode="value", value=v, unit=u)
        else:
            v, u = pre.split()
            b.update(mode="value", value=float(v), unit=u)
    elif mode.startswith("length"):
        kf, kp = f"{k}_lfrac", f"{k}_lpx"
        st.session_state.setdefault(kf, 0.22); st.session_state.setdefault(kp, round(0.22 * width_px))
        c2.slider("Bar length (fraction of image width)", 0.02, 0.9, step=0.01, key=kf, on_change=_sync_len, args=(kf, kp, width_px, False))
        c2.number_input("Bar length (pixels)", 1.0, float(width_px), step=1.0, key=kp, on_change=_sync_len, args=(kp, kf, width_px, True))
        u = c2.selectbox("Label unit", ["auto", "nm", "µm", "mm", "Å"], key=f"{k}_lunit")
        b.update(mode="length", length_px=float(st.session_state[kp]), unit=u)
    else:
        b["mode"] = "auto"
    if shared is not None:
        b.update(shared)
    else:
        b.update(bar_style_controls(c3, k))
    g = P.scalebar_geometry(pix, width_px, b) if pix else None
    if g:
        st.caption(f"→ bar = {g[0]:.0f} px = {g[1]:.4g} nm, label '{g[2]}' (pixel size {pix:.4g} nm/px)")
    return b


def bar_style_controls(col, k):
    s = {}
    c1, c2 = col.columns(2)
    s["position"] = c1.selectbox("Position", ["lower right", "lower left", "upper right", "upper left"], key=f"{k}_pos")
    s["color"] = c2.selectbox("Colour", ["white", "black"], key=f"{k}_col")
    s["thickness"] = c1.slider("Thickness (% of height)", 0.3, 6.0, 1.8, 0.1, key=f"{k}_th") / 100
    s["font_size"] = c2.number_input("Font size (pt, 0 = figure font)", 0.0, 30.0, 0.0, 0.5, key=f"{k}_fs")
    s["box"] = c1.checkbox("Background box", False, key=f"{k}_box")
    s["show_label"] = c2.checkbox("Show number", True, key=f"{k}_lab")
    s["label"] = col.text_input("Custom label text (optional, overrides the number)", "", key=f"{k}_ltxt")
    return s


SEM_IMG_EXT = ["tif", "tiff", "png", "jpg", "jpeg", "bmp", "dm3", "dm4", "ser", "emd", "mrc"]
CAL_UNITS = ["nm", "µm", "mm", "Å", "pm"]


@st.cache_data(show_spinner="Reading the scale bar and magnification from the info strip …", max_entries=64)
def _sem_calibrate(img8, meta_pixel, meta_mag, meta_src):
    from mfs import semcal
    return semcal.calibrate(img8, meta_pixel, meta_mag, meta_src)


def _cal_preview(img8, cal, top, marks=None, x_range=None, zoom=1.0, max_w=1100):
    """RGB preview of the info-strip region with the detected bar (cyan), label (yellow), magnification (magenta)
    and manual marks (red). Returns (PIL image, top, x0, scale)."""
    from PIL import Image, ImageDraw
    h, w = img8.shape
    xa, xb = x_range or (0, w)
    reg = img8[top:, xa:xb]
    s = max(0.25, min(zoom, 8.0)) if x_range else min(1.0, max_w / max(1, xb - xa))
    im = Image.fromarray(reg).convert("RGB")
    im = im.resize((max(1, int(im.width * s)), max(1, int(im.height * s))), Image.NEAREST if s > 1 else Image.LANCZOS)
    d = ImageDraw.Draw(im)
    T = lambda x, y: ((x - xa) * s, (y - top) * s)
    lw = max(1, int(round(s)))
    b = cal.get("bar")
    if b and not marks:
        d.rectangle([*T(b["x0"], b["y0"] - 2), *T(b["x1"], b["y1"] + 2)], outline=(0, 229, 255), width=lw + 1)
    for key, col in (("label", (255, 214, 0)), ("mag_ocr", (255, 64, 200))):
        if cal.get(key) and cal[key].get("box"):
            x0, y0, x1, y1 = cal[key]["box"]
            d.rectangle([*T(x0 - 2, y0 - 2), *T(x1 + 2, y1 + 2)], outline=col, width=lw)
    if marks:
        for x in marks:
            if x is not None:
                X, _ = T(x, top); d.line([(X, 0), (X, im.height)], fill=(255, 40, 40), width=max(2, lw))
        if len([x for x in marks if x is not None]) == 2:
            yb = im.height * 0.5
            d.line([(T(marks[0], top)[0], yb), (T(marks[1], top)[0], yb)], fill=(255, 40, 40), width=max(2, lw))
    return im, s


def _fmt_mag(v):
    return f"×{v:,.0f}".replace(",", " ") if v else ""


def _auto_calibrate(hk, best, bar, lab, crop0):
    """⚡ Auto-calibrate: select the best automatic source and reset every override to the detected values."""
    ss = st.session_state
    ss[f"{hk}_src"] = best
    if bar:
        ss[f"{hk}_blen"] = float(bar["length_px"])
        ss[f"{hk}_m0"], ss[f"{hk}_m1"] = float(bar["x0"]), float(bar["x1"])
    if lab:
        for k in ("bval", "mval"):
            ss[f"{hk}_{k}"] = float(lab["value"])
        if lab["unit"] in CAL_UNITS:
            ss[f"{hk}_bunit"] = ss[f"{hk}_munit"] = lab["unit"]
    ss[f"{hk}_c"] = int(crop0)
    ss[f"{hk}_auto_done"] = True


def calibration_panel(hk, m, fname):
    """Calibration from the ORIGINAL image: metadata, the detected scale bar + OCR label, magnification, or two clicked
    bar ends. Returns (nm_per_px, crop_row, mag_overlay dict|None)."""
    from mfs import semcal
    img = m["image"]; h, w = img.shape
    meta_mag = m["meta"].get("magnification")
    cal = _sem_calibrate(img, m["pixel_nm"], meta_mag, m["meta"].get("calibration"))
    bar, lab, mo = cal.get("bar"), cal.get("label"), cal.get("mag_ocr")
    box = st.container(border=True, key=f"cal_{hk}")
    box.markdown('<div class="rt-sec">📐 Calibration from the original image</div>', unsafe_allow_html=True)
    c = box.columns([1.15, 1, 1, 1.1, 1.3])
    c[0].metric("Metadata (nm/px)", f"{m['pixel_nm']:.4g}" if m["pixel_nm"] else "none",
                help=m["meta"].get("calibration") or "No pixel size stored in the file / sidecar")
    c[1].metric("Scale bar", f"{bar['length_px']:.0f} px" if bar else "not found",
                help=(f"{bar['polarity']} bar at x {bar['x0']}–{bar['x1']}, y {bar['y0']}–{bar['y1']} (thickness {bar['thickness_px']} px)" if bar else None))
    c[2].metric("Label (OCR)", lab["shown"] if lab else ("—" if cal["ocr_available"] else "OCR n/a"),
                help=(f"read as '{lab['text']}'" if lab else None))
    bar_nmpx = (lab["nm"] / bar["length_px"]) if (bar and lab) else None
    c[3].metric("Scale bar → nm/px", f"{bar_nmpx:.4g}" if bar_nmpx else "—")
    mval = meta_mag or (mo["value"] if mo else None)
    msrc = ("file + OCR" if mo else "file") if meta_mag else ("OCR" if mo else "")
    c[4].metric(f"Magnification ({msrc})" if msrc else "Magnification", _fmt_mag(mval) or "—",
                help=(f"printed: '{mo['text']}' → {mo['value']:g}×" if mo else "Not found in the info strip")
                + (f"; metadata: {meta_mag:g}×" if meta_mag else ""))
    for lvl, msg in cal.get("notes", []):
        (box.success if lvl == "ok" else box.warning if lvl == "warn" else box.info)(msg)
    if not cal["ocr_available"]:
        box.info("OCR engine not found – the bar is still measured; type its value below (or install Tesseract).")
    # ---- choose the source
    opts = (["metadata"] if m["pixel_nm"] else []) + (["scale bar (detected)"] if bar else []) + ["scale bar (mark ends)"] \
        + ["magnification", "manual (nm/px)"]
    default = "metadata" if m["pixel_nm"] else ("scale bar (detected)" if bar and lab else ("magnification" if (meta_mag or mo) else
                                                                                           "scale bar (mark ends)"))
    strip_top = cal.get("strip_row")
    crop0 = int(m["databar_row"] or strip_top or h)
    if strip_top is None and m["databar_row"] is None and bar and bar["y0"] > 0.88 * h:
        # bar burned into the image itself (no strip, e.g. Olympus analySIS): crop just above it so the old bar and
        # its label don't appear next to the new one
        crop0 = int(max(0, min(bar["y0"], (lab["box"][1] if lab else bar["y0"])) - 4))
    ab1, ab2 = box.columns([1.2, 2.8])
    ab1.button("⚡ Auto-calibrate", key=f"{hk}_auto", type="primary", on_click=_auto_calibrate,
               args=(hk, default, bar, lab, crop0), help="Use the most reliable automatic result (metadata › detected scale bar + "
               "OCR label › magnification), reset any manual overrides and crop the instrument strip / old bar.")
    if st.session_state.get(f"{hk}_auto_done"):
        ab2.caption(f"⚡ Auto-calibrated → **{default}**. Overrides below are still possible.")
    cc = box.columns([1.2, 2.8])
    if st.session_state.get(f"{hk}_src") not in opts:
        st.session_state[f"{hk}_src"] = default
    src = cc[0].radio("Pixel size from", opts, key=f"{hk}_src",
                      help="Order of reliability: metadata › scale bar in the image › magnification. You can always override.")
    region_top = strip_top if strip_top is not None else (max(0, bar["y0"] - 80) if bar else int(h * 0.75))
    right = cc[1]
    if src == "metadata":
        pix = float(m["pixel_nm"])
        right.caption(f"**{pix:.4g} nm/px** from {m['meta'].get('calibration') or 'the file'}"
                      + (f" · scale bar says {bar_nmpx:.4g} nm/px ({100 * (bar_nmpx / pix - 1):+.1f} %)" if bar_nmpx else ""))
        im, _ = _cal_preview(img, cal, region_top)
        right.image(im, caption="Info strip – cyan: detected bar, yellow: label, magenta: magnification", width="stretch")
    elif src == "scale bar (detected)":
        d1, d2, d3 = right.columns(3)
        st.session_state.setdefault(f"{hk}_blen", float(bar["length_px"]))
        st.session_state.setdefault(f"{hk}_bval", float(lab["value"]) if lab else 1.0)
        st.session_state.setdefault(f"{hk}_bunit", lab["unit"] if lab and lab["unit"] in CAL_UNITS else "µm")
        L = d1.number_input("Bar length (px)", 1.0, float(w), step=1.0, key=f"{hk}_blen",
                            help="Measured from the bar's first to last pixel column. Override if needed.")
        v = d2.number_input("Printed value", 0.001, 1e6, format="%g", key=f"{hk}_bval")
        u = d3.selectbox("Unit", CAL_UNITS, key=f"{hk}_bunit")
        pix = semcal.nm_per_px_from_points(0, L, v, u)
        right.markdown(f"→ **{v:g} {u} / {L:.0f} px = {pix:.4g} nm/px**" + ("" if lab else "  \n*(label not read – check the value)*"))
        im, _ = _cal_preview(img, cal, region_top)
        right.image(im, caption="cyan: detected bar, yellow: label read by OCR, magenta: magnification", width="stretch")
    elif src == "scale bar (mark ends)":
        pix = _manual_bar(right, hk, img, cal, region_top, lab)
    elif src == "magnification":
        d1, d2 = right.columns(2)
        mag0 = float(meta_mag or (mo["value"] if mo else 20000.0))
        mag = d1.number_input("Magnification (×)", 1.0, 1e8, mag0, format="%g", key=f"{hk}_mag")
        ref = d2.number_input("Reference width (mm)", 1.0, 1000.0, 127.0, key=f"{hk}_ref",
                              help="Philips/FEI XL series: magnification is relative to a 127 mm wide display. Check once against the databar.")
        pix = R.pixel_from_magnification(mag, w, ref)
        right.caption(f"→ {pix:.4g} nm/px" + (f" · scale bar says {bar_nmpx:.4g} nm/px" if bar_nmpx else ""))
    else:
        pix = right.number_input("Pixel size (nm/px)", 1e-5, 1e5, float(m["pixel_nm"] or cal.get("nm_per_px") or 10.0), format="%.5g", key=f"{hk}_pix")
    # ---- confirm / crop / magnification label
    e1, e2, e3 = box.columns([1.2, 1.4, 1.4])
    ok = e1.checkbox("Calibration checked", False, key=f"{hk}_ok", help="Tick once you've compared the highlighted bar with the image.")
    e1.caption(("✅ " if ok else "⚠️ unconfirmed – ") + f"{src}: {pix:.4g} nm/px")
    st.session_state.setdefault(f"{hk}_c", crop0)
    crop = e2.number_input("Crop rows below (info strip)", 0, h, key=f"{hk}_c",
                           help="Auto-detected start of the instrument's info strip; set to the image height to keep it.")
    box.image(_full_preview(img, cal, crop), width="stretch",
              caption="Original image – cyan: detected scale bar, yellow: its label, red line: crop (everything below is removed)")
    mag_auto = meta_mag or (mo["value"] if mo else None)
    mon = e3.checkbox("Magnification label on figure", False, key=f"{hk}_magon")
    mag = None
    if mon:
        t = e3.text_input("Magnification text", _fmt_mag(mag_auto) or "×10 000", key=f"{hk}_magtxt")
        pos = e3.selectbox("Position", ["lower left", "upper left", "lower right", "upper right"], key=f"{hk}_magpos")
        mag = {"text": t, "position": pos}
    st.session_state.setdefault("_sem_cal_info", {})[fname] = {"meta": m["pixel_nm"], "bar": bar_nmpx,
                                                               "label": lab["shown"] if lab else None,
                                                               "bar_px": bar["length_px"] if bar else None}
    return pix, crop, mag, src


def _full_preview(img, cal, crop, max_w=900):
    from PIL import Image, ImageDraw
    h, w = img.shape
    s = min(1.0, max_w / w)
    im = Image.fromarray(img).convert("RGB").resize((int(w * s), int(h * s)), Image.LANCZOS)
    d = ImageDraw.Draw(im)
    b = cal.get("bar")
    if b:
        d.rectangle([b["x0"] * s - 3, b["y0"] * s - 3, b["x1"] * s + 3, b["y1"] * s + 3], outline=(0, 229, 255), width=3)
    for key, col in (("label", (255, 214, 0)), ("mag_ocr", (255, 64, 200))):
        if cal.get(key) and cal[key].get("box"):
            x0, y0, x1, y1 = cal[key]["box"]
            d.rectangle([x0 * s - 2, y0 * s - 2, x1 * s + 2, y1 * s + 2], outline=col, width=2)
    if crop and crop < h:
        d.line([(0, crop * s), (im.width, crop * s)], fill=(255, 40, 40), width=2)
    return im


def _stale_bar_notes(names):
    from mfs import semcal
    info = {n: st.session_state.get("_sem_cal_info", {}).get(n) for n in names}
    return semcal.stale_bar_notes(info)


def _manual_bar(col, hk, img, cal, top, lab):
    """Two clicks on the zoomed info strip mark the bar ends; fine-tune with the number boxes."""
    from mfs import semcal
    h, w = img.shape
    col.info("📏 **Select the scale bar length** – click the **left end**, then the **right end** of the scale bar in the "
             "zoomed strip below. A third click starts again. Use zoom + pan for pixel accuracy.")
    k0, k1, kt = f"{hk}_m0", f"{hk}_m1", f"{hk}_mt"
    bar = cal.get("bar")
    st.session_state.setdefault(k0, float(bar["x0"]) if bar else round(w * 0.6))
    st.session_state.setdefault(k1, float(bar["x1"]) if bar else round(w * 0.85))
    z1, z2 = col.columns(2)
    zoom = z1.select_slider("Zoom", [1, 2, 3, 4, 6, 8], 1, key=f"{hk}_zoom")
    span = w / zoom
    if zoom == 1:
        cx = w / 2
        z2.caption("Zoom in to pan along the strip.")
    else:
        lo, hi = int(span / 2), int(w - span / 2)
        kp = f"{hk}_pan"
        mid = (st.session_state[k0] + st.session_state[k1]) / 2
        if kp not in st.session_state or not lo <= st.session_state[kp] <= hi or st.session_state.get(f"{kp}_z") != zoom:
            st.session_state[kp] = int(min(max(mid, lo), hi))      # centre on the bar when the zoom changes
        st.session_state[f"{kp}_z"] = zoom
        cx = z2.slider("Pan (centre x, px)", lo, hi, key=kp)
    xa = int(max(0, cx - span / 2)); xb = int(min(w, xa + span))
    disp = 1000
    im, s = _cal_preview(img, cal, top, marks=[st.session_state[k0], st.session_state[k1]], x_range=(xa, xb), zoom=disp / (xb - xa))
    try:
        from streamlit_image_coordinates import streamlit_image_coordinates as sic
        with col:
            ev = sic(im, key=f"{hk}_click", width="stretch")
    except Exception as e:                       # component missing: number boxes still work
        col.image(im, width="stretch"); ev = None
        col.caption(f"(click selection unavailable: {e})")
    if ev and ev.get("unix_time") != st.session_state.get(kt):
        st.session_state[kt] = ev.get("unix_time")
        x = xa + ev["x"] * (im.width / max(ev.get("width") or im.width, 1)) / s
        n = st.session_state.get(f"{hk}_mn", 0)
        if n % 2 == 0:
            st.session_state[k0] = round(x, 1)
        else:
            st.session_state[k1] = round(x, 1)
            if st.session_state[k1] < st.session_state[k0]:
                st.session_state[k0], st.session_state[k1] = st.session_state[k1], st.session_state[k0]
        st.session_state[f"{hk}_mn"] = n + 1
        st.rerun()
    n1, n2, n3, n4 = col.columns(4)
    x0 = n1.number_input("Left end x (px)", 0.0, float(w), step=0.5, key=k0)
    x1 = n2.number_input("Right end x (px)", 0.0, float(w), step=0.5, key=k1)
    st.session_state.setdefault(f"{hk}_mval", float(lab["value"]) if lab else 1.0)
    st.session_state.setdefault(f"{hk}_munit", lab["unit"] if lab and lab["unit"] in CAL_UNITS else "µm")
    v = n3.number_input("Bar value", 0.001, 1e6, format="%g", key=f"{hk}_mval")
    if not lab:
        col.warning("The bar's label could not be read – type the value printed on the bar (default 1 µm is only a placeholder).")
    u = n4.selectbox("Unit", CAL_UNITS, key=f"{hk}_munit")
    L = abs(x1 - x0)
    pix = semcal.nm_per_px_from_points(x0, x1, v, u) if L > 0 else (cal.get("nm_per_px") or 10.0)
    col.markdown(f"→ **{v:g} {u} / {L:.1f} px = {pix:.4g} nm/px**"
                 + (f" · auto-detected bar {bar['length_px']:.0f} px" if bar else ""))
    return pix


def page_sem():
    st.header("SEM, TEM and EDX maps")
    st.caption("SEM/TEM TIFF, PNG/JPG, Gatan .dm3/.dm4, FEI .ser, Velox .emd, .mrc. Calibration comes from the ORIGINAL image: "
               "file metadata (FEI/Thermo PixelWidth, Zeiss 'Image Pixel Size', ImageJ, DM/EMD, JEOL/Hitachi .txt or TESCAN .hdr "
               "sidecar – upload it with the image), or the scale bar burned into the info strip (measured + label read by OCR), "
               "or the magnification, or two clicks on the bar ends. The info strip is detected and cropped.")
    files = files_input("SEM / TEM images (+ optional .txt / .hdr sidecars)", SEM_IMG_EXT + ["txt", "hdr"], "SEM/*.tif", "sem")
    if st.checkbox("➕ Add the built-in calibration example (no metadata, JEOL-style info strip)", False, key="sem_calex"):
        files = list(files) + [_F(p) for p in sorted(glob.glob(os.path.join(EX, "SEM_calibration", "*.tif")))]
    side = {os.path.splitext(f.name)[0].lower(): f for f in files if f.name.lower().endswith((".txt", ".hdr"))}
    files = [f for f in files if not f.name.lower().endswith((".txt", ".hdr"))]
    if not files:
        return
    same = st.checkbox("Same scale-bar style (position, colour, thickness, font, box) for all images", True, key="sem_same")
    shared = None
    if same:
        with st.container(border=True):
            st.markdown("**Scale-bar style (all images)**")
            shared = bar_style_controls(st, "sem_all")
    with st.container(border=True):
        st.markdown("**Batch options**")
        b1, b2 = st.columns(2)
        match = b1.selectbox("Matched scale bars: one bar value on every panel", ["per image"] + BAR_PRESETS[1:-1], key="sem_match",
                             help="Every panel gets a bar of the same physical length (its pixel length follows each image's calibration).")
        onecal = b2.radio("Calibration", ["per image", "apply the first image's calibration to all", "one typed value for all"],
                          key="sem_onecal", help="Same microscope, same magnification: calibrate once and apply it to the whole batch.")
        cal_all = b2.number_input("nm/px for all images", 1e-4, 1e6, 1.0, format="%g", key="sem_calall") if onecal.startswith("one") else None
    match_bar = None if match == "per image" else {"mode": "value", "value": float(match.split()[0]), "unit": match.split()[1]}
    ims, allims, pix0 = [], [], None
    _seen = set()
    for k, f in enumerate(files):
        sc = side.get(os.path.splitext(f.name)[0].lower())
        try:
            sct = None
            if sc is not None:
                raw = sc.getvalue() if hasattr(sc, "getvalue") else open(sc.path, "rb").read()
                sct = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else str(raw)
            m = R.read_sem(f, f.name, sidecar=sct) if os.path.splitext(f.name)[1].lower() not in R.EM_EXT else R.read_sem(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        hk = "sem_" + hashlib.md5(f.name.encode()).hexdigest()[:8]
        if hk in _seen:                       # same file name twice -> separate widget keys
            hk += f"_{k}"
        _seen.add(hk)
        with st.expander(f"{f.name}" + (" (TEM)" if m.get("kind") == "TEM" else "") + (f"  ·  sidecar {sc.name}" if sc is not None else ""), expanded=k == 0):
            c1, c2 = st.columns([1, 3])
            title = c1.text_input("Label (top right)", guess_name(f.name).replace("_lowmag", ""), key=f"{hk}_t")
            inc = c1.checkbox("include", True, key=f"{hk}_i")
            pix, crop, mag, src = calibration_panel(hk, m, f.name)
            if cal_all:
                pix, src = float(cal_all), "typed value for all images"
            elif onecal.startswith("apply") and pix0:
                pix, src = pix0, "first image's calibration (batch)"
            if pix0 is None and pix:
                pix0 = pix
            bar = scalebar_controls(hk, pix, m["image"].shape[1], shared)
            if match_bar:
                bar.update(match_bar)
            g = P.scalebar_geometry(pix, m["image"].shape[1], bar) if pix else None
            if g:
                st.success(f"✅ **Applied to figure: {pix:.4g} nm/px** ({src}) → new scale bar '{g[2]}' = {g[0]:.1f} px of the "
                           f"original image (width {m['image'].shape[1]} px); image cropped at row {crop}"
                           + (f"; magnification label '{mag['text']}'" if mag else "") + ".")
            e_ = {"image": m["image"], "pixel_nm": pix, "crop_row": crop, "bar": bar, "title": title, "file": f.name,
                  "mag_text": mag["text"] if mag else None, "mag_pos": mag["position"] if mag else None}
            allims.append(e_)
            if inc:
                ims.append(e_)
    if onecal != "per image" and allims:
        st.info(f"Batch calibration: {float(cal_all or pix0 or 0):.4g} nm/px applied to all {len(allims)} images"
                + (f"; matched scale bars: {match} on every panel." if match_bar else "."))
    for msg in _stale_bar_notes([f.name for f in files]):
        st.warning("⚠️ " + msg + " Metadata is used by default; choose “scale bar (detected)” in that image's panel to reproduce the printed bar instead.")
    if ims:
        sem_panel_fig(ims)
    sem_particles(allims)
    sem_edx(allims, shared)


def sem_panel_fig(ims):
    c1, c2, c3 = st.columns(3)
    ncol = c1.slider("Columns", 1, 4, min(2, len(ims)))
    lab = c2.selectbox("Panel labels", ["(a)", "a", "A", "none"])
    afs = c3.number_input("Annotation font (pt, 0 = auto per panel size)", 0.0, 14.0, 0.0, 0.5, key="sem_panel_fs",
                          help="Auto: one font for all panels, sized to the panel (4.5–9 pt) and shrunk until every scale-bar "
                               "label lies inside its own panel.")
    labels = None if lab == "(a)" else ([""] * len(ims) if lab == "none" else [chr(97 + i) if lab == "a" else chr(65 + i) for i in range(len(ims))])
    fig = P.fig_sem_panel(ims, ncol, max(WIDTH, 3.25), labels, font_size=afs or None)
    if getattr(fig, "_mfs_panel_fs", None):
        st.caption(f"Annotation font {fig._mfs_panel_fs:.1f} pt (uniform across panels, inside each panel).")
    show(fig, "SEM_panel", series=[])


def _img8(im):
    a = np.asarray(im, float)
    if a.ndim == 3:
        a = a.mean(axis=2)
    lo, hi = np.percentile(a, (0.2, 99.8))
    return (np.clip((a - lo) / ((hi - lo) or 1), 0, 1) * 255).astype(np.uint8)


PARTICLE_POLARITY = ["dark", "bright"]     # v2.9: dark particles (TEM) is the default everywhere
ALG_ADAPTIVE = "adaptive v2.10 (background flattening + particle-scale markers) – recommended"


@st.cache_data(show_spinner="Segmenting particles …", max_entries=32)
def _particles_auto(img, pix, dark, thr_pct, smooth, open_px, split, border, max_px):
    from mfs import particles as PA
    a = _img8(img)
    thr = None if thr_pct is None else float(np.percentile(a, thr_pct))
    lab, f, t = PA.segment(a, dark_particles=dark, threshold=thr, smooth=smooth, open_px=open_px, split=split,
                           clear_border=border, max_px=max_px)
    return lab, f, t, PA.measure(lab, pix, f)


@st.cache_data(show_spinner="Watershed segmentation …", max_entries=32)
def _particles_ws(img, dark, thr_pct, smooth, open_px, mdf, hf, meth, holes, border, sol, areaf, max_px):
    from mfs import particles as PA
    a = _img8(img)
    thr = None if thr_pct is None else float(np.percentile(a, thr_pct))
    lab, f, t, info = PA.segment_ws(a, dark_particles=dark, threshold=thr, smooth=smooth, open_px=open_px, min_dist_factor=mdf,
                                    h_factor=hf, method=meth, exclude_holes=holes, clear_border=border, solidity_min=sol,
                                    area_factor=areaf, max_px=max_px)
    info = {k: v for k, v in info.items() if k != "_hole_mask"}
    return lab, f, t, info


@st.cache_data(show_spinner="Adaptive segmentation …", max_entries=32)
def _particles_adaptive(img, dark, smooth, low_k, size_px, mks, sol, areaf, border, max_px, prior, pstr, pk):
    from mfs import particles as PA
    a = _img8(img)
    lab, f, t, info = PA.segment_adaptive(a, dark_particles=dark, smooth=smooth, low_k=low_k, size_px=size_px or None,
                                          marker_sigma=mks, solidity_min=sol, area_factor=areaf, clear_border=border,
                                          max_px=max_px, shape_prior=prior, prior_strength=pstr, prior_k=pk)
    info = {k: v for k, v in info.items() if k != "_contrast"}
    return lab, f, t, info


def _template_png(t, px=56):
    """Small PNG of a learned template outline r(θ) (teal fill)."""
    t = np.asarray(t, float); th = np.linspace(0, 2 * np.pi, len(t), endpoint=False)
    fig = plt.figure(figsize=(px / 100, px / 100), dpi=100); ax = fig.add_axes([0, 0, 1, 1])
    xs, ys = np.r_[t * np.cos(th), t[0]], np.r_[t * np.sin(th), 0]
    ax.fill(xs, ys, color="#14B8A6", alpha=0.3, lw=0); ax.plot(xs, ys, color="#0F766E", lw=1.2)
    m = 1.15 * float(np.max(np.abs(t))); ax.set_xlim(-m, m); ax.set_ylim(-m, m); ax.set_aspect("equal"); ax.set_axis_off()
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=100, transparent=True); plt.close(fig)
    return b.getvalue()


def shape_hint_counts(info):
    """(n similar = clean single particles with solidity ≥ 0.9, n aggregates) from a segmentation's info dict."""
    ps = [v for k, v in info.items() if isinstance(k, (int, np.integer)) and isinstance(v, dict)]
    n_sim = sum(1 for v in ps if v.get("kind") == "particle" and float(v.get("solidity", 0) or 0) >= 0.9)
    n_agg = sum(1 for v in ps if v.get("kind") == "aggregate")
    return n_sim, n_agg


def _shape_model_hint(ph, info, prior_on):
    """Fill the slot next to the Shape-model toggle: template thumbnail + what it did, or a hint when it would help."""
    pr = info.get("_prior") or {}
    n_sim, n_agg = shape_hint_counts(info)
    with ph.container():
        if prior_on and pr.get("n", 0) >= 3 and pr.get("template") is not None:
            st.session_state["_ps_tpl"] = list(np.asarray(pr["template"], float))
            c_a, c_b = st.columns([1, 6], vertical_alignment="center")
            c_a.image(_template_png(pr["template"]), width=56)
            c_b.markdown(f"Learned from **{pr['n']}** isolated particles ({pr.get('corners', '?')} corners, weight {pr.get('weight', 0):.2f}) · "
                         f"split **{pr.get('regions_split', 0)}** clump(s) into **{pr.get('new_particles', 0)}** particles (dotted).")
        elif prior_on:
            st.caption(pr.get("reason") or "Too few clean isolated particles in this image to learn a shape.")
        elif n_sim >= 8 and n_agg >= 1:
            c_a, c_b = st.columns([1, 6], vertical_alignment="center")
            if st.session_state.get("_ps_tpl"):
                c_a.image(_template_png(st.session_state["_ps_tpl"]), width=56)
            c_b.info(f"💡 **{n_sim} similar particles found** — try the **shape model** to split the {n_agg} clump(s).")
        else:
            st.caption("Shape model: learns the outline of the clean particles and uses it to split clumps "
                       "(useful when many particles look alike).")


def _prior_panel(pr, pix_nm, f):
    """Learned shape prior: template outline, N example particles, size histogram and the weight."""
    from mfs import shapeprior as SP
    if not pr:
        return
    with st.expander(f"🔷 Learned shape prior – {pr.get('n', 0)} isolated example(s), weight {pr.get('weight', 0):.2f}"
                     + (" (applied)" if pr.get("applied") else " (not applied)"), expanded=bool(pr.get("applied"))):
        if pr.get("n", 0) < 3:
            st.info(pr.get("reason", "Too few isolated single particles to learn a shape.")); return
        st.caption(f"w = n/(n+k)·consistency·strength = {pr['n']}/({pr['n']}+{pr['k']:g}) · {pr['consistency']:.2f} · "
                   f"{pr['strength']:.2f} = **{pr['weight']:.2f}**.  Outline: {pr['corners']} corners "
                   f"(mean IoU of the examples with the template {pr['iou_mean']:.2f}); aspect ratio {pr['aspect_mean']:.2f} ± "
                   f"{pr['aspect_sd']:.2f}; orientation spread ±{pr['orientation_sd_deg']:.0f}°; median ECD "
                   f"{2 * np.sqrt(pr['area_p'][1] / np.pi) * f * pix_nm:.0f} nm.  Aggregates tried {pr['regions_tried']}, "
                   f"split {pr['regions_split']} → {pr['new_particles']} particle(s) (dotted outline; click to exclude).")
        crops = pr.get("crops") or []
        nc = min(6, max(1, len(crops)))
        fig = plt.figure(figsize=(6.4, 2.4))
        gs = fig.add_gridspec(2, nc + 3)
        a0 = fig.add_subplot(gs[:, 0:2]); t = pr["template"]; th = np.linspace(0, 2 * np.pi, len(t), endpoint=False)
        a0.fill(np.r_[t * np.cos(th), t[0]], np.r_[t * np.sin(th), 0], color="#14B8A6", alpha=0.25, lw=0)
        a0.plot(np.r_[t * np.cos(th), t[0]], np.r_[t * np.sin(th), 0], color="#0F766E", lw=1)
        a0.set_aspect("equal"); a0.set_axis_off(); a0.set_title("template", fontsize=7)
        for j, cimg in enumerate(crops[: 2 * nc]):
            ax = fig.add_subplot(gs[j // nc, 2 + j % nc]); ax.imshow(-np.asarray(cimg), cmap="gray"); ax.set_axis_off()
        ah = fig.add_subplot(gs[:, nc + 2]); e = 2 * np.sqrt(np.asarray(pr["area_px"]) / np.pi) * f * pix_nm
        ah.hist(e, bins=min(12, max(4, len(e) // 3)), color="#14B8A6"); ah.set_title("ECD (nm)", fontsize=7); ah.tick_params(labelsize=6)
        fig.tight_layout(); st.pyplot(fig, width="content"); plt.close(fig)


def _particle_editor(im, img, lab, f, info):
    """Manual corrections on the label image: delete / merge (2 clicks) / split line (2 clicks) / add disc; coloured
    labels per particle. Edits are kept per image in the session and re-applied after every automatic run."""
    from mfs import particles as PA
    from PIL import Image as _I
    k = "pse_" + hashlib.md5(im["file"].encode()).hexdigest()[:8]
    st.session_state.setdefault(k, []); st.session_state.setdefault(k + "_p", None)
    PS.track_undo(k)
    e1, e2, e3 = st.columns([2, 1, 1])
    tool = e1.radio("Edit tool", ["delete", "merge", "split", "add"], horizontal=True, key=k + "_tool",
                    help="delete: click a particle · merge: click two particles · split: click two points of a cut line · "
                         "add: click the centre of a missed particle (disc of the typical size)")
    if e2.button("Undo last edit", key=k + "_undo") and st.session_state[k]:
        st.session_state[k].pop(); st.session_state[k + "_p"] = None
    if e3.button("Clear edits", key=k + "_clr"):
        st.session_state[k] = []; st.session_state[k + "_p"] = None
    lab2 = PA.apply_edits(lab, st.session_state[k])
    fig, ax = plt.subplots(figsize=(6, 6 * img.shape[0] / img.shape[1]))
    PA.overlay_labels(ax, _img8(img), lab2, f, info, numbers=True, fs=5)
    fig.subplots_adjust(0, 0, 1, 1)
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=120); plt.close(fig)
    pv = _I.open(io.BytesIO(b.getvalue())).convert("RGB")
    s_ = pv.size[0] / (lab.shape[1])                       # preview px per working-grid px
    st.caption(f"{len(st.session_state[k])} manual edit(s). Particles are numbered; aggregates have a red dashed outline."
               + (" First point set – click the second." if st.session_state[k + "_p"] else ""))
    try:
        from streamlit_image_coordinates import streamlit_image_coordinates as sic
        ev = sic(pv, key=k + "_click")
    except Exception as e:
        st.image(pv); ev = None; st.caption(f"(click editing unavailable: {e})")
    if ev and ev.get("unix_time") != st.session_state.get(k + "_t"):
        st.session_state[k + "_t"] = ev.get("unix_time")
        x, y = ev["x"] / s_, ev["y"] / s_
        if tool == "delete":
            st.session_state[k].append(("delete", (x, y)))
        elif tool == "add":
            st.session_state[k].append(("add", (x, y), float(info.get("_r_est_px", 5)) * 1.05))
        elif st.session_state[k + "_p"] is None:
            st.session_state[k + "_p"] = (x, y)
        else:
            st.session_state[k].append((tool, st.session_state[k + "_p"], (x, y))); st.session_state[k + "_p"] = None
        st.rerun()
    return lab2


def _size_summary(d, PA):
    try:
        fit = PA.lognormal_fit(d)
    except Exception:
        return None, {}
    return fit, {"n": fit["n"], "mean (nm)": fit["mean"], "SD (nm)": fit["sd"], "median (nm)": fit["median"],
                 "log-normal median (nm)": fit["geo_mean"], "σg": fit["geo_sd"], "mode (nm)": fit["mode"], "KS p": fit["ks_p"]}


def SS_get(sample, section):
    from mfs import samplestore as SS
    try:
        r = SS.get(st.session_state, sample, section)
        return r if isinstance(r, dict) else None
    except Exception:
        return None


@st.fragment
def sem_particles(allims):
    """Particle sizing: automatic threshold + split of touching particles, or click-to-measure; histogram + log-normal fit."""
    from mfs import particles as PA
    from mfs import particle_picker as PK
    st.subheader("Particle sizing")
    cal = [im for im in allims if im.get("pixel_nm")]
    if not cal:
        st.info("Calibrate an image first (nm/px) to size particles."); return
    names = [im["file"] for im in cal]
    c1, c2 = st.columns([2, 2])
    _cur = st.session_state.get("ps_imgs")
    if _cur is not None and [n_ for n_ in _cur if n_ in names] != list(_cur):   # v2.10: new upload replaced the images
        st.session_state["ps_imgs"] = [n_ for n_ in _cur if n_ in names] or names[:1]
    sel = c1.multiselect("Images to size (same sample → one combined histogram)", names, default=names[:1], key="ps_imgs")
    meth = c2.radio("Method", ["automatic (threshold + split touching)", "click-to-measure (two clicks per particle)"], key="ps_meth")
    if not sel:
        return
    ims = [cal[names.index(n_)] for n_ in sel]
    metric = st.selectbox("Size reported", ["ECD (nm)", "Feret max (nm)", "Feret min (nm)"], key="ps_metric",
                          help="ECD = equivalent circle diameter from the area; Feret = caliper diameters (36 directions).")
    rows = []
    if meth.startswith("automatic"):
        # v2.11: the shape model is a prominent toggle next to the detection (one line of explanation, learned template
        # thumbnail, contextual hint); its strength / k stay in the adaptive settings below.
        _ad_now = str(st.session_state.get("ps_algo", ALG_ADAPTIVE)).startswith("adaptive")
        _sbox = st.container(border=True, key="ps_shape_box")
        d1, d2, d3 = _sbox.columns([1.25, 1.35, 3.4], vertical_alignment="center")
        d1.markdown("**🔍 Detect particles**<br><span style='font-size:0.8em;color:#64748B'>runs automatically on every change</span>",
                    unsafe_allow_html=True)
        d2.toggle("🔷 Shape model", False, key="ps_prior", disabled=not _ad_now,
                  help="Learns the typical outline and size from the clean, isolated particles of this image and uses it to "
                       "split clumps into single particles. Splits made by the model get a dotted outline and can be clicked away.")
        _shape_hint = d3.empty()
        k1, k2, k3, k4 = st.columns(4)
        dark = k1.radio("Particles are", PARTICLE_POLARITY, index=0, horizontal=True, key="ps_dark",
                        format_func=lambda v: {"dark": "dark (on a light background, TEM)", "bright": "bright (on a dark background, SEM)"}[v],
                        help="Default: dark particles on a lighter background (TEM). The bundled SEM example images have bright particles.") == "dark"
        tm = k1.radio("Threshold", ["Otsu (auto)", "manual"], horizontal=True, key="ps_tm")
        tp = k1.slider("Threshold (intensity percentile)", 1.0, 99.0, 50.0, 0.5, key="ps_tp") if tm == "manual" else None
        smooth = k2.slider("Smoothing σ (px)", 0.0, 4.0, 1.0, 0.25, key="ps_sm")
        opx = k2.slider("Opening (px)", 0, 5, 1, key="ps_op")
        split = k3.checkbox("Split touching particles (distance-transform markers)", True, key="ps_split")
        border = k3.checkbox("Exclude particles cut by the image edge", True, key="ps_border")
        mres = k3.select_slider("Working resolution (px)", [512, 768, 1024, 1536, 2048], 1536, key="ps_res")
        dmin = k4.number_input("Min. size (nm)", 0.0, 1e5, 0.0, format="%g", key="ps_dmin")
        dmax = k4.number_input("Max. size (nm, 0 = none)", 0.0, 1e6, 0.0, format="%g", key="ps_dmax")
        algo = st.radio("Splitting", [ALG_ADAPTIVE, "watershed v2.8 (h-maxima markers, aggregates, support holes)",
                                       "legacy v2.6 (distance markers + nearest-marker split)"], key="ps_algo")
        ad = algo.startswith("adaptive")
        ws = algo.startswith("watershed") or ad
        if ad:
            w1, w2, w3, w4 = st.columns(4)
            lowk = w1.slider("Sensitivity: lower threshold (× background noise σ)", 2.0, 8.0, 4.0, 0.5, key="ps_lowk",
                             help="After background flattening, pixels darker than the background by this many noise σ grow "
                                  "out from the clear particles (hysteresis): lower = also fainter / lighter-grey particles.")
            expnm = w1.number_input("Expected particle width (nm, 0 = auto)", 0.0, 1e5, 0.0, format="%g", key="ps_expnm",
                                    help="Auto: estimated from the image. Sets the marker scale and the aggregate size limit.")
            mks = w2.slider("Marker scale (× particle width)", 0.15, 0.4, 0.25, 0.01, key="ps_mks",
                            help="σ of the blob filter that finds one marker per particle inside clusters.")
            prior = bool(st.session_state.get("ps_prior", False))
            pstr = w2.slider("Prior strength", 0.0, 1.0, 0.7, 0.05, key="ps_pstr", disabled=not prior)
            pk_ = w2.number_input("Prior k (examples for half weight)", 1.0, 500.0, 10.0, 1.0, key="ps_pk", disabled=not prior)
            sol = w3.slider("Solidity limit", 0.6, 0.98, 0.80, 0.01, key="ps_sol2",
                            help="Area / convex-hull area below this → aggregate.")
            areaf = w3.slider("Aggregate size limit (× particle area)", 1.3, 5.0, 2.0, 0.1, key="ps_areaf2")
            pmode = w4.radio("On the image", ["include / exclude (click)", "edit outlines (merge / split / delete / add)"],
                             key="ps_pmode", help="Click mode: one click toggles a particle in/out of the histogram, drag a box "
                             "to exclude many (Shift = include). Aggregates start excluded. Edit mode changes the outlines themselves.")
            nums = w4.checkbox("Show particle numbers", False, key="ps_nums", help="Off = clean overlay; hover a particle for its number.")
            edit = pmode.startswith("edit")
            mdf = hf = meth = holes = None
            smp_guess = st.session_state.get("ps_sample") or ims[0].get("title") or guess_name(ims[0]["file"])
            _sem_saved = (SS_get(smp_guess, "sem") or {})
        elif ws:
            w1, w2, w3, w4 = st.columns(4)
            mdf = w1.slider("Min. marker distance (× typical radius)", 0.3, 1.5, 0.7, 0.05, key="ps_mdf",
                            help="The typical radius is estimated automatically from the distance map; lower = more splitting.")
            hf = w1.slider("Marker depth h (× radius)", 0.02, 0.5, 0.10, 0.01, key="ps_hf", help="h-maxima: a neck must be this much "
                           "shallower than both particles to split there; lower = more splitting.")
            meth = w2.selectbox("Flood on", ["distance", "combined", "gradient"], key="ps_wsm",
                                format_func={"distance": "distance map (shape)", "combined": "distance + image gradient",
                                             "gradient": "image gradient (facet edges)"}.get)
            holes = w2.checkbox("Exclude support holes (bright round areas)", True, key="ps_holes", disabled=not dark)
            sol = w3.slider("Solidity limit", 0.6, 0.98, 0.85, 0.01, key="ps_sol", help="Area / convex-hull area below this → "
                            "re-split; if it cannot be split it is classed as an aggregate.")
            areaf = w3.slider("Oversize limit (× median area)", 1.3, 5.0, 2.0, 0.1, key="ps_areaf")
            pmode = w4.radio("On the image", ["include / exclude (click)", "edit outlines (merge / split / delete / add)"],
                             key="ps_pmode", help="Click mode: one click toggles a particle in/out of the histogram, drag a box "
                             "to exclude many (Shift = include). Aggregates start excluded. Edit mode changes the outlines themselves.")
            nums = w4.checkbox("Show particle numbers", False, key="ps_nums", help="Off = clean overlay; hover a particle for its number.")
            edit = pmode.startswith("edit")
            smp_guess = st.session_state.get("ps_sample") or ims[0].get("title") or guess_name(ims[0]["file"])
            _sem_saved = (SS_get(smp_guess, "sem") or {})
        segs = []
        for im in ims:
            img = im["image"][:im.get("crop_row") or None]
            if ws:
                if ad:
                    lab, f, t, info = _particles_adaptive(img, dark, smooth, lowk, (expnm / float(im["pixel_nm"])) if expnm else 0.0,
                                                          mks, sol, areaf, border, mres, prior, pstr, pk_)
                    if im is ims[0]:
                        _prior_panel(info.get("_prior"), float(im["pixel_nm"]), f)
                        _shape_model_hint(_shape_hint, info, prior)
                    algo_par = (lowk, expnm, mks, prior, pstr, pk_)
                else:
                    lab, f, t, info = _particles_ws(img, dark, tp, smooth, opx, mdf, hf, meth, holes, border, sol, areaf, mres)
                    algo_par = ()
                if edit and im is ims[0]:
                    lab = _particle_editor(im, img, lab, f, info)
                r_ = PA.measure(lab, float(im["pixel_nm"]), f)
                for x in r_:
                    x["kind"] = info.get(x["id"], {}).get("kind", "particle") if isinstance(info.get(x["id"]), dict) else "particle"
                    x["solidity"] = info.get(x["id"], {}).get("solidity", np.nan) if isinstance(info.get(x["id"]), dict) else np.nan
                    x["prior"] = bool(isinstance(info.get(x["id"]), dict) and info[x["id"]].get("split") == "prior")
                nag = sum(1 for x in r_ if x["kind"] == "aggregate")
                sig = hashlib.md5(repr((im["file"], dark, tp, smooth, opx, mdf, hf, meth, holes, border, sol, areaf, mres, algo, algo_par,
                                        len(st.session_state.get("pse_" + hashlib.md5(im["file"].encode()).hexdigest()[:8], [])))).encode()).hexdigest()[:10]
                xk = f"psx_{sig}"
                if xk not in st.session_state:          # restore from the sample store if the segmentation is the same
                    sv = _sem_saved.get(f"excluded::{im['file']}")
                    st.session_state[xk] = [int(i) for i in sv] if (sv is not None and _sem_saved.get(f"excl_sig::{im['file']}") == sig) \
                        else PK.auto_excluded(r_)
                if not edit:
                    r_m = [dict(x, image=im["file"]) for x in r_ if x[metric] >= dmin and (not dmax or x[metric] <= dmax)]
                    src_, sc_, (pw_, ph_) = PK.preview_src(_img8(img))
                    parts = PK.particles_payload(r_m, lab, f, sc_, float(im["pixel_nm"]), metric)
                    ev = PK.picker(src_, pw_, ph_, parts, st.session_state[xk], sig, metric.split(" (")[0], "nm",
                                   snap=max(20, 3 * float(np.median([p_["r"] for p_ in parts])) if parts else 40), component_key="pk_" + sig,
                                   numbers=nums, rev=st.session_state.get("_ps_rev", 0))
                    if ev and ev.get("key") == sig and ev.get("t") != st.session_state.get(xk + "_evt") \
                            and ev.get("excluded") != st.session_state[xk]:
                        st.session_state[xk] = [int(i) for i in ev["excluded"]]
                    if ev:                                # v2.11: each browser change is applied once (page undo can revert it)
                        st.session_state[xk + "_evt"] = ev.get("t")
                PS.track_undo(xk)
                excl_now = st.session_state[xk]
                st.session_state.setdefault("_ps_excl", {})[im["file"]] = (sig, excl_now); PS.track("_ps_excl")
                st.session_state.setdefault("_ps_ids", {})[im["file"]] = [int(x["id"]) for x in r_
                                                                         if x[metric] >= dmin and (not dmax or x[metric] <= dmax)]
                r_ = PK.apply(r_, excl_now)
            else:
                lab, f, t, r_ = _particles_auto(img, float(im["pixel_nm"]), dark, tp, smooth, opx, split, border, mres)
                info, nag = {}, 0
            r_ = [dict(x, image=im["file"]) for x in r_ if x[metric] >= dmin and (not dmax or x[metric] <= dmax)]
            keep = np.isin(lab, [x["id"] for x in r_])
            segs.append((im, img, np.where(keep, lab, 0), f, info)); rows += r_
            st.caption(f"{im['file']}: threshold {t:.0f} (8-bit), working grid ÷{f}, {len(r_)} particles kept"
                       + (f"; typical radius {info.get('_r_est_px', 0) * f * float(im['pixel_nm']):.0f} nm, {info.get('_holes', 0)} support "
                          f"hole(s) excluded, {nag} aggregate(s) flagged, {len(excl_now)} particle(s) excluded by click/auto" if ws else "")
                       + (f", {sum(1 for x in r_ if x.get('prior'))} split by the shape prior" if ws and ad and prior else "") + ".")
    else:
        rows = _click_measure(ims[0]); segs = []
    if not rows:
        st.warning("No particles measured yet."); return
    df = pd.DataFrame(rows)
    d = df[metric].values
    fit, summ = _size_summary(d, PA)
    if summ:
        st.dataframe(pd.DataFrame([summ]).style.format(precision=3), hide_index=True, width="stretch")
        st.markdown(f"**{metric.split(' (')[0]} = {summ['mean (nm)']:.0f} ± {summ['SD (nm)']:.0f} nm** (mean ± SD, n = {summ['n']}); "
                    f"log-normal median {summ['log-normal median (nm)']:.0f} nm, σg = {summ['σg']:.2f}.")
        smp = st.text_input("Store as sample", ims[0].get("title") or guess_name(ims[0]["file"]), key="ps_sample")
        _ex = {}
        for fn_, (sg_, ex_) in st.session_state.get("_ps_excl", {}).items():
            if fn_ in sel:
                _ex[f"excluded::{fn_}"] = list(ex_); _ex[f"excl_sig::{fn_}"] = sg_
        store_put(smp, "sem", {"metric": metric, **{k: float(v) for k, v in summ.items()}, **_ex}, "SEM particle sizing")
    nb = st.slider("Histogram bins (0 = auto)", 0, 60, 0, key="ps_bins")
    col = st.color_picker("Histogram colour", "#0072B2", key="ps_col")

    def draw_h(ax, d=d, fit=fit):
        PA.draw_hist(ax, d, fit, bins=nb or "auto", color=col)
        ax.set_xlabel(f"{metric.split(' (')[0]} (nm)"); ax.set_xlim(left=0)
    if segs:
        S.apply_style(FONT)
        im0, img0, lab0, f0, info0 = segs[0]
        ar = img0.shape[0] / img0.shape[1]
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(max(WIDTH, 6.5), max(WIDTH, 6.5) / 2 * max(ar, 0.75)),
                                     gridspec_kw={"width_ratios": [1, 1]})
        if info0:
            PA.overlay_labels(a1, _img8(img0), lab0, f0, info0, numbers=bool(st.session_state.get("ps_nums")), fs=3.5)
        else:
            PA.overlay_outlines(a1, _img8(img0), lab0, f0)
        P.draw_scalebar(a1, *img0.shape[:2], im0["pixel_nm"], {"mode": "auto"})
        a1.set_title(f"{im0['file']} – {sum(r_['image'] == im0['file'] for r_ in rows)} particles", fontsize=7)
        draw_h(a2); a2.set_box_aspect(max(ar, 0.75))
        show(fig, "SEM_particle_outlines_histogram", df, series=[])
    fig, ax = newfig(); draw_h(ax)
    show(fig, "particle_size_histogram", df, draw_h, "Particle size histogram")


def _click_measure(im):
    """Two clicks per particle (diameter ends) on a down-scaled copy; returns rows with the size in nm."""
    from PIL import Image as _I, ImageDraw as _D
    img = _img8(im["image"][:im.get("crop_row") or None]); h, w = img.shape
    k = "psc_" + hashlib.md5(im["file"].encode()).hexdigest()[:8]
    st.session_state.setdefault(k, []); st.session_state.setdefault(k + "_p", None)
    PS.track_undo(k)
    s = min(1.0, 1000 / w)
    pv = _I.fromarray(img).convert("RGB").resize((int(w * s), int(h * s)))
    dr = _D.Draw(pv)
    for (x0, y0, x1, y1) in st.session_state[k]:
        dr.line([(x0 * s, y0 * s), (x1 * s, y1 * s)], fill=(255, 59, 48), width=2)
    if st.session_state[k + "_p"]:
        x0, y0 = st.session_state[k + "_p"]; dr.ellipse([x0 * s - 3, y0 * s - 3, x0 * s + 3, y0 * s + 3], outline=(255, 214, 10), width=2)
    st.caption("Click both ends of a particle's diameter; each pair adds one particle. Rows can be deleted in the table.")
    try:
        from streamlit_image_coordinates import streamlit_image_coordinates as sic
        ev = sic(pv, key=k + "_click")
    except Exception as e:
        st.image(pv); ev = None; st.caption(f"(click selection unavailable: {e})")
    if ev and ev.get("unix_time") != st.session_state.get(k + "_t"):
        st.session_state[k + "_t"] = ev.get("unix_time")
        x, y = ev["x"] / s, ev["y"] / s
        if st.session_state[k + "_p"] is None:
            st.session_state[k + "_p"] = (x, y)
        else:
            x0, y0 = st.session_state[k + "_p"]; st.session_state[k].append((x0, y0, x, y)); st.session_state[k + "_p"] = None
        st.rerun()
    if st.button("Clear clicks", key=k + "_clr"):
        st.session_state[k] = []; st.session_state[k + "_p"] = None; st.rerun()
    pix = float(im["pixel_nm"])
    tab = pd.DataFrame([{"x0": a, "y0": b, "x1": c, "y1": e, "d (nm)": float(np.hypot(c - a, e - b) * pix)} for a, b, c, e in st.session_state[k]],
                       columns=["x0", "y0", "x1", "y1", "d (nm)"])
    ed = st.data_editor(tab, num_rows="dynamic", key=k + "_ed", hide_index=True, disabled=["d (nm)"])
    ed = ed.dropna()
    return [{"id": i + 1, "image": im["file"], "ECD (nm)": float(np.hypot(r.x1 - r.x0, r.y1 - r.y0) * pix),
             "Feret max (nm)": float(np.hypot(r.x1 - r.x0, r.y1 - r.y0) * pix), "Feret min (nm)": float(np.hypot(r.x1 - r.x0, r.y1 - r.y0) * pix)}
            for i, r in enumerate(ed.itertuples())]


EDX_PALETTE = {"Hf": "#56B4E9", "Zr": "#E69F00", "S": "#F0E442", "O": "#D55E00", "C": "#999999", "Pt": "#CC79A7", "N": "#009E73", "Ti": "#0072B2"}


def sem_edx(allims, shared):
    from PIL import Image
    from mfs import edx as EX_
    from scipy import ndimage as ndi
    st.subheader("EDX element maps")
    ef = files_input("EDX maps (greyscale PNG/TIFF, one per element; optional SE/BSE image named *_SE*)", ["png", "tif", "tiff", "jpg"],
                     "EDX/*.png", "edx")
    if not ef:
        return
    se = [f for f in ef if re.search(r"_(SE|BSE|SEM)\.\w+$", f.name, re.I)]
    ef = [f for f in ef if f not in se]
    if not ef:
        return
    maps = []
    for f in ef:
        el = re.search(r"_(\w{1,2})\.\w+$", f.name); el = el.group(1) if el else f.name
        maps.append((el, np.array(Image.open(f).convert("L"))))
    c1, c2, c3 = st.columns(3)
    epix = c1.number_input("EDX map pixel size (nm/px, 0 = no scale bar)", 0.0, 1e5, 0.0, format="%g", key="edx_pix")
    mode = c2.radio("Show", ["separate maps", "overlay (one image)", "both"], index=2, key="edx_mode")
    bases = ["none"] + [f.name for f in se] + [im["file"] for im in allims]
    bsel = c3.selectbox("Overlay base image", bases, index=1 if se else 0, key="edx_base")
    ba = c3.slider("Base image weight", 0.0, 1.0, 0.6, 0.05, key="edx_ba")
    esm = c1.slider("Smooth maps (σ, px; overlay and line scan)", 0.0, 3.0, 1.0, 0.25, key="edx_sm",
                    help="EDX maps are count-noisy; a small Gaussian makes overlays and profiles readable. 0 = raw counts.")
    st.markdown("**Element colours**")
    cc = st.columns(max(1, len(maps)))
    colors, gains = {}, {}
    for i, (el, _) in enumerate(maps):
        colors[el] = cc[i].color_picker(el, EDX_PALETTE.get(el, "#009E73"), key=f"edx_c_{el}")
        gains[el] = cc[i].slider(f"{el} gain", 0.2, 3.0, 1.0, 0.1, key=f"edx_g_{el}")
    ebar = scalebar_controls("edx", epix, maps[0][1].shape[1], shared) if epix else None
    if mode != "overlay (one image)":
        fns = [(lambda ax, el=el, im=im: P.draw_edx_map(ax, im, colors[el], pixel_nm=epix or None, element=el, bar=ebar)) for el, im in maps]
        fig, axs = C.compose(fns, ncols=len(fns), width=max(WIDTH, 1.6 * len(fns)), panel_aspect=0.8, wspace=0.05, label_x=0.02, label_y=1.02)
        show(fig, "EDX_maps", series=[])
    base = None
    if bsel != "none":
        if bsel in [f.name for f in se]:
            base = np.array(Image.open(se[[f.name for f in se].index(bsel)]).convert("L"))
        else:
            b_ = allims[[im["file"] for im in allims].index(bsel)]; base = b_["image"][:b_.get("crop_row") or None]
    rgb = EX_.overlay(dict(maps), colors, base=base, base_alpha=ba, gains=gains, smooth=esm)

    def draw_ov(ax, rgb=rgb):
        from matplotlib.patches import Patch
        ax.imshow(rgb, interpolation="nearest"); ax.set_xticks([]); ax.set_yticks([])
        if epix:
            P.draw_scalebar(ax, *rgb.shape[:2], epix, ebar)
        lg = ax.legend(handles=[Patch(color=colors[el], label=el) for el, _ in maps], loc="upper left", fontsize=6, frameon=True,
                       facecolor="black", edgecolor="none", framealpha=0.45, handlelength=0.9, borderaxespad=0.3)
        for t_ in lg.get_texts():
            t_.set_color("white")
    if mode != "separate maps":
        fig, ax = newfig(maps[0][1].shape[0] / maps[0][1].shape[1]); draw_ov(ax)
        show(fig, "EDX_overlay", None, draw_ov, "EDX overlay")
    if st.checkbox("Line scan across the maps", False, key="edx_ls"):
        h, w = maps[0][1].shape
        l1, l2, l3 = st.columns(3)
        y0 = l1.slider("Start y (px)", 0, h - 1, h // 2, key="edx_y0"); x0 = l1.slider("Start x (px)", 0, w - 1, w // 8, key="edx_x0")
        y1 = l2.slider("End y (px)", 0, h - 1, h // 2, key="edx_y1"); x1 = l2.slider("End x (px)", 0, w - 1, 7 * w // 8, key="edx_x1")
        bw = l3.slider("Band width (px)", 1, 25, 5, key="edx_bw")
        norm = l3.checkbox("Normalise each profile to its max", True, key="edx_lsn")
        prof = {}
        for el, im in maps:
            dd, pp = EX_.line_scan(ndi.gaussian_filter(im.astype(float), esm) if esm else im, (x0, y0), (x1, y1), width=bw)
            prof[el] = pp / (pp.max() or 1) if norm else pp
        dist = dd * (epix if epix else 1); du = "nm" if epix else "px"
        if epix and dist.max() > 2000:
            dist, du = dist / 1000, "µm"
        tab = pd.DataFrame({f"distance ({du})": dist, **{el: v for el, v in prof.items()}})

        def draw_ls(ax, tab=tab):
            for el, _ in maps:
                ax.plot(tab.iloc[:, 0], tab[el], color=colors[el] if colors[el].lower() not in ("#ffffff", "#f0e442") else "#B8A000", lw=1, label=el)
            ax.set_xlabel(f"Distance ({du})"); ax.set_ylabel("Normalised counts" if norm else "Counts"); ax.legend(frameon=False, fontsize=6)
        S.apply_style(FONT)
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(max(WIDTH, 6.5), max(WIDTH, 6.5) * 0.4))
        a1.imshow(rgb, interpolation="nearest"); a1.set_axis_off()
        a1.plot([x0, x1], [y0, y1], color="w", lw=1.2); a1.plot([x0], [y0], "o", color="w", ms=3)
        draw_ls(a2)
        show(fig, "EDX_linescan", tab, series=[])


def page_sds():
    st.header("SDS-PAGE – Bio-Rad Gel Doc EZ / Image Lab")
    st.caption("Image Lab: Analysis Table → 'Lane and Bands' → Export → CSV/Excel. Add a time column (or map lane → time below).")
    files = files_input("Image Lab export", ["csv", "txt", "xlsx", "xls"], "SDS-PAGE/*.csv", "sds", multiple=False)
    if not files:
        return
    b = R.read_imagelab(files[0], files[0].name)
    if "time" not in b or b["time"].isna().all():
        lanes = sorted(b.lane.unique())
        lt = st.data_editor(pd.DataFrame({"lane": lanes, "time": [float(i) for i in range(len(lanes))]}), key="lt")
        b = b.merge(lt, on="lane")
    renorm = st.checkbox("Re-normalise band % to 100 % per lane", True)
    if renorm:
        b["band_pct"] = 100 * b.band_pct / b.groupby("lane").band_pct.transform("sum")
    st.dataframe(b, width="stretch")
    tunit = st.selectbox("time unit", ["h", "min", "d"])
    def draw(ax, b=b):
        P.draw_bands(ax, b, time_label=f"Time ({tunit})")
    fig, ax = newfig(); draw(ax)
    show(fig, "SDS-PAGE_bands", b, draw, "SDS-PAGE", units={"x": "time_h"})
    st.caption("Tip: fit the intact-protein band with the Kinetics tab (first-order) to report k_obs for HEWL.")


def page_zeta():
    st.header("Zeta potential / DLS – Malvern Zetasizer")
    st.caption("Zetasizer software: File → Export (choose parameters Sample Name, Z-Ave, PdI, Zeta Potential, Zeta Deviation and Sizes[]/Intensities[]) → CSV.")
    files = files_input("Zetasizer CSV", ["csv", "txt"], "Zeta-DLS/*.csv", "zeta", multiple=False)
    if not files:
        return
    z = dict(R.read_zetasizer(files[0], files[0].name))
    _sel = sample_names_list(list(pd.unique(z["summary"]["sample"].astype(str))), "zeta_samples")
    z["summary"] = _apply_sample_list(z["summary"], "sample", None, _sel)
    if z.get("dist") is not None:
        z["dist"] = _apply_sample_list(z["dist"], "sample", None, _sel)
    zs = z["summary"].groupby("sample", sort=False).agg(zeta=("zeta", "mean"), zeta_sd=("zeta", "std"), z_avg=("z_avg", "mean"), pdi=("pdi", "mean")).reset_index()
    st.dataframe(zs, width="stretch")
    zs["label"] = zs["sample"].str.replace("ZrO2", "ZrO$_2$")
    def draw_zeta(ax, zs=zs):
        cols = S.colours_for(zs["label"]); xs = np.arange(len(zs))
        ax.bar(xs, zs.zeta, 0.65, yerr=zs.zeta_sd, color=cols, edgecolor="k", lw=0.6, error_kw=dict(lw=0.7, capsize=2.5))
        ax.axhline(0, color="k", lw=0.6); ax.set_xticks(xs, zs.label, rotation=30, ha="right")
        ax.tick_params(axis="x", which="minor", bottom=False, top=False); ax.set_ylabel(r"$\zeta$-potential (mV)")
    fig, ax = newfig(); draw_zeta(ax)
    show(fig, "zeta", zs, draw_zeta, "Zeta potential")
    if z["dist"] is not None:
        dd = z["dist"].groupby(["sample", "size_nm"], sort=False).intensity.mean().reset_index()
        tr = [(s_.replace("ZrO2", "ZrO$_2$"), g_.size_nm.values, g_.intensity.values) for s_, g_ in dd.groupby("sample", sort=False)]
        def draw_dls(ax, tr=tr):
            P.draw_overlay(ax, tr, "Hydrodynamic diameter (nm)", "Intensity (%)", logx=True); ax.set_ylim(bottom=0)
        fig, ax = newfig(); draw_dls(ax)
        show(fig, "DLS", dd, draw_dls, "DLS")


def page_xps():
    from mfs import xps as XP
    st.header("XPS")
    st.caption("CasaXPS / Avantage exported region (BE, counts) as CSV/TXT; simple VAMAS (.vms) regions. Shirley or Tougaard "
               "background, doublet fitting (Zr 3d, Hf 4f, S 2p with fixed spin–orbit splitting and area ratio), atomic-% with RSFs.")
    multi = st.checkbox("➕ Use the multi-region EXAMPLE instead (ZrO₂@Sulfo-MOF-801: Zr 3d, Hf 4f, S 2p, O 1s, C 1s – synthetic)", False, key="xps_multi")
    files = files_input("Spectra", ["csv", "txt", "vms", "dat"], "XPS/*Hf4f*.csv", "xps")
    if multi:
        files = [_F(p_) for p_ in sorted(glob.glob(os.path.join(EX, "XPS", "ZrO2@Sulfo-MOF801_*.csv")))]
    if not files:
        return
    samples = sample_table(files, "xps_tab" + ("_m" if multi else ""))
    spec = []
    for nm, f in samples:
        try:
            d = R.read_xps(f, f.name)
        except Exception as e:
            st.error(f"{f.name}: {e}"); continue
        reg = XP.guess_region(d.x.values)
        grp = re.sub(r"\s*\(\d+\)$", "", nm)
        grp = re.sub(r"[ _-]*(Zr ?3d|Hf ?4[fd]|S ?2p|O ?1s|C ?1s|N ?1s|F ?1s|Ti ?2p)\b", "", grp, flags=re.I).strip(" _-") or nm
        spec.append({"name": nm, "group": grp, "region": reg, "x": d.x.values, "y": d.y.values})
    st.dataframe(pd.DataFrame([{"spectrum": d_["name"], "sample": d_["group"], "region (from BE range)": d_["region"],
                                "BE range (eV)": f"{d_['x'].min():.1f}–{d_['x'].max():.1f}"} for d_ in spec]), hide_index=True, width="stretch")
    regions = sorted({d_["region"] for d_ in spec if d_["region"]})
    c1, c2, c3 = st.columns(3)
    show_reg = c1.selectbox("Overview: region", regions or ["(all)"], key="xps_ovreg")
    bgk = c2.radio("Background", ["Shirley", "Tougaard", "linear", "none"], horizontal=True, key="xps_bg")
    region = c3.text_input("Region label", show_reg if regions else "Hf 4f", key="xps_lab_" + str(show_reg))
    tr = [(d_["name"], d_["x"], d_["y"]) for d_ in spec if not regions or d_["region"] == show_reg]
    lo_, hi_ = (min(t[1].min() for t in tr), max(t[1].max() for t in tr)) if tr else (13.0, 23.0)
    lo, hi = st.slider("Binding energy (eV)", 0.0, 1200.0, (float(np.floor(lo_)), float(np.ceil(hi_))), 0.5, key="xps_rng_" + str(show_reg))

    def draw(ax, tr=tr):
        P.draw_xps(ax, tr, (lo, hi), bgk == "Shirley"); ax.text(0.03, 0.96, region, transform=ax.transAxes, va="top")
    fig, ax = newfig(); draw(ax)
    show(fig, "XPS", None, draw, "XPS")

    st.subheader("Peak fitting")
    names = [d_["name"] for d_ in spec]
    c1, c2, c3 = st.columns(3)
    pick = c1.selectbox("Spectrum", names, key="xps_fit_pick", format_func=lambda n_: re.sub(r"\$_(\d)\$", lambda m_: chr(0x2080 + int(m_.group(1))), n_))
    d_ = spec[names.index(pick)]
    regs = list(XP.DOUBLETS) + list(XP.SINGLETS) + ["other (singlets)"]
    rsel = c2.selectbox("Region / doublet constraint", regs, index=regs.index(d_["region"]) if d_["region"] in regs else len(regs) - 1,
                        key="xps_fit_reg_" + pick)
    eta = c3.slider("Lorentzian fraction (GL mix)", 0.0, 1.0, 0.3, 0.05, key="xps_eta")
    dd = XP.DOUBLETS.get(rsel)
    if dd:
        st.caption(f"{rsel}: {dd['major']} + {dd['minor']} at +{dd['split']} eV, area {dd['minor']}/{dd['major']} = {dd['ratio']:.3f}, same FWHM.")
    x, y = d_["x"], d_["y"]
    bg0 = XP.shirley(x, y)
    xpk = float(x[np.argmax(y - bg0)])
    defaults = {"Zr 3d": [("Zr–O (MOF / ZrO₂)", xpk)], "Hf 4f": [("Hf–O", xpk)], "S 2p": [("sulfonate / sulfate", xpk)],
                "O 1s": [("O lattice (M–O)", 530.2), ("O carboxylate / OH", 531.8)], "C 1s": [("C–C / C–H", 284.8), ("O–C=O", 288.6)]}
    comp0 = pd.DataFrame([{"name": n_, "BE guess (eV)": round(c_, 2), "doublet": bool(dd), "fix BE": False}
                          for n_, c_ in defaults.get(rsel, [("component 1", xpk)])])
    comps = st.data_editor(comp0, num_rows="dynamic", key="xps_comp_" + pick + rsel, hide_index=True)
    try:
        fr = XP.fit(x, y, [{"name": r_["name"], "center": float(r_["BE guess (eV)"]), "region": rsel if (r_["doublet"] and dd) else None,
                            "fix_center": bool(r_["fix BE"])} for _, r_ in comps.dropna(subset=["BE guess (eV)"]).iterrows()],
                    background=bgk.lower(), eta=eta)
        tab = pd.DataFrame(fr["rows"]); st.dataframe(tab.style.format(precision=2), hide_index=True, width="stretch")
        st.caption(f"Background: {bgk}{' (universal cross-section B·T/(C+T²)², C = 1643 eV², B scaled to the high-BE end)' if bgk == 'Tougaard' else ''}. "
                   f"Residual RMS = {100 * fr['rms_rel']:.1f} % of the peak height.")
        st.session_state.setdefault("xps_areas", {}).setdefault(d_["group"], {})[rsel] = float(tab["area (total)"].sum())

        def draw_f(ax, fr=fr):
            XP.draw_fit(ax, fr); ax.legend(frameon=False, fontsize=6); ax.text(0.03, 0.96, f"{pick}\n{rsel}", transform=ax.transAxes, va="top", fontsize=7)
        fig, ax = newfig(0.85); draw_f(ax)
        show(fig, f"XPS_fit_{safe_name(pick)}", pd.DataFrame({"BE": fr["x"], "counts": fr["y"], "background": fr["bg"], "envelope": fr["fit"]}),
             draw_f, f"XPS fit {pick}")
    except Exception as e:
        st.warning(f"Fit failed: {e}")

    st.subheader("Atomic % (RSF-corrected)")
    st.caption("Areas above the chosen background, per region, per sample (whole doublets). Default RSFs: Scofield cross-sections "
               "(Al Kα, C 1s = 1) without transmission/IMFP correction – replace them with your instrument library's RSFs.")
    rsf = st.data_editor(pd.DataFrame({"region": list(XP.RSF_SCOFIELD), "RSF": list(XP.RSF_SCOFIELD.values())}), key="xps_rsf", hide_index=True)
    rsfd = dict(zip(rsf.region, rsf.RSF))
    areas = {}
    for s_ in spec:
        if not s_["region"]:
            continue
        b_ = XP.shirley(s_["x"], s_["y"]) if bgk == "Shirley" else (XP.tougaard(s_["x"], s_["y"])[0] if bgk == "Tougaard" else
                                                                      np.interp(s_["x"], [s_["x"][0], s_["x"][-1]], [s_["y"][:3].mean(), s_["y"][-3:].mean()]))
        o = np.argsort(s_["x"]); A_ = float(np.trapezoid(np.clip((s_["y"] - b_)[o], 0, None), s_["x"][o]))
        areas.setdefault(s_["group"], {})[s_["region"]] = A_
    if areas:
        at = []
        for g_, ar_ in areas.items():
            pc = XP.atomic_percent(ar_, rsfd)
            row = {"sample": g_}; row.update({f"{k} at.%": v for k, v in pc.items()})
            M_ = {k: pc.get(k, 0) for k in ("Zr 3d", "Hf 4f")}
            if sum(M_.values()):
                row["Hf/(Zr+Hf)"] = M_["Hf 4f"] / sum(M_.values())
                if "S 2p" in pc:
                    row["S/(Zr+Hf)"] = pc["S 2p"] / sum(M_.values())
            at.append(row)
            store_put(g_, "xps", {k: v for k, v in row.items() if k != "sample"}, "XPS page")
        att = pd.DataFrame(at); st.dataframe(att.style.format(precision=2), hide_index=True, width="stretch")
        dl(st, "CSV atomic %", att.to_csv(index=False).encode("utf-8-sig"), "XPS_atomic_percent.csv", "xps_at_csv")


def page_composer():
    st.header("Figure composer")
    st.caption("Panels added from the technique tabs are drawn again into one figure (vector output), with (a), (b), (c) labels.")
    comp = st.session_state.composer
    if not comp:
        st.info("No panels yet – use '➕ Add this panel…' under any figure."); return
    order = st.data_editor(pd.DataFrame({"panel": list(comp), "use": True, "position": range(len(comp))}), key="comp_ed", hide_index=True)
    sel = [r.panel for r in order.sort_values("position").itertuples() if r.use]
    from mfs import templates as TP
    t1, t2, t3 = st.columns([1.2, 1.2, 2.6])
    tpl = t1.selectbox("Journal template", ["custom"] + list(TP.TEMPLATES), key="comp_tpl")
    colw = t2.radio("Column", ["single", "1.5", "double"], index=2, horizontal=True, key="comp_col", disabled=tpl == "custom")
    if tpl != "custom":
        t3.caption(TP.TEMPLATES[tpl]["note"])
    c1, c2, c3, c4, c5 = st.columns(5)
    ncols = c1.slider("Columns", 1, 4, 2)
    if tpl == "custom":
        width = S.WIDTHS[c2.selectbox("Width", list(S.WIDTHS), index=1)]
        lab = c4.selectbox("Labels", ["(a)", "a", "A", "a)"])
        fsz = c5.slider("Font (pt)", 5.0, 10.0, float(FONT), 0.5, key="comp_fs")
    else:
        T_ = TP.TEMPLATES[tpl]; width = TP.width(tpl, colw)
        labs = ["(a)", "a", "A", "a)"]
        lab = c4.selectbox("Labels", labs, index=labs.index(T_["label_style"]), key=f"comp_lab_{tpl}")
        fsz = c5.slider("Font (pt)", float(T_["min_font"]), 10.0, float(T_["font_size"]), 0.5, key=f"comp_fs_{tpl}")
        c2.metric("Width", f"{width:.2f} in / {width * 2.54:.1f} cm")
    asp = c3.slider("Panel aspect", 0.4, 1.5, 0.8, 0.05)
    if st.button("Clear composer"):
        st.session_state.composer = {}; st.rerun()
    fig, axs = C.compose([comp[k] for k in sel], ncols=ncols, width=width, panel_aspect=asp, label_style=lab, font_size=fsz)
    if tpl != "custom":
        for lvl, msg in TP.check(fig, tpl, colw):
            (st.success if lvl == "ok" else st.warning)(f"{tpl} check: {msg}")
    show(fig, f"composed_figure{'_' + tpl if tpl != 'custom' else ''}")


# ---------------------------------------------------------------- Batch E: Structure (CIF)
@st.cache_data(show_spinner=False, max_entries=8)
def _cif_text(raw: bytes):
    return raw.decode("utf-8", "replace")


@st.cache_data(show_spinner=False, max_entries=8)
def _cif_pxrd(text, lo, hi, fw):
    from mfs import cifpxrd as CPX
    return CPX.simulate(text.encode(), (lo, hi), fwhm=fw)


def _cif_struct(text, name):
    from pathlib import Path  # noqa
    from mfs import cif as CX
    key = ("cifS", hash(text))
    if key not in st.session_state:
        st.session_state[key] = CX.read(text.encode(), name)
    return st.session_state[key]


def page_structure():
    from pathlib import Path
    from mfs import cif as CX, cifplot as CP, cifsi as CS, cifpxrd as CPX
    st.header("Structure (CIF)")
    st.caption("Single-crystal or MOF CIF → ORTEP / ball-and-stick / stick / space-filling figures, H-bond tables with "
               "symmetry codes, SI tables (crystal data, bonds, angles, torsions with ESDs) as DOCX/XLSX, simulated PXRD "
               "with hkl labels (and overlay on a measured pattern) and, for MOFs, geometric pore analysis. Pure matplotlib.")
    refs = {p.name: p for p in sorted(Path(REF).glob("*.cif"))}
    c1, c2 = st.columns([2, 1])
    up = c1.file_uploader("CIF file", type=["cif"], key="cx_up")
    ex = c2.selectbox("…or a bundled example", ["—"] + list(refs), key="cx_ex")
    if up is not None:
        text, name = _cif_text(up.getvalue()), up.name
    elif ex != "—":
        text, name = refs[ex].read_text(errors="replace"), ex
    else:
        st.info("Upload a CIF (SHELXL/Olex2 output or a COD/CSD file) or pick an example."); return
    try:
        blocks = _cif_struct(text, Path(name).stem)
    except Exception as e:
        st.error(f"Could not read the CIF: {e}"); return
    if not blocks:
        st.error("No data block with atoms found."); return
    S = blocks[st.selectbox("Data block", range(len(blocks)), format_func=lambda i: blocks[i].name, key="cx_blk")] if len(blocks) > 1 else blocks[0]
    kn = CS.key_numbers(S)
    st.success(" · ".join(f"{k} {v}" for k, v in kn.items() if v not in ("", None)))
    is_mof = any(a["el"] in CX.METALS for a in S.atoms)
    tabs = st.tabs(["Figure", "H-bonds", "SI tables", "Simulated PXRD"] + (["Pores (MOF)"] if is_mof else []) + ["Hirshfeld (experimental)"])
    # ---- figure
    with tabs[0]:
        a1, a2, a3, a4 = st.columns(4)
        style = a1.selectbox("Style", ["ortep", "ball", "stick", "space"], key="cx_style",
                             format_func={"ortep": "ORTEP ellipsoids", "ball": "Ball-and-stick", "stick": "Stick", "space": "Space-filling"}.get)
        prob = a2.select_slider("Ellipsoid probability", [0.3, 0.5, 0.7], 0.5, key="cx_prob", format_func=lambda p: f"{p:.0%}")
        mode = a3.selectbox("Grow", ["asu", "molecules", "packing", "cell", "supercell"], index=1 if not is_mof else 3, key="cx_mode",
                            format_func={"asu": "Asymmetric unit", "molecules": "Complete molecules", "packing": "Packing (unit cell)",
                                         "cell": "Atoms in the cell", "supercell": "Supercell"}.get)
        labels = a4.selectbox("Labels", ["non-H", "all", "none"], key="cx_lab")
        b1, b2, b3, b4 = st.columns(4)
        hmode = b1.selectbox("Hydrogen atoms", ["all", "polar", "none"], key="cx_h",
                             format_func={"all": "All H", "polar": "Polar H only (O–H, N–H)", "none": "No H"}.get)
        showh = hmode != "none"
        has_dis = any(str(a["dis_group"]).strip(".?-") for a in S.atoms)
        part = b2.selectbox("Disorder", ["major", "minor", "all"] if has_dis else ["all"], key="cx_part",
                            format_func={"major": "Major part (1)", "minor": "Minor part (2)", "all": "All parts"}.get)
        cell = b3.checkbox("Cell edges", True, key="cx_cell"); triad = b3.checkbox("Axis triad", True, key="cx_triad")
        hbo = b4.checkbox("H-bonds (dashed)", True, key="cx_hb")
        poly = b4.checkbox("Metal polyhedra", is_mof, key="cx_poly", disabled=not is_mof)
        e1, e2, e3 = st.columns([1, 1, 2])
        hsc = e1.slider("H atom size", 0.3, 1.5, 0.7, 0.1, key="cx_hsc", help="Scale of the H spheres (ORTEP / ball-and-stick)")
        ttl = e2.checkbox("Title on figure", False, key="cx_ttl", help="Off for publication exports (caption goes in the manuscript)")
        offs_txt = e3.text_input("Manual label offsets (Å, screen x y)", "", key="cx_offs", placeholder="C10A 0 -0.8; O5A 0.5 0",
                                 help="Shift single labels after the automatic placement, e.g. 'C10A 0 -0.8; O5A 0.5 0'")
        offs = {}
        for chunk in offs_txt.split(";"):
            tok = chunk.split()
            if len(tok) == 3:
                try:
                    offs[tok[0]] = (float(tok[1]), float(tok[2]))
                except ValueError:
                    st.warning(f"Offset '{chunk.strip()}' not understood (label dx dy)")
        sc = (1, 1, 1)
        if mode == "supercell":
            s1, s2, s3 = st.columns(3)
            sc = (s1.number_input("na", 1, 4, 2, key="cx_na"), s2.number_input("nb", 1, 4, 1, key="cx_nb"), s3.number_input("nc", 1, 4, 1, key="cx_nc"))
        v1, v2, v3 = st.columns([1, 1, 2])
        view = v1.selectbox("View along", ["a", "b", "c", "uvw", "hkl", "best"], key="cx_view",
                            format_func={"a": "a axis", "b": "b axis", "c": "c axis", "uvw": "vector [uvw]", "hkl": "plane normal (hkl)", "best": "best view (mean plane)"}.get)
        vec = v2.text_input("u v w / h k l", "1 1 0", key="cx_vec", disabled=view not in ("uvw", "hkl"))
        with v3:
            r1, r2, r3 = st.columns(3)
            rot = (r1.slider("rot x (°)", -180, 180, 0, 5, key="cx_rx"), r2.slider("rot y (°)", -180, 180, 0, 5, key="cx_ry"),
                   r3.slider("rot z (°)", -180, 180, 0, 5, key="cx_rz"))
        try:
            vv = [float(t) for t in vec.replace(",", " ").split()][:3]
        except ValueError:
            vv = [1, 1, 0]
        gmode = "packing" if (mode == "supercell" and not is_mof) else mode   # molecular crystals: whole molecules per cell
        sites, bonds = S.grow(gmode, part=part, supercell=sc, include_h=showh)
        if len(sites) >= 6000:
            st.warning("Atom limit (6000) reached; drawing a partial structure.")
        if hmode == "polar":
            sites, bonds = CP.polar_h_only(sites, bonds)
        V = CP.view_matrix(S, view, uvw=vv if view == "uvw" else None, hkl=vv if view == "hkl" else None, rot=rot, sites=sites)
        hb = CX.hbond_pairs(S, sites) if (hbo and showh) else []
        n_mol = len(CX._components(sites, bonds)) if sites else 0
        st.caption(f"{len(sites)} atoms, {len(bonds)} bonds, {n_mol} fragment(s), {len(hb)} H-bond(s) drawn"
                   + (" · polymeric framework (molecules = connected network)" if getattr(S, "polymeric", False) else ""))
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        CP.draw(ax, S, sites, bonds, style=style, prob=prob, labels=labels, hbonds=hb, show_cell=cell, triad=triad, view=V,
                polyhedra=poly and is_mof, cells=sc, h_scale=hsc, label_offsets=offs)
        if ttl:
            ax.set_title(f"{S.name}: {dict(asu='asymmetric unit', molecules='molecules', packing='packing', cell='unit cell', supercell='supercell')[mode]}"
                         f", view {view}", fontsize=7)
        if style == "ortep":
            ax.text(0.01, 0.01, f"{prob:.0%} probability ellipsoids", transform=ax.transAxes, fontsize=6, color="#555")
        xyz = pd.DataFrame([{"label": s["label"], "el": s["el"], "x": s["xyz"][0], "y": s["xyz"][1], "z": s["xyz"][2],
                             "symmetry": S.symcode(s["op"], s["T"])} for s in sites])
        show(fig, f"{S.name}_{style}_{mode}_{view}", data=xyz, controls=False)
    # ---- H-bonds
    hbr = CX.hbonds(S)
    with tabs[1]:
        if hbr:
            d = pd.DataFrame(hbr).drop(columns=["_A"])
            d["D···A"] = [CX.fmt_esd(v, e) for v, e in zip(d["D···A (Å)"], d["σ D···A"])]
            st.dataframe(d.drop(columns=["σ D···A"]).round(3), hide_index=True, width="stretch")
            st.caption("D–H···A with H···A < Σr_vdW (Bondi) and ∠DHA ≥ 120°; symmetry code of the acceptor "
                       "(SHELX style: operator + lattice translation). D···A ESD from the coordinate s.u. only.")
            dl(st, "H-bonds (CSV)", lambda: d.to_csv(index=False).encode(), f"{S.name}_hbonds.csv", "cx_hb_csv")
        else:
            st.info("No D–H···A hydrogen bonds found (no H on N/O donors, or no contacts within the vdW criterion).")
    # ---- SI tables
    with tabs[2]:
        ih = st.checkbox("Include bonds to H", False, key="cx_si_h")
        po = st.checkbox("Only rows flagged for publication (_geom_*_publ_flag yes)", False, key="cx_si_publ")
        geom = CX.geom_tables(S, include_h=ih, publ_only=po)
        n0 = st.number_input("First table number (S#)", 1, 99, 1, key="cx_si_n")
        secs = CS.si_sections(S, geom, hbr, int(n0))
        for title, df, cols, note in secs:
            with st.expander(f"{title} ({len(df)} rows)", expanded=title.startswith(f"Table S{int(n0)}.")):
                st.dataframe((df[cols] if cols else df).astype(str), hide_index=True, width="stretch")
                if note:
                    st.caption(note)
        c1, c2 = st.columns(2)
        dl(c1, "SI tables (DOCX)", lambda: CS.docx(secs), f"{S.name}_SI_tables.docx", "cx_si_docx")
        dl(c2, "SI tables (XLSX)", lambda: CS.xlsx(secs), f"{S.name}_SI_tables.xlsx", "cx_si_xlsx")
        _si_put = st.session_state.setdefault("si_tables_extra", {}); _si_put[S.name] = secs
    # ---- PXRD
    with tabs[3]:
        p1, p2, p3 = st.columns(3)
        lo, hi = p1.slider("2θ range (°)", 2.0, 90.0, (3.0, 50.0), 0.5, key="cx_px_rng")
        fw = p2.number_input("FWHM (°)", 0.02, 1.0, 0.10, 0.01, key="cx_px_fw")
        nl = p3.number_input("Label strongest n", 0, 40, 12, key="cx_px_nl")
        with st.spinner("Simulating…"):
            try:
                sim, sticks = _cif_pxrd(text, lo, hi, fw)
            except Exception as e:
                st.error(f"Simulation failed: {e}"); sim = None
        if sim is not None:
            mf = st.file_uploader("Measured pattern for overlay (.xy / .txt / .csv)", key="cx_px_meas")
            meas = fitr = None
            if mf is not None:
                try:
                    from mfs import readers as RD
                    m = RD.read_xy(mf); mx, my = np.asarray(m.iloc[:, 0], float), np.asarray(m.iloc[:, 1], float)
                    sel = (mx >= lo) & (mx <= hi); meas = (mx[sel], my[sel]); fitr = CPX.fit_to(meas[0], meas[1], sim)
                    st.caption(f"Least-squares scale {fitr['scale']:.4g} + linear background; Rwp = {fitr['Rwp']:.3f} "
                               "(profile agreement only — not a Rietveld refinement).")
                except Exception as e:
                    st.error(f"Measured pattern not read: {e}")
            if meas is not None:
                fig, (ax, rax) = plt.subplots(2, 1, figsize=(5.5, 3.6), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
            else:
                fig, ax = plt.subplots(figsize=(5.5, 3.0)); rax = None
            CPX.draw(ax, sim, sticks, label_top=int(nl), meas=meas, fitres=fitr, rax=rax, name=S.name)
            show(fig, f"{S.name}_simulated_PXRD", data=pd.DataFrame({"2theta": sim.x, "I_sim": sim.y}))
            st.dataframe(CPX.table(sticks, 30), hide_index=True, width="stretch")
    # ---- pores
    if is_mof:
        with tabs[4]:
            from mfs import poretools as PT
            q1, q2, q3 = st.columns(3)
            spc = q1.select_slider("Grid spacing (Å)", [0.25, 0.3, 0.4, 0.5, 0.7], 0.4, key="cx_pt_sp")
            pr = q2.number_input("Probe radius (Å)", 0.5, 3.0, 1.82, 0.01, key="cx_pt_pr", help="1.82 Å ≈ N₂ (kinetic diameter 3.64 Å)")
            if q3.button("Analyse pores", key="cx_pt_go", type="primary"):
                with st.spinner("Grid + percolation + Monte Carlo…"):
                    st.session_state["cx_pt_res"] = (S.name, PT.analyse(S, spc, pr, part="major"))
            r = st.session_state.get("cx_pt_res")
            if r and r[0] == S.name:
                res = r[1]
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Di (largest included sphere)", f"{res['Di (Å)']:.2f} Å"); m2.metric("Df (largest free sphere)", f"{res['Df (Å)']:.2f} Å")
                m3.metric("Void fraction", f"{res['void fraction (probe 0)']:.3f}"); m4.metric("Accessible SA", f"{res['ASA (m² g⁻¹)']:.0f} m² g⁻¹")
                tbl = pd.DataFrame([(kk, (f"{v:.4g}" if isinstance(v, (float, np.floating)) else v)) for kk, v in res.items() if kk != "method"],
                                   columns=["quantity", "value"]).astype(str)
                st.dataframe(tbl, hide_index=True, width="stretch")
                st.caption("Method: " + res["method"] + " Compare with Zeo++/RASPA for publication; values depend on radii and grid.")
                dl(st, "Pore analysis (CSV)", lambda: pd.concat([tbl, pd.DataFrame([["method", res["method"]]], columns=tbl.columns)]).to_csv(index=False).encode(),
                   f"{S.name}_pores.csv", "cx_pt_csv")
    # ---- Hirshfeld
    with tabs[-1]:
        st.warning("EXPERIMENTAL: approximate promolecule densities (one exponential per atom, not the Spackman tables of "
                   "CrystalExplorer). Shapes and contact percentages are indicative; use CrystalExplorer for publication.")
        if is_mof or getattr(S, "polymeric", False):
            st.info("Hirshfeld surfaces need a discrete molecule; not available for frameworks.")
        elif st.button("Compute Hirshfeld surface (first molecule)", key="cx_hs_go"):
            from mfs import hirshfeld as HS
            with st.spinner("Ray search for w = 0.5…"):
                mol, mb = S.grow("molecules", part="major")
                comps = CX._components(mol, mb); m = [mol[i] for i in comps[0]]
                env, _ = S.grow("packing", part="major", supercell=(3, 3, 3))
                env = [dict(e, frac=e["frac"] - 1, xyz=S.cart(e["frac"] - 1)) for e in env]
                st.session_state["cx_hs"] = (S.name, HS.surface(S, m, env, n_dir=160))
        r = st.session_state.get("cx_hs")
        if r and r[0] == S.name:
            from mfs import hirshfeld as HS
            hs = r[1]
            fig, (a1_, a2_) = plt.subplots(1, 2, figsize=(6.5, 3.2))
            sc_ = HS.draw_dnorm(a1_, hs, CP.view_matrix(S, "a")); fig.colorbar(sc_, ax=a1_, shrink=0.7, label=r"$d_\mathrm{norm}$")
            HS.draw_fingerprint(a2_, hs, highlight=["O···H", "H···O"])
            show(fig, f"{S.name}_Hirshfeld_experimental", data=pd.DataFrame({"di": hs["di"], "de": hs["de"], "d_norm": hs["d_norm"], "pair": hs["pair"]}), controls=False)
            st.dataframe(HS.contacts(hs).round(1), hide_index=True)


def page_guide():
    st.header("Per-instrument export guide")
    st.markdown(open(os.path.join(HERE, "docs", "EXPORT_GUIDE.md"), encoding="utf-8").read())


# ============================================================== Batch D: batch processing, SI bundle, sample database
def page_batch():
    from mfs import batchproc as BP, templates as TP
    st.header("Batch processing")
    st.caption("All files of one folder (or an upload) → the standard figure of the technique with ONE preset → PNG/TIFF/PDF/SVG + CSV, "
               "zipped and optionally written to a folder. Preset = a figure-settings JSON saved under any figure (🎛️ → Presets) "
               "or the current house style.")
    c1, c2 = st.columns(2)
    kind = c1.selectbox("Technique", list(BP.KINDS), key="bp_kind")
    src = c2.radio("Files from", ["a folder on this computer", "upload", "the bundled EXAMPLE folder"], index=2, key="bp_src")
    exdir = {"PXRD": "PXRD", "FTIR": "FTIR", "Raman": "Raman", "XPS": "XPS", "TGA": "TGA", "DSC": "DSC", "N2 isotherm": "Sorption",
             "SEM/TEM particles": "SEM"}[kind]
    files = []
    if src.startswith("a folder"):
        folder = st.text_input("Folder (sub-folders included)", os.path.expanduser("~"), key="bp_folder")
        if os.path.isdir(folder):
            fl = BP.list_folder(folder, kind); files = [(p_, p_) for p_ in fl]
            st.caption(f"{len(fl)} {kind} file(s) found" + (": " + ", ".join(os.path.basename(p_) for p_ in fl[:8]) + (" …" if len(fl) > 8 else "") if fl else ""))
        else:
            st.warning("Folder not found.")
    elif src == "upload":
        up = st.file_uploader("Files", accept_multiple_files=True, key="bp_up")
        files = [(u, u.name) for u in up or []]
    else:
        fl = BP.list_folder(os.path.join(EX, exdir), kind); files = [(p_, p_) for p_ in fl]
        st.caption(f"{len(fl)} example file(s): " + ", ".join(os.path.basename(p_) for p_ in fl))
    popts = None
    if kind == "SEM/TEM particles":       # v2.11: the SEM page's adaptive detection + optional shape model
        with st.container(border=True, key="bp_particles"):
            q1, q2, q3, q4 = st.columns(4)
            pdark = q1.radio("Particles are", ["dark (TEM)", "bright (SEM)"], index=1 if src.startswith("the bundled") else 0,
                             key="bp_pdark", horizontal=True).startswith("dark")
            pprior = q2.toggle("🔷 Shape model", False, key="bp_prior",
                               help="Learn the outline of the clean isolated particles of each image and use it to split clumps.")
            pstr_ = q2.slider("Shape-model strength", 0.0, 1.0, 0.7, 0.05, key="bp_pstr", disabled=not pprior)
            pnm = q3.number_input("Fallback nm/px (files without metadata)", 0.0, 1e4, 0.0, format="%g", key="bp_pnm")
            plk = q3.slider("Sensitivity (× noise σ)", 2.0, 8.0, 4.0, 0.5, key="bp_plk")
            pagg = q4.checkbox("Exclude aggregates", True, key="bp_pagg")
            popts = dict(dark=pdark, shape_prior=pprior, prior_strength=pstr_, nm_per_px=pnm or None, low_k=plk, exclude_aggregates=pagg)
            st.caption("One figure per file = outlines on the image (aggregates red dashed, shape-model splits dotted); "
                       "'all files overlaid' = ECD histograms. CSV = every particle (ECD, Feret, kind, included).")
    k1, k2, k3, k4 = st.columns(4)
    pre = k1.radio("Preset", ["none (defaults)", "house style", "settings JSON"], key="bp_pre")
    cfg = None
    if pre == "house style":
        cfg = FU.house() or None
        if not cfg:
            k1.caption("No house style loaded – defaults are used.")
    elif pre == "settings JSON":
        pj = k1.file_uploader("Figure settings (JSON)", type=["json"], key="bp_json")
        if pj is not None:
            try:
                cfg = F.merge(F.default_settings(), json.loads(pj.getvalue().decode("utf-8")))
            except Exception as e:
                k1.error(f"Could not read the preset: {e}")
    layout = k2.radio("Layout", ["one figure per file", "all files overlaid", "stacked"], key="bp_lay")
    tpl = k3.selectbox("Journal template", ["sidebar width / font"] + list(TP.TEMPLATES), key="bp_tpl")
    colw = k3.radio("Column", ["single", "1.5", "double"], horizontal=True, key="bp_col") if tpl in TP.TEMPLATES else None
    fmts = k4.multiselect("Formats", ["png", "tiff", "pdf", "svg"], ["png", "pdf"], key="bp_fmt")
    dpi = k4.selectbox("dpi", [300, 600, 1200], index=1, key="bp_dpi")
    width = TP.width(tpl, colw) if tpl in TP.TEMPLATES else WIDTH
    fsz = TP.TEMPLATES[tpl]["font_size"] if tpl in TP.TEMPLATES else FONT
    out = st.text_input("Also write to folder (empty = ZIP download only)", "", key="bp_out",
                        placeholder=os.path.join(EXPORT_DIR, f"batch_{kind.split()[0]}"))
    if not files or not fmts:
        return
    if st.button(f"▶ Process {len(files)} file(s)", type="primary", key="bp_run"):
        pb = st.progress(0.0)
        z, log = BP.run(kind, files, cfg, layout, tuple(fmts), dpi, width, ASPECT, fsz, out or None, progress=lambda v: pb.progress(min(1.0, v)),
                        particle_opts=popts)
        st.session_state["bp_result"] = (z, log, kind)
    if st.session_state.get("bp_result"):
        z, log, k_ = st.session_state["bp_result"]
        lg = pd.DataFrame(log)
        ok = int((lg.status == "ok").sum()) if len(lg) else 0
        (st.success if ok == len(lg) else st.warning)(f"Batch done: {ok} / {len(lg)} file(s) read" + (f", written to {out}" if out else "") + ".")
        st.dataframe(lg, hide_index=True, width="stretch")
        dl(st, "⬇️ ZIP (figures + data + log + preset)", z, f"batch_{safe_name(k_)}.zip", "bp_zip")


def page_si():
    from mfs import si_bundle as SB
    st.header("SI export bundle")
    st.caption("Every figure shown while “Collect figures for the SI bundle” (sidebar) is on is kept here with its data and the "
               "settings used. One ZIP: figures (PDF vector, PNG at the sidebar dpi, SVG), data CSVs, SI_parameters.xlsx "
               "(figure list, analysis + figure settings, results per sample) and the defect SI table.")
    figs = st.session_state.get("si_figs", {})
    if not figs:
        st.info("No figures collected yet – open the technique pages (collection is on by default)."); return
    tab = pd.DataFrame([{"include": True, "figure": k, "page": v.get("page", ""), "captured": v.get("time", "")} for k, v in figs.items()])
    tab = st.data_editor(tab, hide_index=True, disabled=["figure", "page", "captured"], key="si_tab")
    sel = [r.figure for r in tab.itertuples() if r.include]
    inc = st.multiselect("Include", ["pdf", "png", "svg", "csv"], ["pdf", "png", "csv"], key="si_inc")
    c1, c2 = st.columns(2)
    if c2.button("🗑️ Clear the collection", key="si_clear"):
        st.session_state["si_figs"] = {}; st.rerun()
    if not sel:
        return

    def _bundle(sel=tuple(sel), inc=tuple(inc)):
        out = {}
        for k in sel:
            f = figs[k]; fig = f["fig"]
            d = {kk: f.get(kk) for kk in ("page", "time", "params", "cfg", "title", "n_rows")}
            with _RENDER_LOCK:
                for ext in ("pdf", "png", "svg"):
                    if ext in inc:
                        d[ext] = S.export(fig, ext, DPI)
            if "csv" in inc and f.get("data") is not None:
                d["csv"] = f["data"].to_csv(index=False).encode("utf-8-sig")
            out[k] = d
        return SB.build(out, st.session_state.get("sample_store", {}), __version__, inc)
    dl(c1, f"⬇️ SI bundle ({len(sel)} figures)", _bundle, "Reticula_SI_bundle.zip", "si_zip")
    figt, part, rest = SB.tables({k: figs[k] for k in sel}, st.session_state.get("sample_store", {}))
    st.markdown("**Parameters table (preview)**")
    st.dataframe(part, hide_index=True, width="stretch", height=260)


def _db_path():
    return os.environ.get("MFS_DB") or os.path.join(os.path.dirname(STORE_PATH), "reticula.db")


def _pipeline_candidates():
    h = os.path.expanduser("~")
    c = [os.path.join(h, d, "Synthesis_Pipeline") for d in ("Documents", "Downloads", "Desktop", "OneDrive\\Documents", "")]
    c += [os.path.join(os.path.dirname(STORE_PATH), "Synthesis_Pipeline"), os.path.join(HERE, "references", "Synthesis_Pipeline")]
    return [p_ for p_ in c if os.path.exists(os.path.join(p_, "Master_sample_overview.xlsx"))]


def page_db():
    from mfs import sampledb as DB
    st.header("Sample database")
    path = _db_path()
    st.caption(f"Local SQLite database `{path}`: sample IDs from the synthesis pipeline (protocol, family, parent, source), "
               "every result stored by the technique pages, and data files linked by the sample ID in their names "
               "(e.g. `Z801-F25-B1_PXRD.xy`).")
    con = DB.connect(path)
    cands = _pipeline_candidates()
    with st.expander("1 · Import sample IDs from the synthesis pipeline", expanded=not con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]):
        c1, c2 = st.columns([3, 1])
        root = c1.text_input("Synthesis_Pipeline folder (contains Master_sample_overview.xlsx and PaperN/Protocols)",
                             cands[0] if cands else "", key="db_root")
        if cands and os.path.join(HERE, "references", "Synthesis_Pipeline") == cands[0]:
            c1.caption("Using the copy of the master overview bundled with Reticula (protocol files not included) – point this to your "
                       "own Synthesis_Pipeline folder to link the protocol PDFs/DOCX.")
        if c2.button("Import / update", key="db_imp"):
            try:
                n = DB.import_master(con, os.path.join(root, "Master_sample_overview.xlsx")); st.success(f"{n} sample IDs imported from {root}.")
            except Exception as e:
                st.error(f"Import failed: {e}")
    nres = DB.sync_store(con, st.session_state.get("sample_store", {}))
    with st.expander("2 · Link data files"):
        c1, c2 = st.columns([3, 1])
        dfold = c1.text_input("Data folder (sub-folders included)", os.path.join(os.path.expanduser("~"), "Documents"), key="db_data")
        if c2.button("Scan and link", key="db_scan"):
            if os.path.isdir(dfold):
                n, un = DB.link_folder(con, dfold)
                st.success(f"{n} file(s) linked; {len(un)} without a known sample ID.")
            else:
                st.error("Folder not found.")
    ov = DB.overview(con)
    st.caption(f"{len(ov)} samples · {nres} result values synced from the sample store this run.")
    q = st.text_input("🔎 Search (ID, name, family, paper)", "", key="db_q")
    if q:
        m_ = ov.astype(str).apply(lambda c: c.str.contains(q, case=False, regex=False)).any(axis=1); ov = ov[m_]
    st.dataframe(ov, hide_index=True, width="stretch", height=300)
    if not len(ov):
        return
    sid = st.selectbox("Sample", list(ov["Sample ID"]), key="db_pick")
    s_ = DB.sample(con, sid)
    st.subheader(f"{sid} – {s_.get('name') or ''}")
    a, b = st.columns(2)
    a.markdown(f"**Paper:** {s_.get('paper') or '–'}  \n**Family:** {s_.get('family') or '–'}  \n**Parent:** {s_.get('parent') or '–'}  \n"
               f"**Children:** {', '.join(s_['children']) or '–'}  \n**Basis:** {s_.get('basis') or '–'}"
               + (f"  \n**Source:** [{s_['doi']}]({s_['doi']})" if s_.get("doi") and str(s_["doi"]).startswith("http") else ""))
    for lab, k_ in (("protocol PDF", "protocol_pdf"), ("protocol DOCX", "protocol_docx"), ("paper matrix", "matrix")):
        p_ = s_.get(k_)
        if p_ and os.path.exists(p_):
            dl(b, f"⬇️ {lab}: {os.path.basename(p_)}", open(p_, "rb").read(), os.path.basename(p_), f"db_{k_}_{sid}")
            if b.button(f"Open {lab}", key=f"db_open_{k_}_{sid}"):
                open_with_default_app(p_)
        elif p_:
            b.caption(f"{lab}: `{p_}` (not found on this computer)")
    if s_["results"]:
        st.markdown("**Results**")
        st.dataframe(pd.DataFrame(s_["results"]).astype({"value": str}), hide_index=True, width="stretch")
    if s_["files"]:
        st.markdown("**Data files**")
        st.dataframe(pd.DataFrame(s_["files"]), hide_index=True, width="stretch")
    nt = st.text_area("Notes", s_.get("notes") or "", key=f"db_notes_{sid}")
    if st.button("Save notes", key=f"db_savenotes_{sid}"):
        DB.set_notes(con, sid, nt); st.success("Saved.")
    with open(path, "rb") as fh:
        dl(st, "⬇️ Database file (reticula.db)", fh.read(), "reticula.db", "db_dl")


def page_uvvis():
    from mfs import uvvis_page as UVP
    UVP.render(globals())


UI.header(HERE, page)
try:                                   # v2.11: opt-in update check (asks once; background thread; silent offline)
    from mfs import update as _UPD
    _UPD.banner(os.path.dirname(STORE_PATH), __version__)
except Exception:
    pass
_store_boot()
PS.start_page(page)
_PAGEFN = {"Home": page_home, "PXRD": page_pxrd, "FTIR / Raman": page_ftir, "1H NMR (digest)": page_nmr, "N2 sorption (BET/PSD)": page_sorption,
 "TGA": page_tga, "DSC": page_dsc, "Defect calculator": page_defects, "Kinetics & adsorption": page_kinetics, "ICP-OES / XRF": page_icp, "SEM / TEM / EDX": page_sem,
 "SDS-PAGE (Image Lab)": page_sds, "Zeta / DLS": page_zeta, "XPS": page_xps, "Figure composer": page_composer, "Batch processing": page_batch,
 "SI export bundle": page_si, "Sample database": page_db, "Structure (CIF)": page_structure, "UV/Vis": page_uvvis, "Export guide": page_guide}
try:
    if os.environ.get("MFS_PROFILE"):        # dev: cProfile of the page function -> <MFS_PROFILE>_<page>.prof
        import cProfile as _cp
        _pr = _cp.Profile(); _pr.enable()
        try:
            _PAGEFN[page]()
        finally:
            _pr.disable(); _pr.dump_stats(f"{os.environ['MFS_PROFILE']}_{PS.slug(page)}.prof")
    else:
        _PAGEFN[page]()
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
    with st.sidebar.expander("💾 Saved page work", expanded=False):
        sp = PS.saved_pages()
        st.caption("Saved with the dialog or 💾 Save; restored on the next start. " + (", ".join(f"{k} ({v[5:16].replace('T', ' ')})"
                   for k, v in sp.items() if v) if sp else "Nothing saved yet."))
        if sp and st.button("Forget all saved page work", key="_ps_forget_all"):
            PS.forget_all(); st.rerun()

