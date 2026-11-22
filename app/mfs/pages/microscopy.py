"""Reticula page module: SEM / TEM / EDX page (scale bars, calibration, particle sizing, EDX maps).

Moved verbatim out of app.py in v2.13. Names that the code takes from the app (helpers such as show(), dl(), files_input(),
the sidebar values WIDTH / ASPECT / FONT, the cached reader wrapper R …) are listed in INJECTED and filled in by
mfs.pages.bind(globals()) on every run of app.py, so the code behaves exactly as it did inside app.py."""
# ruff: noqa: F821
from __future__ import annotations
import io, os, re, json, glob, hashlib
import numpy as np
import pandas as pd
import streamlit as st
from mfs import numentry as NE  # every slider takes a typed number too (2.19.6)
import matplotlib.pyplot as plt

from mfs import plots as P, style as S, composer as C
from mfs import pagestate as PS

INJECTED = (
    'EX',
    'FONT',
    'R',
    'WIDTH',
    '_F',
    'files_input',
    'guess_name',
    'newfig',
    'show',
    'store_put',
)


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


def _particles_auto(img, pix, *args):
    from mfs import segpool as SPL                 # v2.14: own cache (+ parallel prefetch); same results as before
    return SPL.get("auto", img, (pix,) + args)


def _particles_ws(img, *args):
    from mfs import segpool as SPL
    return SPL.get("watershed", img, args)


def _particles_adaptive(img, *args):
    from mfs import segpool as SPL
    return SPL.get("adaptive", img, args)


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


_TPL_SKIP = ("ps_imgs", "ps_sample", "ps_meth", "ps_tpl_name", "ps_tpl_pick", "ps_nums")


def _tpl_path():
    base = os.path.dirname(os.environ.get("MFS_STORE") or os.path.join(os.path.expanduser("~"), "Documents", "Reticula", "sample_store.json"))
    return os.path.join(base, "particle_templates.json")


def _tpl_load():
    ss = st.session_state
    if "ps_templates" not in ss:
        try:
            with open(_tpl_path(), encoding="utf-8") as fh:
                ss["ps_templates"] = json.load(fh)
        except Exception:
            ss["ps_templates"] = {}
    return ss["ps_templates"]


def _tpl_save_file():
    try:
        os.makedirs(os.path.dirname(_tpl_path()), exist_ok=True)
        with open(_tpl_path(), "w", encoding="utf-8") as fh:
            json.dump(st.session_state.get("ps_templates", {}), fh, ensure_ascii=False, indent=1)
    except Exception as e:
        st.session_state["_ps_tpl_err"] = str(e)


def _tpl_cb(action):
    """Callbacks (run before the widgets are drawn, so the stored values can be written into them)."""
    ss = st.session_state
    tpls = _tpl_load()
    name = (ss.get("ps_tpl_name") or "").strip()
    if action == "save" and name:
        tpls[name] = {k: v for k, v in ss.items() if isinstance(k, str) and k.startswith("ps_") and not k.startswith("ps_tpl") and k not in _TPL_SKIP
                      and isinstance(v, (bool, int, float, str))}
        _tpl_save_file()
        ss["ps_tpl_msg"] = f"Template '{name}' saved ({len(tpls[name])} settings)."
    elif action == "apply" and ss.get("ps_tpl_pick") in tpls:
        for k, v in tpls[ss["ps_tpl_pick"]].items():
            ss[k] = v
        ss["ps_tpl_msg"] = f"Template '{ss['ps_tpl_pick']}' applied."
    elif action == "delete" and ss.get("ps_tpl_pick") in tpls:
        tpls.pop(ss["ps_tpl_pick"], None)
        _tpl_save_file()
        ss["ps_tpl_msg"] = "Template deleted."


def _size_templates_ui():
    """v2.15: the detection settings that work for one kind of sample (e.g. 'MOF-801 SEM', 'ZrO2 TEM') are saved under a name and applied
    with one click - no need to tune the same sliders again. Kept across sessions in a file next to the sample store and in project files."""
    tpls = _tpl_load()
    with st.expander(f"💾 Sizing templates ({len(tpls)} saved) – keep the settings per sample type", expanded=False):
        c1, c2, c3 = st.columns([2, 1, 1], vertical_alignment="bottom")
        c1.text_input("Save the current settings as", "", key="ps_tpl_name", placeholder="e.g. MOF-801 SEM")
        c2.button("💾 Save", key="ps_tpl_save", on_click=_tpl_cb, args=("save",), width="stretch")
        if tpls:
            d1, d2, d3 = st.columns([2, 1, 1], vertical_alignment="bottom")
            d1.selectbox("Saved templates", list(tpls), key="ps_tpl_pick")
            d2.button("↩ Apply", key="ps_tpl_apply", on_click=_tpl_cb, args=("apply",), type="primary", width="stretch")
            d3.button("🗑 Delete", key="ps_tpl_del", on_click=_tpl_cb, args=("delete",), width="stretch")
        if st.session_state.get("ps_tpl_msg"):
            st.caption("✅ " + st.session_state.pop("ps_tpl_msg"))


def _confidence_bar(info):
    """v2.15: how sure the sphere model is about each particle: one stacked bar (high ≥ 0.8, medium, low = excluded / orange)."""
    cs = [float(v["conf"]) for k, v in info.items() if isinstance(k, int) and isinstance(v, dict) and "conf" in v]
    if not cs:
        return
    cm = float((info.get("_sphere") or {}).get("conf_min", 0.5))
    hi = sum(c >= 0.8 for c in cs)
    lo = sum(c < cm for c in cs)
    mid = len(cs) - hi - lo
    w = lambda n: 100.0 * n / len(cs)
    st.markdown(
        "<div style='font-size:0.8em;color:#64748B;margin:2px 0'>Confidence per particle</div>"
        "<div style='display:flex;height:10px;border-radius:5px;overflow:hidden;background:#E2E8F0'>"
        f"<div style='width:{w(hi):.1f}%;background:#10B981' title='high'></div>"
        f"<div style='width:{w(mid):.1f}%;background:#FBBF24' title='medium'></div>"
        f"<div style='width:{w(lo):.1f}%;background:#F97316' title='low (excluded)'></div></div>"
        f"<div style='font-size:0.78em;color:#64748B'>🟢 {hi} high · 🟡 {mid} medium · 🟠 {lo} low (excluded; threshold {cm:.2f})</div>",
        unsafe_allow_html=True)


def _sphere_hint(ph, info):
    """v2.12: one line next to the sphere controls: what the circle model decided and did."""
    s_ = info.get("_sphere") or {}
    rd = s_.get("roundness") or {}
    with ph.container():
        if s_.get("active"):
            st.markdown(f"⚪ **Spheres{' detected' if s_.get('mode') == 'auto' else ''}**: {s_.get('n', 0)} circle fits · "
                        f"{s_.get('n_split', 0)} in touching / overlapping groups · {s_.get('n_occluded', 0)} occluded · "
                        f"**{s_.get('n_lowconf', 0)} low-confidence** (orange, excluded).")
            _confidence_bar(info)
        elif s_.get("mode") == "auto":
            st.caption("Shape: " + (f"faceted / irregular (roundness {rd['roundness']:.3f}) – outline model used."
                                    if rd.get("roundness") is not None else rd.get("reason", "not decided – outline model used.")))
        else:
            st.caption("Sphere model off – outlines of the adaptive detection.")


def _prior_panel(pr, pix_nm, f):
    """Learned shape prior: template outline, N example particles, size histogram and the weight."""
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
    if meth.startswith("automatic"):
        _size_templates_ui()
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
        # v2.12: sphere model (circle fit on the real edge, splits touching / overlapping spheres, confidence per particle)
        e1, e2, e3 = _sbox.columns([1.6, 1.6, 2.8], vertical_alignment="center")
        sph_mode = e1.selectbox("⚪ Particle shape", ["auto", "spheres", "off"], key="ps_sphere", disabled=not _ad_now,
                                format_func={"auto": "auto (detect round / faceted)", "spheres": "spheres – circle fit",
                                             "off": "outlines only (v2.11)"}.get,
                                help="auto: when the clean particles are round, every particle is fitted with a circle on its "
                                     "real edge (sub-pixel), touching pairs / triples are split and overlapping spheres are "
                                     "fitted on their visible arcs. Faceted particles (cubes, octahedra) keep the outline model.")
        conf_min = e2.slider("Low-confidence threshold", 0.0, 1.0, 0.5, 0.05, key="ps_conf", disabled=not _ad_now,
                             help="Sphere fits below this confidence (little visible edge, poor circle fit, partly outside the "
                                  "particle) are excluded automatically and drawn orange – click one to include it.")
        _sph_hint = e3.empty()
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
        mres = k3.select_slider("Working resolution (px)", [512, 768, 1024, 1536, 2048], 1536, key="ps_res",
                                help="Keep 1536 or more for sizing: below ~1000 px the sphere model switches off and the sizes can be 20 % off.")
        if mres < 1024:
            k3.caption("⚠️ A low working resolution is only for a quick look: sizes can be off by up to 20 %.")
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
        from mfs import segpool as SPL                             # v2.14: all images at once on all CPU cores, with progress
        _jobs = []
        for im in ims:
            _img = im["image"][:im.get("crop_row") or None]
            if ws and ad:
                _jobs.append(("adaptive", _img, (dark, smooth, lowk, (expnm / float(im["pixel_nm"])) if expnm else 0.0,
                                                 mks, sol, areaf, border, mres, prior, pstr, pk_, sph_mode, conf_min)))
            elif ws:
                _jobs.append(("watershed", _img, (dark, tp, smooth, opx, mdf, hf, meth, holes, border, sol, areaf, mres)))
            else:
                _jobs.append(("auto", _img, (float(im["pixel_nm"]), dark, tp, smooth, opx, split, border, mres)))
        _todo = [j for j in _jobs if not SPL.cached(*j)]
        if len(_todo) > 1:
            _bar = st.empty(); _st0 = {}

            def _cb(done, total, sec):
                base = _st0.setdefault("base", done)
                new = done - base
                left = (sec / new * (total - done)) if new else None
                _bar.progress(min(0.99, done / total), text=f"Segmenting particles · {done} of {total} images"
                              + (f" · {sec:.0f} s elapsed · about {left:.0f} s left" if left is not None else " · working …"))
            SPL.prefetch(_jobs, progress=_cb)
            _bar.empty()
            if SPL.LAST.get("mode") == "parallel":
                st.caption(f"⚡ {len(_todo)} images segmented on {SPL.LAST['workers']} CPU cores in {SPL.LAST['seconds']:.0f} s.")
        elif _todo:
            with st.spinner("Segmenting particles …"):
                SPL.prefetch(_jobs, parallel=False)
        for im in ims:
            img = im["image"][:im.get("crop_row") or None]
            if ws:
                if ad:
                    lab, f, t, info = _particles_adaptive(img, dark, smooth, lowk, (expnm / float(im["pixel_nm"])) if expnm else 0.0,
                                                          mks, sol, areaf, border, mres, prior, pstr, pk_, sph_mode, conf_min)
                    if im is ims[0]:
                        _prior_panel(info.get("_prior"), float(im["pixel_nm"]), f)
                        _shape_model_hint(_shape_hint, info, prior)
                        _sphere_hint(_sph_hint, info)
                    algo_par = (lowk, expnm, mks, prior, pstr, pk_, sph_mode, conf_min)
                else:
                    lab, f, t, info = _particles_ws(img, dark, tp, smooth, opx, mdf, hf, meth, holes, border, sol, areaf, mres)
                    algo_par = ()
                if edit and im is ims[0]:
                    lab = _particle_editor(im, img, lab, f, info)
                r_ = PA.measure(lab, float(im["pixel_nm"]), f, info=info)
                for x in r_:
                    x["kind"] = info.get(x["id"], {}).get("kind", "particle") if isinstance(info.get(x["id"]), dict) else "particle"
                    x["solidity"] = info.get(x["id"], {}).get("solidity", np.nan) if isinstance(info.get(x["id"]), dict) else np.nan
                    x["prior"] = bool(isinstance(info.get(x["id"]), dict) and info[x["id"]].get("split") in ("prior", "circle"))
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
    nb = NE.slider("Histogram bins (0 = auto)", 0, 60, 0, key="ps_bins")
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
