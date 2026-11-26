"""SEM/TEM calibration from the ORIGINAL image: instrument metadata, the burned-in scale bar (detected + OCR of its
label) and the magnification printed in the info strip. Everything returns plain dicts so the app can show which
method was used and let the user confirm or override it.

Order of reliability used by calibrate():
  1. instrument metadata (FEI/Thermo PixelWidth, Zeiss 'Image Pixel Size', ImageJ, JEOL/Hitachi/TESCAN sidecars, DM/EMD)
  2. the scale bar burned into the image: longest solid horizontal bar in the info strip, length measured in pixels,
     label read by OCR (Tesseract; bundled with the Windows app) -> nm/px = value / length
  3. magnification (metadata or OCR) with the instrument's reference display width (127 mm for Philips/FEI XL)
  4. manual: the user marks the two bar ends and types the value."""
from __future__ import annotations
import math, os, re, shutil, sys
import numpy as np

UNIT_NM = {"pm": 1e-3, "å": 0.1, "a": 0.1, "nm": 1.0, "µm": 1e3, "um": 1e3, "μm": 1e3, "mm": 1e6, "m": 1e9}

# ---------------------------------------------------------------------------------------------- text parsing
_LABEL = re.compile(r"(?<![\w.,])(\d+(?:[.,]\d+)?)\s*(pm|nm|[µμu]m|mm|Å)(?![a-zA-Z])", re.I)
_MAG = [
    re.compile(r"mag(?:nification)?\s*[=:]?\s*(\d+(?:[.,]\d+)?)\s*([kKmM])?\s*[xX×]", re.I),          # Mag = 10.00 K X
    re.compile(r"(?<![\w.])[xX×]\s*(\d{1,3}(?:[ ,.\u2009]\d{3})+|\d+(?:[.,]\d+)?)\s*([kKmM])?(?![\w])"),  # x10,000  ×5.00k
    re.compile(r"(?<![\w.])(\d{1,3}(?:[ ,.\u2009]\d{3})+|\d+(?:[.,]\d+)?)\s*([kKmM])?\s*[xX×](?![a-wyzA-WYZ])"),  # 10.0kx  50 000 x
]


def parse_scale_label(text: str):
    """'1 µm', '500nm', '1um', '2 μm', '0.5 µm' -> (value, unit, nm). Returns None if nothing found.
    OCR often reads µ as 'u' or 'p'; 'pm' next to a bar in an SEM databar is almost always µm, so a bare
    'pm' with a value >= 1 is NOT re-interpreted (PM-scale bars do exist on HR-TEM)."""
    if not text:
        return None
    t = text.replace("μ", "µ").replace("\u00b5", "µ")
    m = _LABEL.search(t)
    if not m:
        return None
    v = float(m.group(1).replace(",", "."))
    u = m.group(2).lower().replace("μ", "µ")
    u = "µm" if u in ("um", "µm") else ("Å" if u == "å" else u)
    f = UNIT_NM.get(u.lower() if u != "Å" else "å")
    return (v, u, v * f) if f and v > 0 else None


def _num(s: str) -> float:
    s = s.strip()
    if re.fullmatch(r"\d{1,3}(?:[ ,.\u2009]\d{3})+", s):        # thousands separators: 10,000 / 50 000 / 10.000
        return float(re.sub(r"[ ,.\u2009]", "", s))
    return float(s.replace(",", "."))


def parse_magnification(text: str):
    r = match_magnification(text)
    return r[0] if r else None


def match_magnification(text: str):
    """'x10,000', '10.0kx', '50 000 x', 'Mag = 10.00 K X', '×5.00k', '10000x', 'MAG: 25.0 kx' -> (10000.0, matched text)."""
    if not text:
        return None
    t = text.replace("×", "x")
    for rx in _MAG:
        m = rx.search(t)
        if m:
            try:
                v = _num(m.group(1))
            except ValueError:
                continue
            k = (m.group(2) or "").lower()
            v *= 1e3 if k == "k" else (1e6 if k == "m" else 1)
            if 1 <= v <= 1e8:
                return v, m.group(0).strip()
    return None


# ---------------------------------------------------------------------------------------------- image analysis
def find_info_strip(img8: np.ndarray):
    """First row of the instrument's info/data strip at the bottom of the image, or None.
    A strip is a band of rows whose median is (near) black or white or whose texture differs abruptly from the
    image above it (low row-to-row correlation boundary)."""
    h, w = img8.shape
    med = np.median(img8, axis=1)
    flat = (med <= 12) | (med >= 243)
    if flat[-1]:
        for r in range(int(h * 0.55), h - 3):
            if flat[r:r + 3].all() and flat[r:].mean() > 0.6:
                return r
    # abrupt boundary: a row whose fraction of near-black/white pixels jumps
    ext = ((img8 <= 12) | (img8 >= 243)).mean(axis=1)
    lo = int(h * 0.6)
    if h - lo > 10:
        d = np.diff(ext[lo:])
        k = int(np.argmax(d))
        if d[k] > 0.45 and ext[lo + k + 1:].mean() > 0.5:
            return lo + k + 1
    return None


def _runs(mask_row):
    """(start, stop) of True runs in a 1-D bool array."""
    x = np.concatenate([[0], mask_row.astype(np.int8), [0]])
    d = np.diff(x)
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def detect_bar(img8: np.ndarray, strip_row=None, min_frac=0.03, max_frac=0.9):
    """Longest solid horizontal bar (filled line/rectangle, white or black) in the info strip (or the bottom 25 %).
    Returns dict(x0, x1, y0, y1, length_px, polarity, region_top) or None. length_px = x1 - x0 (pixel edges)."""
    h, w = img8.shape
    top = strip_row if strip_row is not None else int(h * 0.75)
    reg = img8[top:]
    best = None
    pols = [("white", reg >= 200), ("black", reg <= 55)]
    bg = np.median(reg)
    if strip_row is not None:          # a bar contrasts with the strip background
        pols = [pols[0]] if bg < 128 else [pols[1]]
    for pol, mask in pols:
        # bars: runs that are not the whole strip background
        rows = {}
        for y in range(mask.shape[0]):
            for a, b in _runs(mask[y]):
                L = b - a
                if min_frac * w <= L <= max_frac * w:
                    rows.setdefault((a, b), []).append(y)
        # merge runs whose ends agree within 2 px into one bar; require >= 2 rows (thickness) unless very long
        segs = []
        for (a, b), ys in rows.items():
            for s in segs:
                if abs(s["x0"] - a) <= 2 and abs(s["x1"] - b) <= 2 and (min(abs(y - s["y1"]) for y in ys) <= 2 or min(abs(y - s["y0"]) for y in ys) <= 2):
                    s["ys"] += ys; s["x0s"].append(a); s["x1s"].append(b)
                    s["y0"], s["y1"] = min(s["ys"]), max(s["ys"]); break
            else:
                segs.append({"x0": a, "x1": b, "ys": list(ys), "x0s": [a], "x1s": [b], "y0": min(ys), "y1": max(ys)})
        for s in segs:
            thick = len(set(s["ys"]))
            if thick < 2 or thick > 0.35 * reg.shape[0]:
                continue
            # a bar is solid: the region it covers is (almost) all on, and it is not the strip background
            x0, x1 = int(np.median(s["x0s"])), int(np.median(s["x1s"]))
            blk = mask[s["y0"]:s["y1"] + 1, x0:x1]
            if blk.mean() < 0.85:
                continue
            # contrast to the rows just above/below (the bar must stand out)
            ya, yb = max(s["y0"] - 2, 0), min(s["y1"] + 3, mask.shape[0])
            ring = np.concatenate([mask[ya:s["y0"], x0:x1].ravel(), mask[s["y1"] + 1:yb, x0:x1].ravel()])
            if ring.size and ring.mean() > 0.6:
                continue
            if x0 <= 0 or x1 >= w:      # touches the image edge: strip border / background, not a bar
                continue
            L = x1 - x0
            score = L * (1 + 0.02 * min(thick, 10))
            if best is None or score > best["score"]:
                best = {"x0": x0, "x1": x1, "y0": s["y0"] + top, "y1": s["y1"] + top, "length_px": float(L),
                        "thickness_px": thick, "polarity": pol, "region_top": top, "score": score}
    if best:
        best.pop("score")
    return best


# ---------------------------------------------------------------------------------------------- OCR
def tesseract_cmd():
    """Bundled Tesseract (Windows app: <app>/_internal/tesseract/tesseract.exe), $MFS_TESSERACT, or the one on PATH."""
    cands = [os.environ.get("MFS_TESSERACT")]
    base = getattr(sys, "_MEIPASS", None)
    if base:
        cands.append(os.path.join(base, "tesseract", "tesseract.exe" if os.name == "nt" else "tesseract"))
    cands.append(shutil.which("tesseract"))
    if os.name == "nt":
        cands.append(r"C:\Program Files\Tesseract-OCR\tesseract.exe")
    else:                                    # Homebrew (Apple silicon / Intel), MacPorts
        cands += ["/opt/homebrew/bin/tesseract", "/usr/local/bin/tesseract", "/opt/local/bin/tesseract"]
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def ocr_available():
    try:
        import pytesseract  # noqa: F401
    except Exception:
        return False
    return tesseract_cmd() is not None


def _prep(crop: np.ndarray, polarity="auto", scale=3, ink=None):
    scale = max(1, int(scale))
    from PIL import Image
    a = crop.astype(float)
    if ink is not None:
        # burned-in annotation: the text has exactly the bar's grey value (e.g. Olympus analySIS overlays drawn in
        # palette index 0 with a white outline on top of the image) -> keep only those pixels
        a = np.where(np.abs(a - ink) <= 6, 0.0, 255.0)
        polarity = "black"
    if polarity == "auto":
        polarity = "white" if np.median(a) < 128 else "black"
    if polarity == "white":      # white text on dark strip -> invert to dark-on-white for Tesseract
        a = 255 - a
    im = Image.fromarray(np.clip(a, 0, 255).astype(np.uint8))
    im = im.resize((im.width * scale, im.height * scale), Image.LANCZOS)
    thr = (np.array(im) > 150).astype(np.uint8) * 255
    pad = 12 * scale                     # Tesseract needs a quiet margin around the text
    thr = np.pad(thr, pad, constant_values=255)
    return Image.fromarray(thr), scale, pad


LABEL_CHARS = "0123456789.,µumnkKxXÅ"      # no "p": Tesseract reads µ as p


def ocr_words(region: np.ndarray, x_off=0, y_off=0, psm=11, scale=2, whitelist=None, ink=None):
    """Words with boxes (in original image px) from a region. [] if OCR is unavailable."""
    if not ocr_available():
        return []
    import pytesseract
    pytesseract.pytesseract.tesseract_cmd = tesseract_cmd()
    im, s, pad = _prep(region, scale=scale, ink=ink)
    tdir = os.path.join(os.path.dirname(tesseract_cmd()), "tessdata")
    if os.path.isdir(tdir):                 # bundled copy: point Tesseract at its own tessdata (no quoting issues on Windows)
        os.environ["TESSDATA_PREFIX"] = tdir
    cfg = f"--psm {psm}"
    if whitelist:
        cfg += f" -c tessedit_char_whitelist={whitelist}"
    try:
        d = pytesseract.image_to_data(im, config=cfg, output_type=pytesseract.Output.DICT, timeout=30)
    except Exception:
        return []
    out = []
    for i, t in enumerate(d["text"]):
        t = (t or "").strip()
        if not t:
            continue
        L, T = d["left"][i] - pad, d["top"][i] - pad
        out.append({"text": t, "x0": L / s + x_off, "y0": T / s + y_off,
                    "x1": (L + d["width"][i]) / s + x_off, "y1": (T + d["height"][i]) / s + y_off,
                    "conf": float(d["conf"][i]), "line": (d["block_num"][i], d["par_num"][i], d["line_num"][i])})
    return out


def _lines(words):
    """Join words that sit on the same OCR line, keep their union box."""
    L = {}
    for w in words:
        L.setdefault(w["line"], []).append(w)
    out = []
    for ws in L.values():
        ws = sorted(ws, key=lambda w: w["x0"])
        # split a line where the gap is large (databar fields are spaced far apart)
        grp = [ws[0]]
        for w in ws[1:]:
            hgt = max(1.0, w["y1"] - w["y0"])
            if w["x0"] - grp[-1]["x1"] > 2.5 * hgt:
                out.append(grp); grp = [w]
            else:
                grp.append(w)
        out.append(grp)
    res = []
    for g in out:
        res.append({"text": " ".join(w["text"] for w in g), "x0": min(w["x0"] for w in g), "x1": max(w["x1"] for w in g),
                    "y0": min(w["y0"] for w in g), "y1": max(w["y1"] for w in g)})
    return res


_EXCL = re.compile(r"(WD|W\.D\.|HFW|FWD|FOV|Width|Spot|Dwell|HV|EHT|kV|Det|Probe|Tilt|Pressure)\s*[=:]?\s*$", re.I)


def _fix_ocr(t: str) -> str:
    """Common OCR slips in a scale label: l/I/| for 1, O for 0, 'urn' for 'um'."""
    t = re.sub(r"(?<![A-Za-z])[lI|](?=\s*(?:\d|[.,]?\d*\s*(?:[µu]m|nm)))", "1", t)
    for _ in range(3):      # every O of a run after a digit: "1OO" -> "10O" -> "100" (one pass only fixed the first O)
        t2 = re.sub(r"(?<=\d)[oO]|[oO](?=[.,]?\d)", "0", t)
        if t2 == t:
            break
        t = t2
    t = t.replace("urn", "um").replace("pum", "µm")
    return re.sub(r"(?<=\d)(\s*)(?:[yjqµ]u?|u)m(?![a-zA-Z])", r"\1µm", t)      # µ read as y / j / q


def _excluded(text: str, p) -> bool:
    """'WD 8.0mm', 'HFW 2.98 µm' etc. are not scale-bar labels."""
    m = _LABEL.search(_fix_ocr(text).replace("μ", "µ"))
    return bool(m and _EXCL.search(_fix_ocr(text)[:m.start()]))


def _near_bar(box, bar, h):
    if bar is None:
        return True
    L = bar["x1"] - bar["x0"]
    th = max(box["y1"] - box["y0"], bar.get("thickness_px", 4), 8)
    xin = box["x1"] >= bar["x0"] - 0.6 * L - 60 and box["x0"] <= bar["x1"] + 0.6 * L + 60
    cy = (box["y0"] + box["y1"]) / 2
    yin = (bar["y0"] - 4 * th - 10) <= cy <= (bar["y1"] + 4 * th + 10)
    return xin and yin


def _is_standard(nm) -> bool:
    """A scale-bar length of 1, 2 or 5 × 10^n (nm): the lengths instruments draw bars with (400 nm is not one)."""
    if not nm or nm <= 0:
        return False
    m = nm / 10 ** math.floor(math.log10(nm))
    return any(abs(m - s) < 1e-9 for s in (1.0, 2.0, 5.0))


def _pick_vote(votes: dict):
    """The scale label with the most votes over the OCR passes. A tie goes to a standard bar length, then to the reading
    that came first. votes: {(value, unit, nm): [count, text, boxes]}. Returns (label, its entry, the other readings)."""
    ranked = sorted(votes.items(), key=lambda kv: (-kv[1][0], not _is_standard(kv[0][2])))
    return ranked[0][0], ranked[0][1], [k for k, _ in ranked[1:]]


def read_strip(img8: np.ndarray, strip_row=None, bar=None):
    """OCR the info strip: the scale-bar label nearest to the bar and the magnification. Returns dict."""
    h, w = img8.shape
    top = strip_row if strip_row is not None else (max(0, bar["y0"] - 60) if bar else int(h * 0.75))
    region = img8[top:]
    words = ocr_words(region, 0, top)
    lines = _lines(words) if words else []
    res = {"ocr": bool(words) or ocr_available(), "text": " | ".join(l["text"] for l in lines), "label": None, "mag": None}
    # scale label: parse every line (and every word pair), prefer the one closest to the bar centre
    cands = []
    for l in lines:
        p = parse_scale_label(_fix_ocr(l["text"]))
        if p and not _excluded(l["text"], p) and _near_bar(l, bar, h):
            cx, cy = (l["x0"] + l["x1"]) / 2, (l["y0"] + l["y1"]) / 2
            d = 0.0
            if bar:
                bx, by = (bar["x0"] + bar["x1"]) / 2, (bar["y0"] + bar["y1"]) / 2
                d = np.hypot((cx - bx) / w, (cy - by) / h)
            cands.append((d, p, l))
    if cands:
        d, p, l = min(cands, key=lambda c: c[0])
        res["label"] = {"value": p[0], "unit": p[1], "nm": p[2], "text": l["text"], "box": (l["x0"], l["y0"], l["x1"], l["y1"])}
    first = res["label"]
    if bar is not None:
        # targeted second pass: single text line just below / above / right of the bar
        L, th = bar["x1"] - bar["x0"], max(bar.get("thickness_px", 4), 3)
        hh = int(max(14, 6 * th, 0.035 * h))
        xa, xb = int(max(0, bar["x0"] - 0.25 * L)), int(min(w, bar["x1"] + 0.25 * L))
        boxes = [(bar["y1"] + 1, bar["y1"] + 1 + int(1.5 * hh), xa, xb),       # room for descenders (µ) (bar["y0"] - hh, bar["y0"], xa, xb),
                 (bar["y0"] - hh // 2, bar["y1"] + hh // 2, bar["x1"] + 1, int(min(w, bar["x1"] + 0.6 * L + 40))),
                 (bar["y0"] - hh // 2, bar["y1"] + hh // 2, int(max(0, bar["x0"] - 0.6 * L - 40)), bar["x0"])]
        for ya, yb, x0, x1 in boxes:
            ya, yb = max(0, ya), min(h, yb)
            if yb - ya < 6 or x1 - x0 < 10:
                continue
            votes = {}
            if first and _near_bar({"x0": first["box"][0], "y0": first["box"][1], "x1": first["box"][2], "y1": first["box"][3]}, bar, h) \
                    and x0 <= (first["box"][0] + first["box"][2]) / 2 <= x1 and ya <= (first["box"][1] + first["box"][3]) / 2 <= yb:
                votes[(first["value"], first["unit"], first["nm"])] = [1, first["text"], [{"x0": first["box"][0], "y0": first["box"][1], "x1": first["box"][2], "y1": first["box"][3]}]]
            ink = float(np.median(img8[bar["y0"]:bar["y1"] + 1, bar["x0"]:bar["x1"]]))
            for sc, psm, wl, ik in [(2, 7, LABEL_CHARS, None), (2, 11, None, None), (2, 7, LABEL_CHARS, ink), (2, 11, None, ink),
                                    (1, 7, LABEL_CHARS, None), (1, 11, None, None), (3, 7, LABEL_CHARS, ink),
                                    (3, 7, LABEL_CHARS, None), (3, 11, None, None)]:
                ws = ocr_words(img8[ya:yb, x0:x1], x0, ya, psm=psm, scale=sc, whitelist=wl, ink=ik)
                for i in range(len(ws)):          # the label is one or two words ('1', 'µm')
                    txt = " ".join(t["text"] for t in ws[i:i + 2])
                    q = parse_scale_label(_fix_ocr(txt))
                    if q and not _excluded(" ".join(t["text"] for t in ws[:i + 1]), q):
                        sel = ws[i:i + 2] if parse_scale_label(_fix_ocr(ws[i]["text"])) is None else ws[i:i + 1]
                        votes.setdefault(q, [0, " ".join(t["text"] for t in sel), sel])[0] += 1
                        break
                if votes and max(v[0] for v in votes.values()) >= 2:
                    break
            p = None
            if votes:
                p, (n, txt, ws), others = _pick_vote(votes)
                if others:       # the passes disagree (e.g. 100 vs 400 nm): kept for the notes in calibrate()
                    res["label_conflict"] = [f"{k[0]:g} {k[1]} ({votes[k][0]}×)" for k in [p] + others]
            if p:
                res["label"] = {"value": p[0], "unit": p[1], "nm": p[2], "text": txt,
                                "box": (min(t["x0"] for t in ws), min(t["y0"] for t in ws), max(t["x1"] for t in ws), max(t["y1"] for t in ws))}
                res["text"] += " | " + txt
                break
    mags = [(match_magnification(l["text"]), l) for l in lines]
    mags = [(m, l) for m, l in mags if m]
    if mags:
        # 'just above the scale bar' on the user's instrument: prefer the one closest to the bar, else the first
        if bar:
            bx = (bar["x0"] + bar["x1"]) / 2
            mags.sort(key=lambda ml: abs((ml[1]["x0"] + ml[1]["x1"]) / 2 - bx) / w + abs(ml[1]["y1"] - bar["y0"]) / h)
        (m, mt), l = mags[0]
        res["mag"] = {"value": m, "text": mt, "line": l["text"], "box": (l["x0"], l["y0"], l["x1"], l["y1"])}
    return res


# ---------------------------------------------------------------------------------------------- one call
def calibrate(img8: np.ndarray, meta_pixel_nm=None, meta_mag=None, meta_source=None, ref_width_mm=127.0):
    """Run every automatic method and return
    dict(method, nm_per_px, strip_row, bar, label, mag_ocr, mag_meta, candidates=[...], notes=[...])."""
    h, w = img8.shape
    strip = find_info_strip(img8)
    bar = detect_bar(img8, strip)
    if bar is None and strip is not None:
        bar = detect_bar(img8, None)
    ocr = read_strip(img8, strip, bar)
    cands, notes = [], []
    if meta_pixel_nm:
        cands.append({"method": f"metadata ({meta_source or 'file'})", "nm_per_px": float(meta_pixel_nm)})
    mag = meta_mag or (ocr["mag"]["value"] if ocr["mag"] else None)
    if bar and ocr["label"] and ocr["label"]["unit"] == "pm":
        # OCR reads 'µm' as 'pm' very often; SEM info strips never carry pm bars. Decide with the magnification if any.
        v = ocr["label"]["value"]
        as_um, as_pm = v * 1e3 / bar["length_px"], v * 1e-3 / bar["length_px"]
        ref = ref_width_mm * 1e6 / (mag * w) if mag else None
        if ref is None or abs(np.log(as_um / ref)) < abs(np.log(as_pm / ref)):
            ocr["label"].update(unit="µm", nm=v * 1e3)
            notes.append(("info", f"OCR read '{ocr['label']['text']}' – interpreted as {v:g} µm (µ is often read as p)."))
    if ocr["label"]:
        ocr["label"]["shown"] = f"{ocr['label']['value']:g} {ocr['label']['unit']}"
    if ocr["label"] and not _is_standard(ocr["label"]["nm"]) and not meta_pixel_nm:
        # the value is used as read; a non-standard bar length (not 1, 2 or 5 x 10^n nm) is flagged for a check by the user
        notes.append(("warn", f"The OCR label '{ocr['label']['shown']}' is not a standard bar length (1, 2 or 5 x 10^n). "
                      "It was used as read and there is no file metadata to cross-check it. Check it against the printed bar."))
    if ocr.get("label_conflict"):
        notes.append(("warn", "The OCR passes read different scale labels: " + ", ".join(ocr["label_conflict"])
                      + f". Used {ocr['label']['shown']} (most votes; a tie goes to a standard bar length). "
                      "Check the printed label against the bar, or give the scale by hand."))
    if bar and ocr["label"]:
        cands.append({"method": "scale bar (detected) + OCR label", "nm_per_px": ocr["label"]["nm"] / bar["length_px"],
                      "bar_px": bar["length_px"], "value": ocr["label"]["value"], "unit": ocr["label"]["unit"]})
    if mag:
        cands.append({"method": f"magnification {'(metadata)' if meta_mag else '(OCR)'} × {ref_width_mm:g} mm display",
                      "nm_per_px": ref_width_mm * 1e6 / (mag * w)})
    if meta_mag and ocr["mag"]:
        r = ocr["mag"]["value"] / meta_mag
        notes.append(("ok" if abs(r - 1) < 0.02 else "warn",
                      f"Magnification: metadata {meta_mag:g}× vs printed '{ocr['mag']['text']}' → {ocr['mag']['value']:g}× "
                      + ("(agree)" if abs(r - 1) < 0.02 else f"(differ by {100 * (r - 1):+.1f} %)")))
    for c in cands[1:]:
        a, b = cands[0]["nm_per_px"], c["nm_per_px"]
        dv = 100 * (b - a) / a
        if c["method"].startswith("scale bar") and abs(dv) >= 2:
            # the OCR label disagrees with the metadata pixel size: say which numbers disagree; nothing is changed automatically
            notes.append(("warn", f"The OCR scale label '{ocr['label']['shown']}' does not match the metadata pixel size "
                          f"{a:.4g} nm/px. The detected scale bar is {c['bar_px']:.0f} px, which gives {b:.4g} nm/px "
                          f"({dv:+.1f} %). The label and the metadata disagree; neither was changed automatically. "
                          "Check the printed label against the bar."))
            continue
        is_mag = c["method"].startswith("magnification")
        lvl = "ok" if abs(dv) < 2 else ("info" if is_mag else "warn")
        notes.append((lvl, f"Cross-check – {c['method']}: {b:.4g} nm/px ({dv:+.1f} % vs {cands[0]['method']})"
                      + (" – the 127 mm display width is the Philips/FEI XL reference; other instruments differ, so only the "
                         "scale bar or metadata is authoritative." if is_mag and abs(dv) >= 2 else "")))
    best = cands[0] if cands else None
    return {"method": best["method"] if best else None, "nm_per_px": best["nm_per_px"] if best else None,
            "strip_row": strip, "bar": bar, "label": ocr["label"], "mag_ocr": ocr["mag"], "mag_meta": meta_mag,
            "ocr_text": ocr["text"], "ocr_available": ocr_available(), "candidates": cands, "notes": notes}


def nm_per_px_from_points(x0, x1, value, unit, y0=None, y1=None):
    """Manual calibration: two clicked bar ends (image px) + the printed value -> nm/px."""
    L = float(np.hypot(x1 - x0, (y1 or 0) - (y0 or 0)))
    f = UNIT_NM.get(unit.lower().replace("μ", "µ") if unit != "Å" else "å")
    return (value * f / L) if L > 0 and f else None


def stale_bar_notes(info: dict):
    """info: {image name: dict(meta=nm/px from metadata, bar=nm/px from the burned-in bar, label, bar_px)}.
    A burned-in bar that disagrees with its own metadata but equals ANOTHER image's calibration was drawn with a stale
    magnification setting in the acquisition software (seen with Olympus analySIS on a JEOL JEM-1010)."""
    out = []
    for a, ia in info.items():
        if not ia or not ia.get("meta") or not ia.get("bar") or abs(ia["bar"] / ia["meta"] - 1) < 0.03:
            continue
        hits = [b for b, ib in info.items() if b != a and ib and ib.get("meta") and abs(ia["bar"] / ib["meta"] - 1) < 0.005]
        msg = (f"**{a}**: the burned-in {ia.get('label')} bar ({ia.get('bar_px') or 0:.0f} px → {ia['bar']:.4g} nm/px) disagrees with "
               f"the file's calibration ({ia['meta']:.4g} nm/px, {100 * (ia['bar'] / ia['meta'] - 1):+.0f} %).")
        if hits:
            msg += f" It matches the calibration of **{', '.join(hits)}** exactly, so the bar was drawn with a stale magnification setting."
        out.append(msg)
    return out


def ocr_missing_message():
    """v2.19.2 (item 10): the text shown when Tesseract is not found, with the install command for each system.
    Returns None when OCR can be used. Scale-bar OCR itself is not changed (excluded from the gate)."""
    if tesseract_cmd() is not None:
        return None
    return ("Tesseract OCR was not found, so the scale bar cannot be read automatically. Enter the scale by hand, or install "
            "Tesseract: Windows - the app ships a copy, or install from github.com/UB-Mannheim/tesseract; Debian/Ubuntu - "
            "sudo apt install tesseract-ocr; macOS - brew install tesseract. You can also point $MFS_TESSERACT at the executable.")
