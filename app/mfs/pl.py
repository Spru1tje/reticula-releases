"""v2.14: photoluminescence - steady-state spectra and temperature-dependent PL (thermal quenching, ratiometric thermometry).

Formulas (the usual ones of the luminescence-thermometry literature, e.g. Brites, Balabhadra, Carlos, Adv. Opt. Mater. 2019):
  thermal quenching     I(T) = I0 / (1 + A exp(-Ea / kB T))            Ea activation energy, T50 = Ea / (kB ln A) (half of I0)
  Arrhenius plot        ln(I0/I - 1) against 1/T is a straight line of slope -Ea/kB
  ratiometric readout   D = I1 / I2 = alpha exp(-dE / kB T)            two thermally coupled levels (e.g. Er3+ 2H11/2 and 4S3/2)
  relative sensitivity  Sr = 100 % x |(1/D) dD/dT| = 100 % x dE / (kB T^2)
  temperature uncertainty (photon counting)   sigma_T = (1/Sr) sqrt(1/I1 + 1/I2)   (I1, I2 in counts)
Only numpy / scipy; the drawing functions take an axes so that every panel can also go into the combined figure."""
from __future__ import annotations

import io
import re

import numpy as np

KB_EV = 8.617333262e-5          # eV / K
EV_TO_CM = 8065.54429           # cm-1 per eV
HC = 1239.84198                 # eV nm

# ----------------------------------------------------------------------------------------------- reading
_NUM = re.compile(r"^[-+]?(\d+([.,]\d*)?|[.,]\d+)([eE][-+]?\d+)?$")


def _is_num(s):
    return bool(_NUM.match(str(s).strip()))


def _f(s):
    return float(str(s).strip().replace(",", "."))


def parse_temperature(text):
    """(value, unit) from a header / file name such as '293.15 K', '25C', '50 °C', 'T=300K', or a plain number (unit None).
    Returns (None, None) when there is no temperature in the text."""
    t = str(text)
    m = re.search(r"(?<![A-Za-z0-9.])([-+]?\d+(?:[.,]\d+)?)\s*(?:°|º|deg\s*)?\s*C(?![A-Za-z])", t)
    if m:
        return _f(m.group(1)), "C"
    m = re.search(r"(?<![A-Za-z0-9.])(\d+(?:[.,]\d+)?)\s*K(?![A-Za-z])", t)
    if m:
        return _f(m.group(1)), "K"
    if _is_num(t.strip()):
        return _f(t), None
    return None, None


def to_kelvin(values_units, default="auto"):
    """[(value, unit)] -> list of K (None where unknown). Plain numbers follow default: 'K', 'C' or 'auto' (all > 150 -> K, else degrees C)."""
    plain = [v for v, u in values_units if v is not None and u is None]
    if default == "auto":
        default = "K" if (plain and max(plain) > 150) else "C"
    out = []
    for v, u in values_units:
        if v is None:
            out.append(None)
            continue
        u = u or default
        out.append(v + 273.15 if u == "C" else v)
    return out


def _rows_from_bytes(raw):
    if raw[:2] == b"PK" or raw[:4] == b"\xd0\xcf\x11\xe0":
        import pandas as pd
        df = pd.read_excel(io.BytesIO(raw), sheet_name=0, header=None)
        return [["" if (isinstance(c, float) and np.isnan(c)) else str(c) for c in row] for row in df.values.tolist()]
    for enc in ("utf-8-sig", "utf-16", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except Exception:
            continue
    lines = [l for l in text.replace("\r", "").split("\n") if l.strip()]
    cnt = {d: sum(l.count(d) for l in lines[-60:]) for d in ("\t", ";", ",")}
    d = max(cnt, key=cnt.get) if max(cnt.values()) else None
    return [[c.strip() for c in (l.split(d) if d else l.split())] for l in lines]


def read_pl(src, name=""):
    """Reads a spectrum file: two columns (wavelength, intensity) or a table (first column wavelength, one column per spectrum /
    temperature; the header gives the names). CSV / TXT / DAT / Excel. Returns a list of dicts {name, T_text, wl, I}."""
    if isinstance(src, (bytes, bytearray)):
        raw = bytes(src)
    elif hasattr(src, "read"):
        raw = src.read()
        try:
            src.seek(0)
        except Exception:
            pass
        name = name or getattr(src, "name", "")
    else:
        with open(src, "rb") as fh:
            raw = fh.read()
        name = name or str(src)
    rows = _rows_from_bytes(raw)
    data_idx = [i for i, r in enumerate(rows) if len(r) >= 2 and _is_num(r[0]) and _is_num(r[1])]
    if len(data_idx) < 5:
        raise ValueError("No spectrum found: expected a column of wavelengths (nm) followed by one or more intensity columns.")
    lens = [len(rows[i]) for i in data_idx]
    ncol = max(set(lens), key=lens.count)
    use = [i for i in data_idx if len(rows[i]) == ncol and all(_is_num(c) for c in rows[i][:ncol])]
    arr = np.array([[_f(c) for c in rows[i][:ncol]] for i in use], float)
    header = None
    for j in range(use[0] - 1, max(-1, use[0] - 8), -1):                 # the closest earlier row with text in it
        r = rows[j]
        if len([c for c in r if c and not _is_num(c)]) >= max(1, min(2, ncol - 1)):
            header = r
            break
    stem = re.sub(r"\.[A-Za-z0-9]{2,5}$", "", re.split(r"[\\/]", name)[-1]) if name else "spectrum"
    wl = arr[:, 0]
    order = np.argsort(wl)
    out = []
    for k in range(1, ncol):
        if header is not None and len(header) >= ncol and header[k]:
            label = header[k]
        elif ncol == 2:
            label = stem
        else:
            label = f"{stem} #{k}"
        out.append({"name": label, "T_text": label if (header is not None and len(header) >= ncol and header[k]) else stem,
                    "wl": wl[order], "I": arr[order, k]})
    return out


# ----------------------------------------------------------------------------------------------- processing
def smooth(y, window):
    y = np.asarray(y, float)
    if window and window >= 5 and len(y) > window:
        from scipy.signal import savgol_filter
        w = int(window) | 1
        return savgol_filter(y, w, 3)
    return y


def baseline_subtract(wl, y, mode="none", region=None):
    y = np.asarray(y, float)
    if mode == "minimum":
        return y - np.min(y)
    if mode == "region" and region:
        m = (wl >= min(region)) & (wl <= max(region))
        return y - (np.mean(y[m]) if m.any() else 0.0)
    return y


def to_axis(wl, y, axis="nm"):
    """Spectrum on the wavelength (nm), energy (eV) or wavenumber (cm-1) axis. Energy / wavenumber spectra get the Jacobian
    (I(E) = I(lambda) lambda^2 / hc) so that the integrated intensity is conserved. Returns x ascending, y."""
    wl = np.asarray(wl, float); y = np.asarray(y, float)
    if axis == "nm":
        return wl, y
    ok = wl > 0
    wl, y = wl[ok], y[ok]
    if axis == "eV":
        x, yy = HC / wl, y * wl ** 2 / HC
    elif axis == "cm-1":
        x, yy = 1e7 / wl, y * wl ** 2 / 1e7
    else:
        raise ValueError(axis)
    o = np.argsort(x)
    return x[o], yy[o]


def normalise(y, mode="none", ref=None):
    y = np.asarray(y, float)
    if mode == "max":
        return y / (np.max(y) or 1.0)
    if mode == "area":
        return y / (abs(np.trapezoid(y)) or 1.0)
    if mode == "ref" and ref:
        return y / ref
    return y


def band_integral(wl, y, lo, hi, baseline="linear"):
    """Area between lo and hi (nm) - optionally above the straight line through the window ends."""
    wl = np.asarray(wl, float); y = np.asarray(y, float)
    m = (wl >= min(lo, hi)) & (wl <= max(lo, hi))
    if m.sum() < 3:
        return 0.0
    x, v = wl[m], y[m]
    if baseline == "linear":
        v = v - np.interp(x, [x[0], x[-1]], [v[0], v[-1]])
    return float(np.trapezoid(v, x))


def peak_stats(wl, y, lo=None, hi=None):
    """Position (parabolic), height, FWHM (half-maximum crossings above the window minimum) and area of the strongest peak in lo..hi."""
    wl = np.asarray(wl, float); y = np.asarray(y, float)
    m = np.ones_like(wl, bool) if lo is None else (wl >= min(lo, hi)) & (wl <= max(lo, hi))
    if m.sum() < 5:
        return {"pos": np.nan, "height": np.nan, "fwhm": np.nan, "area": 0.0}
    x, v = wl[m], y[m]
    v = smooth(v, max(5, (len(v) // 40) | 1))                 # a broad, noisy band: locate it on a lightly smoothed copy
    i = int(np.argmax(v))
    pos = float(x[i])
    if 0 < i < len(v) - 1:
        a, b, c = v[i - 1], v[i], v[i + 1]
        den = a - 2 * b + c
        if den:
            pos = float(x[i] + 0.5 * (a - c) / den * (x[i + 1] - x[i]))
    base = float(np.min(v))
    h = float(v[i] - base)
    half = base + h / 2

    def cross(idx_range):
        prev = i
        for j in idx_range:
            if v[j] <= half:
                t = (v[prev] - half) / ((v[prev] - v[j]) or 1e-12)
                return float(x[prev] + t * (x[j] - x[prev]))
            prev = j
        return None
    xl, xr = cross(range(i - 1, -1, -1)), cross(range(i + 1, len(v)))
    fwhm = float(xr - xl) if (xl is not None and xr is not None) else np.nan
    return {"pos": pos, "height": float(v[i]), "fwhm": fwhm, "area": band_integral(x, v, x[0], x[-1], "linear")}


# ----------------------------------------------------------------------------------------------- fits
def _r2(y, f):
    y = np.asarray(y, float); f = np.asarray(f, float)
    ss = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - float(np.sum((y - f) ** 2)) / ss if ss > 0 else 1.0


def quench_model(T, I0, lnA, Ea):
    return I0 / (1.0 + np.exp(lnA - Ea / (KB_EV * np.asarray(T, float))))


def fit_quench(T, I):
    """Single-channel thermal quenching. Returns dict(I0, A, Ea_eV, Ea_meV, T50, R2, status, f) (f = callable) or status only."""
    from scipy.optimize import curve_fit
    T = np.asarray(T, float); I = np.asarray(I, float)
    o = np.argsort(T); T, I = T[o], I[o]
    if len(T) < 4 or np.max(I) <= 0:
        return {"status": "needs at least 4 temperatures"}
    Imax = float(np.max(I))
    drop = float(I[-1] / Imax) if I[np.argmax(I)] != I[-1] else 1.0
    if min(I[-3:].mean() / Imax, 1.0) > 0.97:
        return {"status": "the intensity does not quench in this temperature range"}
    sel = (I < 0.97 * Imax) & (I > 0.01 * Imax)
    p0 = [Imax * 1.02, 8.0, 0.15]
    if sel.sum() >= 3:
        z = np.log(np.clip(Imax * 1.02 / I[sel] - 1, 1e-9, None))
        slope, icpt = np.polyfit(1.0 / T[sel], z, 1)
        p0 = [Imax * 1.02, float(np.clip(icpt, -10, 40)), float(np.clip(-slope * KB_EV, 0.005, 2.5))]
    try:
        p, _ = curve_fit(quench_model, T, I, p0=p0, bounds=([0.5 * Imax, -15, 0.0], [3 * Imax, 45, 3.0]), maxfev=20000)
    except Exception as e:
        return {"status": f"fit failed: {e}"}
    I0, lnA, Ea = [float(v) for v in p]
    f = lambda t: quench_model(t, I0, lnA, Ea)
    return {"status": "ok", "I0": I0, "A": float(np.exp(lnA)), "lnA": lnA, "Ea_eV": Ea, "Ea_meV": 1000 * Ea,
            "T50": (Ea / (KB_EV * lnA) if lnA > 0 else None), "R2": _r2(I, f(T)), "f": f}


def fit_boltzmann(T, D, offset=False):
    """D = alpha exp(-dE / kB T) (+ C). Returns dict(alpha, dE_eV, dE_cm1, C, R2, f)."""
    from scipy.optimize import curve_fit
    T = np.asarray(T, float); D = np.asarray(D, float)
    ok = np.isfinite(T) & np.isfinite(D) & (D > 0)
    T, D = T[ok], D[ok]
    if len(T) < 3:
        return {"status": "needs at least 3 temperatures"}
    slope, icpt = np.polyfit(1.0 / T, np.log(D), 1)
    p0 = [float(np.exp(icpt)), float(-slope * KB_EV)] + ([0.0] if offset else [])
    mdl = (lambda t, a, de, c: a * np.exp(-de / (KB_EV * t)) + c) if offset else (lambda t, a, de: a * np.exp(-de / (KB_EV * t)))
    try:
        p, _ = curve_fit(mdl, T, D, p0=p0, maxfev=20000)
    except Exception as e:
        return {"status": f"fit failed: {e}"}
    f = (lambda t: mdl(np.asarray(t, float), *p))
    return {"status": "ok", "alpha": float(p[0]), "dE_eV": float(p[1]), "dE_cm1": float(p[1] * EV_TO_CM),
            "C": float(p[2]) if offset else 0.0, "R2": _r2(D, f(T)), "f": f}


def relative_sensitivity(T, dE_eV):
    """Sr in % per K of a Boltzmann thermometer."""
    T = np.asarray(T, float)
    return 100.0 * dE_eV / (KB_EV * T ** 2)


def absolute_sensitivity(T, fit):
    """Sa = |dD/dT| of the fitted ratio (per K); for D = alpha exp(-dE/kT) + C it is (D - C) dE / (kB T^2)."""
    T = np.asarray(T, float)
    return (fit["f"](T) - fit.get("C", 0.0)) * fit["dE_eV"] / (KB_EV * T ** 2)


def temperature_uncertainty(sr_pct, I1, I2):
    """sigma_T in K for photon-counting statistics (I1, I2 in counts)."""
    sr = np.asarray(sr_pct, float) / 100.0
    return np.sqrt(1.0 / np.clip(I1, 1e-12, None) + 1.0 / np.clip(I2, 1e-12, None)) / sr


# ----------------------------------------------------------------------------------------------- drawing
CMAPS = ["viridis", "plasma", "cividis", "magma", "coolwarm", "turbo"]


def _tcolors(T, cmap):
    import matplotlib as mpl
    from matplotlib.colors import Normalize
    T = np.asarray(T, float)
    cm = mpl.colormaps[cmap]
    lo, hi = (float(np.nanmin(T)), float(np.nanmax(T))) if np.isfinite(T).any() else (0.0, 1.0)
    norm = Normalize(lo, hi if hi > lo else lo + 1)
    return [cm(norm(t)) if np.isfinite(t) else (0.3, 0.3, 0.3, 1) for t in T], cm, norm


def draw_spectra(ax, spectra, xlabel="Wavelength (nm)", ylabel="Intensity (a.u.)", cmap="viridis", colorbar=True, lw=0.9, labels=None,
                 colours=None):
    """Overlaid spectra; coloured by temperature (with a colour bar in K) when every spectrum has one, else by the sample palette."""
    from . import style as S
    T = [s.get("T") for s in spectra]
    have_t = all(t is not None for t in T) and len(spectra) > 1
    if have_t:
        cols, cm, norm = _tcolors(T, cmap)
    else:
        cols = colours or S.colours_for([s["name"] for s in spectra])
    for c, s in zip(cols, spectra):
        ax.plot(s["x"], s["y"], color=c, lw=lw, label=None if have_t else s["name"])
    if have_t and colorbar:
        import matplotlib as mpl
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cm)
        cb = ax.figure.colorbar(sm, ax=ax, pad=0.02, fraction=0.05, aspect=25)
        cb.set_label("Temperature (K)")
        cb.ax.tick_params(direction="out")
    elif not have_t and len(spectra) > 1:
        ax.legend(loc="best")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if xlabel.startswith("Wavelength"):
        ax.set_xlim(min(s["x"].min() for s in spectra), max(s["x"].max() for s in spectra))
    return ax


def draw_ridge(ax, spectra, xlabel="Wavelength (nm)", ylabel="Intensity (a.u.) + offset", cmap="viridis", step=0.35, lw=0.9, label_every=1):
    """Stacked (waterfall) spectra, the lowest temperature at the bottom; each one hides what lies behind it. The temperature is
    written at the right end of every trace."""
    T = [s.get("T") for s in spectra]
    have_t = all(t is not None for t in T)
    order = np.argsort(T) if have_t else np.arange(len(spectra))
    cols, _, _ = _tcolors(T, cmap) if have_t else ([(0.1, 0.1, 0.1, 1)] * len(spectra), None, None)
    ymax = max(float(np.max(s["y"])) for s in spectra) or 1.0
    for rank, k in enumerate(order):
        s = spectra[k]
        off = rank * step * ymax
        z = 2 + (len(order) - rank) * 0.01                    # the lowest temperature is in FRONT and hides the traces behind it
        ax.fill_between(s["x"], off, s["y"] + off, color="white", zorder=z, lw=0)
        ax.plot(s["x"], s["y"] + off, color=cols[k], lw=lw, zorder=z + 0.005)
        if rank % max(1, label_every) == 0:
            lab = f"{s['T']:.0f} K" if have_t else s["name"]
            ax.text(s["x"].max(), off + float(s["y"][-1]), " " + lab, color="#2B3440", fontsize=plt_fs(ax) - 1.5, va="bottom", ha="left", clip_on=False)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_yticks([])
    return ax


def plt_fs(ax):
    import matplotlib as mpl
    return float(mpl.rcParams.get("font.size", 8))


def draw_map(ax, spectra, xlabel="Wavelength (nm)", cmap="magma", zlabel="Intensity (a.u.)", ylabel="Temperature (K)"):
    """Emission map: colour = intensity, x = wavelength, y = temperature."""
    T = np.array([s["T"] for s in spectra], float)
    o = np.argsort(T)
    xs = spectra[0]["x"]
    Z = np.array([np.interp(xs, spectra[k]["x"], spectra[k]["y"]) for k in o])
    if len(T) >= 2:
        mesh = ax.pcolormesh(xs, T[o], Z, shading="auto", cmap=cmap, rasterized=True)
    else:
        mesh = ax.pcolormesh(xs, [T[0] - 0.5, T[0] + 0.5], np.vstack([Z, Z]), shading="auto", cmap=cmap)
    cb = ax.figure.colorbar(mesh, ax=ax, pad=0.02, fraction=0.05, aspect=25)
    cb.set_label(zlabel); cb.ax.tick_params(direction="out")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    return ax


def draw_trend(ax, T, y, ylabel, yerr=None, fit_T=None, fit_y=None, colour="#0072B2", fit_label=None, marker="o", legend=True):
    ax.errorbar(T, y, yerr=yerr, fmt=marker, ms=3.2, color=colour, ecolor=colour, elinewidth=0.6, capsize=1.8, label="measured", zorder=3)
    if fit_T is not None:
        ax.plot(fit_T, fit_y, color="#D55E00", lw=1.0, label=fit_label or "fit", zorder=2)
    ax.set_xlabel("Temperature (K)"); ax.set_ylabel(ylabel)
    if legend and fit_T is not None:
        up = len(y) > 1 and float(np.asarray(y, float)[-1]) > float(np.asarray(y, float)[0])
        ax.legend(loc="lower right" if up else "upper right")
    return ax


def draw_arrhenius(ax, T, I, fit, colour="#0072B2"):
    """ln(I0/I - 1) against 1000/T: a straight line with slope -Ea/kB for a single quenching channel."""
    T = np.asarray(T, float); I = np.asarray(I, float)
    z = fit["I0"] / np.clip(I, 1e-12, None) - 1
    ok = z > 0.05                                            # where the intensity has dropped, ln(I0/I - 1) is meaningful
    if ok.sum() < 3:
        ok = z > 1e-6
    ax.plot(1000.0 / T[ok], np.log(z[ok]), "o", ms=3.2, color=colour, label="measured")
    tt = np.linspace(T[ok].min(), T[ok].max(), 100)
    ax.plot(1000.0 / tt, fit["lnA"] - fit["Ea_eV"] / (KB_EV * tt), color="#D55E00", lw=1.0,
            label=f"$E_a$ = {fit['Ea_meV']:.0f} meV")
    ax.set_xlabel("1000 / T (K$^{-1}$)"); ax.set_ylabel(r"ln($I_0$/$I$ $-$ 1)")
    ax.legend(loc="best")
    return ax
