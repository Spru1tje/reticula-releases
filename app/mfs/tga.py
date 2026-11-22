"""TGA helpers: automatic mass-loss steps from DTG, step labels, TGA–MS (ion current) reader and overlay."""
from __future__ import annotations
import io, re
import numpy as np
import pandas as pd


def steps(a, min_loss=0.5, min_prom=0.03):
    """Mass-loss steps from an analysis.tga_analysis() dict: every significant DTG peak is a step; the boundaries are
    the DTG minima between neighbouring peaks (data start / end outside). Steps < min_loss %-points are merged into
    their larger neighbour. Returns list of dicts T_start, T_peak, T_end, w_start, w_end, loss (%-points)."""
    from scipy.signal import find_peaks
    T, w, d = np.asarray(a["T"]), np.asarray(a["w"]), np.asarray(a["dtg"])
    pk, _ = find_peaks(d, prominence=max(np.nanmax(d) * min_prom, 1e-5))
    pk = list(pk)
    while True:
        b = [0] + [p0 + int(np.argmin(d[p0:p1 + 1])) for p0, p1 in zip(pk[:-1], pk[1:])] + [len(T) - 1]
        losses = [w[b[i]] - w[b[i + 1]] for i in range(len(pk))]
        small = [i for i, l in enumerate(losses) if l < min_loss]
        if not small or len(pk) <= 1:
            break
        pk.pop(min(small, key=lambda i: losses[i]))
    return [{"T_start": float(T[b[i]]), "T_peak": float(T[p]), "T_end": float(T[b[i + 1]]), "w_start": float(w[b[i]]),
             "w_end": float(w[b[i + 1]]), "loss": float(w[b[i]] - w[b[i + 1]])} for i, p in enumerate(pk)
            if w[b[i]] - w[b[i + 1]] >= min_loss]


def label_steps(ax, a, st, color="0.2", fs=6, normalise_end=False, fmt="−{:.1f} %"):
    """Bracket + text for every step on the TGA axis (gid mfs_steplabel)."""
    f = (100 / a["W_end"]) if normalise_end else 1.0
    out = []
    for s in st:
        x = s["T_end"]; y1, y2 = s["w_start"] * f, s["w_end"] * f
        ax.plot([s["T_peak"], s["T_peak"]], [y1, y2], color=color, lw=0.5, ls=":")[0].set_gid("mfs_steplabel")
        t = ax.annotate(fmt.format(s["loss"] * f) + f"\n{s['T_peak']:.0f} °C", (s["T_peak"], 0.5 * (y1 + y2)), xytext=(4, 0),
                        textcoords="offset points", fontsize=fs, color=color, va="center")
        t.set_gid("mfs_steplabel"); out.append(t)
    return out


def read_ms(src, name=""):
    """Ion-current table: a temperature (or time) column + one column per m/z (headers like 'm/z 18', '18', 'QMID 44',
    'Ion current 44 (A)'). Returns DataFrame with T or t and columns 'm/z N'."""
    raw = src.read() if hasattr(src, "read") else open(src, "rb").read()
    txt = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw
    lines = [l for l in txt.splitlines() if l.strip() and not l.lstrip().startswith("#")]
    sep = ";" if lines[0].count(";") > lines[0].count(",") else ("\t" if "\t" in lines[0] else ",")
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, decimal="," if sep == ";" else ".")
    cols = {}
    for c in df.columns:
        s = str(c).lower()
        if re.search(r"temp|^t\s*\(°?c\)|°c", s):
            cols[c] = "T"
        elif re.search(r"time|^t\s*\(min|min\)", s):
            cols[c] = "t"
        else:
            m = re.search(r"(\d+(?:\.\d+)?)", s)
            if m:
                cols[c] = f"m/z {float(m.group(1)):g}"
    df = df.rename(columns=cols)[[v for v in cols.values()]]
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if "T" not in df and "t" not in df:
        raise ValueError("no temperature or time column found")
    return df


MS_LABELS = {"18": "H$_2$O", "44": "CO$_2$", "28": "CO / N$_2$", "64": "SO$_2$", "48": "SO", "30": "NO", "46": "HCOOH / NO$_2$",
             "32": "O$_2$", "17": "OH / NH$_3$", "12": "C", "45": "HCOO"}


def draw_ms(ax, ms, channels, colors=None, log=False, x="T", norm=True):
    cols = colors or ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
    for k, ch in enumerate(channels):
        y = ms[ch].values.astype(float)
        if norm:
            y = (y - np.nanmin(y)) / ((np.nanmax(y) - np.nanmin(y)) or 1)
        lab = ch + (f" ({MS_LABELS[ch.split()[-1]]})" if ch.split()[-1] in MS_LABELS else "")
        ax.plot(ms[x], y + (0 if log else 1.1 * k * norm), color=cols[k % len(cols)], lw=0.8, label=lab)
    if log:
        ax.set_yscale("log")
    ax.set_ylabel("Ion current (norm.)" if norm else "Ion current (A)")
    ax.legend(frameon=False, fontsize=6, loc="upper left")


def step_uncertainty(m0_mg: float, resolution_mg: float = 0.001) -> dict:
    """v2.19.2: uncertainty of a mass-loss step from the balance resolution (v2.19.2 item 8).

    A step is the difference of two readings (before and after the loss), so the mass-loss uncertainty is sqrt(2) times the
    resolution. Expressed in % of the initial mass m0 (mg). Known answer: m0 = 10 mg, resolution 0.001 mg -> sqrt(2) × 0.01 % = 0.0141 %.
    General practice: the balance resolution is only one part of the error; buoyancy, purge-gas flow and the sample weighing
    contribute too. Quote this value as the lower limit of the uncertainty."""
    m0 = float(m0_mg)
    res = float(resolution_mg)
    if not (m0 > 0 and res >= 0):
        return {"status": "needs a positive initial mass and a non-negative resolution"}
    sigma_mg = (2.0 ** 0.5) * res
    return {"status": "ok", "sigma_mg": sigma_mg, "sigma_pct": 100.0 * sigma_mg / m0,
            "note": "balance resolution only (lower limit); see General practice in the help document"}


def step_with_uncertainty(step_pct: float, m0_mg: float, resolution_mg: float = 0.001) -> dict:
    """A step in % with its uncertainty from the balance resolution, and the mass lost in mg (step_pct × m0 / 100)."""
    u = step_uncertainty(m0_mg, resolution_mg)
    if u.get("status") != "ok":
        return u
    return {"status": "ok", "step_pct": float(step_pct), "step_err_pct": u["sigma_pct"],
            "mass_loss_mg": float(step_pct) * float(m0_mg) / 100.0, "mass_loss_err_mg": u["sigma_mg"]}
