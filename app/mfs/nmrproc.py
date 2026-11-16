"""v2.13 Bruker NMR raw import + processing (1H / 13C / DEPT; 1D only).

Reads Bruker experiment folders written by XWIN-NMR / TopSpin 1.x-4.x (Avance I/II/III/NEO):
fid (int32 DTYPA=0 or float64 DTYPA=2, either byte order), acqus, optional pdata/<n>/procs + 1r/1i.
Processing: digital-filter removal (GRPDLY or DECIM/DSPFVS table), exponential apodisation, zero-filling, FFT,
automatic phase (ACME-like entropy + negative-penalty, coarse grid then simplex), automatic baseline (signal-free point
detection + Whittaker smoother) and automatic referencing (TMS / DSS or the solvent residual signal).
Pure numpy/scipy + nmrglue; pathlib-safe; no OS-specific calls."""
from __future__ import annotations

import io
import os
import re
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# residual solvent signals (Fulmer et al., Organometallics 2010, 29, 2176) — 1H and 13C reference values (ppm)
SOLVENTS = {
    "CDCl3": {"1H": 7.26, "13C": 77.16, "alias": ["cdcl3", "chloroform"]},
    "DMSO-d6": {"1H": 2.50, "13C": 39.52, "alias": ["dmso", "dmso-d6", "dmsod6"]},
    "D2O": {"1H": 4.79, "13C": None, "alias": ["d2o", "h2o+d2o", "h2o_d2o"]},
    "CD3OD": {"1H": 3.31, "13C": 49.00, "alias": ["meod", "cd3od", "methanol-d4", "meoh"]},
    "Acetone-d6": {"1H": 2.05, "13C": 29.84, "alias": ["acetone", "acetone-d6", "aceton"]},
    "C6D6": {"1H": 7.16, "13C": 128.06, "alias": ["c6d6", "benzene"]},
    "CD3CN": {"1H": 1.94, "13C": 1.32, "alias": ["cd3cn", "acetonitrile", "acn"]},
    "CD2Cl2": {"1H": 5.32, "13C": 53.84, "alias": ["cd2cl2", "dcm"]},
    "THF-d8": {"1H": 3.58, "13C": 67.21, "alias": ["thf"]},
    "Toluene-d8": {"1H": 2.08, "13C": 20.43, "alias": ["toluene", "tol"]},
    "TFA-d": {"1H": 11.50, "13C": 164.2, "alias": ["tfa"]},
}


def solvent_key(name: str | None):
    """Bruker SOLVENT string (e.g. 'CDCl3', 'DMSO', 'MeOD', 'H2O+D2O') → key of SOLVENTS or None."""
    s = (name or "").strip().lower().replace(" ", "")
    if not s:
        return None
    for k, v in SOLVENTS.items():
        if s == k.lower() or any(s.startswith(a) for a in v["alias"]):
            return k
    return None


# DSPFVS/DECIM → group delay (points) for digitally filtered data recorded without GRPDLY (Westler & Abildgaard table)
_DSP = {10: {2: 44.75, 3: 33.5, 4: 66.625, 6: 59.0833, 8: 68.5625, 12: 60.375, 16: 69.5313, 24: 61.0208, 32: 70.0156,
             48: 61.3438, 64: 70.2578, 96: 61.5052, 128: 70.3789, 192: 61.5859, 256: 70.4395, 384: 61.6263, 512: 70.4697,
             768: 61.6465, 1024: 70.4849, 1536: 61.6566, 2048: 70.4924},
        11: {2: 46.0, 3: 36.5, 4: 48.0, 6: 50.1667, 8: 53.25, 12: 69.5, 16: 72.25, 24: 70.1667, 32: 72.75, 48: 70.5,
             64: 73.0, 96: 70.6667, 128: 72.5, 192: 71.3333, 256: 72.25, 384: 71.6667, 512: 72.125, 768: 71.8333,
             1024: 72.0625, 1536: 71.9167, 2048: 72.0313},
        12: {2: 46.311, 3: 36.53, 4: 47.87, 6: 50.229, 8: 53.289, 12: 69.551, 16: 71.6, 24: 70.184, 32: 72.138,
             48: 70.528, 64: 72.348, 96: 70.7, 128: 72.524}}


def parse_jcamp_params(text: str) -> dict:
    """Bruker acqus/procs (JCAMP-DX style) → dict; arrays (##$P= (0..63)) become lists of floats/strings."""
    out, key, buf = {}, None, []

    def _val(v):
        v = v.strip()
        if v.startswith("<") and v.endswith(">"):
            return v[1:-1]
        try:
            return float(v) if any(c in v for c in ".eE") else int(v)
        except ValueError:
            return v

    def _flush():
        if key is None:
            return
        raw = " ".join(buf).strip()
        m = re.match(r"^\((\d+)\.\.(\d+)\)\s*(.*)$", raw, re.S)
        if m:
            items = re.findall(r"<[^>]*>|\S+", m.group(3))
            out[key] = [_val(x) for x in items]
        else:
            out[key] = _val(raw)

    for line in text.splitlines():
        if line.startswith("$$"):
            continue
        if line.startswith("##"):
            _flush()
            m = re.match(r"^##\$?([^=]+)=\s?(.*)$", line)
            if not m:
                key, buf = None, []
                continue
            key, buf = m.group(1).strip(), [m.group(2)]
        elif key is not None:
            buf.append(line)
    _flush()
    return out


@dataclass
class Experiment:
    path: str
    acqus: dict
    procs: dict = field(default_factory=dict)
    fid: np.ndarray | None = None           # complex
    proc_1r: np.ndarray | None = None       # TopSpin's processed real spectrum (if present)
    name: str = ""

    @property
    def nucleus(self):
        n = str(self.acqus.get("NUC1", "1H"))
        return n if n not in ("off", "") else "1H"

    @property
    def pulprog(self):
        return str(self.acqus.get("PULPROG", ""))

    @property
    def kind(self):
        """'1H', '13C', 'DEPT135', 'DEPT90', 'DEPT45', 'DEPTQ', 'APT' or the nucleus string."""
        pp = self.pulprog.lower()
        if self.nucleus == "13C":
            for t in ("deptq", "dept135", "dept90", "dept45"):
                if t in pp:
                    return t.upper().replace("DEPTQ", "DEPTQ")
            if pp.startswith("dept"):
                ang = self.acqus.get("P", [0] * 13)
                return "DEPT135"
            if "apt" in pp or "jmod" in pp:
                return "APT"
        return self.nucleus

    def info(self) -> dict:
        a = self.acqus
        title = str(a.get("TITLE", ""))
        ver = re.sub(r"\s+", " ", title.replace("Parameter file,", "")).strip()
        return {"nucleus": self.nucleus, "experiment": self.kind, "pulse program": self.pulprog,
                "spectrometer (MHz)": round(float(a.get("SFO1", a.get("BF1", 0))), 2), "solvent": a.get("SOLVENT", ""),
                "scans": a.get("NS"), "TD": a.get("TD"), "SW (ppm)": round(float(a.get("SW", 0)), 3),
                "D1 (s)": (a.get("D") or [0, None])[1] if isinstance(a.get("D"), list) else None,
                "temperature (K)": a.get("TE"), "software": ver, "DSPFVS": a.get("DSPFVS"),
                "instrument": a.get("INSTRUM", ""), "fid dtype": "float64" if a.get("DTYPA") == 2 else "int32",
                "probe": a.get("PROBHD", ""), "date": _date(a.get("DATE"))}


# ------------------------------------------------------------------------------------------------ reading
def _date(v):
    try:
        import datetime
        return datetime.datetime.fromtimestamp(float(v)).strftime("%Y-%m-%d %H:%M")
    except Exception:
        return v


def _read_fid(path: Path, acqus: dict) -> np.ndarray:
    raw = path.read_bytes()
    endian = ">" if int(acqus.get("BYTORDA", 0)) == 1 else "<"
    dt = np.dtype(endian + ("f8" if int(acqus.get("DTYPA", 0)) == 2 else "i4"))
    d = np.frombuffer(raw[: len(raw) // dt.itemsize * dt.itemsize], dtype=dt).astype(float)
    td = int(acqus.get("TD", len(d)))
    d = d[:td] if 0 < td <= len(d) else d
    d = d[: len(d) // 2 * 2]
    return d[0::2] + 1j * d[1::2]


def _read_1r(pdir: Path, procs: dict):
    f = pdir / "1r"
    if not f.exists():
        return None
    endian = ">" if int(procs.get("BYTORDP", 0)) == 1 else "<"
    dt = np.dtype(endian + ("f8" if int(procs.get("DTYPP", 0)) == 2 else "i4"))
    d = np.frombuffer(f.read_bytes(), dtype=dt).astype(float)
    return d * 2.0 ** float(procs.get("NC_proc", 0))


def find_experiments(root) -> list[Path]:
    """All experiment folders (acqus + fid) below `root` (folder path)."""
    root = Path(root)
    res = [p.parent for p in root.rglob("acqus") if (p.parent / "fid").exists()]
    res += [p.parent for p in root.rglob("acqus") if not (p.parent / "fid").exists() and (p.parent / "pdata").exists()
            and any((p.parent / "pdata").rglob("1r"))]
    return sorted(set(res), key=lambda p: (str(p.parent), _num(p.name)))


def _num(s):
    try:
        return int(s)
    except ValueError:
        return 10 ** 9


def read_experiment(d) -> Experiment:
    d = Path(d)
    acqus = parse_jcamp_params((d / "acqus").read_text(errors="replace"))
    acqus["TITLE"] = acqus.get("TITLE", "")
    fid = _read_fid(d / "fid", acqus) if (d / "fid").exists() else None
    procs, r1 = {}, None
    pds = sorted((d / "pdata").glob("*/procs"), key=lambda p: _num(p.parent.name)) if (d / "pdata").exists() else []
    if pds:
        procs = parse_jcamp_params(pds[0].read_text(errors="replace"))
        r1 = _read_1r(pds[0].parent, procs)
    nm = d.parent.name + "/" + d.name if d.parent.name else d.name
    return Experiment(str(d), acqus, procs, fid, r1, nm)


def load(src, name: str = "") -> list[Experiment]:
    """Upload (bytes / file-like .zip) or a folder path → list of experiments (1D only; 2D 'ser' data are skipped)."""
    if isinstance(src, (str, os.PathLike)) and Path(src).is_dir():
        return [e for e in (read_experiment(p) for p in find_experiments(src)) if e.fid is not None or e.proc_1r is not None]
    raw = src.read() if hasattr(src, "read") else (Path(src).read_bytes() if isinstance(src, (str, os.PathLike)) else src)
    if raw[:2] != b"PK":
        raise ValueError("Bruker data: upload a .zip of the experiment (or dataset) folder")
    tmp = Path(tempfile.mkdtemp(prefix="nmr_"))
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for m in z.infolist():
            p = (tmp / m.filename).resolve()
            if str(p).startswith(str(tmp.resolve())):     # zip-slip guard
                z.extract(m, tmp)
    exps = load(tmp)
    for e in exps:
        e.name = (Path(name).stem + ":" if name else "") + os.path.relpath(e.path, tmp).replace(os.sep, "/")
    return exps


# ------------------------------------------------------------------------------------------------ processing
def group_delay(acqus: dict) -> float:
    g = acqus.get("GRPDLY", -1)
    try:
        g = float(g)
    except (TypeError, ValueError):
        g = -1
    if g > 0:
        return g
    dsp = int(acqus.get("DSPFVS", 0) or 0); dec = acqus.get("DECIM", 1)
    try:
        dec = int(round(float(dec)))
    except (TypeError, ValueError):
        dec = 1
    return float(_DSP.get(dsp, {}).get(dec, 0.0))


def whittaker(y, w, lam):
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    n = len(y)
    D = sparse.diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(n, n - 2))
    W = sparse.diags(w)
    return spsolve((W + lam * D.dot(D.T)).tocsc(), w * y)


def noise_sigma(y):
    d = np.diff(y)
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) + 1e-12


def baseline_points(y, k=3.0, grow=None):
    """Signal-free points: small first derivative AND small smoothed amplitude, dilated around peaks."""
    from scipy import ndimage as ndi
    n = len(y); grow = grow or max(8, n // 2000)
    s = noise_sigma(y)
    dy = np.abs(np.gradient(ndi.uniform_filter1d(y, 5)))
    sig = dy > k * s
    sig = ndi.binary_dilation(sig, iterations=grow)
    return ~sig


def auto_baseline(y, lam=None):
    n = len(y)
    if n < 64:
        return y, np.zeros_like(y)
    f = max(1, n // 16384)                              # solve on a decimated grid (fast), interpolate back
    ys = y[: n // f * f].reshape(-1, f).mean(1) if f > 1 else y
    mask = baseline_points(ys)
    if mask.sum() < 20:
        return y, np.zeros_like(y)
    lam = lam or (len(ys) / 50.0) ** 2 * 10
    b = whittaker(ys, mask.astype(float) + 1e-6, lam)
    for _ in range(2):                                   # refit only below-signal points (robust to small missed peaks)
        r = ys - b; s = 1.4826 * np.median(np.abs(r[mask])) + 1e-12
        m2 = mask & (np.abs(r) < 3 * s)
        b = whittaker(ys, m2.astype(float) + 1e-6, lam)
    bl = np.interp(np.arange(n), np.arange(len(ys)) * f + (f - 1) / 2, b)
    return y - bl, bl


def phase(spec, p0, p1, pivot=0.0):
    n = len(spec); x = np.arange(n) / n - pivot
    return spec * np.exp(1j * np.deg2rad(p0 + p1 * x))


def _phase_cost(spec, p0, p1, gamma):
    r = phase(spec, p0, p1).real
    d = np.abs(np.diff(r)); h = d / (d.sum() + 1e-30)
    ent = -np.sum(h * np.log(h + 1e-30))
    neg = r[r < 0]
    return ent + gamma * np.sum(neg ** 2) / (np.sum(r ** 2) + 1e-30)


def _true_lines(mag, pk, heights, n, ratio=8.0, win=None):
    """Drop truncation side-lobes / wiggles: a maximum within `win` points of a line `ratio`× stronger is not a line."""
    win = win or max(20, n // 1500)
    order = np.argsort(heights)[::-1]
    keep = []
    for j in order:
        i = pk[j]
        if any(abs(i - pk[k]) < win and heights[k] > ratio * heights[j] for k in keep):
            continue
        keep.append(j)
    return pk[np.sort(np.array(keep, int))] if keep else pk[:0]


def _apex_phase_fit(spec, p1_range=180.0, n_peaks=40, axial=False, prior=90.0):
    """Phase from the line apices: after correction every line's phase must be 0 (or 0/180 for DEPT: axial=True),
    so p0 + p1·x_i ≈ −φ_i. Weighted circular regression over a p1 grid with a Gaussian prior on p1 (width `prior`):
    after digital-filter removal Bruker data need only a small first-order correction."""
    from scipy.signal import find_peaks
    mag = np.abs(spec); n = len(spec)
    s = noise_sigma(mag)
    pk, pr = find_peaks(mag, height=10 * s, distance=max(3, n // 4000))
    if len(pk) < 1:
        return None
    pk = _true_lines(mag, pk, pr["peak_heights"], n)
    pk = pk[np.argsort(mag[pk])[::-1][:n_peaks]]
    # line phase = phase of the complex sum over a symmetric window (dispersion is odd and cancels) — much less noisy
    # than the single apex point
    hw = 3
    zs = np.array([spec[max(0, i - hw): i + hw + 1].sum() for i in pk])
    x = pk / n; phi = np.angle(zs); w = mag[pk] ** 0.75
    k = 2.0 if axial else 1.0
    best = None
    for p1 in np.arange(-p1_range, p1_range + 0.1, 0.5):
        th = k * (-phi - np.deg2rad(p1 * x))
        z = np.sum(w * np.exp(1j * th)); R = np.abs(z) / w.sum()
        sc = R - 0.15 * (p1 / prior) ** 2
        if best is None or sc > best[0]:
            best = (sc, float(np.rad2deg(np.angle(z)) / k), float(p1), R)
    return best[1], best[2], best[3], len(pk)


def apex_p0(spec, p1, kind="13C", ppm=None):
    """Zero-order phase only (p1 given): weighted circular mean of the line phases (axial for DEPT), then the sign."""
    from scipy.signal import find_peaks
    mag = np.abs(spec); n = len(spec); s = noise_sigma(mag)
    pk, pr = find_peaks(mag, height=8 * s, distance=max(3, n // 4000))
    if not len(pk):
        return 0.0
    pk = _true_lines(mag, pk, pr["peak_heights"], n)
    pk = pk[np.argsort(mag[pk])[::-1][:40]]
    zs = np.array([spec[max(0, i - 3): i + 4].sum() for i in pk])
    dept = str(kind).upper().startswith(("DEPT", "APT")); k = 2.0 if dept else 1.0
    th = k * (-np.angle(zs) - np.deg2rad(p1 * pk / n))
    p0 = float(np.rad2deg(np.angle(np.sum(np.sqrt(mag[pk]) * np.exp(1j * th)))) / k)
    if _sign_score(phase(spec, p0, p1).real, kind, ppm) < 0:
        p0 += 180.0
    return (p0 + 180) % 360 - 180


def _same_receiver(a, b):
    ka = ("SW_h", "DECIM", "DSPFVS", "GRPDLY", "TD", "O1", "NUC1")
    return all(str(a.acqus.get(k)) == str(b.acqus.get(k)) for k in ka if k != "TD") and a.acqus.get("NUC1") == b.acqus.get("NUC1")


def process_all(exps, **kw):
    """Process a dataset: 13C first; DEPT/APT experiments recorded with the same receiver settings reuse its first-order
    phase (DEPT spectra have few, weak and partly negative lines — p1 alone is poorly defined there)."""
    out = {}
    c13 = [e for e in exps if e.kind == "13C" and e.fid is not None]
    for e in exps:
        if e.fid is None and e.proc_1r is None:
            continue
        if False and e.kind.startswith(("DEPT", "APT")) and e.fid is not None:   # p1 differs (pulse timing) — measured
            ref = next((c for c in c13 if _same_receiver(c, e)), None)
            if ref is not None:
                if ref.name not in out:
                    out[ref.name] = process(ref, **kw)
                out[e.name] = process(e, p1=out[ref.name].p1, **kw)
                out[e.name].ref_note += f" · p1 from {ref.name}"
                continue
        if e.name not in out:
            out[e.name] = process(e, **kw)
    res = [out[e.name] for e in exps if e.name in out]
    _dept_sign_from_dept90(res)
    return res


def line_list(sp, k=12.0, max_n=60):
    "[(ppm, height)] of the true lines (side-lobes removed), strongest first."
    from scipy.signal import find_peaks
    y = sp.y; n = len(y)
    s = 1.4826 * np.median(np.abs(y - np.median(y))) + 1e-12
    pk, pr = find_peaks(np.abs(y), height=k * s, distance=max(3, n // 4000))
    if not len(pk):
        return []
    pk = _true_lines(np.abs(y), pk, pr["peak_heights"], n)
    pk = pk[np.argsort(np.abs(y[pk]))[::-1][:max_n]]
    return [(float(sp.ppm[i]), float(y[i])) for i in pk]


def _dept_sign_from_dept90(sps):
    "DEPT-135 sign check: lines also present in the DEPT-90 of the same dataset are CH and must point up."
    d90 = [s for s in sps if s.kind == "DEPT90"]
    for s in sps:
        if s.kind != "DEPT135" or not d90:
            continue
        ch = [p for p, h in line_list(d90[0]) if h > 0]
        if not ch:
            continue
        sc = sum(np.sign(h) * np.sqrt(abs(h)) for p, h in line_list(s) if min(abs(p - c) for c in ch) < 0.08)
        if sc < 0:
            s.y = -s.y; s.p0 = (s.p0 + 360) % 360 - 180
            s.ref_note += " · sign from DEPT-90 (CH up)"


def _acme_cost(spec, p0, p1, gamma=1000.0, pen=0.0):
    r = phase(spec, p0, p1).real
    d = np.abs(np.diff(r)); h = d / (d.sum() + 1e-30)
    ent = -np.sum(h * np.log(h + 1e-30))
    neg = r[r < 0]
    return ent + gamma * np.sum(neg ** 2) / (np.sum(r ** 2) + 1e-30) + pen * (p1 / 90.0) ** 2


def _sign_score(y, kind, ppm=None):
    """>0 if the phase sign is right. 13C/DEPT90: lines positive. DEPT135: CH (aromatic / olefinic, > 100 ppm) and
    typical CH3 (< 25 ppm) up, CH2 down — weighted by region since the sign of a single line is not known a priori.
    Works on a roughly flattened copy (running median) so a baseline offset cannot decide the sign."""
    from scipy.signal import find_peaks
    from scipy import ndimage as ndi
    n = len(y); f = max(1, n // 8192)
    ys = y[: n // f * f].reshape(-1, f).mean(1) if f > 1 else y
    bl = ndi.median_filter(ys, size=max(51, len(ys) // 40) | 1)
    y0 = y - np.interp(np.arange(n), np.arange(len(ys)) * f + (f - 1) / 2, bl)
    s = 1.4826 * np.median(np.abs(y0 - np.median(y0))) + 1e-12
    pk, pr = find_peaks(np.abs(y0), height=12 * s, distance=max(3, n // 4000))
    if not len(pk):
        return float(np.sum(y0))
    pk = _true_lines(np.abs(y0), pk, pr["peak_heights"], n)
    pk = pk[np.argsort(np.abs(y0[pk]))[::-1][:30]]
    a = y0[pk]
    k = str(kind).upper()
    if k == "DEPT135" and ppm is not None:
        w = np.where(ppm[pk] > 100, 1.0, np.where(ppm[pk] < 25, 0.4, 0.05))
        return float(np.sum(w * np.sign(a) * np.sqrt(np.abs(a))))
    if k == "APT" and ppm is not None:
        w = np.where(ppm[pk] > 100, -1.0, 0.1)
        return float(np.sum(w * np.sign(a) * np.sqrt(np.abs(a))))
    return float(np.sum(np.sign(a) * np.abs(a)))


def auto_phase(spec, nucleus="1H", kind=None, p1_range=180.0, ppm=None):
    """Automatic phase (p0, p1 in degrees; x from 0 at the left edge).
    1H: ACME (entropy of the derivative + negative-area penalty; Chen et al. 2002) on a coarse grid, then simplex,
        with a weak penalty on large p1.
    13C / DEPT: weighted circular regression of the line-apex phases (robust to noise; axial for DEPT where lines may
        be negative), then the sign is chosen so that the strongest lines are positive."""
    from scipy.optimize import minimize
    kind = str(kind or nucleus)
    n = len(spec)
    dept = kind.upper().startswith(("DEPT", "APT"))
    if nucleus == "1H" and not dept:
        f = max(1, n // 8192)
        sd = spec[::f]
        best = None
        for p1 in np.linspace(-p1_range, p1_range, 13):
            for p0 in np.arange(-180, 180, 15):
                c = _acme_cost(sd, p0, p1, pen=0.3)
                if best is None or c < best[0]:
                    best = (c, p0, p1)
        r = minimize(lambda q: _acme_cost(sd, q[0], q[1], pen=0.3), [best[1], best[2]], method="Nelder-Mead",
                     options={"xatol": 0.05, "fatol": 1e-9, "maxiter": 600})
        r = minimize(lambda q: _acme_cost(spec, q[0], q[1], pen=0.3), r.x, method="Nelder-Mead",
                     options={"xatol": 0.02, "fatol": 1e-10, "maxiter": 300})
        p0, p1 = float(r.x[0]), float(r.x[1])
        if phase(spec, p0, p1).real.sum() < 0:
            p0 += 180.0
    else:
        ap = _apex_phase_fit(spec, p1_range=max(p1_range, 720.0), axial=dept, prior=400.0)
        p0, p1 = (0.0, 0.0) if ap is None else (ap[0], ap[1])
        y = phase(spec, p0, p1).real
        if _sign_score(y, kind, ppm) < 0:
            p0 += 180.0
    return (p0 + 180) % 360 - 180, p1


@dataclass
class Spectrum:
    ppm: np.ndarray
    y: np.ndarray                 # processed real (phased, baseline-corrected), arbitrary units
    nucleus: str
    kind: str
    sf: float                     # MHz used for the ppm axis
    info: dict
    p0: float = 0.0
    p1: float = 0.0
    ref_shift: float = 0.0        # ppm added by referencing
    ref_note: str = ""
    name: str = ""
    complex_spec: np.ndarray | None = None

    def hz_per_pt(self):
        return abs(self.ppm[1] - self.ppm[0]) * self.sf


def process(e: Experiment, lb=None, zf=2, p0=None, p1=None, baseline=True, reference="auto", use_topspin=False) -> Spectrum:
    """Full automatic processing of one experiment. `reference`: 'auto' | 'topspin' | 'none' | float (ppm shift)."""
    a, pr = e.acqus, e.procs
    nuc = e.nucleus
    sfo1 = float(a.get("SFO1", a.get("BF1", 400.0))); bf1 = float(a.get("BF1", sfo1))
    sw_hz = float(a.get("SW_h", float(a.get("SW", 10)) * sfo1))
    o1 = float(a.get("O1", (sfo1 - bf1) * 1e6))
    sf = float(pr.get("SF", bf1)) if pr else bf1
    if use_topspin and e.proc_1r is not None and pr:
        n = len(e.proc_1r)
        ppm = float(pr["OFFSET"]) - np.arange(n) * (float(pr["SW_p"]) / float(pr["SF"])) / n
        y = e.proc_1r[::-1][::-1].astype(float)
        sp = Spectrum(ppm, y, nuc, e.kind, float(pr["SF"]), e.info(), ref_note="TopSpin processed (1r)", name=e.name)
        if baseline:
            sp.y, _ = auto_baseline(sp.y)
        return sp
    fid = e.fid.copy()
    gd = group_delay(a)
    n0 = len(fid)
    lb = (0.3 if nuc == "1H" else 1.0) if lb is None else lb
    t = np.arange(n0) / sw_hz
    fid = fid * np.exp(-np.pi * lb * t)
    fid[0] *= 0.5
    si = 1 << int(np.ceil(np.log2(n0 * max(1, zf))))
    si = max(si, n0)
    spec = np.fft.fftshift(np.fft.fft(fid, si))
    if gd:                                              # digital filter: first-order phase of 360°·GRPDLY over the window
        k = np.arange(si) - si / 2
        spec = spec * np.exp(2j * np.pi * gd * k / si)
    spec = spec[::-1]                                   # left (high ppm) → right
    # ppm axis (TopSpin): centre at O1, referenced against SF (procs) = BF1 + SR
    hz = (np.arange(si)[::-1] - si / 2) * sw_hz / si + o1 + (bf1 - sf) * 1e6
    ppm = hz / sf
    if p0 is None and p1 is not None:             # first-order phase known (e.g. from the 13C of the same dataset)
        p0 = apex_p0(spec, p1, e.kind, ppm)
    elif p0 is None or p1 is None:
        q0, q1 = auto_phase(spec, nuc, e.kind, ppm=ppm)
        p0 = q0 if p0 is None else p0; p1 = q1 if p1 is None else p1
    y = phase(spec, p0, p1).real
    if baseline:
        y, _ = auto_baseline(y)
    sp = Spectrum(ppm, y, nuc, e.kind, sf, e.info(), p0, p1, name=e.name, complex_spec=spec)
    if reference == "auto":
        auto_reference(sp)
    elif isinstance(reference, (int, float)) and not isinstance(reference, bool):
        sp.ppm = sp.ppm + float(reference); sp.ref_shift = float(reference); sp.ref_note = f"manual {reference:+.3f} ppm"
    else:
        sp.ref_note = "TopSpin SR" if pr else "spectrometer frequency"
    return sp


def _local_max_near(ppm, y, target, win):
    m = (ppm > target - win) & (ppm < target + win)
    if m.sum() < 3:
        return None
    idx = np.where(m)[0]
    i = idx[np.argmax(y[idx])]
    s = noise_sigma(y)
    if y[i] < 8 * s:
        return None
    if 0 < i < len(y) - 1:                               # parabolic interpolation
        a_, b_, c_ = y[i - 1], y[i], y[i + 1]
        den = a_ - 2 * b_ + c_
        d = 0.5 * (a_ - c_) / den if den else 0.0
        return float(ppm[i] + d * (ppm[i + 1] - ppm[i]))
    return float(ppm[i])


def auto_reference(sp: Spectrum, solvent=None, win=None):
    """TMS (CDCl3, organic) / DSS (D2O) at 0 ppm if present, else the residual solvent signal.
    For 13C the solvent carbon (e.g. CDCl3 triplet centre 77.16). Only small corrections (≤ win) are applied."""
    key = solvent_key(solvent or sp.info.get("solvent"))
    nuc = sp.nucleus
    win = win or (0.25 if nuc == "1H" else 2.0)
    cands = []
    if nuc == "1H":
        cands.append(("TMS/DSS", 0.0, 0.08))
    if key and SOLVENTS[key].get(nuc) is not None:
        cands.append((key + " residual", SOLVENTS[key][nuc], win))
    if nuc == "13C":
        cands.append(("TMS", 0.0, 1.0))
    for nm, target, w in cands:
        if nm.startswith("TMS") and nuc == "1H":
            pk = _local_max_near(sp.ppm, sp.y, target, w)
            if pk is None:
                continue
            # TMS must be a sharp singlet among the bigger signals (avoid grease/noise)
            if np.max(sp.y[(sp.ppm > -0.1) & (sp.ppm < 0.1)]) < 0.01 * np.max(sp.y):
                continue
        else:
            pk = _local_max_near(sp.ppm, sp.y, target, w)
        if pk is None:
            continue
        if nuc == "13C" and key == "CDCl3" and nm.endswith("residual"):
            pk = _centre_of_triplet(sp, pk)
        d = target - pk
        sp.ppm = sp.ppm + d; sp.ref_shift += d; sp.ref_note = f"{nm} → {target:g} ppm (shift {d:+.3f} ppm)"
        return sp
    sp.ref_note = "TopSpin SR (no reference signal found)"
    return sp


def _centre_of_triplet(sp, pk):
    """CDCl3 13C 1:1:1 triplet (J ≈ 32 Hz): take the centre line."""
    jppm = 32.0 / sp.sf
    m = (sp.ppm > pk - 1.5 * jppm - 0.3) & (sp.ppm < pk + 1.5 * jppm + 0.3)
    if m.sum() < 5:
        return pk
    from scipy.signal import find_peaks
    idx = np.where(m)[0]
    pi, _ = find_peaks(sp.y[idx], height=0.3 * sp.y[idx].max())
    if len(pi) == 3:
        return float(sp.ppm[idx[pi[1]]])
    return pk


def to_frame(sp: Spectrum):
    import pandas as pd
    df = pd.DataFrame({"x": sp.ppm, "y": sp.y})
    df.attrs.update({"source": sp.name, "SF": sp.sf, "nucleus": sp.nucleus, "kind": sp.kind, "ref": sp.ref_note})
    return df


# ------------------------------------------------------------------------------------------------ JCAMP-DX (processed)
def _f(dic, k, default=None):
    v = dic.get(k)
    if v is None:
        return default
    v = v[0] if isinstance(v, list) else v
    try:
        return float(str(v).split(",")[0])
    except ValueError:
        return default


def load_jcamp(src, name: str = "") -> Spectrum:
    """JCAMP-DX 1D NMR spectrum (Bruker/TopSpin NTUPLES real+imag, or XYDATA in Hz/ppm) → Spectrum.
    Real + imaginary are kept so the phase can still be corrected."""
    import nmrglue as ng
    import warnings
    raw = src.read() if hasattr(src, "read") else (Path(src).read_bytes() if isinstance(src, (str, os.PathLike)) else src)
    tmp = Path(tempfile.mkdtemp(prefix="jdx_")) / "s.jdx"
    tmp.write_bytes(raw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dic, data = ng.jcampdx.read(str(tmp))
    nuc = str((dic.get(".OBSERVENUCLEUS") or ["^1H"])[0]).replace("^", "").strip() or "1H"
    sf = _f(dic, "$SF") or _f(dic, ".OBSERVEFREQUENCY", 400.0)
    if isinstance(data, list):
        re_, im_ = np.asarray(data[0], float), (np.asarray(data[1], float) if len(data) > 1 else None)
    else:
        re_, im_ = np.asarray(data, float), None
    n = len(re_)
    if "FIRST" in dic:                                    # NTUPLES: FIRST = Hz of the first point; $OFFSET = its ppm
        first = _f(dic, "FIRST"); last = float(str(dic["LAST"][0]).split(",")[0])
        off = _f(dic, "$OFFSET", first / sf)
        hz = np.linspace(first, last, n)
        ppm = off - (first - hz) / sf
    else:
        fx, lx = _f(dic, "FIRSTX"), _f(dic, "LASTX")
        unit = str((dic.get("XUNITS") or ["PPM"])[0]).upper()
        x = np.linspace(fx, lx, n)
        ppm = x / sf if "HZ" in unit else x
    if ppm[0] < ppm[-1]:
        ppm = ppm[::-1]; re_ = re_[::-1]; im_ = im_[::-1] if im_ is not None else None
    pp = str((dic.get(".PULSESEQUENCE") or dic.get("$PULPROG") or [""])[0]).strip("<> ")
    kind = nuc
    if nuc == "13C":
        for t in ("deptq", "dept135", "dept90", "dept45"):
            if t in pp.lower():
                kind = t.upper()
    solv = str((dic.get(".SOLVENTNAME") or dic.get("$SOLVENT") or [""])[0]).strip("<> ")
    info = {"nucleus": nuc, "experiment": kind, "pulse program": pp, "spectrometer (MHz)": round(sf, 2), "solvent": solv,
            "software": "JCAMP-DX " + str((dic.get("JCAMPDX") or [""])[0]).split("$$")[0].strip(), "fid dtype": "processed"}
    cs = re_ + 1j * im_ if im_ is not None else None
    y, _ = auto_baseline(re_)
    sp = Spectrum(ppm, y, nuc, kind, sf, info, name=name or "jcamp", complex_spec=cs,
                  ref_note="as exported (" + str((dic.get(".SHIFTREFERENCE") or ["file"])[0]) + ")")
    return sp


def rephase(sp: Spectrum, p0: float, p1: float, baseline=True) -> Spectrum:
    """Apply an additional manual phase to a processed spectrum (needs the complex spectrum)."""
    if sp.complex_spec is None:
        return sp
    y = phase(sp.complex_spec, p0, p1).real
    if baseline:
        y, _ = auto_baseline(y)
    sp.y = y; sp.p0, sp.p1 = p0, p1
    return sp
