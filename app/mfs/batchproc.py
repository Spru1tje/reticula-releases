"""Batch processing: every file of a folder (or an upload) read with the page's reader, drawn with the standard plot,
one saved figure-settings preset (or the house style) applied, exported in the chosen formats + CSV, zipped."""
from __future__ import annotations
import io, os, re, zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from . import readers as R, plots as P, style as S, figctl as F, analysis as A

KINDS = {
    "PXRD": dict(ext=(".xy", ".xye", ".txt", ".csv", ".dat", ".raw", ".brml", ".asc", ".rasx", ".ras"), read=lambda f, n: R.read_pxrd(f, n), kind="pxrd",
                 xl=r"2$\theta$ (°)", yl="Intensity (a.u.)"),
    "FTIR": dict(ext=(".dpt", ".csv", ".txt", ".spa", ".0", ".1", ".dat"), read=lambda f, n: R.read_ftir(f, n), kind="ftir",
                 xl=r"Wavenumber (cm$^{-1}$)", yl="Transmittance (a.u.)"),
    "Raman": dict(ext=(".dpt", ".csv", ".txt", ".dat"), read=lambda f, n: R.read_ftir(f, n), kind="raman",
                  xl=r"Raman shift (cm$^{-1}$)", yl="Intensity (a.u.)"),
    "XPS": dict(ext=(".csv", ".txt", ".vms", ".dat"), read=lambda f, n: R.read_xps(f, n), kind="xps", xl="Binding energy (eV)", yl="Intensity (a.u.)"),
    "TGA": dict(ext=(".txt", ".csv", ".xlsx", ".xls", ".dat"), read=lambda f, n: R.read_tga(f, n), kind="tga", xl="Temperature (°C)", yl="Mass (%)"),
    "DSC": dict(ext=(".txt", ".csv", ".xlsx", ".dat"), read=None, kind="dsc", xl="Temperature (°C)", yl="Heat flow (W g$^{-1}$), exo up"),
    "N2 isotherm": dict(ext=(".xls", ".xlsx", ".csv", ".txt", ".rep", ".qps"), read=lambda f, n: R.read_sorption(f, n), kind="iso",
                        xl=r"$p/p_0$", yl=r"$V_\mathrm{ads}$ (cm$^3$ g$^{-1}$ STP)"),
    "SEM/TEM particles": dict(ext=(".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp", ".dm3", ".dm4", ".ser", ".emd", ".mrc"),
                              read=None, kind="particles", xl="ECD (nm)", yl="Count"),
}
# v2.11: particle sizing in batch = the adaptive detection (+ optional learned shape model) of the SEM page
PARTICLE_DEFAULTS = dict(dark=True, shape_prior=False, prior_strength=0.7, prior_k=10.0, nm_per_px=None, low_k=4.0,
                         solidity_min=0.8, area_factor=2.0, max_px=1536, exclude_aggregates=True,
                         sphere="auto", conf_min=0.5)          # v2.12: circle model for round particles, low-confidence cut


SKIP = {"TGA": r"qms|ion.?current|_ms[_.]|tga.?ms", "DSC": r"cp_blank|cp_sapphire", "N2 isotherm": r"kernel"}


def list_folder(folder, kind):
    """Files of one technique below `folder`; README files and companion files (TGA–MS ion currents …) are skipped."""
    ext = KINDS[kind]["ext"]; skip = SKIP.get(kind)
    return sorted(os.path.join(dp, f) for dp, _, fs in os.walk(folder) for f in fs
                  if f.lower().endswith(ext) and not f.startswith((".", "~$")) and not f.upper().startswith("README")
                  and not (skip and re.search(skip, f, re.I)))


def _name(n):
    return re.sub(r"\.[^.]+$", "", os.path.basename(str(n)))


def load(kind, src, name):
    """Returns a list of (name, x, y) traces (TGA: (name, analysis dict); isotherms: (name, iso df))."""
    k = KINDS[kind]; nm = _name(name)
    if kind == "DSC":
        from . import dsc as D
        d = D.read_dsc(src, name)
        return [(f"{nm} {c['name']}" if len(d["curves"]) > 1 else nm, np.asarray(c["T"]), np.asarray(c["q"])) for c in d["curves"]]
    if kind == "SEM/TEM particles":
        return [(nm, particles_one(src, name, PARTICLE_OPTS))]
    d = k["read"](src, name)
    if kind == "TGA":
        return [(nm, A.tga_analysis(d["T"].values, d["w"].values, 21))]
    if kind == "N2 isotherm":
        return [(nm, d["iso"])]
    return [(nm, d.x.values.astype(float), d.y.values.astype(float))]


PARTICLE_OPTS = dict(PARTICLE_DEFAULTS)


def particles_one(src, name, opts=None):
    """One image -> adaptive segmentation (+ shape model) -> dict(image, lab, f, info, rows, pixel_nm, summary)."""
    from . import particles as PA
    o = dict(PARTICLE_DEFAULTS, **(opts or {}))
    m = R.read_sem(src, name)
    img = np.asarray(m["image"])
    if m.get("databar_row"):
        img = img[: int(m["databar_row"])]
    pix = m.get("pixel_nm") or o.get("nm_per_px")
    if not pix:
        raise ValueError("no pixel size in the file metadata – enter a fallback nm/px")
    if img.dtype != np.uint8:
        lo_, hi_ = np.percentile(img, (0.5, 99.5)); img = np.clip((img - lo_) / max(hi_ - lo_, 1e-9) * 255, 0, 255).astype(np.uint8)
    lab, f, _t, info = PA.segment_adaptive(img, dark_particles=o["dark"], low_k=o["low_k"], solidity_min=o["solidity_min"],
                                           area_factor=o["area_factor"], max_px=o["max_px"], shape_prior=o["shape_prior"],
                                           prior_strength=o["prior_strength"], prior_k=o["prior_k"],
                                           sphere=o["sphere"], conf_min=o["conf_min"])
    rows = PA.measure(lab, float(pix), f, info=info)
    for r in rows:
        v = info.get(r["id"]) if isinstance(info.get(r["id"]), dict) else {}
        r["kind"] = v.get("kind", "particle"); r["prior split"] = v.get("split") == "prior"
        r["included"] = not ((o["exclude_aggregates"] and r["kind"] == "aggregate") or r["kind"] == "lowconf")
    e = np.array([r["ECD (nm)"] for r in rows if r["included"]])
    pr = info.get("_prior") or {}
    summ = {"n": int(len(e)), "mean": float(e.mean()) if len(e) else np.nan, "sd": float(e.std(ddof=1)) if len(e) > 1 else np.nan,
            "aggregates": sum(1 for r in rows if r["kind"] == "aggregate"), "prior_splits": sum(1 for r in rows if r["prior split"]),
            "prior_n": pr.get("n"), "pixel_nm": float(pix),
            "spheres": bool((info.get("_sphere") or {}).get("active")), "low_confidence": sum(1 for r in rows if r["kind"] == "lowconf")}
    return {"image": img, "lab": lab, "f": f, "info": {k: v for k, v in info.items() if k != "_contrast"}, "rows": rows,
            "pixel_nm": float(pix), "summary": summ}


def draw(ax, kind, traces, layout):
    k = KINDS[kind]
    if kind == "SEM/TEM particles":
        if len(traces) == 1 and layout != "stacked":
            from . import particles as PA
            nm, d = traces[0]
            ex = [r["id"] for r in d["rows"] if not r["included"]]
            PA.overlay_labels(ax, d["image"], d["lab"], d["f"], d["info"], excluded=ex)
            s_ = d["summary"]
            ax.set_title(f"{nm}: N = {s_['n']}, ECD {s_['mean']:.0f} ± {s_['sd']:.0f} nm"
                         + (f", {s_['prior_splits']} shape-model splits" if s_["prior_splits"] else ""), fontsize="small")
            return
        for nm, d in traces:                       # several images: ECD histograms (included particles)
            e = [r["ECD (nm)"] for r in d["rows"] if r["included"]]
            ax.hist(e, bins=20, histtype="step", lw=1.2, label=f"{nm} (N = {len(e)})")
        ax.set_xlabel(k["xl"]); ax.set_ylabel(k["yl"]); ax.legend(frameon=False, fontsize="small")
        return
    if kind == "TGA":
        P.draw_tga(ax, traces); return
    if kind == "N2 isotherm":
        P.draw_isotherms(ax, traces); return
    if layout == "stacked" and k["kind"] in ("pxrd", "ftir", "raman", "xps"):
        P.draw_stacked(ax, traces, kind=k["kind"]); ax.set_xlabel(k["xl"]); return
    P.draw_overlay(ax, traces, k["xl"], k["yl"], invert=k["kind"] in ("ftir", "xps"), legend=len(traces) > 1)


def _csv(kind, traces):
    if kind == "SEM/TEM particles":
        return pd.concat([pd.DataFrame(d["rows"]).assign(sample=n) for n, d in traces]).to_csv(index=False).encode()
    if kind == "TGA":
        return pd.concat([pd.DataFrame({"sample": n, "T": a["T"], "w": a["w"]}) for n, a in traces]).to_csv(index=False).encode()
    if kind == "N2 isotherm":
        return pd.concat([d.assign(sample=n) for n, d in traces]).to_csv(index=False).encode()
    return pd.concat([pd.DataFrame({"sample": n, "x": x, "y": y}) for n, x, y in traces]).to_csv(index=False).encode()


def run(kind, files, cfg=None, layout="one figure per file", formats=("png", "pdf"), dpi=600, width=S.SINGLE, aspect=0.8,
        font_size=8.0, out_dir=None, progress=None, particle_opts=None):
    """files: list of (src, name). Returns (zip bytes, log rows). With out_dir the files are also written there."""
    global PARTICLE_OPTS
    PARTICLE_OPTS = dict(PARTICLE_DEFAULTS, **(particle_opts or {}))
    S.apply_style(font_size)
    data, log = [], []
    for i, (src, name) in enumerate(files):
        try:
            tr = load(kind, src, name); data.append((name, tr))
            row = {"file": os.path.basename(str(name)), "status": "ok", "traces": len(tr)}
            if kind == "SEM/TEM particles":
                sm = tr[0][1]["summary"]
                row.update({"N": sm["n"], "ECD mean (nm)": round(sm["mean"], 1), "ECD SD (nm)": round(sm["sd"], 1),
                            "aggregates": sm["aggregates"], "shape-model splits": sm["prior_splits"], "nm/px": sm["pixel_nm"]})
            log.append(row)
        except Exception as e:
            log.append({"file": os.path.basename(str(name)), "status": f"error: {e}"[:160], "traces": 0})
        if progress:
            progress((i + 1) / max(len(files), 1) * 0.5)
    groups = [(_name(n), tr) for n, tr in data] if layout == "one figure per file" else \
             [(f"{kind}_all_{'stacked' if layout == 'stacked' else 'overlay'}", [t for _, tr in data for t in tr])] if data else []
    buf = io.BytesIO(); zf = zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED)
    for j, (stem, tr) in enumerate(groups):
        S.apply_style(font_size)
        fig, ax = plt.subplots(figsize=(width, width * aspect))
        draw(ax, kind, tr, "stacked" if layout == "stacked" else "overlay")
        if cfg:
            F.apply(fig, cfg)
        stem = re.sub(r"[^A-Za-z0-9._@+-]+", "_", stem)
        for fmt in formats:
            b = S.export(fig, fmt, dpi); zf.writestr(f"{stem}.{fmt}", b)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True); open(os.path.join(out_dir, f"{stem}.{fmt}"), "wb").write(b)
        c = _csv(kind, tr); zf.writestr(f"data/{stem}.csv", c)
        if out_dir:
            os.makedirs(os.path.join(out_dir, "data"), exist_ok=True); open(os.path.join(out_dir, "data", f"{stem}.csv"), "wb").write(c)
        plt.close(fig)
        if progress:
            progress(0.5 + (j + 1) / max(len(groups), 1) * 0.5)
    zf.writestr("batch_log.csv", pd.DataFrame(log).to_csv(index=False))
    if cfg:
        zf.writestr("preset_used.json", F.to_json(cfg))
    zf.close()
    return buf.getvalue(), log
