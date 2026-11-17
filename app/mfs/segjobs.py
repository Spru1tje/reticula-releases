"""v2.14: the three particle segmentations as plain functions (no Streamlit), so that worker processes can run them.

run(kind, img, args) is what mfs.segpool sends to a worker; the Streamlit page calls the same function in-process when
the work is a single image. Results are (labels, working-grid factor, threshold, info | rows)."""
from __future__ import annotations

import numpy as np

KINDS = ("adaptive", "watershed", "auto")


def img8(im):
    a = np.asarray(im, float)
    if a.ndim == 3:
        a = a.mean(axis=2)
    lo, hi = np.percentile(a, (0.2, 99.8))
    return (np.clip((a - lo) / ((hi - lo) or 1), 0, 1) * 255).astype(np.uint8)


def auto(img, pix, dark, thr_pct, smooth, open_px, split, border, max_px):
    from mfs import particles as PA
    a = img8(img)
    thr = None if thr_pct is None else float(np.percentile(a, thr_pct))
    lab, f, t = PA.segment(a, dark_particles=dark, threshold=thr, smooth=smooth, open_px=open_px, split=split,
                           clear_border=border, max_px=max_px)
    return lab, f, t, PA.measure(lab, pix, f)


def watershed(img, dark, thr_pct, smooth, open_px, mdf, hf, meth, holes, border, sol, areaf, max_px):
    from mfs import particles as PA
    a = img8(img)
    thr = None if thr_pct is None else float(np.percentile(a, thr_pct))
    lab, f, t, info = PA.segment_ws(a, dark_particles=dark, threshold=thr, smooth=smooth, open_px=open_px, min_dist_factor=mdf,
                                    h_factor=hf, method=meth, exclude_holes=holes, clear_border=border, solidity_min=sol,
                                    area_factor=areaf, max_px=max_px)
    info = {k: v for k, v in info.items() if k != "_hole_mask"}
    return lab, f, t, info


def adaptive(img, dark, smooth, low_k, size_px, mks, sol, areaf, border, max_px, prior, pstr, pk, sph="auto", cmin=0.5):
    from mfs import particles as PA
    a = img8(img)
    lab, f, t, info = PA.segment_adaptive(a, dark_particles=dark, smooth=smooth, low_k=low_k, size_px=size_px or None,
                                          marker_sigma=mks, solidity_min=sol, area_factor=areaf, clear_border=border,
                                          max_px=max_px, shape_prior=prior, prior_strength=pstr, prior_k=pk,
                                          sphere={"auto": "auto", "spheres": True, "off": False}.get(sph, "auto"), conf_min=cmin)
    info = {k: v for k, v in info.items() if k != "_contrast"}
    return lab, f, t, info


_FUNCS = {"adaptive": adaptive, "watershed": watershed, "auto": auto}


def run(kind, img, args):
    """Entry point for a worker: args is the tuple of parameters after the image."""
    return _FUNCS[kind](img, *args)
