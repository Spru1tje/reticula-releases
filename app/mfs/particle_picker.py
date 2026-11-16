"""Click-to-exclude particle picker: outlines for the browser component and the include/exclude logic (mirrored in
components/particle_picker/index.html, tested here)."""
from __future__ import annotations
import base64, io, os
import numpy as np
from scipy import ndimage as ndi

def component_dir():
    """Folder of the bundled static front-end (also inside the PyInstaller exe, where mfs is in the archive)."""
    import sys
    cands = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "components", "particle_picker")]
    if getattr(sys, "_MEIPASS", None):
        cands.insert(0, os.path.join(sys._MEIPASS, "app", "mfs", "components", "particle_picker"))
    for c in cands:
        if os.path.exists(os.path.join(c, "index.html")):
            return c
    return cands[-1]


def outlines(lab, factor=1, scale=1.0, max_pts=48):
    """{id: polygon [[x, y], ...]} in preview pixels (full image px × scale), the longest contour of each label."""
    import contourpy
    out = {}
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        sub = np.pad((lab[sl] == i).astype(float), 1)
        try:
            lines = contourpy.contour_generator(z=sub).lines(0.5)
        except Exception:
            continue
        if not lines:
            continue
        L = max(lines, key=len)
        step = max(1, len(L) // max_pts)
        L = L[::step]
        x = (L[:, 0] - 1 + sl[1].start + 0.5) * factor * scale; y = (L[:, 1] - 1 + sl[0].start + 0.5) * factor * scale
        out[i] = np.round(np.stack([x, y], 1), 1).tolist()
    return out


def particles_payload(rows, lab, factor, scale, pixel_nm, metric="ECD (nm)"):
    """Particles for the component: id, centroid and outline in preview px, the size shown in the live stats, kind."""
    polys = outlines(lab, factor, scale)
    P = []
    for r in rows:
        i = int(r["id"])
        if i not in polys:
            continue
        P.append({"id": i, "cx": round(float(r["x (px)"]) * scale, 1), "cy": round(float(r["y (px)"]) * scale, 1),
                  "size": float(r[metric]), "kind": r.get("kind", "particle"), "prior": bool(r.get("prior")),
                  "conf": (round(float(r["confidence"]), 2) if r.get("confidence") is not None else None),
                  "r": float(np.sqrt(r["area (nm²)"] / np.pi) / pixel_nm * scale), "poly": polys[i]})
    return P


def preview_src(img8, max_w=900):
    from PIL import Image
    im = Image.fromarray(np.asarray(img8)).convert("L")
    s = min(1.0, max_w / im.size[0])
    if s < 1:
        im = im.resize((int(im.size[0] * s), int(im.size[1] * s)))
    b = io.BytesIO(); im.save(b, format="JPEG", quality=85)
    return "data:image/jpeg;base64," + base64.b64encode(b.getvalue()).decode(), s, im.size


# ---- the logic of the component (same rules as the JS)
def toggle(excluded, pid):
    e = set(excluded); e.symmetric_difference_update({pid}); return sorted(e)


def nearest(particles, x, y, snap=60):
    best, bd = None, 1e18
    for p in particles:
        d = (p["cx"] - x) ** 2 + (p["cy"] - y) ** 2
        if d < bd:
            best, bd = p["id"], d
    return best if best is not None and bd ** 0.5 < snap else None


def box(particles, excluded, x0, y0, x1, y1, include=False):
    xa, xb, ya, yb = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
    ids = {p["id"] for p in particles if xa <= p["cx"] <= xb and ya <= p["cy"] <= yb}
    e = set(excluded)
    return sorted(e - ids if include else e | ids)


def auto_excluded(rows):
    return sorted(int(r["id"]) for r in rows if r.get("kind") in ("aggregate", "lowconf"))


def apply(rows, excluded):
    ex = set(int(i) for i in (excluded or []))
    return [r for r in rows if int(r["id"]) not in ex]


def picker(src, w, h, particles, excluded, state_key, metric="ECD", unit="nm", snap=60, component_key=None, numbers=False, rev=0):
    """Render the component; returns {'key', 'excluded', 't'} after a change (debounced in the browser), else None.
    rev: bump it when `excluded` was changed on the server (page undo/redo) so the browser adopts it."""
    import streamlit.components.v1 as components
    global _COMP
    if _COMP is None:
        _COMP = components.declare_component("reticula_particle_picker", path=component_dir())
    return _COMP(src=src, w=w, h=h, particles=particles, excluded=list(excluded), skey=state_key, metric=metric, unit=unit,
                 snap=snap, debounce=350, numbers=bool(numbers), rev=int(rev), default=None, key=component_key)


_COMP = None
