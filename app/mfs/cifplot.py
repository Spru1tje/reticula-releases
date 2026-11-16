"""Pure-matplotlib structure drawings: ORTEP ellipsoids (principal axes), ball-and-stick, stick, space-filling, labels,
H-bonds, cell edges, axis triad, metal polyhedra. Depth-sorted (painter's algorithm), so PNG/TIFF/SVG/PDF all work."""
from __future__ import annotations
import math
import numpy as np
from scipy import stats
from matplotlib.patches import Ellipse, Circle, Polygon, FancyArrowPatch
from .cif import COV, VDW, COLOR, METALS


def prob_scale(p):
    """Radius (in σ) of the 3D-Gaussian surface enclosing probability p: 30 % 1.097, 50 % 1.538, 70 % 2.027."""
    return float(stats.chi(3).ppf(p))


def view_matrix(S, view="a", uvw=None, hkl=None, rot=(0.0, 0.0, 0.0), sites=None):
    """Rows = screen x, screen y, depth (towards viewer). view: 'a'|'b'|'c'|'uvw'|'hkl'|'best' (normal to the mean plane
    of `sites`, longest extent horizontal: the least-overlap view for one molecule)."""
    M = S.M
    if view == "best" and sites:
        X = np.array([s_["xyz"] for s_ in sites]); X = X - X.mean(0)
        _, _, Vt = np.linalg.svd(X, full_matrices=False)
        V0 = np.vstack([Vt[0], Vt[1], np.cross(Vt[0], Vt[1])])
        heavy = np.array([s_["el"] != "H" for s_ in sites])
        Xh = X[heavy] if heavy.sum() >= 3 else X
        best = (1e9, V0)
        for tx in range(-40, 41, 10):                      # tilt off the mean plane to separate overlapping atoms
            for ty in range(-40, 41, 10):
                cx, sx, cy, sy = math.cos(math.radians(tx)), math.sin(math.radians(tx)), math.cos(math.radians(ty)), math.sin(math.radians(ty))
                R = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]]) @ np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
                Vc = R @ V0; Q = Xh @ Vc[:2].T
                dd = np.sqrt(((Q[:, None] - Q[None]) ** 2).sum(-1))[np.triu_indices(len(Q), 1)]
                score = float(np.sum(np.clip(1.2 - dd, 0, None))) + 0.002 * (abs(tx) + abs(ty))
                if score < best[0]:
                    best = (score, Vc)
        V = best[1]
        rx, ry, rz = map(math.radians, rot)
        Rx = np.array([[1, 0, 0], [0, math.cos(rx), -math.sin(rx)], [0, math.sin(rx), math.cos(rx)]])
        Ry = np.array([[math.cos(ry), 0, math.sin(ry)], [0, 1, 0], [-math.sin(ry), 0, math.cos(ry)]])
        Rz = np.array([[math.cos(rz), -math.sin(rz), 0], [math.sin(rz), math.cos(rz), 0], [0, 0, 1]])
        return Rz @ Ry @ Rx @ V
    if view in ("a", "b", "c"):
        d = M[:, "abc".index(view)]
    elif view == "uvw":
        d = M @ np.asarray(uvw, float)
    elif view == "hkl":
        d = np.linalg.inv(M).T @ np.asarray(hkl, float)          # plane normal (reciprocal vector)
    else:
        d = np.array([0, 0, 1.0])
    z = d / np.linalg.norm(d)
    up = M[:, 2] if view in ("a", "b") else M[:, 1]            # c up when looking along a/b, b up along c
    up = up - (up @ z) * z
    if np.linalg.norm(up) < 1e-6:
        up = np.array([0, 1.0, 0]) - z[1] * z
    y = up / np.linalg.norm(up); x = np.cross(y, z)
    V = np.vstack([x, y, z])
    rx, ry, rz = map(math.radians, rot)
    Rx = np.array([[1, 0, 0], [0, math.cos(rx), -math.sin(rx)], [0, math.sin(rx), math.cos(rx)]])
    Ry = np.array([[math.cos(ry), 0, math.sin(ry)], [0, 1, 0], [-math.sin(ry), 0, math.cos(ry)]])
    Rz = np.array([[math.cos(rz), -math.sin(rz), 0], [math.sin(rz), math.cos(rz), 0], [0, 0, 1]])
    return Rz @ Ry @ Rx @ V


def _shade(c, f):
    c = np.array([int(c[i:i + 2], 16) for i in (1, 3, 5)]) / 255
    return tuple(np.clip(c * f + (1 - f) * (0 if f < 1 else 1) * 0, 0, 1))


def draw(ax, S, sites, bonds, style="ortep", prob=0.5, labels="non-H", hbonds=(), show_cell=True, triad=True, view=None,
         polyhedra=False, h_radius=0.12, label_fs=6, bond_lw=None, cell_origin=(0, 0, 0), cells=(1, 1, 1), principal=True,
         label_filter=None, label_offsets=None, label_repel=True, h_scale=1.0):
    """sites/bonds from Structure.grow; hbonds = [(i_H, i_A)] site indices (dashed)."""
    V = view if view is not None else view_matrix(S, "a")
    h_radius = h_radius * h_scale
    X = np.array([s["xyz"] for s in sites]) if sites else np.zeros((0, 3))
    P = X @ V.T                                               # screen x, y, depth
    k = prob_scale(prob)
    objs = []                                                 # (depth, kind, payload)
    for i, s in enumerate(sites):
        objs.append((P[i, 2], "atom", i))
    for a, b in bonds:
        objs.append(((P[a, 2] + P[b, 2]) / 2 - 0.01, "bond", (a, b)))
    for h, a in hbonds:
        objs.append(((P[h, 2] + P[a, 2]) / 2 - 0.01, "hb", (h, a)))
    if polyhedra:
        from scipy.spatial import ConvexHull
        for i, s in enumerate(sites):
            if s["el"] not in METALS:
                continue
            lig = [b if a == i else a for a, b in bonds if i in (a, b) and sites[b if a == i else a]["el"] in ("O", "N", "F", "Cl")]
            if len(lig) >= 4:
                try:
                    hull = ConvexHull(X[lig])
                except Exception:
                    continue
                for f in hull.simplices:
                    objs.append((P[[lig[j] for j in f], 2].mean(), "face", ([lig[j] for j in f], s["el"])))
    objs.sort(key=lambda o: o[0])
    lw = bond_lw or {"ortep": 1.2, "ball": 2.2, "stick": 4.0, "space": 0}[style]
    for depth, kind, pl in objs:
        if kind == "face":
            idx, el = pl
            ax.add_patch(Polygon(P[idx, :2], closed=True, fc=COLOR.get(el, "#88AACC"), ec="k", lw=0.3, alpha=0.45, zorder=1))
        elif kind == "bond" and style != "space":
            a, b = pl
            pa, pb = P[a, :2], P[b, :2]
            if style == "stick":
                m = (pa + pb) / 2
                for p0, p1, s_ in ((pa, m, sites[a]), (m, pb, sites[b])):
                    hh = s_["el"] == "H"
                    ax.plot(*np.array([p0, p1]).T, color="#C8C8C8" if hh else COLOR.get(s_["el"], "#AAAAAA"), lw=lw * (0.6 if hh else 1),
                            solid_capstyle="butt" if hh else "round", zorder=2)
            else:
                ax.plot(*np.array([pa, pb]).T, color="k" if style == "ortep" else "#666666", lw=lw * (0.6 if "H" in (sites[a]["el"], sites[b]["el"]) else 1),
                        solid_capstyle="round", zorder=2)
        elif kind == "hb":
            h, a = pl
            ax.plot(*np.array([P[h, :2], P[a, :2]]).T, color="#1f77b4", lw=0.8, ls=(0, (2.5, 1.5)), zorder=2)
        elif kind == "atom":
            s = sites[pl]; c = COLOR.get(s["el"], "#FF66CC"); xy = P[pl, :2]
            if style == "space":
                r = VDW.get(s["el"], 1.8)
                ax.add_patch(Circle(xy, r, fc=c, ec="k", lw=0.3, zorder=3))
                ax.add_patch(Circle(xy + 0.3 * r * np.array([-0.4, 0.4]), 0.35 * r, fc="white", ec="none", alpha=0.35, zorder=3))
                continue
            if style == "stick":
                if s["el"] in METALS:
                    ax.add_patch(Circle(xy, 0.35, fc=c, ec="k", lw=0.3, zorder=3))
                continue
            if style == "ball":
                r = 0.18 * h_scale if s["el"] == "H" else 0.25 + 0.25 * COV.get(s["el"], 1.0)
                ax.add_patch(Circle(xy, r, fc=c, ec="k", lw=0.4, zorder=3))
                ax.add_patch(Circle(xy + 0.35 * r * np.array([-0.5, 0.5]), 0.3 * r, fc="white", ec="none", alpha=0.5, zorder=3))
                continue
            # ORTEP
            if s["el"] == "H" or s["U"] is None:
                r = h_radius if s["el"] == "H" else k * math.sqrt(max(s["uiso"], 1e-4))
                ax.add_patch(Circle(xy, r, fc="white" if s["el"] == "H" else c, ec="k", lw=0.5, zorder=3))
                continue
            U2 = V[:2] @ s["U"] @ V[:2].T
            w, ev = np.linalg.eigh(U2); w = np.clip(w, 1e-6, None)
            ang = math.degrees(math.atan2(ev[1, 1], ev[0, 1]))
            ax.add_patch(Ellipse(xy, 2 * k * math.sqrt(w[1]), 2 * k * math.sqrt(w[0]), angle=ang, fc=c, ec="k", lw=0.6, zorder=3))
            if principal:
                w3, e3 = np.linalg.eigh(s["U"])
                for j in range(3):
                    v = (V[:2] @ e3[:, j]) * k * math.sqrt(max(w3[j], 1e-6))
                    ax.plot([xy[0] - v[0], xy[0] + v[0]], [xy[1] - v[1], xy[1] + v[1]], color="k", lw=0.35, zorder=3)
                # front octant shading, the classic ORTEP look
                ax.add_patch(Ellipse(xy + 0.25 * np.array([-k * math.sqrt(w[1]), k * math.sqrt(w[0])]) * 0.5, k * math.sqrt(w[1]) * 0.7,
                                     k * math.sqrt(w[0]) * 0.7, angle=ang, fc="white", ec="none", alpha=0.25, zorder=3))
    if show_cell:
        o = np.asarray(cell_origin, float); n = np.asarray(cells, float)
        corners = np.array([[i, j, k_] for i in (0, n[0]) for j in (0, n[1]) for k_ in (0, n[2])]) + o
        C = S.cart(corners) @ V.T
        for p, q in [(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)]:
            ax.plot(*C[[p, q], :2].T, color="#777777", lw=0.5, zorder=0.5)
    ax.set_aspect("equal"); ax.set_axis_off()
    allp = P[:, :2] if len(P) else np.zeros((1, 2))
    if show_cell:
        allp = np.vstack([allp, C[:, :2]])
    lo, hi = allp.min(0) - 1.2, allp.max(0) + 1.2
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
    # labels on top: placed away from the bonded neighbours, then a repel pass (labels vs labels and vs atoms)
    if labels != "none":
        idx = [i for i, s in enumerate(sites) if not (labels == "non-H" and s["el"] == "H") and not (label_filter and not label_filter(s))]
        texts = [sites[i]["label"] + ("" if sites[i]["op"] == 0 and sites[i]["T"] == (0, 0, 0) else "'") for i in idx]
        pos = place_labels(P[:, :2], idx, texts, bonds, [label_fs * (0.8 if sites[i]["el"] == "H" else 1) for i in idx],
                           _a_per_pt(ax, lo, hi), repel=label_repel,
                           radii=[0.15 * h_scale if s_['el'] == 'H' else 0.42 for s_ in sites])
        for j, (i, t) in enumerate(zip(idx, texts)):
            dx, dy = (label_offsets or {}).get(sites[i]["label"], (0.0, 0.0))
            if np.linalg.norm(pos[j] - P[i, :2]) > 1.6:          # pushed far away: thin leader line
                ax.plot([P[i, 0], pos[j, 0] + dx], [P[i, 1], pos[j, 1] + dy], color="#888888", lw=0.3, zorder=4)
            ax.text(pos[j, 0] + dx, pos[j, 1] + dy, t, fontsize=label_fs if sites[i]["el"] != "H" else label_fs * 0.8, zorder=5,
                    color="k", ha="center", va="center")
    if triad:
        L = 0.09 * max(hi - lo); o2 = lo + 0.06 * (hi - lo)
        for j, (nm, col) in enumerate(zip("abc", ("#D62728", "#2CA02C", "#1F77B4"))):
            v = V[:2] @ (S.M[:, j] / np.linalg.norm(S.M[:, j]))
            if np.linalg.norm(v) < 0.08:
                ax.plot(*o2, marker="o", ms=3, color=col, zorder=6); ax.text(*(o2 + [0.03 * L * 6, 0]), nm, color=col, fontsize=7, zorder=6)
                continue
            ax.add_patch(FancyArrowPatch(o2, o2 + v * L, arrowstyle="-|>", mutation_scale=7, color=col, lw=0.8, zorder=6))
            ax.text(*(o2 + v * L * 1.25), nm, color=col, fontsize=7, ha="center", va="center", zorder=6)
    return P


def polar_h_only(sites, bonds):
    """Drop H atoms bonded to C (keep O–H / N–H for H-bond figures). Returns (sites, bonds) with re-indexed bonds."""
    ch = set()
    for a, b in bonds:
        for h, o in ((a, b), (b, a)):
            if sites[h]["el"] == "H" and sites[o]["el"] == "C":
                ch.add(h)
    keep = [i for i in range(len(sites)) if i not in ch]; new = {o: n for n, o in enumerate(keep)}
    return [sites[i] for i in keep], [(new[a], new[b]) for a, b in bonds if a in new and b in new]


def _a_per_pt(ax, lo, hi):
    """Data units (Å) per typographic point for an equal-aspect axes with limits lo..hi."""
    fw, fh = ax.figure.get_size_inches(); b = ax.get_position()
    w_in, h_in = b.width * fw, b.height * fh
    return max((hi[0] - lo[0]) / w_in, (hi[1] - lo[1]) / h_in) / 72.0


def place_labels(P2, idx, texts, bonds, fss, aa, repel=True, n_iter=150, atom_r=0.42, radii=None):
    """Label centres (screen Å). Start: opposite to the mean direction of the bonded neighbours, ~0.55 Å + half the text
    width from the atom; then iterative repulsion between overlapping label boxes and from atom centres."""
    if not idx:
        return np.zeros((0, 2))
    W = np.array([0.58 * f * len(t) * aa for t, f in zip(texts, fss)]); H = np.array([1.05 * f * aa for f in fss])
    nb = {i: [] for i in idx}
    for a, b in bonds:
        if a in nb: nb[a].append(b)
        if b in nb: nb[b].append(a)
    out = np.zeros((len(idx), 2))
    for j, i in enumerate(idx):
        v = np.zeros(2)
        for m in nb[i]:
            d = P2[m] - P2[i]; n = np.linalg.norm(d)
            if n > 1e-6:
                v += d / n
        d = -v if np.linalg.norm(v) > 0.2 else np.array([0.7, 0.7])
        if np.linalg.norm(v) <= 0.2 and nb[i]:            # neighbours on both sides: perpendicular to the first bond
            b0 = P2[nb[i][0]] - P2[i]; d = np.array([-b0[1], b0[0]])
        d = d / np.linalg.norm(d) if np.linalg.norm(d) > 1e-9 else np.array([0.7071, 0.7071])
        r = atom_r + 0.08 + 0.5 * abs(d[0]) * W[j] + 0.5 * abs(d[1]) * H[j]
        out[j] = P2[i] + d * r
    if not repel:
        return out
    anchor = P2[idx]
    RA = np.full(len(P2), atom_r) if radii is None else np.asarray(radii)
    for _ in range(n_iter):
        moved = False
        for j in range(len(idx)):
            f = np.zeros(2)
            dx = out[j, 0] - out[:, 0]; dy = out[j, 1] - out[:, 1]
            ox = (W[j] + W) / 2 - np.abs(dx); oy = (H[j] + H) / 2 - np.abs(dy)
            hit = (ox > 0) & (oy > 0); hit[j] = False
            for k in np.where(hit)[0]:
                if ox[k] < oy[k]:
                    f[0] += np.sign(dx[k] or 1) * ox[k] * 0.5
                else:
                    f[1] += np.sign(dy[k] or 1) * oy[k] * 0.5
            ax_ = out[j, 0] - P2[:, 0]; ay_ = out[j, 1] - P2[:, 1]       # keep clear of the atoms
            ox = W[j] / 2 + RA - np.abs(ax_); oy = H[j] / 2 + RA - np.abs(ay_)
            for k in np.where((ox > 0) & (oy > 0))[0]:
                if ox[k] < oy[k]:
                    f[0] += np.sign(ax_[k] or 1) * ox[k] * 0.3
                else:
                    f[1] += np.sign(ay_[k] or 1) * oy[k] * 0.3
            f += (anchor[j] - out[j]) * 0.04                 # weak spring back to its atom
            if np.linalg.norm(f) > 1e-3:
                out[j] += np.clip(f, -0.3, 0.3); moved = True
        if not moved:
            break
    return out


def label_overlaps(pos, texts, fss, aa):
    """Number of overlapping label-box pairs (for tests / a warning in the app)."""
    W = np.array([0.58 * f * len(t) * aa for t, f in zip(texts, fss)]); H = np.array([1.05 * f * aa for f in fss])
    n = 0
    for j in range(len(pos)):
        for k in range(j + 1, len(pos)):
            if abs(pos[j, 0] - pos[k, 0]) < (W[j] + W[k]) / 2 * 0.95 and abs(pos[j, 1] - pos[k, 1]) < (H[j] + H[k]) / 2 * 0.95:
                n += 1
    return n
