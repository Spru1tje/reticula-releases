"""v2.21: net matcher. Given indexed reflections (hkl, |F|^2) of a framework, score each net in the own
topology library by how well its node-only structure factor reproduces the observed intensities.

The node model has one scatterer per net node at the node's fractional position (uniform weight). The
scale is fitted by least squares on amplitudes; the score is the weighted R on |F| (lower is better).
A net with no good match is reported as 'low confidence' rather than forced to a winner."""
import numpy as np
from . import topology as TP


def node_types(name):
    """List of node groups (fractional positions). Nets with two node types (tbo: paddlewheel 4-c and ring 3-c)
    return two groups whose relative weight is fitted; the others return one group."""
    if name == "tbo":
        import json, os
        path = os.path.join(os.path.dirname(TP.__file__), "data", "tbo_nodes.json")
        d = json.load(open(path, encoding="utf-8"))
        return [np.array(d["four_c"], float), np.array(d["three_c"], float)]
    return [TP.LIB[name]()[0]]


def _group_sums(name, hkl):
    h = np.asarray(hkl, float)
    return [np.exp(2j * np.pi * (h @ g.T)).sum(axis=1) for g in node_types(name)]


def node_amplitude(name, hkl, weights=None):
    """|F| of the node model of net `name` for the reflections hkl (N x 3); weights per node group."""
    sums = _group_sums(name, hkl)
    w = np.ones(len(sums)) if weights is None else np.asarray(weights, float)
    return np.abs(sum(wi * s for wi, s in zip(w, sums)))


def score_net(name, hkl, Fsq):
    """Return dict with R (weighted, amplitude), fitted scale and number of reflections used.
    All reflections are kept, including the (near-)absent ones: their zeros are what separates
    e.g. fcu (absences for mixed parity) from pcu."""
    Fo = np.sqrt(np.clip(np.asarray(Fsq, float), 0, None))
    n_types = len(node_types(name))
    grid = [None] if n_types == 1 else [np.array([1.0, r]) for r in np.logspace(-1.5, 1.5, 61)]
    best = None
    for w in grid:
        Fc = node_amplitude(name, hkl, w)
        s = float(np.sum(Fo * Fc) / max(np.sum(Fc * Fc), 1e-30))
        R = float(np.sqrt(np.sum((Fo - s * Fc) ** 2) / np.sum(Fo ** 2)))
        if best is None or R < best["R"]:
            best = {"net": name, "R": R, "scale": s, "n": len(Fo),
                    "weight_ratio": None if w is None else float(w[1])}
    return best


def rank_nets(hkl, Fsq, nets=None, margin=0.05):
    """Rank candidate nets (best first). `confident` is True only when the best net beats the
    runner-up by at least `margin` in R; otherwise the match is ambiguous and must be reported."""
    nets = list(nets or TP.LIB.keys())
    res = sorted((score_net(n, hkl, Fsq) for n in nets), key=lambda r: r["R"])
    confident = len(res) < 2 or (res[1]["R"] - res[0]["R"]) >= margin
    return {"ranking": res, "best": res[0]["net"], "R_best": res[0]["R"], "confident": confident}


def metal_pair_fraction(frac, heights, a, K=48, lo=2.3, hi=2.9):
    """Metal-sublattice test from charge-flip density peaks (no CIF).
    Take the K strongest peaks (the heavy metal), and return the fraction whose nearest neighbour lies in
    [lo, hi] A (paddlewheel Cu-Cu ~2.5-2.6 A) and the median nearest-neighbour distance."""
    frac = np.asarray(frac, float)
    h = np.asarray(heights, float)
    o = np.argsort(-h)[:K]
    f = frac[o]
    d = f[:, None, :] - f[None, :, :]
    d -= np.round(d)
    c = np.linalg.norm(d * a, axis=2)
    np.fill_diagonal(c, np.inf)
    nn = c.min(axis=1)
    return {"frac_pair": float(np.mean((nn >= lo) & (nn <= hi))), "nn_median": float(np.median(nn)), "K": int(len(o))}


def classify_metal_sublattice(frac, heights, a, K=48):
    """'paired' when the heavy-metal peaks come in pairs (paddlewheel nets: only tbo in the library),
    otherwise 'single' (fcu, pcu, dia)."""
    m = metal_pair_fraction(frac, heights, a, K=K)
    return ("paired" if m["frac_pair"] >= 0.8 else "single"), m


def rank_nets_pattern(hkl, Fsq, frac, heights, a, nets=None):
    """Pattern-only matcher: keep only nets whose metal sublattice type matches the peak test, then rank by R."""
    kind, meta = classify_metal_sublattice(frac, heights, a)
    allowed = [n for n in (nets or TP.LIB.keys()) if ((n == "tbo") == (kind == "paired"))]
    H, F = complete_reflections(hkl, Fsq, a)
    res = rank_nets(H, F, nets=allowed)
    res.update({"metal_type": kind, "metal_test": meta, "candidates": allowed})
    return res


def complete_reflections(hkl_obs, Fsq_obs, a, dmin=None):
    """Full reflection list up to the resolution of the observed data, with 0 where nothing was observed.
    Absences are the evidence that separates nets (e.g. fcu vs pcu); a list of observed peaks alone hides them."""
    hkl_obs = np.rint(np.asarray(hkl_obs, float)).astype(int)
    Fsq_obs = np.asarray(Fsq_obs, float)
    if dmin is None:
        dmin = a / np.sqrt((np.abs(hkl_obs) ** 2).sum(1)).max()
    hmax = int(np.floor(a / dmin))
    g = np.arange(-hmax, hmax + 1)
    H = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    H = H[(np.abs(H).sum(1) > 0) & (np.sqrt((H ** 2).sum(1)) <= a / dmin)]
    lookup = {tuple(x): v for x, v in zip(hkl_obs, Fsq_obs)}
    Fsq = np.array([lookup.get(tuple(x), 0.0) for x in H])
    return H.astype(float), Fsq
