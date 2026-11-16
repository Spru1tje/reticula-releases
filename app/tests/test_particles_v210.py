"""v2.10 particle recognition: background flattening + hysteresis threshold, particle-scale markers, optional learned
shape prior, clean overlay. Synthetic TEM-like images + the user's real TEM image (hand-checked reference)."""
import os, json, math, time, inspect
import numpy as np
from scipy import ndimage as ndi
from mfs import particles as PA, shapeprior as SP, particle_picker as PK

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")


def _square(H, W, cx, cy, side, ang):
    yy, xx = np.mgrid[:H, :W]
    c, s = math.cos(ang), math.sin(ang)
    u = (xx - cx) * c + (yy - cy) * s; v = -(xx - cx) * s + (yy - cy) * c
    return (np.abs(u) <= side / 2) & (np.abs(v) <= side / 2)


def synth_tem(seed=0, side=26):
    """Uneven background (gradient 230→185) with dark and light-grey rotated squares: isolated ones, touching pairs
    and one 2×2 block. Returns (img uint8, list of (cx, cy, kind)) with kind 'iso' / 'pair' / 'block'."""
    rng = np.random.default_rng(seed)
    H, W = 420, 640
    yy, xx = np.mgrid[:H, :W]
    bg = 230 - 45 * xx / W - 10 * yy / H
    absorb = np.zeros((H, W))
    truth = []
    grid = [(60 + 70 * i, 60 + 70 * j) for i in range(8) for j in range(3)]
    for k, (cx, cy) in enumerate(grid):
        dark = 0.55 if k % 2 == 0 else 0.22                 # every 2nd square is light grey
        a = rng.uniform(0, math.pi / 2)
        absorb += dark * _square(H, W, cx + rng.uniform(-4, 4), cy + rng.uniform(-4, 4), side * rng.uniform(0.92, 1.08), a)
        truth.append((cx, cy, "iso"))
    for k, x0 in enumerate((80, 260, 440)):                  # touching pairs (face to face)
        a = rng.uniform(0, 0.3)
        for dx in (0, side * 1.02):
            absorb += 0.45 * _square(H, W, x0 + dx * math.cos(a), 300 + dx * math.sin(a), side, a)
            truth.append((x0 + dx * math.cos(a), 300 + dx * math.sin(a), "pair"))
    for (dx, dy) in ((0, 0), (1, 0), (0, 1), (1, 1)):        # a 2×2 block
        absorb += 0.4 * _square(H, W, 560 + dx * side * 1.02, 330 + dy * side * 1.02, side, 0.1)
        truth.append((560 + dx * side * 1.02, 330 + dy * side * 1.02, "block"))
    img = bg * np.exp(-absorb) + rng.normal(0, 6, (H, W))
    return np.clip(img, 0, 255).astype(np.uint8), truth


def _recall(lab, info, truth, f=1, kinds=("iso",)):
    inc = {k for k, v in info.items() if isinstance(k, int) and v.get("kind") != "aggregate"}
    pts = [(x, y) for x, y, k in truth if k in kinds]
    owner = [int(lab[int(y / f), int(x / f)]) for x, y in pts]
    cnt = {o: owner.count(o) for o in owner}
    return sum(1 for o in owner if o and o in inc and cnt[o] == 1) / len(pts)


def test_adaptive_finds_light_grey_particles_on_uneven_background():
    img, truth = synth_tem()
    lab, f, t, info = PA.segment_adaptive(img)
    r_new = _recall(lab, info, truth)
    l0, f0, t0, i0 = PA.segment_ws(img, True)
    r_old = _recall(l0, i0, truth, f0)
    light = [(x, y, "iso") for k, (x, y, kk) in enumerate(truth[:24]) if k % 2 == 1]
    assert r_new >= 0.95, (r_new, r_old)
    assert _recall(lab, info, light) >= 0.9                      # the light-grey squares are caught
    assert r_new > r_old
    assert abs(info["_side_px"] - 26) < 6                        # particle scale estimated from the image


def test_adaptive_splits_touching_pairs_and_size():
    img, truth = synth_tem(1)
    lab, f, t, info = PA.segment_adaptive(img)
    assert _recall(lab, info, truth, kinds=("pair",)) >= 0.66
    rows = PA.measure(lab, 1.0, f)
    ecd = np.median([r["ECD (nm)"] for r in rows if info[r["id"]]["kind"] == "particle"])
    assert abs(ecd - 26 * 2 / math.sqrt(math.pi)) < 5            # ECD of a 26 px square ≈ 29.3 px


def test_shape_prior_learns_square_and_weight_formula():
    img, truth = synth_tem(2)
    lab, f, t, info = PA.segment_adaptive(img, clear_border=False)
    c = info["_contrast"]; L = info["_side_px"]
    pr1 = SP.learn(lab, info, c, L, k=10, strength=1.0)
    assert pr1["n"] >= 15 and pr1["corners"] == 4, (pr1["n"], pr1.get("corners"))
    assert pr1["consistency"] > 0.5
    assert abs(pr1["weight"] - pr1["n"] / (pr1["n"] + 10) * pr1["consistency"]) < 1e-9
    pr2 = SP.learn(lab, info, c, L, k=100, strength=1.0)        # larger k → weaker prior
    pr3 = SP.learn(lab, info, c, L, k=10, strength=0.5)
    assert pr2["weight"] < pr1["weight"] and abs(pr3["weight"] - 0.5 * pr1["weight"]) < 1e-9
    assert len(pr1["crops"]) >= 6 and len(pr1["template"]) == SP.NB


def test_shape_prior_off_by_default_and_breaks_block():
    assert inspect.signature(PA.segment_adaptive).parameters["shape_prior"].default in (None, False)
    img, truth = synth_tem(3)
    _, _, _, off = PA.segment_adaptive(img)
    lab, f, t, on = PA.segment_adaptive(img, shape_prior=True, prior_strength=1.0)
    pr = on["_prior"]
    assert pr["applied"] and pr["n"] >= 10
    assert not any(isinstance(v, dict) and v.get("split") == "prior" for v in off.values())
    n_pr = sum(1 for v in on.values() if isinstance(v, dict) and v.get("split") == "prior")
    agg_off = sum(1 for v in off.values() if isinstance(v, dict) and v.get("kind") == "aggregate")
    agg_on = sum(1 for v in on.values() if isinstance(v, dict) and v.get("kind") == "aggregate")
    assert agg_on <= agg_off and n_pr >= 2
    l_off = PA.segment_adaptive(img)[0]
    r_blk_off = _recall(l_off, off, truth, kinds=("block",)); r_blk_on = _recall(lab, on, truth, kinds=("block",))
    assert r_blk_on >= 0.75 and r_blk_on > r_blk_off, (r_blk_off, r_blk_on)   # the fused 2×2 block is split into cubes
    assert _recall(lab, on, truth) >= 0.95                        # isolated particles are untouched


def test_shape_score_and_template():
    t = np.ones(SP.NB)
    sq = _square(80, 80, 40, 40, 30, 0.3)
    p, _ = SP.radial_profile(sq)
    pr = {"template": p / math.sqrt((p ** 2).mean())}
    assert SP.shape_score(_square(80, 80, 40, 40, 30, 0.9), pr) > 0.9       # any rotation of the same square
    yy, xx = np.mgrid[:80, :80]
    assert SP.shape_score((xx - 40) ** 2 / 900 + (yy - 40) ** 2 / 100 < 1, pr) < 0.8   # an elongated ellipse
    m = SP.template_mask(pr["template"], 15, 20)
    assert 0.7 * 900 < m.sum() < 1.3 * 900


def _real():
    from mfs import readers as R
    p = os.path.join(D, "TEM_MOF801_user5_sample1.tif")
    d = R.read_sem(open(p, "rb"), os.path.basename(p))
    return d["image"][:976], d["pixel_nm"]


def _eval(lab, info, f=1):
    ref = json.load(open(os.path.join(D, "TEM_MOF801_user5_reference.json")))
    boxes = json.load(open(os.path.join(D, "TEM_MOF801_user5_dense_boxes.json")))
    inb = lambda x, y: any(a <= x <= c and b <= y <= d for a, b, c, d in boxes)
    ref = [p for p in ref if not inb(*p)]
    inc = {k for k, v in info.items() if isinstance(k, int) and v.get("kind") != "aggregate"}
    own = [int(lab[int(y / f), int(x / f)]) for x, y in ref]
    cnt = {o: own.count(o) for o in own}
    good = {o for o in own if o and o in inc and cnt[o] == 1}
    ids = sorted(inc & set(np.unique(lab).tolist()))
    cm = ndi.center_of_mass(lab > 0, lab, ids)
    outside = [i for i, (cy, cx) in zip(ids, cm) if not inb(cx * f, cy * f)]
    return len(good) / len(ref), sum(1 for i in outside if i in good) / max(1, len(outside))


def test_real_tem_regression_recall_precision():
    img, pix = _real()
    assert abs(pix - 1.1174) < 1e-3                               # Olympus SIS metadata calibration
    t0 = time.time(); lab, f, t, off = PA.segment_adaptive(img); dt = time.time() - t0
    r_off, p_off = _eval(lab, off, f)
    l0, f0, t0_, old = PA.segment_ws(img, True)
    r_old, p_old = _eval(l0, old, f0)
    lab2, f2, t2, on = PA.segment_adaptive(img, shape_prior=True)
    r_on, p_on = _eval(lab2, on, f2)
    assert r_off >= 0.75 and p_off >= 0.85, (r_off, p_off)
    assert r_on >= r_off and r_on >= 0.82 and p_on >= 0.85, (r_on, p_on)
    assert r_off > r_old + 0.1, (r_old, r_off)
    assert 26 < off["_side_px"] < 36                              # MOF-801 cubes ≈ 31 px ≈ 35 nm
    assert on["_prior"]["applied"] and on["_prior"]["corners"] in (4, 6) and on["_prior"]["n"] >= 10
    assert dt < 20


def test_picker_component_v210_style():
    html = open(os.path.join(PK.component_dir(), "index.html"), encoding="utf-8").read()
    assert 'id="nums"' in html and "body.nums .n" in html and ".n{display:none" in html   # numbers off by default
    assert ".p.pr{stroke-dasharray" in html and "stroke-width:1.1" in html
    assert "<title" in html or 'createElementNS(NS,"title")' in html
    lab = np.zeros((40, 40), int); lab[5:15, 5:15] = 1; lab[20:30, 20:30] = 2
    rows = [{"id": 1, "x (px)": 10, "y (px)": 10, "ECD (nm)": 11.0, "area (nm²)": 100.0, "kind": "particle", "prior": True},
            {"id": 2, "x (px)": 25, "y (px)": 25, "ECD (nm)": 11.0, "area (nm²)": 100.0, "kind": "aggregate"}]
    P = PK.particles_payload(rows, lab, 1, 1.0, 1.0)
    assert P[0]["prior"] is True and P[1]["prior"] is False


def test_overlay_numbers_off_by_default():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    assert inspect.signature(PA.overlay_labels).parameters["numbers"].default is False
    lab = np.zeros((60, 60), int); lab[5:25, 5:25] = 1; lab[30:55, 30:55] = 2
    fig, ax = plt.subplots(); PA.overlay_labels(ax, np.full((60, 60), 128, np.uint8), lab, 1, {1: {"kind": "particle"}, 2: {"kind": "aggregate"}})
    assert not ax.texts
    plt.close(fig)
