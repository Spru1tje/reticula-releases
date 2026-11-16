"""v2.12 sphere model on the user's real SiO2-NH2 TEM images (manual sizing: 209 ± 19 nm, N≈100) and synthetic cases."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import numpy as np, pytest
from mfs import readers as R, particles as PA, spheres as SPH, batchproc as BP

D = os.path.join(ROOT, "tests", "data")


def _run(fn, **kw):
    f = os.path.join(D, fn); d = R.read_sem(open(f, "rb"), f); a = d["image"]
    lo, hi = np.percentile(a, (0.2, 99.8)); a8 = (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
    lab, fac, _, info = PA.segment_adaptive(a8, **kw)
    rows = PA.measure(lab, d["pixel_nm"], fac, info=info)
    e = np.array([r["ECD (nm)"] for r in rows if info.get(r["id"], {}).get("kind") == "particle"])
    return e, info, rows


@pytest.mark.parametrize("fn,mean,n_min", [("TEM_SiO2NH2_user6_20k.tif", 209.0, 45), ("TEM_SiO2NH2_user6_60k3.tif", 213.3, 12),
                                           ("TEM_SiO2NH2_user6_100k2.tif", 201.9, 5)])
def test_real_spheres_mean_vs_manual(fn, mean, n_min):
    e, info, rows = _run(fn)
    s = info["_sphere"]
    assert s["active"] and s["roundness"]["round"], s
    assert len(e) >= n_min
    assert abs(e.mean() - mean) < 6, e.mean()                 # regression on the v2.12 numbers
    assert abs(e.mean() - 209) < 15                           # manual sizing 209 ± 19 nm
    assert 5 < e.std(ddof=1) < 30
    assert any("confidence" in r for r in rows)


def test_v210_overestimate_removed_20k():
    e_off, _, _ = _run("TEM_SiO2NH2_user6_20k.tif", sphere=False)
    e_on, info, _ = _run("TEM_SiO2NH2_user6_20k.tif")
    assert e_off.mean() > 235 and e_on.mean() < 215           # v2.10 247 nm → v2.12 209 nm
    s = info["_sphere"]
    assert s["n_split"] >= 30 and s["n_occluded"] >= 2 and s["n_lowconf"] >= 3


def test_lowconf_excluded_and_conf_threshold():
    _, i_def, r_def = _run("TEM_SiO2NH2_user6_100k2.tif")
    _, i_hi, _ = _run("TEM_SiO2NH2_user6_100k2.tif", conf_min=0.9)
    assert i_hi["_sphere"]["n_lowconf"] > i_def["_sphere"]["n_lowconf"]
    from mfs import particle_picker as PK
    for r in r_def:
        r["kind"] = i_def.get(r["id"], {}).get("kind", "particle")
    lc = [r["id"] for r in r_def if r["kind"] == "lowconf"]
    assert lc and set(lc) <= set(PK.auto_excluded(r_def))


def test_synthetic_touching_and_overlapping_spheres():
    """Blurred absorbing spheres: isolated, a touching pair, a touching triple and an overlapping pair (darker lens)."""
    H, W = 420, 640; yy, xx = np.mgrid[:H, :W]; R0 = 30.0
    cs = [(80, 80), (220, 80), (280, 80), (420, 80), (480, 80), (450, 132), (120, 300), (160, 310), (400, 300)]
    t = np.zeros((H, W))
    for cx, cy in cs:
        t += 2 * np.sqrt(np.clip(R0 ** 2 - (xx - cx) ** 2 - (yy - cy) ** 2, 0, None))
    from scipy import ndimage as ndi
    im = ndi.gaussian_filter(220 * np.exp(-t / 25.0), 1.2) + np.random.default_rng(1).normal(0, 3, (H, W))
    lab, f, _, info = PA.segment_adaptive(np.clip(im, 0, 255).astype(np.uint8), sphere=True)
    circ = [v for k, v in info.items() if isinstance(k, (int, np.integer)) and v.get("circle")]
    found = [v for v in circ if v["kind"] == "particle"]
    assert len(circ) == len(cs), len(circ)
    for cx, cy in cs:
        best = min(circ, key=lambda v: np.hypot(v["circle"][0] - cx, v["circle"][1] - cy))
        assert np.hypot(best["circle"][0] - cx, best["circle"][1] - cy) < 4
    d = np.array([v["d_px"] for v in found])
    assert len(found) == len(cs) and np.ptp(d) < 1.5, d        # touching/overlapping fits = isolated fit (no split bias)
    assert abs(np.median(d) - 2 * R0) < 0.12 * 2 * R0, d       # absolute edge level is calibrated on real TEM, not on this toy model


def test_cubes_unchanged_by_auto():
    f = os.path.join(D, "TEM_MOF801_user5_sample1.tif")
    if not os.path.exists(f):
        pytest.skip("cube image missing")
    d = R.read_sem(open(f, "rb"), f); a = d["image"]
    lo, hi = np.percentile(a, (0.2, 99.8)); a8 = (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
    l1, _, _, i1 = PA.segment_adaptive(a8, sphere="auto"); l2, _, _, _ = PA.segment_adaptive(a8, sphere=False)
    assert not i1["_sphere"]["active"] and np.array_equal(l1, l2)


def test_batch_particles_spheres():
    o = BP.particles_one(os.path.join(D, "TEM_SiO2NH2_user6_60k3.tif"), "x.tif")
    assert o["summary"]["spheres"] and abs(o["summary"]["mean"] - 213) < 8 and o["summary"]["low_confidence"] >= 1
