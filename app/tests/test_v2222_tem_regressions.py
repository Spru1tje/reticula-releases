"""v2.22.2 regression tests for the TEM particle module: exclusion circles, picture-edge and film-rim exclusion, the
threshold cap (99.5th percentile) and the per-sample agglomerate tags. Synthetic images only; the real micrographs are
covered by test_particles_v210.py and test_real_olympus_tem.py."""
import numpy as np
import pytest

from mfs import rods as R
from mfs import particles as PA
from mfs import agglo_learn as AL
from mfs import semcal as SC


def _obj(cx, cy, **kw):
    o = {"id": 0, "cx": float(cx), "cy": float(cy), "area_px": 40, "length_nm": 50.0, "width_nm": 5.0, "aspect": 10.0,
         "solidity": 0.95, "film_edge": False, "cut_by_edge": False, "agglomerate": False}
    o.update(kw)
    return o


def _dark_rod(im, y0, x0, y1, x1, w=4, val=120):
    """Paint a dark straight rod (width w px) into im, in place."""
    H, W = im.shape
    n = int(max(abs(y1 - y0), abs(x1 - x0))) * 2 + 2
    for t in np.linspace(0, 1, n):
        y = int(round(y0 + t * (y1 - y0))); x = int(round(x0 + t * (x1 - x0)))
        im[max(0, y - w // 2): min(H, y + w // 2 + 1), max(0, x - w // 2): min(W, x + w // 2 + 1)] = val
    return im


def _film(seed, H=300, W=400, noise=3.0):
    rng = np.random.default_rng(seed)
    return 200 + noise * rng.standard_normal((H, W))


# ------------------------------------------------------------------ exclusion circles
def test_exclusion_circle_flags_centres_inside_only_and_boundary_counts():
    objs = [_obj(10, 10), _obj(40, 10), _obj(13, 14), _obj(16, 10)]     # centre; far; on the rim (d=5); just outside (d=6)
    R.in_circles(objs, [(10.0, 10.0, 5.0)])
    assert [o["manual_excluded"] for o in objs] == [True, False, True, False]
    s = R.summary(objs)
    assert s["n"] == 2 and s["n_manual_excluded"] == 2


def test_no_circles_excludes_nothing_and_summary_ignores_them():
    objs = R.in_circles([_obj(10, 10), _obj(20, 20)], [])
    assert not any(o["manual_excluded"] for o in objs)
    assert R.summary(objs)["n"] == 2
    assert R.summary([])["n"] == 0


# ------------------------------------------------------------------ randuitsluiting: picture edge (rods)
def test_rod_cut_by_picture_edge_is_not_counted():
    im = _dark_rod(_film(0), 150, 120, 150, 200)                       # interior rod
    _dark_rod(im, 60, 0, 62, 40)                                       # rod cut by the left picture edge
    _dark_rod(im, 230, 300, 260, 330)                                  # interior rod
    lab, objs = R.segment_rods(im, 1.0)
    assert len(objs) == 3
    assert sorted(o["cut_by_edge"] for o in objs) == [False, False, True]
    s = R.summary(objs)
    assert s["n"] == 2 and s["n_cut_by_edge"] == 1


# ------------------------------------------------------------------ randuitsluiting: film rim (rods)
def test_rod_in_the_rim_of_a_carbon_hole_is_flagged_not_counted():
    im = _film(1)
    yy, xx = np.mgrid[:im.shape[0], :im.shape[1]]
    im[(yy - 100) ** 2 + (xx - 300) ** 2 < 30 ** 2] = 255                # bright carbon-film hole, radius 30 px
    _dark_rod(im, 85, 335, 115, 335, w=4)                                # rod 33-37 px from the centre: inside the 10 px rim
    _dark_rod(im, 200, 60, 200, 100, w=4)                                # rod far from any hole
    lab, objs = R.segment_rods(im, 1.0)
    rim = [o for o in objs if o["film_edge"]]
    free = [o for o in objs if not o["film_edge"]]
    assert len(rim) == 1 and abs(rim[0]["cx"] - 335) < 3
    assert len(free) == 1 and abs(free[0]["cx"] - 80) < 3
    s = R.summary(objs)
    assert s["n"] == 1 and s["n_film_edge"] == 1


def test_minimum_area_drops_small_objects():
    im = _film(2)
    im[150:155, 150:155] = 120                                           # 5 x 5 = 25 px square
    _, small = R.segment_rods(im, 1.0, min_area_px=12)
    _, big = R.segment_rods(im, 1.0, min_area_px=60)
    assert any(abs(o["cx"] - 152) < 3 for o in small)
    assert not any(abs(o["cx"] - 152) < 3 for o in big)


# ------------------------------------------------------------------ randuitsluiting: particles cut by the image edge
def test_particle_cut_by_image_edge_is_dropped_with_clear_border():
    H, W = 300, 400
    yy, xx = np.mgrid[:H, :W]
    im = 170 - 0.05 * xx + 4 * np.random.default_rng(3).standard_normal((H, W))
    im[(np.abs(xx - 200) <= 13) & (np.abs(yy - 150) <= 13)] -= 70        # interior square
    im[(np.abs(xx - 2) <= 13) & (np.abs(yy - 60) <= 13)] -= 70           # square cut by the left edge
    im = np.clip(im, 0, 255).astype(np.uint8)
    lab, f, hi, info = PA.segment_adaptive(im, sphere=False, clear_border=True)
    assert not (lab[:, 0].any() or lab[:, -1].any() or lab[0].any() or lab[-1].any())
    assert info["_edge_excluded"] == 1
    assert lab[150, 200] > 0                                             # the interior particle is kept
    lab2, _, _, info2 = PA.segment_adaptive(im, sphere=False, clear_border=False)
    assert len([k for k in info2 if isinstance(k, int)]) == len([k for k in info if isinstance(k, int)]) + 1


# ------------------------------------------------------------------ drempelcap (99.5th percentile)
def test_threshold_is_capped_at_99_5th_percentile_on_textured_film():
    # pure noise film (no particles): the noise estimate puts 8·σ above every contrast value. v2.22.1 then found 0 objects;
    # v2.22.2 caps the seed level at the 99.5th percentile of the contrast, so the picture is segmented again.
    img = np.clip(_film(4, noise=12.0), 0, 255).astype(np.uint8)
    lab, f, hi, info = PA.segment_adaptive(img, sphere=False)
    c = info["_contrast"]
    p995 = float(np.percentile(c, 99.5))
    assert 8 * info["_noise"] > p995                                     # the condition that used to break the segmentation
    assert hi <= p995 + 1e-12 and np.isclose(hi, p995)
    assert info["_lo"] <= 0.7 * hi + 1e-12
    assert len([k for k in info if isinstance(k, int)]) > 0


# ------------------------------------------------------------------ per-sample agglomerate tags
def test_tag_store_is_one_file_per_sample(tmp_path, monkeypatch):
    monkeypatch.setenv("RETICULA_DATA", str(tmp_path))
    pa, pb = AL.default_store("Sample A"), AL.default_store("Sample B")
    assert pa != pb and pa.parent == pb.parent == tmp_path / "tags"
    assert AL.default_store("a/b c").name == "a_b_c.json"
    assert AL.default_store(None).name == "default.json"


def test_tags_do_not_leak_between_samples(tmp_path, monkeypatch):
    monkeypatch.setenv("RETICULA_DATA", str(tmp_path))
    pa, pb = AL.default_store("Sample A"), AL.default_store("Sample B")
    o = {"area_px": 50, "length_nm": 5, "width_nm": 2, "aspect": 2.5, "solidity": 0.9, "cx": 10.0, "cy": 20.0}
    AL.save_tags(AL.add_tag([], "img_a.tif", o, 0.274, 50.0, 0), pa)
    assert len(AL.load_tags(pa)) == 1
    assert AL.load_tags(pb) == []                                        # a new sample starts without the previous tags


def test_deleting_one_samples_tags_keeps_the_other(tmp_path, monkeypatch):
    monkeypatch.setenv("RETICULA_DATA", str(tmp_path))
    pa, pb = AL.default_store("Sample A"), AL.default_store("Sample B")
    o = {"area_px": 50, "length_nm": 5, "width_nm": 2, "aspect": 2.5, "solidity": 0.9, "cx": 10.0, "cy": 20.0}
    AL.save_tags(AL.add_tag([], "img_a.tif", o, 0.274, 50.0, 0), pa)
    AL.save_tags(AL.add_tag([], "img_b.tif", o, 0.274, 50.0, 1), pb)
    pa.unlink()                                                          # what the "Delete the tags of ..." button does
    assert AL.load_tags(pa) == []
    assert [t["label"] for t in AL.load_tags(pb)] == [1]


# ------------------------------------------------------------------ scale-bar label (OCR reads 100 nm as 400 nm)
def test_ocr_run_of_letter_o_after_digits_becomes_zeros():
    assert SC._fix_ocr("1OOnm") == "100nm"                               # the old single pass gave "10Onm": label lost
    assert SC.parse_scale_label(SC._fix_ocr("1OOnm"))[2] == 100.0
    assert SC.parse_scale_label(SC._fix_ocr("2O0 nm"))[2] == 200.0


def test_standard_scale_values_are_1_2_5_times_power_of_ten():
    assert all(SC._is_standard(v) for v in [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000])
    assert not any(SC._is_standard(v) for v in [400, 300, 250, 4, 0])


def test_vote_tie_between_100_and_400_goes_to_standard_value_whatever_the_order():
    v100 = ((100.0, "nm", 100.0), [1, "100nm", []])
    v400 = ((400.0, "nm", 400.0), [1, "400 nm", []])
    for votes in ({v400[0]: v400[1], v100[0]: v100[1]}, {v100[0]: v100[1], v400[0]: v400[1]}):
        label, entry, others = SC._pick_vote(votes)
        assert label == (100.0, "nm", 100.0) and others == [(400.0, "nm", 400.0)]


def test_vote_majority_still_wins_over_standard_value():
    votes = {(400.0, "nm", 400.0): [2, "400 nm", []], (100.0, "nm", 100.0): [1, "100nm", []]}
    assert SC._pick_vote(votes)[0] == (400.0, "nm", 400.0)


# ------------------------------------------------------------------ OCR passes that disagree (synthetic, no private file)
_BAR = {"x0": 993, "x1": 1365, "y0": 980, "y1": 981, "length_px": 372.0, "thickness_px": 2, "polarity": "black", "region_top": 774}


def _fake_ocr(monkeypatch, *, scale1_text="400 nm", scale3_text="100 nm"):
    """Synthetic OCR: only the single-line label passes return text. Pass order (scale, psm, ink) is fixed in read_strip:
    (1,7) first, then (3,7,ink), (3,7) and (3,11); so one reading (400, scale 1) comes first and one (100, scale 3) second."""
    def fake(region, x_off=0, y_off=0, psm=11, scale=2, whitelist=None, ink=None):
        text = None
        if scale == 1 and psm == 7 and whitelist and ink is None:
            text = scale1_text
        elif scale == 3 and psm == 7 and whitelist and ink is not None:
            text = scale3_text
        if text is None:
            return []
        return [{"text": text, "x0": 1130.0, "y0": 985.0, "x1": 1200.0, "y1": 1000.0, "conf": 90.0, "line": (1, 1, 1)}]
    monkeypatch.setattr(SC, "ocr_words", fake)
    monkeypatch.setattr(SC, "ocr_available", lambda: True)


def test_read_strip_tie_between_100_and_400_takes_standard_value(monkeypatch):
    _fake_ocr(monkeypatch)                                              # one vote each: 400 nm first, 100 nm second
    img = np.zeros((1032, 1376), np.uint8)
    r = SC.read_strip(img, strip_row=None, bar=_BAR)
    assert r["label"]["value"] == 100.0 and r["label"]["unit"] == "nm"
    assert r["label_conflict"] == ["100 nm (1×)", "400 nm (1×)"]


def test_read_strip_no_conflict_when_passes_agree(monkeypatch):
    _fake_ocr(monkeypatch, scale1_text="100 nm")                        # both readings say 100 nm
    img = np.zeros((1032, 1376), np.uint8)
    r = SC.read_strip(img, strip_row=None, bar=_BAR)
    assert r["label"]["value"] == 100.0 and "label_conflict" not in r


def test_calibrate_warns_when_ocr_readings_conflict(monkeypatch):
    _fake_ocr(monkeypatch)
    monkeypatch.setattr(SC, "find_info_strip", lambda img8: None)
    monkeypatch.setattr(SC, "detect_bar", lambda img8, strip=None: dict(_BAR))
    img = np.zeros((1032, 1376), np.uint8)
    cal = SC.calibrate(img)                                             # no metadata: the OCR label sets nm/px
    assert cal["label"]["value"] == 100.0
    assert cal["nm_per_px"] == pytest.approx(100.0 / 372.0)
    warn = [t for lvl, t in cal["notes"] if lvl == "warn"]
    assert len(warn) == 1 and "400 nm (1×)" in warn[0] and "100 nm (1×)" in warn[0] and "Used 100 nm" in warn[0]

