"""v2.21: learned agglomerate model. Synthetic objects only: long thin rods (label 0) and compact clumps (label 1).
Checks the model separates them, the cross-validation is honest (no leakage), tags are stored once per object and the
rule is used when there are too few tags."""
import numpy as np

from mfs import agglo_learn as AL
from mfs import rods as R

NM = 0.274


def _rod(rng, cx, cy):
    L = rng.uniform(12, 30) / NM / 4.0 * 4.0
    w = rng.uniform(1.8, 3.5)
    area = int(max(12, L * w * 0.8))
    return {"id": 0, "area_px": area, "length_nm": L * NM, "width_nm": w, "aspect": L * NM / w, "solidity": rng.uniform(0.85, 1.0),
            "cx": cx, "cy": cy, "film_edge": False}


def _clump(rng, cx, cy):
    d = rng.uniform(12, 28)
    area = int(np.pi * (d / 2) ** 2 / NM ** 2 * rng.uniform(0.9, 1.1))
    return {"id": 0, "area_px": area, "length_nm": d * rng.uniform(1.0, 1.4), "width_nm": d * rng.uniform(0.7, 0.95),
            "aspect": rng.uniform(1.0, 1.6), "solidity": rng.uniform(0.6, 0.85), "cx": cx, "cy": cy, "film_edge": False}


def _set(seed, n=60):
    rng = np.random.default_rng(seed)
    tags = []
    for k in range(n):
        tags = AL.add_tag(tags, "img_a", _rod(rng, 10.0 * k, 0.0), NM, 80.0, 0)
        tags = AL.add_tag(tags, "img_b", _clump(rng, 10.0 * k, 5.0), NM, 80.0, 1)
    return tags


def test_model_separates_rods_from_clumps():
    tags = _set(0)
    cv = AL.cross_validate(tags)
    assert cv is not None and cv["accuracy"] >= 0.9, cv
    assert cv["balanced_accuracy"] >= 0.9, cv
    m = AL.fit(tags)
    assert m is not None and m["n"] == len(tags)
    rod = AL.predict_proba(m, tags[0]["features"]); clump = AL.predict_proba(m, tags[1]["features"])
    assert rod < 0.5 < clump


def test_too_few_tags_gives_no_model_and_rule_is_used():
    tags = _set(1, n=5)
    assert AL.fit(tags) is None and AL.cross_validate(tags) is None
    objs = [{"id": 1, "area_px": 9000, "length_nm": 40, "width_nm": 30, "aspect": 1.3, "solidity": 0.7, "film_edge": False,
             "cx": 0, "cy": 0}]
    summ = AL.apply(objs, None, NM, fallback=R.flag_agglomerates)
    assert summ["method"] == "rule"


def test_one_tag_per_object_and_replace():
    tags = []
    o = {"area_px": 50, "length_nm": 5, "width_nm": 2, "aspect": 2.5, "solidity": 0.9, "cx": 10.0, "cy": 20.0}
    tags = AL.add_tag(tags, "x", o, NM, 50.0, 0)
    tags = AL.add_tag(tags, "x", dict(o, cx=11.0), NM, 50.0, 1)        # same object (within 3 px): replaces the tag
    assert len(tags) == 1 and tags[0]["label"] == 1
    tags = AL.add_tag(tags, "y", o, NM, 50.0, 0)                        # other image: a separate tag
    assert len(tags) == 2
    tags = AL.remove_tag(tags, "x", o)
    assert len(tags) == 1 and tags[0]["image"] == "y"


def test_store_roundtrip(tmp_path):
    tags = _set(2, n=3)
    p = AL.save_tags(tags, tmp_path / "t.json")
    back = AL.load_tags(p)
    assert len(back) == len(tags) and back[0]["features"].keys() == tags[0]["features"].keys()


def test_apply_uses_model_when_available():
    tags = _set(3)
    m = AL.fit(tags)
    rng = np.random.default_rng(9)
    objs = [_rod(rng, 0, 0), _clump(rng, 5, 5)]
    summ = AL.apply(objs, m, NM, fallback=R.flag_agglomerates)
    assert summ["method"] == "learned"
    assert objs[0]["agglomerate"] is False and objs[1]["agglomerate"] is True
