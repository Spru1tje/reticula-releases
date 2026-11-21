"""Batch C: particle sizing, EDX overlay/line scan, ICP/XRF ratios with propagated errors, SEM batch options."""
import math
import os
import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _synthetic(n=150, seed=3, mean=55.0, s=0.25, pix=1.0, size=1400):
    rng = np.random.default_rng(seed)
    img = np.full((size, size), 40.0); yy, xx = np.mgrid[:size, :size]; d_true = []
    for _ in range(n * 4):
        if len(d_true) >= n:
            break
        d = rng.lognormal(math.log(mean), s); r = d / 2 / pix
        cx, cy = rng.uniform(r + 2, size - r - 2, 2)
        sub = (slice(int(cy - r - 2), int(cy + r + 3)), slice(int(cx - r - 2), int(cx + r + 3)))
        m = (xx[sub] - cx) ** 2 + (yy[sub] - cy) ** 2 <= r * r
        if (img[sub][m] > 100).any() or (img[sub] > 100).sum() > 0:
            continue
        img[sub][m] = 200; d_true.append(d)
    img += rng.normal(0, 8, img.shape)
    return img, np.array(d_true)


def test_particles_synthetic_lognormal():
    from mfs import particles as PA
    img, d = _synthetic()
    lab, f, t = PA.segment(img, dark_particles=False, split=True, max_px=1400)
    rows = PA.measure(lab, 1.0, f)
    ecd = np.array([r["ECD (nm)"] for r in rows])
    assert abs(len(ecd) - len(d)) <= 0.1 * len(d)
    assert abs(ecd.mean() - d.mean()) / d.mean() < 0.06
    fit = PA.lognormal_fit(ecd)
    assert abs(fit["geo_sd"] - math.exp(0.25)) < 0.08 and fit["n"] == len(ecd)


def test_particles_split_touching_pair():
    from mfs import particles as PA
    img = np.zeros((200, 300)); yy, xx = np.mgrid[:200, :300]
    for cx in (110, 178):
        img[(xx - cx) ** 2 + (yy - 100) ** 2 <= 36 ** 2] = 1
    lab, f, _ = PA.segment(img, dark_particles=False, split=True, smooth=0, open_px=0, max_px=400)
    assert lab.max() == 2
    lab1, _, _ = PA.segment(img, dark_particles=False, split=False, smooth=0, open_px=0, max_px=400)
    assert lab1.max() == 1


def test_edx_overlay_and_linescan():
    from mfs import edx as E
    a = np.zeros((50, 80)); a[:, :40] = 10; b = np.zeros((50, 80)); b[:, 40:] = 10
    rgb = E.overlay({"Zr": a, "S": b}, {"Zr": "#FF0000", "S": "#0000FF"})
    assert rgb.shape == (50, 80, 3) and rgb[25, 10, 0] > 0.9 and rgb[25, 70, 2] > 0.9 and rgb[25, 10, 2] < 0.05
    d, p = E.line_scan(a, (0, 25), (79, 25), width=3)
    assert p[5] > 9 and p[-5] < 1 and abs(d[-1] - 79) < 1e-6


def test_ratios_propagation():
    from mfs import ratios as RT
    r = RT.zr_hf_s({"Zr": 20.0, "Hf": 10.0, "S": 2.0}, {"Zr": 0.4, "Hf": 0.3, "S": 0.1})
    nz, nh, ns = 20 / 91.224, 10 / 178.49, 2 / 32.06
    assert abs(r["Hf/(Zr+Hf)"][0] - nh / (nz + nh)) < 1e-12
    assert abs(r["S per M6 node"][0] - 6 * ns / (nz + nh)) < 1e-12
    # Monte-Carlo check of σ(S per node)
    rng = np.random.default_rng(0); N = 200000
    z = rng.normal(20, 0.4, N) / 91.224; h = rng.normal(10, 0.3, N) / 178.49; s = rng.normal(2, 0.1, N) / 32.06
    assert abs(r["S per M6 node"][1] - np.std(6 * s / (z + h))) / r["S per M6 node"][1] < 0.03
    assert abs(r["Hf/(Zr+Hf)"][1] - np.std(h / (z + h))) / r["Hf/(Zr+Hf)"][1] < 0.03
    zr, hf, su = r["Zr:Hf:S (per M6)"]; assert abs(zr + hf - 6) < 1e-9


# ---------------- app ----------------
def _goto(at, page):
    at.button(key="nv_" + page).click().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    return at


@pytest.fixture
def at(tmp_path, monkeypatch):
    monkeypatch.setenv("MFS_STORE", str(tmp_path / "store.json"))
    from streamlit.testing.v1 import AppTest
    a = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=400); a.run()
    return a


def test_app_sem_particles_edx(at):
    _goto(at, "SEM / TEM / EDX")
    at.radio(key="ps_dark").set_value("bright").run(); assert not at.exception   # v2.9 default is dark; example = bright
    assert at.multiselect(key="ps_imgs").value
    assert any("± " in m.value and "n = " in m.value for m in at.markdown), "no size summary"
    secs = {s for v in at.session_state["sample_store"].values() for s in v["sections"]}
    assert "sem" in secs
    at.selectbox(key="sem_match").set_value("1 µm").run(); assert not at.exception
    at.radio(key="sem_onecal").set_value("apply the first image's calibration to all").run(); assert not at.exception
    assert any("batch" in s.value for s in at.success)
    at.radio(key="edx_mode").set_value("both").run(); assert not at.exception
    at.checkbox(key="edx_ls").check().run(); assert not at.exception, [e.value for e in at.exception][:2]
    at.radio(key="ps_meth").set_value("click-to-measure (two clicks per particle)").run()
    assert not at.exception, [e.value for e in at.exception][:2]


def test_app_icp_xrf_ratios_to_defects(at):
    _goto(at, "ICP-OES / XRF")
    assert any("S per M" in str(h.value) or "ratios" in str(h.value) for h in at.subheader)
    at.radio(key="icp_kind").set_value("XRF / any table (sample, element, value, sd)").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    store = at.session_state["sample_store"]
    x = [v for v in store.values() if "xrf" in v["sections"] and v["sections"]["xrf"]["values"].get("S_per_M6_sd") is not None]
    assert x, "no S per M6 ± σ stored"
    _goto(at, "Defect calculator")
    at.selectbox(key="def_pull").set_value([k for k, v in store.items() if v in x][0]).run()
    assert not at.exception
    assert any("S/M = " in s.value and "±" in s.value for s in at.success), ([s.value for s in at.success], x[0])


def test_molar_ratios_generic_matches_zr_hf_s_and_works_for_other_metals():
    from mfs import ratios as RT
    wt, sd = {"Zr": 20.0, "Hf": 10.0, "S": 2.0}, {"Zr": 0.4, "Hf": 0.3, "S": 0.1}
    old = RT.zr_hf_s(wt, sd); new = RT.molar_ratios(wt, sd, ["Zr", "Hf"], per=6)
    assert abs(new["S per 6 (Zr+Hf)"][0] - old["S per M6 node"][0]) < 1e-12 and abs(new["S per 6 (Zr+Hf)"][1] - old["S per M6 node"][1]) < 1e-12
    assert abs(new["Hf/(Zr+Hf)"][0] - old["Hf/(Zr+Hf)"][0]) < 1e-12 and abs(new["Hf/(Zr+Hf)"][1] - old["Hf/(Zr+Hf)"][1]) < 1e-12
    cu = RT.molar_ratios({"Cu": 30.0, "Zn": 3.0}, {"Cu": 0.5, "Zn": 0.1}, ["Cu"])      # any other metal pair
    assert abs(cu["Zn/(Cu)"][0] - (3.0 / 65.38) / (30.0 / 63.546)) < 1e-3 and cu["Zn/(Cu)"][1] > 0
    assert RT.molar_ratios({"Cu": 1.0}, {}, ["Fe"]) == {}
