"""v2.18.4: Edinburgh .FS reader and the standard background subtraction of luminescence spectra (mfs.plbackground)."""
import os
import struct
import sys
import warnings

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import pl as PL, plbackground as BG, readers as R  # noqa: E402


# ------------------------------------------------------------------------------------------------------------------------------ a synthetic .FS file (same layout as the real ones)
def _block(name, x0, x1, y, kind=b"EAI Emission Scan"):
    ax = kind + b"\x00\r\n" + struct.pack("<iii", 150, 0, 0) + struct.pack("<ff", x0, x1) + b"Wavelength\x00nm\x00\x00" + struct.pack("<ff", 12000.0, 70000.0) + b"Counts\x00" + b"2.16.2.1\x00" + b"e\x00\x00\x00\x01\x00\x00\x00"
    return (ax + name.encode() + b"\x00" + b"\x00" * 20 + b"CNI \x00" + b"\x00" * 14 + b"NIR PMT\x00" + b"\x00" * 10 + b"Reference Det\x00" + b"\x00" * 6 + b"External Laser (MCS)\x00" + b"\x00\x00"
            + np.asarray(y, "<f4").tobytes())


def make_fs(scans, multi=False):
    """scans = [(name, x0, x1, y)]."""
    if multi:
        head = b"EAI Multiple Scans\x00\r\n" + struct.pack("<iii", 1000, 0, 0) + struct.pack("<ff", scans[0][1], scans[0][2]) + b"Wavelength\x00nm\x00\x00" + b"Counts\x00" + b"2.16.2.1\x00" + b"\x00" * 8
        return head + b"".join(_block(n, a, b, y) for n, a, b, y in scans)
    n, a, b, y = scans[0]
    head = b"EAI Emission Scan\x00\r\n" + struct.pack("<iii", 1000, 0, 0) + struct.pack("<ff", a, b) + b"Wavelength\x00nm\x00\x00" + b"Counts\x00" + b"2.16.2.1\x00" + b"\x00" * 8
    return head + _block(n, a, b, y)


def _series(n_spec=7, pedestal=12500.0, slope=-0.6, seed=0):
    rng = np.random.default_rng(seed)
    x = np.linspace(900.0, 1250.0, 351)
    g = lambda c, w, a: a * np.exp(-0.5 * ((x - c) / w) ** 2)
    ys, T = [], []
    for k in range(n_spec):
        T.append(20 + 5 * k)
        em = g(995, 20, 50000 - 2500 * k) + g(960, 8, 22000 - 1500 * k) + g(1020, 14, 45000 - 2200 * k) + g(1190, 5, 2800) + g(1150, 6, 700)
        ys.append(pedestal + slope * (x - 900) + 40 * k + em + rng.normal(0, 90, len(x)))
    return x, ys, T


def test_single_fs_roundtrip():
    x, ys, _ = _series(1)
    sc = R.read_edinburgh_fs(make_fs([("sample_20C", 900.0, 1250.0, ys[0])]), "a.FS")
    assert len(sc) == 1 and sc[0]["name"] == "sample_20C" and sc[0]["n"] == 351
    assert np.allclose(sc[0]["x"], x) and np.allclose(sc[0]["y"], ys[0], rtol=1e-6)
    assert sc[0]["detector"] == "NIR PMT" and "Laser" in sc[0]["excitation"] and sc[0]["step"] == pytest.approx(1.0)


def test_multi_fs_gives_every_scan_and_short_names():
    x, ys, T = _series(4)
    names = [f"Sample_{t}C(Multi Scans)(Multi Scans)" if k < 3 else f"Sample_{t}C" for k, t in enumerate(T)]
    raw = make_fs([(n, 900.0, 1250.0, y) for n, y in zip(names, ys)], multi=True)
    sp = PL.read_pl(raw, "multi.FS")
    assert [s["name"] for s in sp] == [f"{t}C" for t in T]
    assert [PL.parse_temperature(s["T_text"])[0] for s in sp] == T
    assert all(np.allclose(s["I"], y, rtol=1e-6) for s, y in zip(sp, ys))


def test_not_an_fs_is_refused():
    with pytest.raises(ValueError):
        R.read_edinburgh_fs(b"nothing", "x.FS")


# ------------------------------------------------------------------------------------------------------------------------------ background
def test_automatic_plan_finds_the_pedestal_in_the_signal_free_windows():
    x, ys, _ = _series(7)
    plan = BG.autoselect(x, ys)
    assert plan["method"] in ("offset", "poly") and (plan["order"] or 0) <= 2 and plan["kind"] == "NIR"
    wins = plan["windows"]
    assert wins and all(not (940 < (a + b) / 2 < 1060) for a, b in wins)            # none in the emission band
    res = BG.apply_plan([{"wl": x, "I": y} for y in ys], plan)
    for r in res:
        free = BG.mask_from_windows(x, wins)
        assert abs(r["y"][free].mean()) < 0.25 * plan["sigma"]                         # pedestal gone in the free windows
        assert r["y"].max() > 30000                                                     # the emission is not touched
    # the corrected emission equals the true emission within the noise
    g = lambda c, w, a: a * np.exp(-0.5 * ((x - c) / w) ** 2)
    true0 = g(995, 20, 50000) + g(960, 8, 22000) + g(1020, 14, 45000) + g(1190, 5, 2800) + g(1150, 6, 700)
    i0 = np.argmin(abs(x - 995))
    assert abs(res[0]["y"][i0] - true0[i0]) < 4 * plan["sigma"] + 1500


def test_constant_pedestal_gets_order_zero():
    x, ys, _ = _series(5, slope=0.0)
    plan = BG.autoselect(x, ys)
    assert plan["method"] == "offset" and plan["order"] == 0


def test_manual_windows_and_dark_file():
    x, ys, _ = _series(3)
    plan = BG.autoselect(x, ys, windows=[(1100, 1140), (1205, 1250)])
    assert plan["windows"] == [(1100, 1140), (1205, 1250)]
    dark = (x, 12500.0 - 0.6 * (x - 900) + np.random.default_rng(1).normal(0, 20, len(x)))
    pd_ = BG.autoselect(x, ys, dark=dark)
    assert pd_["method"] == "dark"
    r = BG.apply_plan([{"wl": x, "I": ys[0]}], pd_, dark)[0]
    assert abs(r["y"][(x > 1100) & (x < 1140)].mean()) < 60


def test_no_free_window_falls_back_to_a_flexible_baseline():
    x = np.linspace(400, 700, 301)
    y = 1000 + 3 * (x - 400) + 5000 * np.exp(-0.5 * ((x - 520) / 60) ** 2) + 4000 * np.exp(-0.5 * ((x - 640) / 50) ** 2) + np.random.default_rng(3).normal(0, 20, len(x))
    plan = BG.autoselect(x, [y])
    assert plan["method"] == "rubber" and not plan["windows"]
    r = BG.apply_plan([{"wl": x, "I": y}], plan)[0]
    assert r["y"].min() > -150 and r["y"].max() > 3500                                  # the floor goes, the broad bands stay; the rubber band never rises above the data
    assert {c["method"] for c in plan["candidates"]} >= {"rubber", "arpls", "snip", "modpoly"}


def test_despike_removes_a_single_point_spike_only():
    x, ys, _ = _series(1)
    y = ys[0].copy(); y[100] += 9000
    out, n = BG.despike(y)
    assert n == 1 and abs(out[100] - ys[0][100]) < 600 and np.allclose(np.delete(out, 100), np.delete(y, 100))


def test_ratio_reports_the_spread_over_background_alternatives():
    """The ratio of two bands on a pedestal + slope: right within a few %, and the background choice appears as its own part of the systematic error."""
    from mfs import plratio as PR
    x = np.linspace(900.0, 1250.0, 351)
    g = lambda c, w, a: a * np.exp(-0.5 * ((x - c) / w) ** 2)
    rng = np.random.default_rng(5)
    T = [293.15 + 10 * k for k in range(6)]
    true = [0.5 + 0.06 * k for k in range(6)]
    ys = [12500.0 - 0.6 * (x - 900) + 30 * k + g(980, 10, r * 20000) + g(1030, 10, 20000) + rng.normal(0, 100, len(x)) for k, r in enumerate(true)]
    plan = BG.autoselect(x, ys)
    assert plan["method"] in ("offset", "poly")
    alts = BG.alternatives(x, plan)
    assert alts and all(isinstance(lab, str) and p["method"] in ("offset", "poly") for lab, p in alts)
    mk = lambda out: [{"T": t, "x": x, "y": o["y"], "var": y} for t, o, y in zip(T, out, ys)]
    main = mk(BG.apply_plan([{"wl": x, "I": y} for y in ys], plan))
    alt_sets = [(lab, mk(BG.apply_plan([{"wl": x, "I": y} for y in ys], p))) for lab, p in alts]
    fit = PR.auto_fit(main, True, 900.0, 1250.0, 1, 0, bands=[(980.0, 23.5), (1030.0, 23.5)])
    assert fit is not None
    fit = PR.robust_ratio(main, fit["bands"], bg_order=1, counts=True, main=fit, bg_alternatives=alt_sets)
    D = np.array([r["ratio"] for r in fit["rows"]])
    assert np.allclose(D, true, rtol=0.04)
    assert fit["bg_alt_labels"] and np.all(np.isfinite(fit["ratio_sys_bg"])) and np.all(fit["ratio_sys_bg"] >= 0)
    assert np.all(fit["ratio_total"] >= fit["ratio_sys"] - 1e-12)
    assert np.allclose(fit["ratio_sys"] ** 2, fit["ratio_sys_fit"] ** 2 + fit["ratio_sys_bg"] ** 2)


def test_range_kind():
    assert BG.range_kind(np.linspace(900, 1250, 5)) == "NIR" and BG.range_kind(np.linspace(450, 700, 5)) == "UV/VIS" and BG.range_kind(np.linspace(500, 1100, 5)) == "VIS-NIR"


# ------------------------------------------------------------------------------------------------------------------------------ the real files (private lab data): set MFS_TEST_FS to a folder
FS_DIR = os.environ.get("MFS_TEST_FS")


@pytest.mark.skipif(not FS_DIR, reason="set MFS_TEST_FS to the folder with the Edinburgh .FS test files")
def test_real_edinburgh_files():
    import glob
    for fn in glob.glob(os.path.join(FS_DIR, "*.FS")):
        sc = R.read_edinburgh_fs(open(fn, "rb").read(), fn)
        assert all(s["n"] == 351 and s["x_start"] == 900.0 and s["x_end"] == 1250.0 for s in sc)
    multi = [f for f in glob.glob(os.path.join(FS_DIR, "*20-50C*Multi*.FS"))][0]
    sp = PL.read_pl(open(multi, "rb").read(), multi)
    assert [round(PL.parse_temperature(s["T_text"])[0]) for s in sp] == [20, 25, 30, 35, 40, 45, 50]
    x = sp[0]["wl"]
    plan = BG.autoselect(x, [s["I"] for s in sp])
    assert plan["method"] in ("offset", "poly") and plan["kind"] == "NIR" and len(plan["windows"]) >= 2


def test_page_reads_fs_and_shows_the_background_block(monkeypatch):
    """The PL page with an .FS series (the file reader is patched in, because AppTest cannot upload files)."""
    from streamlit.testing.v1 import AppTest
    warnings.filterwarnings("ignore")
    x, ys, T = _series(7)
    names = [f"S_{t}C" for t in T]
    raw = make_fs([(n, 900.0, 1250.0, y) for n, y in zip(names, ys)], multi=True)
    orig = PL.read_pl
    monkeypatch.setattr(PL, "read_pl", lambda src, name="": orig(raw, "multi.FS"))
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300)
    at.session_state["_lv"] = {"global": "Standard"}
    at.run()
    at.sidebar.radio[0].set_value("Photoluminescence").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert any("Background subtraction" in e.label for e in at.expander), [e.label for e in at.expander]
    assert any(m.label == "Signal-free windows" for m in at.metric)
    assert any("Thermometry" in t.label for t in at.tabs)
