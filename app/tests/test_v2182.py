"""v2.18.2: Zetasizer 8.x .dts (Z-average and PDI as stored), FTIR with several spectra per file, TEM instruments in the Methods text,
Reticon full screen / editor on start-up (static checks)."""
import io
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import readers as R, summary as SM  # noqa: E402

# Real 8.02 .dts files are lab data and are not in the repo: point MFS_TEST_DTS_MUA / MFS_TEST_DTS_ZRO2 at them to run the value checks.
DTS_MUA = os.environ.get("MFS_TEST_DTS_MUA")
DTS_ZRO2 = os.environ.get("MFS_TEST_DTS_ZRO2")


@pytest.mark.skipif(not DTS_MUA, reason="set MFS_TEST_DTS_MUA to 040626-ll-MUA.dts")
def test_dts_mua_z_average_and_pdi_as_stored():
    z = R.read_zetasizer(open(DTS_MUA, "rb").read(), "040626-ll-MUA.dts")
    s = z["summary"].sort_values("run")
    assert list(s["z_avg"].round(1)) == [167.0, 166.0, 166.3]
    assert list(s["pdi"].round(3)) == [0.139, 0.133, 0.126]
    assert "8.02" in z["note"]


@pytest.mark.skipif(not DTS_ZRO2, reason="set MFS_TEST_DTS_ZRO2 to ZrO2-HfMOF801_XL.dts")
def test_dts_zro2_z_average_and_pdi_as_stored():
    z = R.read_zetasizer(open(DTS_ZRO2, "rb").read(), "ZrO2-HfMOF801_XL.dts")
    s = z["summary"].sort_values("run")
    assert list(s["z_avg"].round(0)) == [2238.0, 2363.0, 2219.0]
    assert list(s["pdi"].round(3)) == [0.155, 0.527, 0.631]


DTS_ZETA = os.environ.get("MFS_TEST_DTS_ZETA")


@pytest.mark.skipif(not DTS_ZETA, reason="set MFS_TEST_DTS_ZETA to 2026-10-08.dts")
def test_dts_zeta_records_equal_the_zeta_reports():
    z = R.read_zetasizer(open(DTS_ZETA, "rb").read(), "2026-10-08.dts")
    s = z["summary"]
    assert len(s) == 11
    r2 = s[(s["sample"] == "LE101-mPEG") & (s["run"] == 2)].iloc[0]          # Zeta Potential Report of record 2
    assert round(r2.zeta, 1) == 11.7 and round(r2.zeta_sd, 1) == 21.0 and round(r2.cond_mS_cm, 4) == 0.0229
    assert [round(r2.peak1_mV, 2), round(r2.peak2_mV, 1), round(r2.peak3_mV, 1)] == [3.54, 44.8, -38.2]
    assert [round(r2.peak1_area_pct, 1), round(r2.peak2_area_pct, 1), round(r2.peak3_area_pct, 1)] == [78.3, 20.8, 0.8]
    r8 = s[(s["sample"] == "LE98-PEG-cooh") & (s["run"] == 2)].iloc[0]
    assert round(r8.zeta, 1) == -17.9 and round(r8.zeta_sd, 1) == 15.0 and round(r8.cond_mS_cm, 4) == 0.0986
    r6 = s[(s["sample"] == "LE100-PEG-NH2") & (s["run"] == 3)].iloc[0]
    assert round(r6.zeta, 1) == 14.7 and round(r6.zeta_sd, 2) == 7.51 and round(r6.count_kcps, 1) == 168.5


def test_dts_zeta_record_none_for_zero_bytes():
    assert R._dts_zeta_record(b"\x00" * 5000) is None


def test_dts_stored_summary_none_for_zero_bytes():
    assert R._dts_stored_summary(b"\x00" * 4096) is None


def _csv_multi(n=3, pts=50):
    x = np.linspace(400, 2000, pts)
    cols = {"Wavenumber (cm-1)": x}
    for i in range(n):
        cols[f"Sample {i + 1}"] = np.sin(x / (100 + i)) + i
    return pd.DataFrame(cols).to_csv(index=False).encode("utf-8")


def test_ftir_multi_csv_gives_one_spectrum_per_column():
    out = R.read_ftir_multi(io.BytesIO(_csv_multi(3)), "many.csv")
    assert [lab for lab, _ in out] == ["Sample 1", "Sample 2", "Sample 3"]
    for _, d in out:
        assert len(d) == 50 and list(d.columns) == ["x", "y"]


def test_ftir_single_csv_still_one_spectrum():
    csv = b"wavenumber,absorbance\n400,0.1\n500,0.2\n600,0.3\n"
    out = R.read_ftir_multi(io.BytesIO(csv), "one.csv")
    assert len(out) == 1 and len(out[0][1]) == 3


def _fake_spa(x0=4000.0, x1=100.0, n=50):
    import struct
    b = bytearray(2000)
    b[0:18] = b"Spectral Data File"
    struct.pack_into("<I", b, 304 + 2, 560); struct.pack_into("<I", b, 304 + 6, 140); b[304] = 2
    struct.pack_into("<I", b, 320 + 2, 1184); struct.pack_into("<I", b, 320 + 6, 4 * n); b[320] = 3
    struct.pack_into("<i", b, 564, n); struct.pack_into("<ff", b, 576, x0, x1)
    y = np.arange(n, dtype="<f4")
    b[1184:1184 + 4 * n] = y.tobytes()
    return bytes(b), y


def test_spa_is_read_ascending_with_the_right_values():
    raw, y = _fake_spa()
    d = R.read_ftir(raw, "a.SPA")
    assert d.x.is_monotonic_increasing and abs(d.x.iloc[0] - 100) < 1e-3 and abs(d.x.iloc[-1] - 4000) < 1e-3
    assert np.allclose(d.y.values, y[::-1])
    assert len(R.read_ftir_multi(raw, "a.SPA")) == 1


def test_not_an_spa_is_refused():
    with pytest.raises(ValueError):
        R.read_spa(b"nothing")


def test_headerless_csv_is_one_spectrum_with_all_points():
    csv = b"4.0e+002,1.0e+001\n4.1e+002,2.0e+001\n4.2e+002,3.0e+001\n"
    out = R.read_ftir_multi(io.BytesIO(csv), "x.CSV")
    assert len(out) == 1 and out[0][0] is None and len(out[0][1]) == 3


def test_methods_text_names_the_jeol_instruments():
    src = open(os.path.join(HERE, "mfs", "summary.py"), encoding="utf-8").read()
    assert "JEOL JEM-1010 at 100 kV" in src and "JEOL JEM-2200FS field-emission microscope at 200 kV" in src


def test_simple_mode_bars_scatter_and_legend():
    import matplotlib.pyplot as plt
    from mfs import figctl as F
    fig, ax = plt.subplots()
    c = np.arange(1, 4)
    ax.bar(c - 0.2, [5, 20, 40], width=0.4, color="C0", label="A")
    ax.bar(c + 0.2, [2, 10, 30], width=0.4, color="C2", label="B")
    ax.scatter([1, 2], [1, 2], color="C0", label="g1")
    ax.scatter([2, 3], [2, 3], color="C3", label="g2")
    ax.legend()
    cfg = F.default_settings(); cfg["simple_black"] = True
    F.apply(fig, cfg)
    hatches = [c_.patches[0].get_hatch() for c_ in ax.containers]
    assert hatches[0] != hatches[1]                                      # grouped bars: one hatch per series
    sc = [col for col in ax.collections]
    assert np.allclose(sc[0].get_facecolor()[0][:3], 0) and np.allclose(sc[1].get_facecolor()[0][:3], 1)   # filled / open
    lg = ax.get_legend()
    by_text = {t.get_text(): h for h, t in zip(lg.legend_handles, lg.get_texts())}
    assert np.allclose(np.asarray(by_text["A"].get_facecolor()).reshape(-1)[:3], [0.86, 0.86, 0.86])   # the existing legend follows the simple look
    assert by_text["A"].get_hatch() != by_text["B"].get_hatch()
    plt.close(fig)


def _si_like_pattern(x, a, prof, bg=100.0):
    """Synthetic line pattern: Si-type reflections (Fd-3m) with the module's own pseudo-Voigt width model (test data only)."""
    lam = 1.5406
    hkl = [(1, 1, 1, 100), (2, 2, 0, 55), (3, 1, 1, 30), (4, 0, 0, 12), (3, 3, 1, 20), (4, 2, 2, 10), (5, 1, 1, 25), (4, 4, 0, 8)]
    y = np.full_like(x, bg)
    for h, k, l, I in hkl:
        tt = 2 * np.degrees(np.arcsin(lam / (2 * a / np.sqrt(h * h + k * k + l * l))))
        th = np.radians(tt / 2); t = np.tan(th)
        from mfs import rietveld as RV
        H, eta = RV._tch(prof["U"] * t * t + prof["V"] * t + prof["W"], prof["X"] / np.cos(th) + prof["Y"] * t)
        dx = x - tt
        y += I * (eta / (1 + 4 * dx * dx / H ** 2) + (1 - eta) * np.exp(-4 * np.log(2) * dx * dx / H ** 2))
    return y


def test_rietveld_standard_profile_is_fixed_and_si_cell_is_the_certificate():
    from mfs import rietveld as RV
    truth = {"U": 0.002, "V": -0.001, "W": 0.006, "X": 0.03, "Y": 0.01}
    x = np.linspace(20, 90, 1500)
    rng_ = np.random.default_rng(3)
    std = RV.instrument_profile_from_standard(x, _si_like_pattern(x, 5.4312, truth) + rng_.normal(0, 0.5, x.size))
    assert abs(std["cell_a"] - 5.4312) < 1e-6                             # the standard's lattice parameter is fixed to the certificate
    sample = _si_like_pattern(x, 5.658, truth) + rng_.normal(0, 0.5, x.size)   # a Ge-type sample measured on the same instrument
    res = RV.quick_refine(x, sample, [RV.Phase.from_cell((5.62, 5.62, 5.62, 90, 90, 90), "Ge")], mode="lebail",
                          instrument={k: std[k] for k in ("U", "V", "W", "X", "Y")} | {"fixed": True})
    assert not any(res["global"][k]["refined"] for k in ("U", "V", "W", "X", "Y"))   # the instrument profile stays fixed
    assert res["stats"]["Rwp"] < 0.03
    assert abs(res["phases"][0]["cell"][0] - 5.658) < 0.01


def test_rietveld_glass_holder_background_automatic_and_measured():
    from mfs import rietveld as RV
    prof = {"U": 0.002, "V": -0.001, "W": 0.006, "X": 0.03, "Y": 0.01}
    x = np.linspace(10, 90, 1600)
    rng_ = np.random.default_rng(5)
    holder = 900 * np.exp(-0.5 * ((x - 23.5) / 5.5) ** 2) + 250 * np.exp(-0.5 * ((x - 37) / 9) ** 2) + 60
    y = _si_like_pattern(x, 5.658, prof, bg=40) * 3 + holder + rng_.normal(0, 1.5, x.size)
    mk = lambda: [RV.Phase.from_cell((5.62, 5.62, 5.62, 90, 90, 90), "Ge")]
    off = RV.quick_refine(x, y, mk(), mode="lebail", glass_bg=False)
    auto = RV.quick_refine(x, y, mk(), mode="lebail")
    meas = RV.quick_refine(x, y, mk(), mode="lebail", glass_bg=False, bg_pattern=(x, holder + rng_.normal(0, 1.5, x.size)))
    assert auto["stats"]["Rwp"] < 0.5 * off["stats"]["Rwp"] and abs(auto["phases"][0]["cell"][0] - 5.658) < 0.01
    assert meas["bg_measured"] and abs(meas["phases"][0]["cell"][0] - 5.658) < 0.01


def test_help_documents_cover_every_module_and_render_as_pdf():
    from mfs import helpdocs as HD, helppdf, ui
    skip = ("Home", "Help and info")
    assert all(p in HD.DOCS for p in ui.PAGES if p not in skip), [p for p in ui.PAGES if p not in skip and p not in HD.DOCS]
    for k, d in HD.DOCS.items():
        assert d["title"] and any(kind == "h" for kind, _ in d["sections"]), k
        assert d["sections"][1][0] == "p" or d["sections"][0][0] == "h", k
        txt = " ".join(t for _, t in d["sections"])
        assert "Terms" in txt, f"{k}: no terms section (abbreviations must be written out first)"
    pdf = helppdf.render_pdf(HD.DOCS["Photoluminescence"])
    assert pdf[:5] == b"%PDF-" and len(pdf) > 20000
    assert HD.DOCS["Photoluminescence"]["title"] == "Info implementatie ratiometric thermometry in Reticula"


def test_help_page_and_info_button_render():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300); at.run()
    at.sidebar.radio[0].set_value("Help and info").run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert any("Download this document" in getattr(d, "label", "") for d in at.get("download_button"))
    at.sidebar.radio[0].set_value("PXRD").run()
    assert not at.exception
    labels = [getattr(b, "label", "") for b in at.get("download_button")]
    assert any("Info about this module" in l for l in labels), labels[:6]


def test_reticon_editor_has_no_opt_in_gate_and_full_screen_is_prominent():
    page = open(os.path.join(HERE, "mfs", "pages", "schemes.py"), encoding="utf-8").read()
    assert "opt_in(\"Open the scheme editor" not in page
    html = open(os.path.join(HERE, "mfs", "components", "scheme_canvas", "index.html"), encoding="utf-8").read()
    assert re.search(r'<div id="bar">\s*<button id="fsbtn"', html)
    assert "body.max{--top:58px}" in html
