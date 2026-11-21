"""v2.9 UV/Vis: readers, conversions, calibration, kinetics, Tauc, page."""
import os, io, glob, zipfile
import numpy as np, pandas as pd, pytest
from mfs import uvvis as U

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXD = os.path.join(HERE, "example_data", "UVVis")
REAL = sorted(glob.glob(os.path.join(HERE, "tests", "data", "uvvis_real*")))


def _rd(p):
    return U.read_any(open(p, "rb").read(), os.path.basename(p))


def test_nanophotometer_excel_example():
    f = _rd(os.path.join(EXD, "NanoPhotometer_BCA_EXAMPLE.xlsx"))
    assert f.instrument == "Implen NanoPhotometer"
    assert len(f.spectra) == 11 and f.spectra[0].x[0] == 200 and f.spectra[0].x[-1] == 900
    assert f.results is not None and "A562" in f.results.columns
    # the spectrum value at 562 nm equals the A562 of the results sheet
    for s in f.spectra:
        r = f.results[f.results["Sample name"] == s.name].iloc[0]
        assert abs(U.value_at(s.x, s.y, 562) - float(r["A562"])) < 6e-5


def test_csv_variants_and_decimal_comma(tmp_path):
    x = np.arange(250, 401, 1.0); y = np.exp(-0.5 * ((x - 320) / 15) ** 2)
    txt = "Wavelength (nm);Abs\n" + "\n".join(f"{a:.0f};{b:.5f}".replace(".", ",") for a, b in zip(x, y))
    f = U.read_any(txt.encode(), "comma.csv")
    assert len(f.spectra) == 1 and abs(U.lambda_max(f.spectra[0].x, f.spectra[0].y)[0] - 320) < 1
    # Cary-style paired columns
    df = pd.DataFrame({"A_nm": x, "A_abs": y, "B_nm": x, "B_abs": 2 * y})
    hdr = "S1,,S2,\nWavelength (nm),Abs,Wavelength (nm),Abs\n"
    f = U.read_any((hdr + df.to_csv(index=False, header=False)).encode(), "cary.csv")
    assert len(f.spectra) == 2 and abs(f.spectra[1].y.max() - 2) < 1e-3
    # row-wise (wavelengths across)
    row = "Sample," + ",".join(f"{v:.0f}" for v in x) + "\nX1," + ",".join(f"{v:.4f}" for v in y) + "\n"
    f = U.read_any(row.encode(), "wide.csv")
    assert len(f.spectra) == 1 and f.spectra[0].name == "X1"


def test_ids_encrypted_and_containers():
    rnd = np.random.default_rng(0).integers(0, 256, 20000, dtype=np.uint8).tobytes()
    with pytest.raises(U.IDSEncrypted, match="Export"):
        U.read_any(rnd, "run.ids")
    x = np.arange(200, 300.0); buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("s1.csv", "nm,A\n" + "\n".join(f"{a},{a/1000}" for a in x))
    f = U.read_any(buf.getvalue(), "run.ids")
    assert len(f.spectra) == 1


def test_conversions_and_km():
    s = U.Spectrum("t", np.arange(200, 300.0), np.full(100, 10.0), "T")
    assert np.allclose(U.to_absorbance(s), 1.0)
    s2 = U.Spectrum("a", s.x, np.full(100, 2.0), "A")
    assert np.allclose(U.to_percent_T(s2), 1.0)
    assert np.isclose(U.kubelka_munk([50.0])[0], 0.25)


def test_calibration_lod_loq_and_unknown():
    c = np.array([0, 10, 20, 40, 60, 80.0]); A = 0.01 + 0.0123 * c + np.array([0.001, -0.002, 0.0015, -0.001, 0.002, -0.0005])
    cal = U.calibration(c, A)
    assert abs(cal["slope"] - 0.0123) < 2e-4 and cal["R2"] > 0.999
    assert np.isclose(cal["LOD"], 3.3 * cal["s_yx"] / cal["slope"]) and np.isclose(cal["LOQ"] / cal["LOD"], 10 / 3.3)
    cu, se = U.predict_conc(cal, [0.5])
    assert abs(cu[0] - (0.5 - cal["intercept"]) / cal["slope"]) < 1e-9 and se[0] > 0
    q = U.calibration(c, 0.02 + 1.25e-2 * c - 4e-5 * c ** 2, "quadratic")
    cq, _ = U.predict_conc(q, [0.02 + 1.25e-2 * 30 - 4e-5 * 900])
    assert abs(cq[0] - 30) < 1e-6
    o = U.calibration(c, 0.0123 * c, "through origin")
    assert o["intercept"] == 0 and abs(U.epsilon(o, 1.0, 1e-6) - 12300) < 1


def test_bca_example_recovers_unknowns():
    f = _rd(os.path.join(EXD, "NanoPhotometer_BCA_EXAMPLE.xlsx"))
    t = f.results
    std = t[t["Standard conc. [µg/mL]"].notna() & (t["Standard conc. [µg/mL]"] != "")]
    cal = U.calibration(std["Standard conc. [µg/mL]"].astype(float), std["A562"].astype(float), "quadratic")
    unk = t[t["Sample name"].str.startswith("Hydrolysate")]
    cu, _ = U.predict_conc(cal, unk["A562"].astype(float).values)
    assert abs(cu[0] - 380) < 15 and abs(cu[1] - 820) < 25


def test_kinetics_first_order_from_names():
    files = sorted(glob.glob(os.path.join(EXD, "OPA_GlyGly_t*min_EXAMPLE.csv")))
    t, a = [], []
    for p in files:
        s = _rd(p).spectra[0]
        t.append(U.time_from_name(os.path.basename(p))); a.append(U.value_at(s.x, s.y, 340))
    fit = U.first_order_fit(t, a)
    assert abs(fit["k"] - 0.012) < 0.0015 and fit["R2"] > 0.99
    assert U.time_from_name("x 2 h", "min") == 120


def test_tauc_band_gaps_and_manual_window():
    f = _rd(os.path.join(EXD, "DRS_reflectance_EXAMPLE.csv"))
    eg = {}
    for s in f.spectra:
        E, T = U.tauc(s.x, s.y, "direct", True)
        eg[s.name] = U.band_gap(E, T)["Eg"]
    assert 4.8 < eg["ZrO2 %R"] < 5.4 and 3.9 < eg["Sample A %R"] < 4.4 and eg["Sample B %R"] < eg["Sample A %R"]
    E = np.linspace(3, 5, 200); T = np.clip(10 * (E - 4.0), 0, None)
    assert abs(U.band_gap(E, T, (4.2, 4.8))["Eg"] - 4.0) < 1e-6


def test_derivative_zero_crossing():
    x = np.arange(200, 400.0); y = np.exp(-0.5 * ((x - 300) / 12) ** 2)
    d = U.derivative(x, y, 1, 11, 3)
    mx = U.maxima_from_derivative(x, y, d)
    assert len(mx) == 1 and abs(mx[0] - 300) < 0.5
    d2 = U.derivative(x, y, 2, 15, 3)
    assert abs(x[np.argmin(d2)] - 300) <= 1


@pytest.mark.skipif(not REAL, reason="real UV/Vis files not available")
@pytest.mark.parametrize("p", REAL)
def test_real_files(p):
    if p.lower().endswith(".pdf"):
        pytest.skip("report")
    try:
        f = _rd(p)
    except U.IDSEncrypted:
        assert p.lower().endswith(".ids"); return
    assert f.spectra or f.results is not None


def test_uvvis_page_apptest(tmp_path):
    from streamlit.testing.v1 import AppTest
    os.environ.update(MFS_STORE=str(tmp_path / "s.json"), MFS_DB=str(tmp_path / "r.db"), MFS_EXPORT_DIR=str(tmp_path), MFS_NAV_GUARD="0")
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600); at.run()
    at.button(key="nv_" + "UV/Vis").click().run()
    assert not at.exception, [e.value for e in at.exception][:2]
    assert len(at.dataframe) >= 2
    labels = [m.label for m in at.metric]
    assert any("LOD" in l for l in labels) and any("k (" in l for l in labels)
    txt = " ".join(str(d.value.to_dict()) for d in at.dataframe)
    assert "Eg (eV)" in txt
    at.button(key="uv_k_send").click().run()
    at.button(key="nv_" + "Kinetics & adsorption").click().run()
    assert not at.exception and at.checkbox(key="kin_from_uv").value
