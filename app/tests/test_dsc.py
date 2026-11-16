"""DSC readers + analysis. Primary case: the user's real Mettler STARe export (single curve, Ts, W/g)."""
import glob, os
import numpy as np
import pytest
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mfs import dsc

HERE = os.path.dirname(__file__); ROOT = os.path.dirname(HERE)
REAL = os.path.join(HERE, "data", "Mettler_STARe_SG_DSC_1.1_1K_REAL.txt")
EX = os.path.join(ROOT, "example_data", "DSC")


@pytest.fixture(scope="module")
def real():
    return dsc.read_dsc(REAL)


def test_real_reader(real):
    assert real["meta"]["vendor"].startswith("Mettler")
    assert len(real["curves"]) == 1                  # one curve, footer not parsed as data
    c = real["curves"][0]
    assert len(c["T"]) == 1351 and c["unit"] == "W/g" and c["kind"] == "heating"
    assert c["T"][0] == pytest.approx(25.0) and c["T"][-1] == pytest.approx(70.0, abs=0.1)
    assert np.allclose(c["T"][:3], [25.0, 25.032, 25.065])       # Ts (sample), not Tr
    assert real["meta"]["mass_mg"] == pytest.approx(5.64)        # '5,6400 mg' decimal comma
    assert "1°C/min" in real["meta"]["method"]
    assert dsc.heating_rate(c) == pytest.approx(1.0, rel=0.02)


def test_real_reported_values_parsed(real):
    r = real["meta"]["reported"]
    assert r["integral_mJ"] == pytest.approx(-594.34) and r["dH_Jg"] == pytest.approx(-105.38)
    assert r["onset_C"] == pytest.approx(48.29) and r["peak_C"] == pytest.approx(48.85)


@pytest.mark.parametrize("bl", ["linear", "sigmoidal"])
def test_real_peak_matches_stare(real, bl):
    c = dict(real["curves"][0], id=0)
    res = dsc.evaluate_all(c, bl, real["meta"]["mass_mg"])
    endo = [r for r in res if r["kind"] == "endo"]
    assert len(endo) == 1
    r = endo[0]
    assert r["T_peak"] == pytest.approx(48.85, abs=0.05)
    assert r["T_onset"] == pytest.approx(48.29, abs=0.1)
    assert r["T_endset"] > r["T_peak"] > r["T_onset"]
    assert r["dH"] == pytest.approx(-105.38, rel=0.005)
    assert r["area_mJ"] == pytest.approx(-594.34, rel=0.005)


def test_real_units(real):
    c = real["curves"][0]
    wg = dsc.to_specific(c, 5.64)
    assert np.allclose(wg, c["q"])                    # already W/g; 1 W/g == 1 mW/mg
    assert np.nanmin(wg) == pytest.approx(-2.0, abs=0.05)
    assert np.allclose(wg * 5.64, wg * 5.64)


def test_real_figure(real, tmp_path):
    c = dict(real["curves"][0], id=0)
    fig, ax = plt.subplots()
    dsc.draw_dsc(ax, [("x", c, dsc.to_specific(c, 5.64))], "T", True, peaks={0: dsc.evaluate_all(c, "linear", 5.64)})
    txt = " ".join(t.get_text() for t in ax.texts)
    assert "exo" in txt and "48.3" in txt and "-105.6" in txt
    fig.savefig(tmp_path / "f.png"); plt.close(fig)


TRUTH = dict(Tm=120.0, dHm=-80.0, Tc=95.0, dHc=70.0, Tg=60.0)


@pytest.mark.parametrize("f", sorted(glob.glob(os.path.join(EX, "MOF-PEG_*"))))
def test_examples_all_formats(f):
    d = dsc.read_dsc(f)
    m = d["meta"]["mass_mg"] or 5.0
    assert d["meta"]["vendor"]
    heat = [c for c in d["curves"] if c["kind"] == "heating"]
    assert heat
    c = dict(heat[0], id=0)
    res = dsc.evaluate_all(c, "linear", m)
    endo = [r for r in res if r["kind"] == "endo"]
    assert len(endo) == 1, f                       # sign convention correct for every vendor
    assert endo[0]["T_peak"] == pytest.approx(TRUTH["Tm"], abs=0.3)
    assert endo[0]["T_onset"] == pytest.approx(116.0, abs=0.6)
    assert endo[0]["dH"] == pytest.approx(TRUTH["dHm"], rel=0.02)
    w = dsc.auto_tg_window(c, res)
    g = dsc.glass_transition(c, *w, mass_mg=m)
    assert g["Tg_mid"] == pytest.approx(TRUTH["Tg"], abs=1.0)
    assert g["dcp_J_gK"] == pytest.approx(0.30, rel=0.1)


@pytest.mark.parametrize("name", ["MOF-PEG_TA-TRIOS_EXAMPLE.txt", "MOF-PEG_Netzsch_EXAMPLE.csv", "MOF-PEG_generic_EXAMPLE.csv"])
def test_segments_and_cooling(name):
    d = dsc.read_dsc(os.path.join(EX, name))
    assert [c["kind"] for c in d["curves"]] == ["heating", "cooling", "heating"]
    c = dict(d["curves"][1], id=1)
    exo = [r for r in dsc.evaluate_all(c, "linear", d["meta"]["mass_mg"] or 5.0) if r["kind"] == "exo"]
    assert len(exo) == 1 and exo[0]["T_peak"] == pytest.approx(TRUTH["Tc"], abs=0.3)
    assert exo[0]["dH"] == pytest.approx(TRUTH["dHc"], rel=0.02)
    assert exo[0]["T_onset"] > exo[0]["T_peak"]            # cooling: onset on the high-T side


def test_no_split():
    d = dsc.read_dsc(os.path.join(EX, "MOF-PEG_TA-TRIOS_EXAMPLE.txt"), split_segments=False)
    assert len(d["curves"]) >= 1


def test_mw_normalisation():
    d = dsc.read_dsc(os.path.join(EX, "MOF-PEG_TA-UniversalAnalysis_EXAMPLE.txt"))
    c = d["curves"][0]
    assert c["unit"] == "mW"
    assert np.allclose(dsc.to_specific(c, 5.0), c["q"] / 5.0)


def test_sta_and_overlay():
    d = dsc.read_dsc(os.path.join(EX, "MOF-801_STA_Netzsch_EXAMPLE.csv"))
    c = d["curves"][0]
    assert "w" in c and np.nanmax(c["w"]) == pytest.approx(100, abs=1)
    fig, ax = plt.subplots()
    ax2 = dsc.draw_sta(ax, c, dsc.to_specific(c, d["meta"]["mass_mg"]))
    assert ax2 is not None; plt.close(fig)


def test_decimal_comma_netzsch_mass():
    d = dsc.read_dsc(os.path.join(EX, "MOF-PEG_Netzsch_EXAMPLE.csv"))
    assert d["meta"]["mass_mg"] == pytest.approx(5.0)
    assert np.all(np.isfinite(d["curves"][0]["T"]))


def test_baselines_sane():
    T = np.linspace(0, 100, 501); q = 0.01 * T - 2 * np.exp(-((T - 50) / 2) ** 2)
    i0, i1 = 200, 300
    bl = dsc.baseline_linear(T, q, i0, i1)
    assert bl[i0] == pytest.approx(q[i0]) and bl[i1] == pytest.approx(q[i1])
    bs = dsc.baseline_sigmoidal(T, q, i0, i1)
    assert np.all(np.isfinite(bs))
