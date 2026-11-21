"""v2.18.11: the new methods and checks, each on data with a known answer."""
import io
import json
import zipfile

import numpy as np

from mfs import zeta_conv as Z
from mfs import isotherm_models as IM
from mfs import avrami as AV
from mfs import uvvis as UV
from mfs import problem_report as PR
from mfs import syscheck as SC


# ---- zeta from mobility
def test_one_micron_cm_per_Vs_is_about_12_8_mV_in_water_at_25C():
    # textbook rule of thumb: 1 µm·cm/(V·s) ≈ 12.8 mV in water at 25 °C (Smoluchowski)
    r = Z.zeta_from_mobility(1.0, ionic_strength_mM=1.0, radius_nm=None, temp_C=25.0)
    assert r["status"] == "ok"
    assert abs(r["zeta_smoluchowski_mV"] - 12.8) < 0.1


def test_large_ka_gives_smoluchowski_and_small_ka_gives_hueckel_factor():
    # ka = 1000: f -> 1.5 (Smoluchowski); ka = 0.001: f -> 1 (Hueckel)
    assert abs(Z.henry_f(1000.0) - 1.5) < 0.5e-2   # the approximation is 0.25 % below 1.5 at ka = 1000
    assert abs(Z.henry_f(0.001) - 1.0) < 2e-3


def test_zeta_from_henry_recovers_smoluchowski_for_a_large_particle():
    # 1 M ionic strength is a very thin layer; with a 5000 nm radius ka is huge, so Henry == Smoluchowski within 0.1 %
    r = Z.zeta_from_mobility(2.0, ionic_strength_mM=100.0, radius_nm=5000.0)
    assert r["ka"] > 300
    assert abs(r["zeta_mV"] - r["zeta_smoluchowski_mV"]) / r["zeta_smoluchowski_mV"] < 1e-3
    assert "valid" in r["note"]


def test_debye_length_of_1_mM_one_to_one_electrolyte_is_about_9_6_nm():
    # standard value for water at 25 °C: 0.304 nm / sqrt(I in M) -> 0.304 / sqrt(0.001) = 9.6 nm
    assert abs(Z.debye_length_m(1.0) * 1e9 - 9.6) < 0.1


# ---- isotherms: Toth and AICc comparison
def _grid():
    return np.linspace(0.02, 1.0, 14)


def test_toth_with_t_equal_one_is_langmuir():
    p = _grid()
    assert np.allclose(IM.toth(p, 3.0, 4.0, 1.0), IM.langmuir(p, 3.0, 4.0))


def test_fit_toth_recovers_known_parameters():
    p = _grid()
    q = IM.toth(p, 5.0, 3.0, 0.6)
    r = IM.fit_toth(p, q)
    assert r["status"] == "ok"
    assert abs(r["qm"] - 5.0) < 1e-3 and abs(r["b"] - 3.0) < 1e-3 and abs(r["t"] - 0.6) < 1e-3


def test_model_comparison_picks_toth_for_toth_data_and_langmuir_for_langmuir_data():
    p = np.linspace(0.02, 1.0, 40)
    rng = np.random.default_rng(0)
    q_toth = IM.toth(p, 5.0, 3.0, 0.5) * (1 + 0.01 * rng.standard_normal(len(p)))
    assert IM.compare_isotherms(p, q_toth)["best"] == "Toth"
    q_lang = IM.langmuir(p, 5.0, 3.0) * (1 + 0.01 * rng.standard_normal(len(p)))
    res = IM.compare_isotherms(p, q_lang)
    assert res["best"] == "Langmuir"
    assert res["table"][0]["delta_AICc"] == 0.0


def test_comparison_needs_six_points():
    assert IM.compare_isotherms([1, 2, 3], [1, 2, 3])["status"].startswith("needs")


# ---- Avrami
def test_avrami_recovers_known_exponent_and_rate():
    t = np.linspace(0.5, 8.0, 60)
    X = AV.avrami_curve(t, k=0.02, n=2.5)
    r = AV.avrami_fit(t, X)
    assert r["status"] == "ok"
    assert abs(r["n"] - 2.5) < 1e-6
    assert abs(r["k"] - 0.02) / 0.02 < 1e-6
    assert abs(r["t_half"] - (np.log(2) / 0.02) ** (1 / 2.5)) < 1e-4


# ---- Tauc: automatic direct or indirect choice
def _edge(eg, kind_n, x_nm):
    E = 1239.841984 / x_nm
    a = np.where(E > eg, np.clip(E - eg, 0, None) ** kind_n, 0.0)
    return a


def test_tauc_auto_chooses_direct_for_a_direct_edge_at_3_eV():
    x = np.linspace(300.0, 700.0, 801)                 # nm, so E from about 1.8 to 4.1 eV
    alpha = _edge(3.0, 0.5, x) / (1239.841984 / x)     # (alpha h nu)^2 = linear in E: direct edge
    r = UV.tauc_auto_kind(x, alpha)
    assert r["kind"] == "direct"
    assert abs(r["Eg"] - 3.0) < 0.02


# ---- problem report: last five errors, in order
def _zip_json(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        return json.loads(z.read("report.json"))


def test_problem_report_keeps_the_last_five_errors_in_order():
    session = {}
    for i in range(7):
        try:
            raise ValueError(f"error {i}")
        except ValueError as exc:
            PR.record_error(session, f"page {i}", exc)
    rep = _zip_json(PR.build(session, "page 6", "2.18.11", [], ""))
    msgs = [e["message"] for e in rep["error_history"]]
    assert msgs == ["error 2", "error 3", "error 4", "error 5", "error 6"]
    assert rep["last_error"]["message"] == "error 6"


# ---- start-up tool check
def test_missing_tesseract_is_reported_with_hint(tmp_path):
    out = SC.missing_tools(which=lambda name: None, tessdata_prefix=str(tmp_path), platform="linux")
    assert [m["tool"] for m in out] == ["tesseract", "eng.traineddata"]
    assert "apt install tesseract-ocr" in out[0]["hint"]


def test_all_tools_present_gives_empty_list(tmp_path):
    (tmp_path / "eng.traineddata").write_bytes(b"x")
    assert SC.missing_tools(which=lambda name: "/usr/bin/tesseract", tessdata_prefix=str(tmp_path), platform="linux") == []
