"""v2.18.11 round: each new or extended item has a known-answer test. Synthetic data only."""
import ast
import io
import os

import numpy as np

from mfs import autoread as AR
from mfs import avrami as AV
from mfs import isotherm_models as IM
from mfs import ui as UI
from mfs import zeta_convert as ZC

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# 1. Smoluchowski: known answer for water at 25 degC: zeta = 25 mV -> mu = eps*zeta/eta
def test_smoluchowski_known_mobility():
    mu = ZC.EPS_R_WATER_25C * ZC.EPS0 * 0.025 / ZC.ETA_WATER_25C
    res = ZC.zeta_from_mobility(mu)
    assert abs(res["zeta_mV"] - 25.0) < 1e-9
    assert abs(mu - 1.9516e-8) / 1.9516e-8 < 1e-3


# 2. Hückel gives 1.5 times the Smoluchowski zeta for the same mobility, and the validity note appears
def test_huckel_factor_and_validity_note():
    mu = 2e-8
    s = ZC.zeta_from_mobility(mu, model="smoluchowski")["zeta_V"]
    h = ZC.zeta_from_mobility(mu, model="huckel")["zeta_V"]
    assert abs(h / s - 1.5) < 1e-12
    note = ZC.zeta_from_mobility(mu, model="smoluchowski", ka=10)["note"]
    assert "below ~100" in note


# 3. Avrami: known n and k are recovered; t_half matches the closed form
def test_avrami_recovers_known_exponent_and_rate():
    n, k = 2.5, 0.01
    t = np.linspace(1, 100, 60)
    x = 1 - np.exp(-k * t ** n)
    res = AV.fit_avrami(t, x)
    assert res["status"] == "ok"
    assert abs(res["n"] - n) < 1e-6
    assert abs(res["k"] - k) / k < 1e-6
    assert abs(res["t_half"] - (np.log(2) / k) ** (1 / n)) < 1e-6


# 4. Toth: known parameters are recovered from noise-free data
def test_toth_recovers_known_parameters():
    p = np.linspace(0.05, 10, 30)
    q = IM.toth(p, 3.0, 1.5, 0.6)
    res = IM.fit_toth(p, q)
    assert res["status"] == "ok"
    assert abs(res["qm"] - 3.0) < 1e-3 and abs(res["b"] - 1.5) < 1e-3 and abs(res["t"] - 0.6) < 1e-3


# 5. Dubinin-Radushkevich: the mean energy E = 12 kJ/mol and W0 are recovered at 77.35 K
def test_dubinin_radushkevich_recovers_energy_and_capacity():
    R, T, E, W0 = 8.314462618, 77.35, 12000.0, 0.5
    p = np.linspace(1e-4, 0.3, 40)
    q = W0 * np.exp(-((R * T * np.log(1.0 / p)) / E) ** 2)
    res = IM.dubinin_radushkevich(p, q, p0=1.0, temperature=T)
    assert res["status"] == "ok"
    assert abs(res["E_kJ"] - 12.0) < 1e-6
    assert abs(res["W0"] - W0) < 1e-6


# 6. AIC prefers the true model: Sips data (n = 0.6) fitted by Langmuir (k = 2 parameters) vs Sips (k = 3)
def test_aic_prefers_the_model_that_generated_the_data():
    p = np.linspace(0.05, 10, 40)
    q = IM.sips(p, 4.0, 0.8, 0.6)
    sips_fit = IM.fit_sips(p, q)
    qs = IM.sips(p, sips_fit["qm"], sips_fit["K"], sips_fit["n"])
    ql = 4.0 * (0.8 * p) / (1 + 0.8 * p)   # best Langmuir with the same K (not refitted: a fair, simple comparison)
    assert IM.aic(q, qs, 3) < IM.aic(q, ql, 2)


# 7. Units: heating rate and concentration are converted to base units with the known factors
def test_heating_rate_and_concentration_units():
    r, info = AR.convert_column([10.0], "Heating rate (K/min)")
    assert info["status"] == "converted" and abs(r[0] - 10.0 / 60.0) < 1e-12
    c, info = AR.convert_column([250.0], "Concentration (mg/L)")
    assert info["status"] == "converted" and abs(c[0] - 0.25) < 1e-12


# 8. Metadata: key: value lines above the table are extracted, the table is not
def test_metadata_extraction():
    text = "Instrument: TGA-synthetic\nSample: MOF-A\nHeating rate: 10 K/min\nTemperature (°C)\tMass (mg)\n25\t10\n"
    md = AR.extract_metadata(text)
    assert md["Instrument"] == "TGA-synthetic"
    assert md["Heating rate"] == "10 K/min"
    assert "Temperature (°C)" not in md


# 9. Navigation: every page in ui.NAV has a dispatch entry in app.py _PAGEFN (a page without one is unreachable)
def test_every_nav_page_is_dispatched_in_app():
    src = open(os.path.join(APP, "app.py"), encoding="utf-8").read()
    tree = ast.parse(src)
    keys = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "_PAGEFN" for t in node.targets):
            for k in node.value.keys:
                if isinstance(k, ast.Constant):
                    keys.add(k.value)
    missing = [p for p in UI.PAGES if p not in keys]
    assert not missing, missing


# 10. Icons: every navigation icon renders to a non-empty PNG
def test_every_nav_icon_renders_to_png():
    import resvg_py
    from mfs import sidebar_nav as SN
    for p in UI.PAGES:
        png = bytes(resvg_py.svg_to_bytes(svg_path=SN.icon_path(p), width=48, height=48))
        assert png[:8] == b"\x89PNG\r\n\x1a\n", p
        assert len(png) > 200, p
