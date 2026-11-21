"""v2.18.11 round: each new or extended item has a known-answer test. Synthetic data only."""
import ast
import io
import os

import numpy as np

from mfs import autoread as AR
from mfs import isotherm_models as IM
from mfs import ui as UI

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# 5. Dubinin-Radushkevich: the mean energy E = 12 kJ/mol and W0 are recovered at 77.35 K
def test_dubinin_radushkevich_recovers_energy_and_capacity():
    R, T, E, W0 = 8.314462618, 77.35, 12000.0, 0.5
    p = np.linspace(1e-4, 0.3, 40)
    q = W0 * np.exp(-((R * T * np.log(1.0 / p)) / E) ** 2)
    res = IM.dubinin_radushkevich(p, q, p0=1.0, temperature=T)
    assert res["status"] == "ok"
    assert abs(res["E_kJ"] - 12.0) < 1e-6
    assert abs(res["W0"] - W0) < 1e-6


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
