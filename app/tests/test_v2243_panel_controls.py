"""v2.24.3 (S-Structuur-1, S-Structuur-2): (1) the grid estimate of x is robust to a small 2θ noise (0.02° sd, fixed seed);
(2) the Vegard panel opens with two HEXAGONAL library end members (the non-cubic branch of the panel).
Synthetic numbers for (1); repo reference entries for (2). pymatgen and gemmi are required (no importorskip)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pytest
from mfs import vegard as VG

WL = 1.540593
CORUNDUM = (4.7606, 4.7606, 12.994, 90.0, 90.0, 120.0)
B_HEX = (4.9, 4.9, 13.5, 90.0, 90.0, 120.0)
HKLS = [(1, 0, 4), (1, 1, 0), (0, 1, 2), (1, 1, 3), (2, 1, 4), (0, 0, 6)]


# Measured (2026-11-27, 40 seeds 0..39, noise sd 0.02°, x_true 0.37): median |dx| 0.0040, p95 0.0091, max 0.0170.
# The bound 0.02 is the measured maximum plus a small margin; it is a measured value, not a guess.
NOISE_SD_DEG = 0.02
MEASURED_MAX_ABS_DX = 0.0170
BOUND_DX = 0.02


def test_grid_estimate_is_robust_to_two_theta_noise_over_many_seeds():
    x_true = 0.37
    cell = VG.cell_from_composition(x_true, CORUNDUM, B_HEX)
    errs = []
    for seed in range(40):
        rng = np.random.default_rng(seed)
        peaks = []
        for h in HKLS:
            tt = VG.two_theta_for_cell(cell, h, WL)
            peaks.append((tt + rng.normal(0.0, NOISE_SD_DEG), h))
        r = VG.estimate_x_grid(peaks, WL, CORUNDUM, B_HEX, tol_deg=0.2)
        assert r["status"] in ("ok", "end member"), (seed, r)
        errs.append(abs(r["x"] - x_true))
    assert max(errs) <= BOUND_DX, max(errs)
    assert max(errs) == pytest.approx(MEASURED_MAX_ABS_DX, abs=0.005)


def test_vegard_panel_opens_with_two_hexagonal_library_members(monkeypatch):
    """The panel with corundum (A) and zincite (B) from the repo library: a hexagonal pair must reach the figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mfs.pages import spectra as SP
    from streamlit.delta_generator import DeltaGenerator
    from mfs import phaseid as PID
    names = [e["name"] for e in PID.library() if VG.system_of(tuple(e["cell"])) == "hexagonal"]
    assert any("corundum" in n.lower() for n in names) and any("zn" in n.lower() for n in names), names

    def fake_selectbox(self, label, options, index=0, key=None, **kw):
        label = str(label)
        if label.startswith("Source of end member"):
            return "Library"
        if label == "End member A":
            return next(o for o in options if "corundum" in o.lower())
        if label == "End member B":
            return next(o for o in options if "zinc" in o.lower() or "zincite" in o.lower())
        return options[index] if options else None

    monkeypatch.setattr(DeltaGenerator, "selectbox", fake_selectbox)
    calls = {}
    monkeypatch.setattr(SP, "newfig", lambda aspect=1.0, **k: plt.subplots(), raising=False)
    monkeypatch.setattr(SP, "show", lambda fig, name, *a, **k: calls.__setitem__("name", name), raising=False)
    x = np.linspace(5, 60, 2000)
    y = 100 + 1000 * np.exp(-((x - 43.4) / 0.1) ** 2)
    SP._vegard_ui({"S1": (x, y)}, 5.0, 60.0)
    assert calls.get("name") == "PXRD_Vegard_S1", "the hexagonal pair did not reach the figure"
