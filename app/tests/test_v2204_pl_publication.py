"""v2.20.4: publication checks for PL ratiometric thermometry: model choice by AIC/BIC, validation against an independent
reference temperature, and the spectral response correction. Synthetic data only (the checks are arithmetic, not a model of a sample)."""
import numpy as np

from mfs import plqa as QA

KB = 8.617333262e-5


def _boltzmann(T, a=5.0, dE_eV=0.1):
    """dE in eV (0.1 eV = 806 cm-1)."""
    return a * np.exp(-dE_eV / (KB * T))


def test_model_choice_prefers_simple_law_when_there_is_no_offset():
    T = np.linspace(100, 400, 12)
    D = _boltzmann(T)
    sD = 0.01 * D
    mc = QA.model_comparison(T, D, sD)
    assert mc["status"] == "ok"
    assert mc["preferred"] == "Boltzmann"
    assert abs(mc["models"]["Boltzmann"]["dE_cm1"] - 0.1 * 8065.544) < 2.0      # the fitted gap is the input gap


def test_model_choice_prefers_offset_when_a_background_is_present():
    T = np.linspace(100, 400, 12)
    D = _boltzmann(T) + 0.8                       # constant offset C = 0.8 added to the ratio
    sD = 0.01 * D
    mc = QA.model_comparison(T, D, sD)
    assert mc["preferred"] == "Boltzmann + offset"
    assert mc["models"]["Boltzmann + offset"]["C"] is not None
    chk = QA.model_check(T, D, sD)
    assert chk["status"] == "warn"


def test_model_check_is_ok_for_clean_boltzmann():
    T = np.linspace(100, 400, 12)
    D = _boltzmann(T)
    assert QA.model_check(T, D, 0.01 * D)["status"] == "ok"


def test_model_comparison_needs_enough_points():
    mc = QA.model_comparison([100, 200, 300], [1, 2, 3])
    assert mc["status"] != "ok"


def test_validation_statuses():
    T = np.array([300.0, 310.0, 320.0])
    assert QA.validation_check(T, T + 0.2)["status"] == "ok"
    assert QA.validation_check(T, T + 2.0)["status"] == "warn"
    assert QA.validation_check(T, T + 5.0)["status"] == "fail"
    assert QA.validation_check(T, None)["status"] == "info"
    assert QA.validation_check(T, [300.0])["status"] == "fail"      # wrong length is a failure, not silently ignored


def test_validation_reports_bias_and_rms():
    T = np.array([300.0, 310.0, 320.0])
    chk = QA.validation_check(T, T - 1.5)
    assert "bias +1.50" in chk["text"]
    assert "RMS 1.50" in chk["text"]


def test_correction_check_states_what_is_missing():
    assert QA.correction_check(True)["status"] == "ok"
    chk = QA.correction_check(False)
    assert chk["status"] == "info"
    assert "energy gap" in chk["detail"]
