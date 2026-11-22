"""v2.19.2 round: known-answer tests for the new methods and checks.

Each test uses synthetic data whose answer is known in advance (noise-free model curves or exact arithmetic)."""
import numpy as np
import pytest

from mfs import gurvich, kinetics_nonisothermal as K, problem_report as PR, saturation as SAT, semcal, tga


# ---- Hill plot (saturation.hill_fit): theta = x^n / (K^n + x^n) with K = 10, n = 2.5
def test_hill_fit_recovers_known_K_and_n():
    x = np.logspace(-1, 2, 25)
    theta = x ** 2.5 / (10.0 ** 2.5 + x ** 2.5)
    out = SAT.hill_fit(x, theta)
    assert out["status"] == "ok"
    assert out["K"] == pytest.approx(10.0, rel=1e-4)
    assert out["n"] == pytest.approx(2.5, rel=1e-4)
    assert out["R2"] > 0.999999


def test_hill_fit_needs_points_inside_the_curve():
    out = SAT.hill_fit([1, 2, 3], [0.0, 1.0, 1.0])
    assert "status" in out and out["status"] != "ok"


# ---- Michaelis-Menten and Lineweaver-Burk: v = Vmax S / (Km + S), Vmax = 2, Km = 5
def test_michaelis_menten_and_lineweaver_burk_recover_known_values():
    S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
    v = 2.0 * S / (5.0 + S)
    out = SAT.michaelis_menten_fit(S, v)
    assert out["status"] == "ok"
    assert out["Vmax"] == pytest.approx(2.0, rel=1e-6)
    assert out["Km"] == pytest.approx(5.0, rel=1e-6)
    assert out["lineweaver_burk"]["Vmax"] == pytest.approx(2.0, rel=1e-6)
    assert out["lineweaver_burk"]["Km"] == pytest.approx(5.0, rel=1e-6)


# ---- Gurvich pore volume: 22.414 cm³ STP/g = 1 mmol/g; liquid N2 molar volume 34.67 cm³/mol -> 0.03467 cm³/g
def test_gurvich_volume_known_answer():
    out = gurvich.gurvich_volume(22.414)
    assert out["status"] == "ok"
    assert out["n_mmol_per_g"] == pytest.approx(1.0, rel=1e-9)
    assert out["V_pore_cm3_per_g"] == pytest.approx(0.03467, rel=1e-9)


# ---- Friedman differential method: first-order model, Ea = 120 kJ/mol, A = 1e7 1/s
def _first_order_curve(beta_k_per_s, Ea=120e3, A=1e7, T0=400.0, T1=950.0):
    from scipy.integrate import solve_ivp
    R = 8.314462618

    def rhs(T, a):
        return [A / beta_k_per_s * np.exp(-Ea / (R * T)) * (1.0 - a[0])]

    T = np.linspace(T0, T1, 4000)
    sol = solve_ivp(rhs, (T0, T1), [0.0], t_eval=T, rtol=1e-10, atol=1e-12)
    return T, sol.y[0]


def test_friedman_recovers_known_activation_energy():
    curves = []
    for beta_k_min in (1.0, 2.0, 4.0):
        T, a = _first_order_curve(beta_k_min / 60.0)
        curves.append((beta_k_min / 60.0, T, a))
    out = K.friedman(curves, 0.5)
    assert out["status"] == "ok"
    assert out["Ea_kJ"] == pytest.approx(120.0, rel=0.01)
    assert out["n_points"] == 3


def test_friedman_needs_three_rates():
    T, a = _first_order_curve(1.0 / 60.0)
    assert K.friedman([(1.0 / 60.0, T, a)], 0.5)["status"] != "ok"


# ---- heating-rate spread check (before Kissinger / KAS)
def test_beta_spread_warning_known_cases():
    assert K.beta_spread_warning([5, 10, 20], [400, 410, 425])["ok"] is True
    narrow = K.beta_spread_warning([5, 5.05, 5.1])
    assert narrow["ok"] is False and any("span only" in w for w in narrow["warnings"])
    falling = K.beta_spread_warning([5, 10, 20], [420, 410, 400])
    assert falling["ok"] is False and any("does not rise" in w for w in falling["warnings"])


# ---- Kissinger against Friedman agreement: 100 vs 110 kJ/mol (9.5 % apart) agrees; 100 vs 130 (26 %) does not
def test_kissinger_friedman_cross_check_known_answer():
    assert K.kissinger_friedman_check(100.0, 110.0)["agree"] is True
    bad = K.kissinger_friedman_check(100.0, 130.0)
    assert bad["agree"] is False
    assert bad["relative_difference"] == pytest.approx(30.0 / 115.0, rel=1e-9)


# ---- TGA step uncertainty from the balance resolution: sqrt(2) × 0.001 mg on 10 mg = 0.014142 %
def test_tga_step_uncertainty_known_answer():
    out = tga.step_uncertainty(10.0, 0.001)
    assert out["sigma_pct"] == pytest.approx(0.0141421356, rel=1e-8)
    assert out["sigma_mg"] == pytest.approx(2 ** 0.5 * 0.001, rel=1e-12)


def test_tga_step_with_uncertainty_mass_loss():
    out = tga.step_with_uncertainty(20.0, 10.0, 0.001)
    assert out["mass_loss_mg"] == pytest.approx(2.0, rel=1e-12)
    assert out["step_err_pct"] == pytest.approx(0.0141421356, rel=1e-8)


# ---- Tesseract message (the OCR itself is excluded from the gate): None when found, install hint when missing
def test_ocr_missing_message_names_the_install_command(monkeypatch):
    monkeypatch.setattr(semcal, "tesseract_cmd", lambda: None)
    msg = semcal.ocr_missing_message()
    assert "sudo apt install tesseract-ocr" in msg and "brew install tesseract" in msg and "$MFS_TESSERACT" in msg


def test_ocr_message_is_none_when_tesseract_is_found(monkeypatch):
    monkeypatch.setattr(semcal, "tesseract_cmd", lambda: "/usr/bin/tesseract")
    assert semcal.ocr_missing_message() is None


# ---- problem report: detected file family of each upload (names and families only, no data)
class _Upload:
    def __init__(self, name, data):
        self.name = name
        self._d = data

    def getvalue(self):
        return self._d


def test_problem_report_records_file_families_without_data():
    jcamp = b"##TITLE=spectrum\n##JCAMP-DX=5.00\n##XYDATA=(X++(Y..Y))\n1 2\n"
    csv = b"Wavelength,Absorbance\n200,0.10\n201,0.12\n202,0.15\n203,0.19\n"
    session = {}
    PR.note_file_families(session, [_Upload("/home/x/run1.jdx", jcamp), _Upload("scan.csv", csv)])
    fams = session["_file_families"]
    assert fams["run1.jdx"] == "jcamp-dx"
    assert fams["scan.csv"] == "delimited table"
    assert "0.12" not in str(fams)          # no data values in the record


def test_problem_report_family_helper_never_raises():
    session = {}
    PR.note_file_families(session, [object()])     # no getvalue(): skipped quietly
    assert session.get("_file_families", {}) == {}
