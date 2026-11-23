"""v2.18.6: the PL analysis report (mfs/plreport.py) builds a valid multi-page PDF with every section."""
import io
import pytest
import numpy as np
import pandas as pd
from mfs import pl as PL
from mfs import plreport as PRPT

KB = 8.617333262e-5


def _fit():
    T = np.arange(100.0, 401.0, 25.0)
    D = 1.5 * np.exp(-700 / 8065.544 / (KB * T)) * (1 + 0.01 * np.sin(np.arange(len(T))))
    Dstat = 0.02 * D
    bf = PL.fit_boltzmann(T, D, sigma=Dstat)
    return T, D, Dstat, bf


def test_report_is_a_pdf_with_all_sections():
    T, D, Dstat, bf = _fit()
    table = pd.DataFrame({"T (K)": T, "Δ": D, "σΔ": Dstat})
    checks = [{"status": "ok", "title": "Temperatures", "detail": "9 temperatures", "advice": ""},
              {"status": "warn", "title": "Boltzmann residuals", "detail": "scatter 1.4x", "advice": "check the ratio errors"}]
    pdf = PRPT.build_pdf(sample="Test sample", version="2.18.6", T=T, D=D, D_stat=Dstat, bf=bf,
                         Sr=PL.relative_sensitivity(T, bf["dE_eV"]), table=table,
                         budget_rows=[{"component": "statistics", "value": "5 cm-1"}],
                         checks=checks, verdict="HIGH", settings={"ratio method": "fit"}, u_dE_cm1=8.0,
                         T_res=np.full(len(T), 0.3), T_cal=np.full(len(T), 0.5))
    assert pdf[:4] == b"%PDF"
    import re as _re
    assert len(_re.findall(rb"/Type\s*/Page[^s]", pdf)) >= 2
    pypdf = pytest.importorskip("pypdf")
    reader = pypdf.PdfReader(io.BytesIO(pdf))
    text = "\n".join(p.extract_text() or "" for p in reader.pages)
    for section in ["1. Verdict", "2. Methods", "3. Model choice and validation", "4. Settings used", "5. Equations",
                    "6. Calibration", "7. Thermometry table", "8. Uncertainty budget", "9. Checks and advice", "10. Limitations"]:
        assert section in text, section


def test_report_builds_without_checks_and_budget():
    T, D, Dstat, bf = _fit()
    pdf = PRPT.build_pdf(sample="Empty", version="2.18.6", T=T, D=D, D_stat=Dstat, bf=bf,
                         Sr=PL.relative_sensitivity(T, bf["dE_eV"]), table=None, budget_rows=[], checks=[],
                         verdict="LOW", settings={}, u_dE_cm1=1.0, T_res=np.ones(len(T)), T_cal=np.ones(len(T)))
    assert pdf[:4] == b"%PDF" and len(pdf) > 5000
