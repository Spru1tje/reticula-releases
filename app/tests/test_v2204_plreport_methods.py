"""v2.20.4: the thermometry PDF has a Methods section with figures, a model section, and a validation section that states clearly
whether an independent reference was given. The text is read back from the PDF (pypdf), so the test checks what the reader sees."""
import io

import numpy as np
import pandas as pd
import pytest

from mfs import pl as PL
from mfs import plqa as QA
from mfs import plreport as R

KB = 8.617333262e-5


def _inputs(with_offset=False):
    T = np.linspace(100, 400, 12)
    D = 0.5 * np.exp(-0.1 / (KB * T)) + (0.02 if with_offset else 0.0)
    sD = 0.01 * D
    bf = PL.fit_boltzmann(T, D, offset=False, sigma=sD)
    mc = QA.model_comparison(T, D, sD)
    T_pred = QA.predict_T(D, bf["alpha"], bf["dE_eV"], bf.get("C", 0.0))
    return T, D, sD, bf, mc, T_pred


def _build(refs=None, with_offset=False, corrected=False):
    T, D, sD, bf, mc, T_pred = _inputs(with_offset)
    return R.build_pdf(sample="test", version="2.20.4", T=T, D=D, D_stat=sD, bf=bf, Sr=PL.relative_sensitivity(T, bf["dE_eV"]),
                       table=pd.DataFrame({"T": T, "D": D}), budget_rows=[], checks=[QA.model_check(T, D, sD)],
                       verdict="moderate", settings={"Boltzmann offset C": False}, u_dE_cm1=20.0,
                       T_res=[0.5] * len(T), T_cal=[0.3] * len(T), model=mc, refs=refs, T_pred=T_pred, corrected=corrected,
                       band_info={"windows": [[500, 520], [540, 560]]})


def _text(pdf):
    pypdf = pytest.importorskip("pypdf")
    reader = pypdf.PdfReader(io.BytesIO(pdf))
    return "\n".join(p.extract_text() or "" for p in reader.pages), len(reader.pages)


def test_methods_section_and_figure_captions_are_in_the_pdf():
    text, pages = _text(_build(refs=None))
    assert "2. Methods" in text
    assert "Figure 1. The analysis chain" in text
    assert "3. Model choice and validation" in text
    assert "Figure 2" in text
    assert pages >= 3


def test_no_reference_says_the_accuracy_is_not_shown():
    text, _ = _text(_build(refs=None))
    assert "no independent reference temperatures were given" in text
    assert "NOT shown" in text


def test_reference_gives_bias_and_rms_and_the_figure():
    T, D, sD, bf, mc, T_pred = _inputs()
    text, _ = _text(_build(refs=T + 0.3))
    assert "Result: bias" in text and "RMS" in text
    assert "Figure 3" in text


def test_correction_status_is_stated_in_the_methods():
    text, _ = _text(_build(refs=None, corrected=False))
    assert "NOT applied" in text
    text2, _ = _text(_build(refs=None, corrected=True))
    assert "applied to every spectrum" in text2


def test_figures_build_and_are_nonempty():
    import matplotlib.pyplot as plt
    T, D, sD, bf, mc, T_pred = _inputs(with_offset=True)
    for fig in (R._flow_figure(), R._model_figure(T, D, sD, mc), R._validation_figure(T_pred, T + 0.2)):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=60)
        assert len(buf.getvalue()) > 2000
        plt.close(fig)
