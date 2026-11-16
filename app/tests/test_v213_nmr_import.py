"""v2.13 Bruker raw import + automatic processing on public test data (no user spectra yet):
cheminfo bruker-data-test (MIT): XWIN-NMR 3.5 / TopSpin 3.5 (DSPFVS 10/12), TopSpin 4.0.6, TopSpin 4.5.0 float64 fid (Avance NEO style);
BMRB bmse000010 (D-fructose, TopSpin 2.1, 1H + 13C + DEPT-90/135)."""
import gzip, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import numpy as np, pytest
from mfs import nmrproc as N

D = os.path.join(ROOT, "tests", "data", "nmr")


def _exps(fn):
    return N.load(open(os.path.join(D, fn), "rb").read(), fn)


def _corr_ts(e, sp):
    pr = e.procs; n = len(e.proc_1r)
    p = float(pr["OFFSET"]) - np.arange(n) * (float(pr["SW_p"]) / float(pr["SF"])) / n
    yt = np.interp(sp.ppm[::-1], p[::-1], e.proc_1r[::-1])[::-1]; yt, _ = N.auto_baseline(yt)
    from scipy import ndimage as ndi
    m = ndi.binary_dilation(np.abs(yt) > 0.05 * np.abs(yt).max(), iterations=20)
    return np.corrcoef(sp.y[m], yt[m])[0, 1]


@pytest.mark.parametrize("fn,ver,dtype", [("aspirin_300_XWIN35.zip", "XWIN-NMR", "int32"),
                                          ("naphtoic_500_TS35_dspfvs12.zip", "TopSpin 3.5", "int32"),
                                          ("strychnine_400_TS406.zip", "TopSpin 4.0.6", "int32"),
                                          ("neo300_TS450_float.zip", "TopSpin 4.5.0", "float64")])
def test_formats_read_and_process(fn, ver, dtype):
    ex = _exps(fn)
    assert len(ex) == 1
    e = ex[0]; inf = e.info()
    assert ver.lower() in inf["software"].lower() and inf["fid dtype"] == dtype and e.nucleus == "1H"
    sp = N.process(e)
    assert np.all(np.diff(sp.ppm) < 0) and sp.y.max() > 50 * N.noise_sigma(sp.y)
    # phased: negligible negative area around the lines
    neg = -sp.y[sp.y < 0].sum() / sp.y[sp.y > 0].sum()
    assert neg < 0.08, neg


def test_strychnine_matches_topspin_processing():
    e = _exps("strychnine_400_TS406.zip")[0]
    sp = N.process(e, reference="none")
    assert _corr_ts(e, sp) > 0.995
    sp2 = N.process(e)
    assert ("TMS" in sp2.ref_note or "CDCl3" in sp2.ref_note) and abs(sp2.ref_shift) < 0.03


def test_reference_solvent_and_tms():
    sp = N.process(_exps("neo300_TS450_float.zip")[0])
    assert "CDCl3" in sp.ref_note
    i = np.argmin(abs(sp.ppm - 7.26)); w = slice(i - 30, i + 30)
    assert abs(sp.ppm[w][np.argmax(sp.y[w])] - 7.26) < 0.003
    sp = N.process(_exps("aspirin_300_XWIN35.zip")[0])
    assert ("TMS" in sp.ref_note or "CDCl3" in sp.ref_note) and abs(sp.ref_shift) < 0.06


def test_group_delay_table_old_avance():
    e = _exps("aspirin_300_XWIN35.zip")[0]
    assert e.acqus.get("DSPFVS") == 10 and abs(N.group_delay(e.acqus) - 61.0208) < 1e-3     # DECIM 24


def test_fructose_13c_dept_tsp21():
    ex = _exps("fructose_bmse000010_TS21_1H_13C_DEPT.zip")
    kinds = sorted(e.kind for e in ex)
    assert kinds == ["13C", "1H", "DEPT135", "DEPT90"], kinds
    sps = N.process_all(ex, reference="none")
    for e, sp in zip(ex, sps):
        if e.kind != "1H":
            assert _corr_ts(e, sp) > 0.95, (e.kind, _corr_ts(e, sp))
    d135 = next(s for s in sps if s.kind == "DEPT135")
    lines = N.line_list(d135)
    assert any(h < 0 for p, h in lines) and any(h > 0 for p, h in lines)   # CH2 down, CH up


def test_jcamp_ntuples():
    sp = N.load_jcamp(gzip.decompress(open(os.path.join(D, "cytisine_1h_600.jdx.gz"), "rb").read()), "cytisine")
    assert sp.nucleus == "1H" and abs(sp.sf - 600.2) < 0.1 and sp.complex_spec is not None
    assert 0.5 < sp.ppm.min() < 1 and 8 < sp.ppm.max() < 8.5


def test_param_parser_arrays():
    d = N.parse_jcamp_params("##$P= (0..3)\n1 2.5 3 4\n##$SOLVENT= <CDCl3>\n##$NS= 16\n")
    assert d["P"] == [1, 2.5, 3, 4] and d["SOLVENT"] == "CDCl3" and d["NS"] == 16
    assert N.solvent_key("DMSO") == "DMSO-d6" and N.solvent_key("MeOD") == "CD3OD" and N.solvent_key("xyz") is None
