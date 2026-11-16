"""Full heating/cooling cycle: segment detection, combining files, multi-curve STARe, peaks per segment, supercooling,
cycle plots. Example = REAL heating (SG_DSC_1.1) + SIMULATED cooling (tools/make_dsc_cycle_example.py)."""
import os, datetime
import numpy as np
import pytest
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mfs import dsc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FC = os.path.join(ROOT, "example_data", "DSC", "full_cycle")
EX = os.path.join(ROOT, "example_data", "DSC")
HEAT, COOL = os.path.join(FC, "SG_DSC_heating_Mettler_REAL.txt"), os.path.join(FC, "SG_DSC_cooling_Mettler_EXAMPLE.txt")
TRUTH_TC_ON, TRUTH_TC, TRUTH_DHC = 40.8, 39.6, 98.0


def _eval(segs, mass=5.64):
    pk = {}
    for c in segs:
        c.setdefault("sample", "SG"); c["uid"] = c["seg"]; c["id"] = c["uid"]
        pk[c["uid"]] = dsc.evaluate_all(c, "linear", mass)
    return pk


def _check_cycle(segs):
    pk = _eval(segs)
    h = next(c for c in segs if c["seg"] == "Heat 1"); cl = next(c for c in segs if c["seg"] == "Cool 1")
    m = dsc.main_peak(pk[h["uid"]], "endo"); x = dsc.main_peak(pk[cl["uid"]], "exo")
    assert m["T_peak"] == pytest.approx(48.85, abs=0.05) and m["dH"] == pytest.approx(-105.38, rel=0.005)
    assert x["T_peak"] == pytest.approx(TRUTH_TC, abs=0.15)
    assert x["T_onset"] == pytest.approx(TRUTH_TC_ON, abs=0.2)
    assert x["T_onset"] > x["T_peak"]                         # cooling: onset on the high-T side
    assert x["dH"] == pytest.approx(TRUTH_DHC, rel=0.06)       # long crystallisation tail: ~5 % below the truth
    s = dsc.cycle_summary(segs, pk)
    r = s.iloc[0]
    assert r["heating"] == "Heat 1" and r["cooling"] == "Cool 1"
    assert r["ΔT = Tm − Tc, peak (K)"] == pytest.approx(48.85 - TRUTH_TC, abs=0.2)
    assert r["ΔT, onset (K)"] == pytest.approx(48.3 - TRUTH_TC_ON, abs=0.3)
    assert 0.8 < r["|ΔHc/ΔHm|"] < 1.0
    return pk


def test_case_a_one_file_all_segments():
    d = dsc.read_dsc(os.path.join(FC, "SG_DSC_full-cycle_Mettler_EXAMPLE.txt"))
    segs = dsc.label_segments(d["curves"])
    assert [c["seg"] for c in segs] == ["Heat 1", "Iso 1", "Cool 1"]
    assert [c["kind"] for c in segs] == ["heating", "isothermal", "cooling"]
    assert segs[0]["T"][-1] > 69.9 and segs[2]["T"][0] < 70.1          # clean boundaries
    assert len(segs[0]["T"]) == 1351                                     # = the real heating export
    _check_cycle(segs)


def test_case_b_separate_files_combined_by_time():
    ds = [dsc.read_dsc(COOL), dsc.read_dsc(HEAT)]                        # uploaded in the 'wrong' order
    assert ds[0]["meta"]["timestamp"] > ds[1]["meta"]["timestamp"]
    segs = dsc.combine_cycle(ds, "time")
    assert [c["seg"] for c in segs] == ["Heat 1", "Cool 1"]
    assert segs[1]["t"][0] > segs[0]["t"][-1]                            # continuous time axis
    _check_cycle(segs)


def test_case_b_user_order():
    segs = dsc.combine_cycle([dsc.read_dsc(COOL), dsc.read_dsc(HEAT)], "user")
    assert [c["seg"] for c in segs] == ["Cool 1", "Heat 1"]
    assert np.all(np.diff(np.concatenate([c["t"] for c in segs])) > 0)


def test_case_c_stare_several_curve_blocks():
    d = dsc.read_dsc(os.path.join(FC, "SG_DSC_3-curves_Mettler_EXAMPLE.txt"))
    assert d["meta"]["n_curves_in_file"] == 3
    segs = dsc.label_segments(d["curves"])
    assert [c["seg"] for c in segs] == ["Heat 1", "Iso 1", "Cool 1"]
    t = np.concatenate([c["t"] for c in segs])
    assert np.all(np.diff(t) > 0)                                        # restarted clocks made continuous
    _check_cycle(segs)


def test_heat_cool_heat_labels_and_pairing():
    d = dsc.read_dsc(os.path.join(EX, "MOF-PEG_TA-TRIOS_EXAMPLE.txt"))
    segs = dsc.label_segments(d["curves"])
    assert [c["seg"] for c in segs] == ["Heat 1", "Cool 1", "Heat 2"]
    pk = _eval(segs, 5.0)
    s = dsc.cycle_summary(segs, pk)
    assert list(s["heating"]) == ["Heat 1", "Heat 2"] and list(s["cooling"]) == ["Cool 1", "–"]
    assert s.iloc[0]["ΔT = Tm − Tc, peak (K)"] == pytest.approx(25.0, abs=0.3)


def test_segment_kind_from_rate():
    t = np.linspace(0, 10, 200)
    assert dsc.segment_kind({"T": 25 + 10 * t, "t": t}) == "heating"
    assert dsc.segment_kind({"T": 200 - 10 * t, "t": t}) == "cooling"
    assert dsc.segment_kind({"T": 70 + 0.001 * np.sin(t), "t": t}) == "isothermal"
    assert dsc.segment_kind({"T": np.linspace(70, 25, 50), "t": None}) == "cooling"


def test_split_one_curve_heat_iso_cool():
    t = np.arange(0, 100, 0.05)
    T = np.where(t < 45, 25 + t, np.where(t < 55, 70, 70 - (t - 55)))
    segs = dsc.label_segments(dsc.split_by_rate({"name": "x", "T": T, "t": t, "q": np.zeros_like(T), "unit": "W/g"}))
    assert [c["seg"] for c in segs] == ["Heat 1", "Iso 1", "Cool 1"]
    assert abs(segs[0]["t"][-1] - 45) < 0.2 and abs(segs[2]["t"][0] - 55) < 0.2


def test_timestamp_and_group_names():
    assert dsc.file_timestamp("Title      03.04.2024 09:18\n") == datetime.datetime(2024, 4, 3, 9, 18)
    assert dsc.file_timestamp("Run date 2024-04-03 10:10:05") == datetime.datetime(2024, 4, 3, 10, 10, 5)
    assert dsc.cycle_group_name("SG_DSC_heating_Mettler_REAL") == dsc.cycle_group_name("SG_DSC_cooling_Mettler_EXAMPLE") == "SG_DSC"
    assert dsc.cycle_group_name("sampleA_H1") == dsc.cycle_group_name("sampleA_C1") == "sampleA"


def test_styles():
    one = dsc.default_styles(["A"])
    assert one[("A", "heating")][0] != one[("A", "cooling")][0]
    two = dsc.default_styles(["A", "B"])
    assert two[("A", "heating")][0] == two[("A", "cooling")][0] and two[("A", "cooling")][1] == "--"
    assert two[("A", "heating")][0] != two[("B", "heating")][0]


@pytest.fixture
def items():
    segs = dsc.combine_cycle([dsc.read_dsc(HEAT), dsc.read_dsc(COOL)], "time")
    pk = _eval(segs)
    return [(c["seg"], c, dsc.to_specific(c, 5.64)) for c in segs], pk


def test_plot_cycle_T(items):
    its, pk = items
    fig, ax = plt.subplots()
    st = {("SG", "heating"): ("#112233", "-"), ("SG", "cooling"): ("#445566", "--")}
    dsc.draw_cycle_T(ax, its, True, st, pk)
    cols = [l.get_color() for l in ax.lines if l.get_label() in ("Heat 1", "Cool 1")]
    assert cols == ["#112233", "#445566"]
    assert [l.get_linestyle() for l in ax.lines if l.get_label() == "Cool 1"] == ["--"]
    arrows = [p for p in ax.patches if p.get_gid() == "mfs_dir_arrow"]
    assert len(arrows) == 4
    # heating arrows point to higher T, cooling arrows to lower T
    for p in arrows:
        (x0, _), (x1, _) = p._posA_posB
        assert x0 != x1
    xs = [(p._posA_posB[1][0] - p._posA_posB[0][0]) for p in arrows]
    assert sum(x > 0 for x in xs) == 2 and sum(x < 0 for x in xs) == 2
    txt = " ".join(t.get_text() for t in ax.texts)
    assert "-105.6" in txt and "39.6" in txt
    plt.close(fig)


def test_plot_cycle_time(items):
    its, pk = items
    fig, ax = plt.subplots()
    tw = dsc.draw_cycle_time(ax, its, True, None, pk, program="Tr")
    assert tw.get_ylabel() == "Temperature (°C)"
    T = tw.lines[0].get_ydata(); assert T.min() > 24 and T.max() < 71
    assert len([p for p in ax.patches if p.get_gid() == "mfs_shade"]) == 2
    assert {"Heat 1", "Cool 1"} <= {t.get_text() for t in ax.texts}
    plt.close(fig)


def test_plot_stacked(items):
    its, pk = items
    fig, axs = plt.subplots(2, 1, sharex=True)
    dsc.draw_cycle_stacked(axs, [i for i in its if i[1]["kind"] == "heating"], [i for i in its if i[1]["kind"] == "cooling"],
                           True, None, pk)
    assert [l.get_label() for l in axs[0].lines][0] == "Heat 1" and [l.get_label() for l in axs[1].lines][0] == "Cool 1"
    plt.close(fig)
