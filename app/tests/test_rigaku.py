"""v2.9: Rigaku SmartLab .rasx (SmartLab Studio II zip) and legacy .ras readers (constructed files – no real one yet)."""
import io, zipfile
import numpy as np, pytest
from mfs import rigaku as RG, readers as R

XML = """<?xml version="1.0" encoding="utf-8"?>
<MeasurementConditions xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <GeneralInformation><SampleName>MOF-801(Hf) batch 7</SampleName><UserName>SG</UserName>
    <StartTime>2026-11-10T10:12:00</StartTime></GeneralInformation>
  <HWConfigurations><XrayGenerator><TargetName>Cu</TargetName><KAlpha1>1.540593</KAlpha1><KAlpha2>1.544414</KAlpha2>
    <KBeta>1.392250</KBeta><Voltage>40</Voltage><Current>30</Current></XrayGenerator></HWConfigurations>
  <ScanInformation><AxisName>{axis}</AxisName><Mode>CONTINUOUS</Mode><Start>{a}</Start><End>{b}</End><Step>0.01</Step>
    <Speed>5</Speed><SpeedUnit>deg/min</SpeedUnit><Unit>deg</Unit></ScanInformation>
</MeasurementConditions>"""


def _scan(a, b, peak):
    x = np.round(np.arange(a, b + 1e-9, 0.01), 4)
    counts = np.round(50 + 4000 * np.exp(-0.5 * ((x - peak) / 0.05) ** 2))
    att = np.where(counts > 3000, 3.2, 1.0)            # attenuator in on the strong peak
    counts = np.round(counts / att)
    return x, counts, att


def _rasx(scans):
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w") as z:
        for i, (axis, a, bb, pk) in enumerate(scans):
            x, c, at = _scan(a, bb, pk)
            z.writestr(f"Data{i}/Profile0.txt", "\n".join(f"{u:.4f}\t{v:.0f}\t{w:.4f}" for u, v, w in zip(x, c, at)))
            z.writestr(f"Data{i}/MesurementConditions0.xml", XML.format(axis=axis, a=a, b=bb))
        z.writestr("root.xml", "<Root/>")
    return b.getvalue()


def test_rasx_single_scan_metadata_and_attenuation():
    data = _rasx([("TwoThetaTheta", 5, 15, 8.6)])
    sc = RG.read_any(data, "a.rasx")
    assert len(sc) == 1
    s = sc[0]
    x, c, at = _scan(5, 15, 8.6)
    assert np.allclose(s["y"], c * at) and s["att"].max() == pytest.approx(3.2)
    assert s["is_2theta_theta"] and s["meta"]["target"] == "Cu" and s["meta"]["wavelength_ka1"] == "1.540593"
    assert s["meta"]["sample"] == "MOF-801(Hf) batch 7" and s["meta"]["step"] == "0.01" and s["meta"]["speed"] == "5"
    v = RG.store_values(s)
    assert v["xrd_wavelength_ka1"] == pytest.approx(1.540593) and v["xrd_target"] == "Cu" and v["scan_is_2theta_theta"]
    assert abs(x[np.argmax(s["y"])] - 8.6) < 0.011


def test_rasx_multi_scan_flags_non_2theta_and_read_pxrd_picks_2theta():
    data = _rasx([("Omega", 0, 10, 4.3), ("TwoThetaTheta", 5, 50, 9.9)])
    sc = RG.read_any(data, "m.rasx")
    assert [s["is_2theta_theta"] for s in sc] == [False, True] and sc[0]["axis"] == "Omega"
    f = io.BytesIO(data); f.name = "m.rasx"
    df = R.read_pxrd(f, "m.rasx")
    assert df.x.min() == pytest.approx(5) and df.x.max() == pytest.approx(50)
    assert len(df.attrs["rigaku_all"]) == 2


def test_rasx_without_profile_raises():
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w") as z:
        z.writestr("x.txt", "nothing")
    with pytest.raises(ValueError, match="Profile"):
        RG.read_any(b.getvalue(), "bad.rasx")


def _ras(blocks):
    out = []
    for axis, a, bb, pk in blocks:
        x, c, at = _scan(a, bb, pk)
        out += ["*RAS_DATA_START", "*RAS_HEADER_START", '*FILE_SAMPLE "ZrO2 ref"', '*HW_XG_TARGET_NAME "Cu"',
                '*HW_XG_WAVE_LENGTH_ALPHA1 "1.540593"', f'*MEAS_SCAN_AXIS_X "{axis}"', f'*MEAS_SCAN_START "{a}"',
                f'*MEAS_SCAN_STOP "{bb}"', '*MEAS_SCAN_STEP "0.01"', '*MEAS_SCAN_SPEED "2"', '*MEAS_SCAN_SPEED_UNIT "deg/min"',
                "*RAS_HEADER_END", "*RAS_INT_START"]
        out += [f"{u:.4f} {v:.0f} {w:.4f}" for u, v, w in zip(x, c, at)]
        out += ["*RAS_INT_END", "*RAS_DATA_END"]
    return "\r\n".join(out).encode("cp932")


def test_legacy_ras_blocks():
    sc = RG.read_any(_ras([("TwoThetaTheta", 20, 40, 30.2), ("Phi", 0, 360, 90)]), "z.ras")
    assert len(sc) == 2 and sc[0]["is_2theta_theta"] and not sc[1]["is_2theta_theta"]
    x, c, at = _scan(20, 40, 30.2)
    assert np.allclose(sc[0]["y"], c * at)
    assert sc[0]["meta"]["sample"] == "ZrO2 ref" and sc[0]["meta"]["speed_unit"] == "deg/min"
    f = io.BytesIO(_ras([("TwoThetaTheta", 20, 40, 30.2)])); f.name = "z.ras"
    assert len(R.read_pxrd(f, "z.ras")) == len(x)


def test_axis_names():
    for a in ("TwoThetaTheta", "2Theta/Theta", "2θ/θ", "Theta/2Theta", ""):
        assert RG._axis_kind(a)[0], a
    for a in ("Omega", "Phi", "TwoThetaChi", "Chi", "Rocking"):
        assert not RG._axis_kind(a)[0], a
