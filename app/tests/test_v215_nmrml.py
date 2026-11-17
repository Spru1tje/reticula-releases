"""v2.15: nmrML import (FID and processed spectrum)."""
import base64, os, sys, zlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from mfs import nmrproc as NP, nmrpeaks as PK

NS = 'xmlns="http://nmrml.org/schema"'


def make_fid(sf=400.13, sw_ppm=20.0, n=16384, lines=((7.26, 1.0), (1.00, 3.0))):
    sw = sw_ppm * sf
    t = np.arange(n) / sw
    fid = sum(a * np.exp(2j * np.pi * (p * sf) * t) for p, a in lines) * np.exp(-t / 0.8)
    return fid, sw


def nmrml_fid(fid, sf=400.13, sw_ppm=20.0, compressed=True):
    a = np.empty(2 * len(fid)); a[0::2] = fid.real; a[1::2] = fid.imag
    raw = a.astype("<f8").tobytes()
    raw = zlib.compress(raw) if compressed else raw
    b64 = base64.b64encode(raw).decode()
    return (f'<?xml version="1.0"?><nmrML {NS} version="1.0.rc1"><acquisition><acquisition1D><acquisitionParameterSet numberOfScans="16">'
            f'<sampleAcquisitionTemperature value="298" unitName="kelvin"/><DirectDimensionParameterSet numberOfDataPoints="{len(fid)}">'
            f'<acquisitionNucleus name="1H"/><irradiationFrequency value="{sf}" unitName="megahertz"/><sweepWidth value="{sw_ppm}" unitName="parts per million"/>'
            f'</DirectDimensionParameterSet></acquisitionParameterSet><fidData byteFormat="Complex128" compressed="{str(compressed).lower()}" '
            f'encodedLength="{len(b64)}">{b64}</fidData></acquisition1D></acquisition></nmrML>').encode()


def test_nmrml_fid_is_processed_like_bruker_data():
    for comp in (True, False):
        fid, _ = make_fid()
        kind, e = NP.load_nmrml(nmrml_fid(fid, compressed=comp), "demo.nmrML")
        assert kind == "fid" and e.nucleus == "1H" and e.acqus["NS"] == 16 and e.acqus["TE"] == 298
        sp = NP.process(e, p0=None, p1=None, reference="none")
        ppm = sp.ppm[np.argmax(sp.y)]
        assert min(abs(ppm - 1.0), abs(ppm + 1.0)) < 0.05 or abs(ppm - 1.0) < 0.05 or abs(abs(ppm) - 1.0) < 0.05, ppm


def test_nmrml_processed_spectrum():
    ppm = np.linspace(10, -1, 4000)
    y = np.exp(-0.5 * ((ppm - 2.1) / 0.01) ** 2).astype("<f4")
    b64 = base64.b64encode(zlib.compress(y.tobytes())).decode()
    xml = (f'<nmrML {NS}><spectrumList><spectrum1D numberOfDataPoints="4000"><spectrumDataArray byteFormat="float32" compressed="true">{b64}</spectrumDataArray>'
           f'<xAxis unitName="parts per million" startValue="10" endValue="-1"/></spectrum1D></spectrumList>'
           f'<acquisition><acquisition1D><acquisitionParameterSet><DirectDimensionParameterSet><irradiationFrequency value="500.0"/></DirectDimensionParameterSet></acquisitionParameterSet></acquisition1D></acquisition></nmrML>').encode()
    kind, sp = NP.load_nmrml(xml, "x.nmrml")
    assert kind == "spectrum" and abs(sp.ppm[np.argmax(sp.y)] - 2.1) < 0.01 and sp.sf == 500.0
