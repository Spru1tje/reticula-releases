"""Rigaku SmartLab readers (v2.9).

* ``.rasx`` (SmartLab Studio II): a ZIP with ``Data0/Profile0.txt`` (2θ, counts, attenuation-factor columns; the
  intensity is counts × attenuation factor) and ``Data0/MesurementConditions0.xml`` (sic – Rigaku's spelling; target
  / wavelength, scan axis, range, step, speed, sample name).  Several ``Data*/Profile*`` scans may be present.
* ``.ras`` (legacy SmartLab Guidance / PDXL text): ``*KEY "value"`` header lines and ``*RAS_INT_START`` …
  ``*RAS_INT_END`` blocks (2θ, counts, attenuation), possibly several blocks per file.

Every scan is returned as a dict ``{x, y, counts, att, meta, name, is_2theta_theta, axis}``.
"""
from __future__ import annotations

import io
import re
import zipfile
import xml.etree.ElementTree as ET

import numpy as np

TWO_THETA_AXES = re.compile(r"^(two\s*theta\s*/?\s*theta|2theta\s*/?\s*theta|2θ\s*/?\s*θ|twothetatheta|theta\s*/\s*2\s*theta|"
                            r"θ\s*/\s*2θ|2-?theta\s*/\s*theta|twotheta\s*theta)$", re.I)
TWO_THETA_ONLY = re.compile(r"^(two\s*theta|2theta|2θ|twothetachi|2theta\s*chi)$", re.I)


def _axis_kind(axis):
    a = (axis or "").strip().replace("_", "")
    if not a:
        return True, "unknown (assumed 2θ/θ)"
    if TWO_THETA_AXES.match(a) or re.fullmatch(r"twothetatheta", a, re.I):
        return True, axis
    return False, axis


def _profile(txt):
    rows = []
    for line in txt.replace("\r", "").split("\n"):
        line = line.strip()
        if not line or line.startswith(("#", "*")):
            continue
        parts = re.split(r"[\s,;]+", line)
        try:
            vals = [float(p) for p in parts[:3]]
        except ValueError:
            continue
        rows.append(vals + [1.0] * (3 - len(vals)) if len(vals) >= 2 else None)
    rows = [r for r in rows if r]
    if not rows:
        raise ValueError("empty profile")
    a = np.asarray(rows, float)
    return a[:, 0], a[:, 1], a[:, 2]


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def _xml_meta(xml_bytes):
    """Flat {path-ish key: text} for every leaf element, plus the fields Reticula uses."""
    root = ET.fromstring(xml_bytes)
    flat = {}

    def walk(el, path):
        kids = list(el)
        name = _local(el.tag)
        p = f"{path}/{name}" if path else name
        if not kids and el.text and el.text.strip():
            flat.setdefault(p, el.text.strip())
        for k, v in el.attrib.items():
            flat.setdefault(f"{p}@{_local(k)}", v)
        for k in kids:
            walk(k, p)
    walk(root, "")

    def find(*names, contains=None):
        for k, v in flat.items():
            leaf = k.rsplit("/", 1)[-1].split("@")[-1]
            if leaf.lower() in [n.lower() for n in names] and (contains is None or contains.lower() in k.lower()):
                return v
        return None
    m = {
        "sample": find("SampleName", "Sample", "Comment"),
        "target": find("TargetName", "Target", "AnodeMaterial"),
        "wavelength_ka1": find("KAlpha1", "WavelengthKAlpha1", "Ka1", "WaveLengthAlpha1"),
        "wavelength_ka2": find("KAlpha2", "WavelengthKAlpha2", "Ka2", "WaveLengthAlpha2"),
        "wavelength_kb": find("KBeta", "WavelengthKBeta"),
        "voltage_kV": find("Voltage"),
        "current_mA": find("Current"),
        "axis": find("AxisName", "ScanAxis", "Axis", contains="Scan") or find("AxisName", "ScanAxis"),
        "mode": find("Mode", "ScanMode", contains="Scan") or find("ScanMode"),
        "start": find("Start", contains="Scan"),
        "end": find("End", "Stop", contains="Scan"),
        "step": find("Step", contains="Scan"),
        "speed": find("Speed", contains="Scan"),
        "speed_unit": find("SpeedUnit", contains="Scan"),
        "unit": find("Unit", contains="Scan"),
        "date": find("StartTime", "Date", "MeasurementDate", "CreateDate"),
        "operator": find("UserName", "Operator"),
        "optics": find("OpticsName", "SlitCondition"),
    }
    return {k: v for k, v in m.items() if v not in (None, "")}, flat


def _num(v):
    try:
        return float(str(v).strip().strip('"'))
    except (TypeError, ValueError):
        return None


def _finish(name, x, counts, att, meta, raw):
    y = counts * att
    ok, axis = _axis_kind(meta.get("axis"))
    o = np.argsort(x)
    return {"name": name, "x": x[o], "y": y[o], "counts": counts[o], "att": att[o], "meta": meta, "raw_meta": raw,
            "is_2theta_theta": ok, "axis": axis}


def read_rasx(data: bytes, name=""):
    z = zipfile.ZipFile(io.BytesIO(data))
    names = z.namelist()
    profs = sorted([n for n in names if re.search(r"(^|/)Profile\d*\.txt$", n, re.I)],
                   key=lambda n: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", n)])
    if not profs:
        raise ValueError(f"{name}: no Data*/Profile*.txt inside the .rasx (members: {', '.join(names[:10])})")
    scans = []
    for p in profs:
        folder = p.rsplit("/", 1)[0] if "/" in p else ""
        idx = re.search(r"Profile(\d*)\.txt$", p, re.I).group(1)
        cands = [n for n in names if n.startswith(folder) and re.search(r"(Mesurement|Measurement)Conditions" + idx + r"\.xml$", n, re.I)] or \
                [n for n in names if n.startswith(folder) and re.search(r"(Mesurement|Measurement)Conditions\d*\.xml$", n, re.I)]
        meta, raw = {}, {}
        if cands:
            try:
                meta, raw = _xml_meta(z.read(cands[0]))
            except ET.ParseError:
                pass
        x, c, a = _profile(z.read(p).decode("utf-8-sig", "replace"))
        scans.append(_finish(f"{folder or 'Data'}/{p.rsplit('/', 1)[-1]}", x, c, a, meta, raw))
    return scans


_RAS_KV = re.compile(r'^\*(\w+)\s+"?(.*?)"?\s*$')


def read_ras(text, name=""):
    scans, header, block, inblock = [], {}, [], False
    for line in text.replace("\r", "").split("\n"):
        s = line.strip()
        if not s:
            continue
        if s.upper().startswith("*RAS_HEADER_START"):
            header = {}
            continue
        if s.upper().startswith("*RAS_INT_START"):
            inblock, block = True, []
            continue
        if s.upper().startswith("*RAS_INT_END"):
            inblock = False
            if block:
                a = np.asarray(block, float)
                att = a[:, 2] if a.shape[1] >= 3 else np.ones(len(a))
                h = header
                meta = {"sample": h.get("FILE_SAMPLE") or h.get("FILE_COMMENT"), "target": h.get("HW_XG_TARGET_NAME"),
                        "wavelength_ka1": h.get("HW_XG_WAVE_LENGTH_ALPHA1"), "wavelength_ka2": h.get("HW_XG_WAVE_LENGTH_ALPHA2"),
                        "wavelength_kb": h.get("HW_XG_WAVE_LENGTH_BETA"), "voltage_kV": h.get("HW_XG_VOLTAGE"),
                        "current_mA": h.get("HW_XG_CURRENT"), "axis": h.get("MEAS_SCAN_AXIS_X"), "mode": h.get("MEAS_SCAN_MODE"),
                        "start": h.get("MEAS_SCAN_START"), "end": h.get("MEAS_SCAN_STOP"), "step": h.get("MEAS_SCAN_STEP"),
                        "speed": h.get("MEAS_SCAN_SPEED"), "speed_unit": h.get("MEAS_SCAN_SPEED_UNIT"),
                        "unit": h.get("MEAS_SCAN_UNIT_X"), "date": h.get("MEAS_SCAN_START_TIME"), "operator": h.get("FILE_OPERATOR")}
                meta = {k: v for k, v in meta.items() if v not in (None, "")}
                scans.append(_finish(f"scan {len(scans) + 1}", a[:, 0], a[:, 1], att, meta, dict(h)))
            continue
        if inblock:
            parts = re.split(r"[\s,]+", s)
            try:
                vals = [float(p) for p in parts[:3]]
                if len(vals) >= 2:
                    block.append(vals + [1.0] * (3 - len(vals)))
            except ValueError:
                pass
            continue
        m = _RAS_KV.match(s)
        if m:
            header[m.group(1).upper()] = m.group(2).strip()
    if not scans:
        raise ValueError(f"{name}: no *RAS_INT_START … *RAS_INT_END data block")
    return scans


def decode_ras(data: bytes):
    for enc in ("utf-8-sig", "cp932", "latin-1"):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1", "replace")


def read_any(data: bytes, name=""):
    if data[:2] == b"PK":
        return read_rasx(data, name)
    return read_ras(decode_ras(data), name)


def store_values(scan):
    """The metadata worth keeping in the sample store (pxrd section)."""
    m = scan["meta"]
    out = {"instrument": "Rigaku SmartLab", "scan_axis": scan["axis"], "scan_is_2theta_theta": bool(scan["is_2theta_theta"])}
    for k in ("target", "wavelength_ka1", "wavelength_ka2", "voltage_kV", "current_mA", "start", "end", "step", "speed",
              "speed_unit", "mode", "date", "sample"):
        if k in m:
            v = _num(m[k])
            out[f"xrd_{k}"] = v if v is not None and k not in ("target", "mode", "date", "sample", "speed_unit") else m[k]
    out["n_points"] = int(len(scan["x"]))
    out["max_attenuation"] = float(np.max(scan["att"])) if len(scan["att"]) else 1.0
    return out
