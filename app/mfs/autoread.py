"""v2.18.8: automatic reading of instrument text exports with unit recognition.

- Finds the header line by itself (the first line that has at least two columns and is followed by numbers), so metadata lines
  above the table are skipped; the delimiter (tab, comma, semicolon) is detected.
- Recognises the unit of each column from its header, e.g. 'Temperature (°C)', 'Power / mW', 'Mass [mg]', 'Time (min)'.
- Converts every recognised column to the base unit of its quantity (K, W, g, s, %, cm⁻¹, Pa) and reports what it did.
Columns whose unit is not recognised are kept unchanged and listed as 'unknown unit', so nothing is silently guessed.
General practice: check the 'units' report before you trust a converted column."""
from __future__ import annotations

import io
import re

import numpy as np

# quantity -> {unit spelling (lower case, no spaces): (base unit, factor to base, offset to base)}
# base value = value * factor + offset
_UNITS = {
    "temperature": {"k": ("K", 1.0, 0.0), "kelvin": ("K", 1.0, 0.0), "°c": ("K", 1.0, 273.15), "c": ("K", 1.0, 273.15),
                    "degc": ("K", 1.0, 273.15), "°f": ("K", 5 / 9, 255.3722), "f": ("K", 5 / 9, 255.3722)},
    "power": {"w": ("W", 1.0, 0.0), "mw": ("W", 1e-3, 0.0), "µw": ("W", 1e-6, 0.0), "uw": ("W", 1e-6, 0.0), "kw": ("W", 1e3, 0.0)},
    "mass": {"g": ("g", 1.0, 0.0), "mg": ("g", 1e-3, 0.0), "µg": ("g", 1e-6, 0.0), "ug": ("g", 1e-6, 0.0), "kg": ("g", 1e3, 0.0)},
    "time": {"s": ("s", 1.0, 0.0), "sec": ("s", 1.0, 0.0), "min": ("s", 60.0, 0.0), "h": ("s", 3600.0, 0.0), "hr": ("s", 3600.0, 0.0),
             "ms": ("s", 1e-3, 0.0)},
    "pressure": {"pa": ("Pa", 1.0, 0.0), "kpa": ("Pa", 1e3, 0.0), "mbar": ("Pa", 100.0, 0.0), "bar": ("Pa", 1e5, 0.0),
                 "mmhg": ("Pa", 133.322, 0.0), "torr": ("Pa", 133.322, 0.0), "atm": ("Pa", 101325.0, 0.0)},
    "percent": {"%": ("%", 1.0, 0.0), "wt%": ("%", 1.0, 0.0), "at%": ("%", 1.0, 0.0), "wt.%": ("%", 1.0, 0.0)},
    "wavenumber": {"cm-1": ("cm-1", 1.0, 0.0), "cm⁻¹": ("cm-1", 1.0, 0.0), "1/cm": ("cm-1", 1.0, 0.0), "cm^-1": ("cm-1", 1.0, 0.0)},
    "rate": {"k/min": ("K/s", 1 / 60.0, 0.0), "k/s": ("K/s", 1.0, 0.0), "°c/min": ("K/s", 1 / 60.0, 0.0)},
    "concentration": {"g/l": ("g/L", 1.0, 0.0), "mg/l": ("g/L", 1e-3, 0.0), "µg/l": ("g/L", 1e-6, 0.0), "ug/l": ("g/L", 1e-6, 0.0)},
}
_QUANTITY_WORDS = {"temperature": "temperature", "temp": "temperature", "t": "temperature", "heat flow": None,
                   "power": "power", "mass": "mass", "weight": "mass", "time": "time", "pressure": "pressure", "p": "pressure",
                   "wavenumber": "wavenumber", "raman shift": "wavenumber", "wavenumbers": "wavenumber",
                   "heating rate": "rate", "rate": "rate", "concentration": "concentration", "conc": "concentration"}
_UNIT_RE = re.compile(r"[\(\[]\s*([^\)\]]+?)\s*[\)\]]|/\s*([^\s/]+)\s*$")


def _split_header(name: str):
    """Returns (quantity or None, unit spelling or None) from a column header such as 'Temperature (°C)' or 'Power / mW'."""
    text = name.strip()
    m = _UNIT_RE.search(text)
    unit = None
    if m:
        unit = (m.group(1) or m.group(2) or "").strip()
    head = text[: m.start()].strip().lower() if m else text.lower()
    head = re.sub(r"[^a-z ]", " ", head).strip()
    quantity = None
    for word, q in _QUANTITY_WORDS.items():
        if q and (head == word or head.endswith(" " + word) or head.startswith(word + " ") or head == word):
            quantity = q
            break
    if unit is None and quantity is None:
        return None, None
    return quantity, (unit.lower().replace(" ", "") if unit else None)


def detect_unit(column_name: str):
    """Returns dict(quantity, unit, base, factor, offset) or dict(quantity=None, status='unknown unit')."""
    quantity, unit = _split_header(column_name)
    if quantity is None and unit is not None:
        for q, table in _UNITS.items():
            if unit in table:
                quantity = q
                break
    if quantity is None or unit is None or unit not in _UNITS.get(quantity, {}):
        return {"quantity": quantity, "unit": unit, "status": "unknown unit"}
    base, factor, offset = _UNITS[quantity][unit]
    return {"quantity": quantity, "unit": unit, "base": base, "factor": factor, "offset": offset, "status": "converted"}


def convert_column(values, column_name: str):
    info = detect_unit(column_name)
    arr = np.asarray(values, float)
    if info.get("status") == "converted":
        arr = arr * info["factor"] + info["offset"]
    return arr, info


def _sniff_delimiter(line: str) -> str:
    for d in ("\t", ";", ","):
        if d in line:
            return d
    return None


def _split(line: str, delim):
    return [c.strip() for c in (line.split(delim) if delim else line.split())]


def _to_float(text: str, decimal: str = ".") -> float:
    """decimal=',' is used with semicolon-separated exports (European style: 298,15)."""
    return float(text.replace(",", ".") if decimal == "," else text)


def _is_number(text: str, decimal: str = ".") -> bool:
    try:
        _to_float(text, decimal)
        return True
    except ValueError:
        return False


def read_text(text: str):
    """Parses an instrument text export. Returns dict(columns, data (dict of np arrays, converted), units (list), skipped_lines)."""
    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    header_idx = None
    delim = None
    for i in range(len(lines) - 1):
        cand = _sniff_delimiter(lines[i])
        cells = _split(lines[i], cand)
        nxt = _split(lines[i + 1], cand)
        dec = "," if cand == ";" else "."
        if len(cells) >= 2 and len(nxt) >= 2 and sum(_is_number(c, dec) for c in nxt) >= 2 and not all(_is_number(c, dec) for c in cells):
            header_idx, delim = i, cand
            break
    if header_idx is None:
        raise ValueError("no table found: need a header line followed by numeric rows")
    names = _split(lines[header_idx], delim)
    dec = "," if delim == ";" else "."
    rows = []
    for ln in lines[header_idx + 1:]:
        cells = _split(ln, delim)
        if len(cells) != len(names) or not all(_is_number(c, dec) for c in cells):
            continue
        rows.append([_to_float(c, dec) for c in cells])
    if not rows:
        raise ValueError("header found, but no numeric rows")
    table = np.array(rows, float)
    data, units = {}, []
    for j, name in enumerate(names):
        col, info = convert_column(table[:, j], name)
        data[name] = col
        units.append({"column": name, **{k: info[k] for k in ("quantity", "unit", "base", "status") if k in info}})
    return {"columns": names, "data": data, "units": units, "skipped_lines": header_idx, "delimiter": delim or "whitespace"}


def read_file(path: str):
    with open(path, "rb") as fh:
        raw = fh.read()
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return read_text(raw.decode(enc))
        except UnicodeDecodeError:
            continue
    raise ValueError("could not decode file")


_META = re.compile(r"^\s*([A-Za-z][A-Za-z0-9 _./()-]{0,60}?)\s*[:=]\s*(.+?)\s*$")


def extract_metadata(text: str) -> dict:
    """Key: value lines above the table (instrument, operator, sample, heating rate, date ...). Keys are kept as written.
    Table lines (tab or semicolon separated) are ignored. Values stay text; use convert_column for numbers with units."""
    out = {}
    for line in text.splitlines():
        if "\t" in line or ";" in line:
            continue
        m = _META.match(line)
        if m:
            out.setdefault(m.group(1).strip(), m.group(2).strip())
    return out
