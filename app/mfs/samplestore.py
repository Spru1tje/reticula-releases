"""Sample-ID store shared by all pages: every analysis page puts its per-sample results here (TGA, NMR digestion,
ICP/XRF, N2, kinetics, DSC, defects), the defect calculator and the kinetics TOF pull them back. Kept in the
Streamlit session and mirrored to a JSON file (Documents\\Reticula\\sample_store.json in the Windows app)."""
from __future__ import annotations
import datetime, json, os, re
import numpy as np

SECTIONS = ("tga", "nmr", "icp", "xrf", "xps", "sem", "n2", "pxrd", "kinetics", "dsc", "defects")


def norm_id(name: str) -> str:
    """'ZrO$_2$@Sulfo-MOF-801 (2)' -> 'zro2sulfomof801(2)': case, spaces, TeX and punctuation are ignored."""
    s = re.sub(r"\$|_\{?|\}|\\mathrm|\\", "", str(name or "")).lower()
    return re.sub(r"[^a-z0-9()@]+", "", s)


def _clean(v):
    if isinstance(v, (np.floating,)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, float):
        return None if not np.isfinite(v) else v
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, dict):
        return {str(k): _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    return v


def _store(state):
    if "sample_store" not in state:
        state["sample_store"] = {}
    return state["sample_store"]


def put(state, sample: str, section: str, values: dict, source: str = "", path: str | None = None):
    """Store/replace one section of one sample. Returns the sample id."""
    st = _store(state)
    sid = norm_id(sample)
    if not sid:
        return None
    rec = st.setdefault(sid, {"name": sample, "sections": {}})
    rec["name"] = rec.get("name") or sample
    rec["sections"][section] = {"values": _clean(dict(values)), "source": source,
                                "updated": datetime.datetime.now().isoformat(timespec="seconds")}
    if path:
        save(state, path)
    return sid


def get(state, sample: str, section: str | None = None):
    rec = _store(state).get(norm_id(sample))
    if rec is None:
        return None
    if section is None:
        return rec
    s = rec["sections"].get(section)
    return s["values"] if s else None


def ids(state):
    return {sid: rec["name"] for sid, rec in _store(state).items()}


def find(state, name: str):
    """Exact id match, else the stored sample whose id contains / is contained in the given one (longest overlap)."""
    sid = norm_id(name); st = _store(state)
    if sid in st:
        return sid
    cand = [k for k in st if k and (k in sid or sid in k)]
    return max(cand, key=len) if cand else None


def summary_rows(state):
    rows = []
    for sid, rec in _store(state).items():
        r = {"sample": rec["name"], "id": sid}
        for sec in SECTIONS:
            r[sec] = "✓" if sec in rec["sections"] else ""
        rows.append(r)
    return rows


def save(state, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(_store(state), fh, indent=1, ensure_ascii=False)
    os.replace(tmp, path)


def load(state, path, merge=True):
    if not os.path.exists(path):
        return 0
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    st = _store(state)
    for sid, rec in data.items():
        if merge and sid in st:
            for sec, v in rec.get("sections", {}).items():
                st[sid]["sections"].setdefault(sec, v)
        else:
            st[sid] = rec
    return len(data)


def clear(state, sample=None):
    st = _store(state)
    if sample is None:
        st.clear()
    else:
        st.pop(norm_id(sample), None)
