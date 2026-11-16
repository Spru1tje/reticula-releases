"""Local SQLite sample database (Documents/Reticula/reticula.db): sample IDs from the synthesis pipeline's master
overview (protocol docx/PDF, family, parent, source DOI, paper matrix), results from the sample store, and data files
linked by sample ID in their names (convention <Sample ID>-B<batch>_…)."""
from __future__ import annotations
import datetime, json, os, re, sqlite3

SCHEMA = """
CREATE TABLE IF NOT EXISTS samples (id TEXT PRIMARY KEY, name TEXT, paper TEXT, family TEXT, parent TEXT, basis TEXT,
    doi TEXT, protocol_docx TEXT, protocol_pdf TEXT, matrix TEXT, notes TEXT, updated TEXT);
CREATE TABLE IF NOT EXISTS results (sample_id TEXT, section TEXT, key TEXT, value TEXT, source TEXT, updated TEXT,
    PRIMARY KEY (sample_id, section, key));
CREATE TABLE IF NOT EXISTS files (sample_id TEXT, path TEXT, technique TEXT, batch TEXT, added TEXT, PRIMARY KEY (sample_id, path));
"""
TECH = [("PXRD", r"pxrd|xrd|\.xy$|\.raw$|\.brml$"), ("FTIR", r"ftir|atr|\.dpt$|\.spa$"), ("Raman", r"raman"),
        ("TGA", r"tga"), ("DSC", r"dsc"), ("N2", r"bet|n2|sorption|3flex|tristar|asap"), ("NMR", r"nmr|\.fid|1h"),
        ("XPS", r"xps|casa|\.vms$"), ("ICP", r"icp"), ("XRF", r"xrf"), ("SEM", r"sem|tem|\.tif+$|\.dm[34]$"),
        ("EDX", r"edx|eds"), ("Kinetics", r"kinet|hydrol|glygly"), ("Zeta/DLS", r"zeta|dls")]


def connect(path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    con = sqlite3.connect(path); con.executescript(SCHEMA); return con


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower().replace("@", "at").replace("+", "plus"))


def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def import_master(con, xlsx_path: str, root: str | None = None):
    """Reads 'All samples' of Master_sample_overview.xlsx. Protocol paths are made absolute (relative to the workbook,
    PaperN/Protocols/<file>), or to `root` = the user's Synthesis_Pipeline folder. Returns the number of samples."""
    import pandas as pd
    raw = pd.read_excel(xlsx_path, sheet_name=0, header=None)
    hr = next(i for i in range(min(10, len(raw))) if "Sample ID" in [str(v).strip() for v in raw.iloc[i].values])
    df = raw.iloc[hr + 1:].copy(); df.columns = [str(c).strip() for c in raw.iloc[hr].values]
    df = df[df["Sample ID"].notna()]
    root = root or os.path.dirname(os.path.abspath(xlsx_path)); n = 0
    for _, r in df.iterrows():
        paper = str(r.get("Paper", "")).strip(); pn = re.sub(r"\D", "", paper)
        pdir = os.path.join(root, f"Paper{pn}", "Protocols") if pn else root
        g = lambda c: (None if c not in r or str(r[c]) in ("nan", "-", "") else str(r[c]).strip())
        con.execute("INSERT OR REPLACE INTO samples VALUES (?,?,?,?,?,?,?,?,?,?,COALESCE((SELECT notes FROM samples WHERE id=?),''),?)",
                    (str(r["Sample ID"]).strip(), g("Name"), paper, g("Protocol family"), g("Parent / starting material"), g("Protocol basis"),
                     g("Main source (DOI)"), os.path.join(pdir, g("Protocol (docx)")) if g("Protocol (docx)") else None,
                     os.path.join(pdir, g("Protocol (PDF)")) if g("Protocol (PDF)") else None,
                     os.path.join(root, f"Paper{pn}", g("Paper matrix")) if g("Paper matrix") and pn else None,
                     str(r["Sample ID"]).strip(), _now()))
        n += 1
    con.commit(); return n


def match_id(con, text: str):
    """Longest sample ID contained in `text` (normalised), e.g. 'Z801-F25-B1_PXRD.xy' → 'Z801-F25'; also matches names."""
    t = norm(text); best = None
    for sid, name in con.execute("SELECT id, name FROM samples"):
        for cand in (sid, name):
            c = norm(cand or "")
            if len(c) >= 3 and c in t and (best is None or len(c) > best[1]):
                best = (sid, len(c))
    return best[0] if best else None


def technique(path: str):
    p = path.lower()
    for name, rx in TECH:
        if re.search(rx, os.path.basename(p)) or re.search(rx, os.path.basename(os.path.dirname(p))):
            return name
    return None


def link_folder(con, folder: str, exts=None):
    """Walks a data folder and links every file whose name contains a sample ID. Returns (linked, unmatched)."""
    linked, unmatched = 0, []
    for dp, _, fs in os.walk(folder):
        for f in fs:
            if f.startswith(".") or (exts and not f.lower().endswith(tuple(exts))):
                continue
            p = os.path.join(dp, f); sid = match_id(con, f) or match_id(con, os.path.basename(dp))
            if not sid:
                unmatched.append(p); continue
            b = re.search(r"-B(\d+)", f, re.I)
            con.execute("INSERT OR REPLACE INTO files VALUES (?,?,?,?,?)", (sid, p, technique(p), b.group(1) if b else None, _now()))
            linked += 1
    con.commit(); return linked, unmatched


def sync_store(con, store: dict):
    """Copies the session sample store (results per technique) into the results table. Store names are matched to
    sample IDs; unknown names become new sample rows. Returns the number of values written."""
    n = 0
    for sid_store, rec in (store or {}).items():
        name = re.sub(r"\$|\\mathrm|[{}]|_(?=\d)", "", str(rec.get("name", sid_store)))     # ZrO$_2$ -> ZrO2
        sid = match_id(con, name) or name
        if not con.execute("SELECT 1 FROM samples WHERE id=?", (sid,)).fetchone():
            con.execute("INSERT INTO samples (id, name, notes, updated) VALUES (?,?,?,?)", (sid, name, "from the sample store", _now()))
        for sec, d in rec.get("sections", {}).items():
            for k, v in (d.get("values") or {}).items():
                con.execute("INSERT OR REPLACE INTO results VALUES (?,?,?,?,?,?)",
                            (sid, sec, k, json.dumps(v, ensure_ascii=False), d.get("source", ""), d.get("updated", _now())))
                n += 1
    con.commit(); return n


def sample(con, sid: str):
    cols = [c[1] for c in con.execute("PRAGMA table_info(samples)")]
    r = con.execute("SELECT * FROM samples WHERE id=?", (sid,)).fetchone()
    if not r:
        return None
    out = dict(zip(cols, r))
    out["results"] = [dict(section=s, key=k, value=json.loads(v), source=src, updated=u)
                      for s, k, v, src, u in con.execute("SELECT section, key, value, source, updated FROM results WHERE sample_id=? ORDER BY section, key", (sid,))]
    out["files"] = [dict(path=p, technique=t, batch=b) for p, t, b in con.execute("SELECT path, technique, batch FROM files WHERE sample_id=? ORDER BY technique, path", (sid,))]
    out["children"] = [c for (c,) in con.execute("SELECT id FROM samples WHERE parent LIKE ?", (f"%{sid}%",))]
    return out


def overview(con):
    import pandas as pd
    return pd.read_sql_query("""SELECT s.id AS "Sample ID", s.name AS Name, s.paper AS Paper, s.family AS Family, s.parent AS Parent,
        (SELECT COUNT(*) FROM files f WHERE f.sample_id=s.id) AS "data files",
        (SELECT GROUP_CONCAT(DISTINCT r.section) FROM results r WHERE r.sample_id=s.id) AS results,
        CASE WHEN s.protocol_pdf IS NOT NULL THEN 'yes' ELSE '' END AS protocol
        FROM samples s ORDER BY s.paper IS NULL, s.paper, s.id""", con)


def set_notes(con, sid, notes):
    con.execute("UPDATE samples SET notes=?, updated=? WHERE id=?", (notes, _now(), sid)); con.commit()
