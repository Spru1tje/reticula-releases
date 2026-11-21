"""Headless smoke test of every Streamlit page with the bundled example data (streamlit.testing AppTest)."""
import os, sys, time
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)
from streamlit.testing.v1 import AppTest
PAGES = ["Home", "Reticon", "PXRD", "FTIR / Raman", "1H NMR (digest)", "N2 sorption (BET/PSD)", "TGA", "DSC", "Defect calculator",
         "Kinetics & adsorption", "ICP-OES / XRF", "SEM / TEM / EDX", "SDS-PAGE (Image Lab)", "Zeta / DLS", "XPS", "Figure composer", "Export guide"]


def run(p):
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=300)
    at.run()
    at.button(key="nv_" + p).click().run()
    return at


def composer_check():
    at = run("PXRD")
    for b in at.button:
        if b.key and b.key.startswith("add_"):
            b.click().run(); break
    at = run_after(at, "N2 sorption (BET/PSD)")
    for b in at.button:
        if b.key and b.key.startswith("add_"):
            b.click().run(); break
    at.button(key="nv_" + "Figure composer").click().run()
    errs = [str(e.value)[:400] for e in at.exception] + [str(e.value)[:400] for e in at.error]
    return len(at.session_state["composer"]), len(at.get("download_button")), errs


def run_after(at, p):
    at.button(key="nv_" + p).click().run(); return at


if __name__ == "__main__":
    bad = 0
    for p in PAGES:
        t = time.time(); at = run(p)
        errs = [str(e.value)[:400] for e in at.exception] + [str(e.value)[:400] for e in at.error]
        bad += bool(errs)
        print(f"{p:26s} errors={len(errs)} t={time.time()-t:5.1f}s {errs if errs else ''}")
    n, ni, e = composer_check(); print("composer panels", n, "errors", e); bad += bool(e)
    sys.exit(1 if bad else 0)

