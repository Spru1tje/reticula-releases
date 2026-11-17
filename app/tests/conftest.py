import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


import pytest

os.environ.setdefault("MFS_NO_UPDATE_CHECK", "1")     # also for the Streamlit servers the browser suites start
os.environ.setdefault("MFS_DEFAULT_LEVEL", "Developer")      # v2.15: also for the servers the browser suites start


@pytest.fixture(autouse=True)
def _nav_guard_off(monkeypatch):
    """AppTest suites navigate freely; the v2.9 unsaved-work dialog is tested explicitly (test_pagestate*.py)."""
    monkeypatch.setenv("MFS_NAV_GUARD", "0")


@pytest.fixture(autouse=True)
def _update_check_off(monkeypatch):
    """No 'check for updates?' banner / network in the suites; tests/test_v211_update.py switches it on explicitly."""
    monkeypatch.setenv("MFS_NO_UPDATE_CHECK", "1")


@pytest.fixture
def example_colour_rules():
    """The name -> colour rules of the bundled MOF-801 example set (not active by default since v2.13)."""
    from mfs import style as S
    old = list(S.SAMPLE_COLOURS); S.SAMPLE_COLOURS[:] = S.EXAMPLE_SAMPLE_COLOURS
    yield
    S.SAMPLE_COLOURS[:] = old


# ---- private fixtures -------------------------------------------------------------------------------------------
# Real measurements (TEM images, a TriStar report) are NOT part of the public repository. Tests that need them are
# skipped when the files are absent; drop the files into tests/data/ to run them locally.
PRIVATE_FIXTURES = ("TEM_MOF801_user5", "TEM_SiO2NH2_user6", "olympus_sis", "TriStar_full_report_REAL")
_PRIVATE_MODULES = {"test_real_tristar_n2": "MOF801SGA2_TriStar_full_report_REAL.xls", "test_real_olympus_tem": "olympus_sis"}


def _private_missing(token):
    import glob
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    return not glob.glob(os.path.join(d, "*" + token + "*"))


def pytest_collection_modifyitems(config, items):
    import inspect
    for item in items:
        mod = getattr(item.module, "__name__", "").split(".")[-1]
        try:
            src = inspect.getsource(item.function)
        except Exception:
            src = ""
        text = src + " " + item.name + " " + _PRIVATE_MODULES.get(mod, "")
        need = [t for t in PRIVATE_FIXTURES if t in text and _private_missing(t)]
        if need:
            item.add_marker(pytest.mark.skip(reason=f"private test data not in the public repository: {', '.join(need)}"))


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """A test that fails only because a private measurement file is not in this checkout is reported as skipped."""
    outcome = yield
    rep = outcome.get_result()
    if call.when == "call" and rep.failed and call.excinfo is not None and call.excinfo.errisinstance(FileNotFoundError):
        if any(t in str(call.excinfo.value) for t in PRIVATE_FIXTURES):
            rep.outcome = "skipped"
            rep.longrepr = (str(item.fspath), (item.location[1] or 0) + 1, "Skipped: private test data not in the public repository")


@pytest.fixture(autouse=True)
def _level_developer(monkeypatch):
    """v2.15: the app starts in 'Standard'; the older suites test every control, so they run at 'Developer'
    (tests/test_v215_b.py::test_levels_* use the real default)."""
    monkeypatch.setenv("MFS_DEFAULT_LEVEL", "Developer")
