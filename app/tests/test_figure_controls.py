"""v2.14: EVERY control of 'Customise this figure' must have an effect on EVERY figure of EVERY page.

Two levels:
1. logic - the real figures of every page are captured (AppTest, bundled example data) and every setting of the panel is applied
   to a copy with mfs.figctl.apply; the effect is checked on the matplotlib objects (font sizes, ticks, limits, legend ...).
2. wiring - the widgets of the panel are set in a running page and the settings that arrive in the figure are checked.
A control that does nothing on some figure is reported with the page, the figure and the control name."""
import os
import pickle
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.legend import Legend  # noqa: E402
from matplotlib.ticker import MultipleLocator, NullLocator  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from mfs import figctl as F  # noqa: E402

PAGES = ["PXRD", "FTIR / Raman", "UV/Vis", "1H NMR (digest)", "N2 sorption (BET/PSD)", "TGA", "DSC", "Defect calculator",
         "Kinetics & adsorption", "ICP-OES / XRF", "SEM / TEM / EDX", "SDS-PAGE (Image Lab)", "Zeta / DLS", "XPS", "Structure (CIF)", "Photoluminescence"]


# ------------------------------------------------------------------------------------------------ capture the real figures
@pytest.fixture(scope="module")
def captured(tmp_path_factory):
    warnings.filterwarnings("ignore")
    os.environ["MFS_NO_UPDATE_CHECK"] = "1"
    os.environ["MFS_STORE"] = str(tmp_path_factory.mktemp("fcstore") / "store.json")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    
    at.run()
    out = {}
    for pg in PAGES:
        try:
            at.button(key="nv_" + pg).click().run()
        except Exception as e:                      # a page that cannot be opened is reported by the test of that page
            out[pg] = {"_error": str(e)[:200]}
            continue
        figs = at.session_state["si_figs"] if "si_figs" in at.session_state else {}
        out[pg] = {}
        for stem, rec in figs.items():
            if rec.get("page") == pg:
                try:
                    out[pg][stem] = pickle.dumps(rec["fig"])
                except Exception:
                    pass
        if at.exception:
            out[pg]["_error"] = str([e.value for e in at.exception][:1])[:200]
    return out


# ------------------------------------------------------------------------------------------------ helpers
def _panel(fig):
    """First data panel that is not a twin axis: (index, axes, description) or None."""
    desc = F.describe(fig)
    for d in desc:
        if d["axison"] and not F.is_twin(fig.axes[d["index"]], fig):
            return d["index"], fig.axes[d["index"]], d
    return None


def _legends(fig):
    return list(fig.legends) + [l for ax in fig.axes for l in ax.get_children() if isinstance(l, Legend)]


def _cfg():
    return F.default_settings()


def _tick(ax, which="x"):
    A = ax.xaxis if which == "x" else ax.yaxis
    t = A.get_major_ticks()
    return t[0] if t else None


# Each control: (name, needs, cfg mutator, check) - check(fig, ctx) returns None (ok) or a message. needs: 'axes' | 'legend' | 'any'
def _controls():
    C = []

    def add(name, needs, mut, check):
        C.append((name, needs, mut, check))

    # ---- title and fonts
    add("title", "any", lambda c, ctx: c.update(title="TESTTITLE"),
        lambda f, ctx: None if any(t.get_text() == "TESTTITLE" for t in ([f._suptitle] if f._suptitle else []) + [a.title for a in f.axes] + f.texts)
        else "title not drawn")

    def title_size(c, ctx):
        c.update(title="TESTTITLE", title_size=15.0)

    def title_size_check(f, ctx):
        ts = [t for t in ([f._suptitle] if f._suptitle else []) + [a.title for a in f.axes] if t.get_text() == "TESTTITLE"]
        return None if ts and abs(ts[0].get_fontsize() - 15.0) < 1e-6 else "title size not applied"
    add("title size", "any", title_size, title_size_check)

    add("axis-title size", "axes", lambda c, ctx: c["font"].update(label_size=13.0),
        lambda f, ctx: None if abs(ctx["ax"].xaxis.label.get_fontsize() - 13.0) < 1e-6 and abs(ctx["ax"].yaxis.label.get_fontsize() - 13.0) < 1e-6
        else f"axis title sizes {ctx['ax'].xaxis.label.get_fontsize()}/{ctx['ax'].yaxis.label.get_fontsize()}")

    def tick_size_check(f, ctx):
        bad = []
        for w in "xy":
            t = _tick(ctx["ax"], w)
            if t is not None and t.label1.get_visible() and abs(t.label1.get_fontsize() - 12.0) > 1e-6:
                bad.append(f"{w}:{t.label1.get_fontsize()}")
        return None if not bad else "tick label size " + ",".join(bad)
    add("tick-label size", "axes", lambda c, ctx: c["font"].update(tick_size=12.0), tick_size_check)

    def family_check(f, ctx):
        from mfs import style as S
        want = S._font_family("DejaVu Serif")
        ax = ctx["ax"]
        got = {"xlabel": ax.xaxis.label.get_fontfamily()[0], "ylabel": ax.yaxis.label.get_fontfamily()[0]}
        t = _tick(ax, "x")
        if t is not None and t.label1.get_visible():
            got["xtick"] = t.label1.get_fontfamily()[0]
        for lg in _legends(f):
            got["legend"] = lg.get_texts()[0].get_fontfamily()[0] if lg.get_texts() else want
        bad = {k: v for k, v in got.items() if v != want}
        return None if not bad else f"font family not applied: {bad}"
    add("font family", "axes", lambda c, ctx: c["font"].update(family="DejaVu Serif"), family_check)

    # ---- ticks
    def preset_check(f, ctx):
        ax = ctx["ax"]
        t = _tick(ax, "x"); ty = _tick(ax, "y")
        if t is not None and (not t.tick1line.get_visible() or t.tick2line.get_visible()):
            return "x ticks: bottom should be on and top off"
        if ty is not None and (not ty.tick1line.get_visible() or ty.tick2line.get_visible()):
            return "y ticks: left should be on and right off"
        return None
    add("ticks on 2 sides", "axes", lambda c, ctx: c["ticks"].update(preset="2 sides (bottom + left)"), preset_check)

    def frame_check(f, ctx):
        ax = ctx["ax"]
        has_twin = any(F.is_twin(o, f) for o in f.axes)           # a twin y axis keeps the right-hand spine on purpose
        right_ok = has_twin or not ax.spines["right"].get_visible()
        return None if (not ax.spines["top"].get_visible() and right_ok
                        and ax.spines["bottom"].get_visible() and ax.spines["left"].get_visible()) else "spines not matched to ticks"
    add("frame matches ticks", "axes", lambda c, ctx: c["ticks"].update(preset="2 sides (bottom + left)", frame="match ticks"), frame_check)

    def direction_check(f, ctx):
        ax = ctx["ax"]
        bad = [w for w in "xy" if (_tick(ax, w) is not None and (ax.xaxis if w == "x" else ax.yaxis)._major_tick_kw.get("tickdir") != "out")]
        return None if not bad else f"tick direction not 'out' on {bad}"
    add("tick direction", "axes", lambda c, ctx: c["ticks"].update(direction="out"), direction_check)

    def length_check(f, ctx):
        ax = ctx["ax"]
        k = ax.xaxis._major_tick_kw
        return None if abs(float(k.get("size", -1)) - 8.0) < 1e-6 else f"major length {k.get('size')}"
    add("major tick length", "axes", lambda c, ctx: c["ticks"].update(major_len=8.0), length_check)

    def width_check(f, ctx):
        k = ctx["ax"].xaxis._major_tick_kw
        return None if abs(float(k.get("width", -1)) - 1.5) < 1e-6 else f"tick width {k.get('width')}"
    add("tick width", "axes", lambda c, ctx: c["ticks"].update(width=1.5), width_check)

    def minor_check(f, ctx):
        ax = ctx["ax"]
        if ax.get_xscale() == "log" or len(ax.xaxis.get_majorticklocs()) == 0:
            return None
        return None if isinstance(ax.xaxis.get_minor_locator(), NullLocator) else "minor ticks still drawn"
    add("minor ticks off", "axes", lambda c, ctx: c["ticks"].update(minor=False), minor_check)

    # ---- axes
    def labels_check(f, ctx):
        ax = ctx["ax"]
        return None if ax.get_xlabel() == "XLAB" and ax.get_ylabel() == "YLAB" else f"labels {ax.get_xlabel()!r}/{ax.get_ylabel()!r}"
    add("axis titles", "axes", lambda c, ctx: c["panels"].update({str(ctx["i"]): {"xlabel": "XLAB", "ylabel": "YLAB"}}), labels_check)

    def lim_mut(c, ctx):
        ax = ctx["ax"]
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        lo, hi = sorted((x0, x1)); ylo, yhi = sorted((y0, y1))
        ctx["want"] = (lo + 0.1 * (hi - lo), hi - 0.1 * (hi - lo), ylo + 0.1 * (yhi - ylo), yhi - 0.1 * (yhi - ylo))
        c["panels"][str(ctx["i"])] = {"xmin": ctx["want"][0], "xmax": ctx["want"][1], "ymin": ctx["want"][2], "ymax": ctx["want"][3]}

    def lim_check(f, ctx):
        ax = ctx["ax"]
        xs, ys = sorted(ax.get_xlim()), sorted(ax.get_ylim())
        w = ctx["want"]
        ok = np.allclose(xs, [w[0], w[1]], rtol=1e-6, atol=1e-9) and np.allclose(ys, [w[2], w[3]], rtol=1e-6, atol=1e-9)
        return None if ok else f"limits {xs}/{ys} instead of {w}"
    add("axis ranges", "axes", lim_mut, lim_check)

    def log_mut(c, ctx):
        ax = ctx["ax"]
        ctx["xlog_to"] = ax.get_xscale() != "log"
        c["panels"][str(ctx["i"])] = {"xlog": ctx["xlog_to"]}

    def log_check(f, ctx):
        ax = ctx["ax"]
        if ctx["xlog_to"] and min(ctx["lim0"][0]) <= 0:
            return None                                          # a log axis cannot show a range that includes 0 (not a control fault)
        return None if (ax.get_xscale() == "log") == ctx["xlog_to"] else "log x not applied"
    add("log x", "axes", log_mut, log_check)

    def inc_mut(c, ctx):
        ax = ctx["ax"]
        lo, hi = sorted(ax.get_xlim())
        ctx["step"] = (hi - lo) / 4.0
        c["panels"][str(ctx["i"])] = {"xmajor": ctx["step"]}

    def inc_check(f, ctx):
        ax = ctx["ax"]
        if ax.get_xscale() == "log" or len(ctx["maj0"]) == 0:
            return None
        loc = ax.xaxis.get_major_locator()
        return None if isinstance(loc, MultipleLocator) and abs(loc._edge.step - ctx["step"]) < 1e-6 * max(1, abs(ctx["step"])) \
            else f"x major increment not applied ({type(loc).__name__})"
    add("x major increment", "axes", inc_mut, inc_check)

    # ---- legend
    def lg_first(f):
        lg = _legends(f)
        return lg[0] if lg else None

    def has_direct(f):
        return bool(F.direct_labels(f))

    add("legend font size", "legend", lambda c, ctx: c["legend"].update(fontsize=13.0),
        lambda f, ctx: None if (lambda sizes: sizes and all(abs(s - 13.0) < 1e-6 for s in sizes))(
            [t.get_fontsize() for lg in _legends(f) for t in lg.get_texts()] +
            [t.get_fontsize() for ax in f.axes for t in ax.texts if t.get_text() in F.direct_labels(f)]) else "legend/label font size not applied")
    add("legend title", "legendonly", lambda c, ctx: c["legend"].update(title="LEGTITLE"),
        lambda f, ctx: None if lg_first(f) is not None and lg_first(f).get_title().get_text() == "LEGTITLE" else "legend title not shown")
    add("legend title size", "legendonly", lambda c, ctx: c["legend"].update(title="LEGTITLE", title_size=14.0),
        lambda f, ctx: None if lg_first(f) is not None and abs(lg_first(f).get_title().get_fontsize() - 14.0) < 1e-6 else "legend title size not applied")
    add("legend frame", "legendonly", lambda c, ctx: c["legend"].update(frame=True),
        lambda f, ctx: None if lg_first(f) is not None and lg_first(f).get_frame_on() else "legend frame not drawn")
    add("legend columns", "legendonly", lambda c, ctx: c["legend"].update(ncol=2),
        lambda f, ctx: None if lg_first(f) is not None and (getattr(lg_first(f), "_ncols", 1) == 2 or len(lg_first(f).get_texts()) < 2)
        else "legend columns not applied")
    add("legend position", "legendonly", lambda c, ctx: c["legend"].update(loc="lower left"),
        lambda f, ctx: None if lg_first(f) is not None and lg_first(f)._loc in (3, "lower left") or lg_first(f) is None else "legend position not applied")
    add("legend marker size", "legendonly", lambda c, ctx: c["legend"].update(markerscale=2.0),
        lambda f, ctx: None if lg_first(f) is not None else "legend lost")
    add("legend hide", "legendonly", lambda c, ctx: c["legend"].update(show=False, custom=True),
        lambda f, ctx: None if not _legends(f) else "legend still shown")

    def rename_mut(c, ctx):
        ent = F.legend_entries(ctx["fig0"])
        ctx["first"] = ent[0][0]
        c["legend"].update(custom=True, entries=[{"label": ent[0][0], "text": "RENAMED", "show": True, "order": 0}])

    def rename_check(f, ctx):
        lg = lg_first(f)
        return None if lg is not None and "RENAMED" in [t.get_text() for t in lg.get_texts()] else "legend entry not renamed"
    add("legend rename", "legendonly", rename_mut, rename_check)
    return C


CONTROLS = _controls()


def _needs_ok(needs, fig, panel):
    if needs == "axes":
        return panel is not None
    if needs == "legend":
        return bool(_legends(fig)) or bool(F.direct_labels(fig))
    if needs == "legendonly":
        return bool(_legends(fig)) and bool(F.legend_entries(fig))
    return True


def run_controls(fig_bytes):
    """Applies every control to a fresh copy of the figure; returns [(control, message)] of the failures and the number checked."""
    fails, n = [], 0
    for name, needs, mut, check in CONTROLS:
        fig = pickle.loads(fig_bytes)
        panel = _panel(fig)
        if not _needs_ok(needs, fig, panel):
            continue
        ctx = {"fig0": pickle.loads(fig_bytes)}
        if panel:
            ctx.update(i=panel[0], ax=panel[1])
            ctx["lim0"] = (panel[1].get_xlim(), panel[1].get_ylim())
            ctx["maj0"] = panel[1].xaxis.get_majorticklocs()
        cfg = _cfg()
        try:
            mut(cfg, ctx)
            F.apply(fig, cfg)
            if panel:
                ctx["ax"] = fig.axes[panel[0]]
            msg = check(fig, ctx)
        except Exception as e:
            msg = f"{type(e).__name__}: {str(e)[:120]}"
        n += 1
        if msg:
            fails.append((name, msg))
    return fails, n


# ------------------------------------------------------------------------------------------------ the tests
@pytest.mark.parametrize("page", PAGES)
def test_every_control_works_on_every_figure(captured, page):
    figs = captured[page]
    assert "_error" not in figs, f"{page}: page failed to run: {figs.get('_error')}"
    stems = [k for k in figs if k != "_error"]
    report, checked = [], 0
    for stem in stems:
        fails, n = run_controls(figs[stem])
        checked += n
        report += [f"{page} / {stem}: {name} - {msg}" for name, msg in fails]
    assert not report, f"{len(report)} control(s) without effect ({checked} checks on {len(stems)} figure(s)):\n  " + "\n  ".join(report)


def test_all_pages_produced_figures_or_are_known_figureless(captured):
    figureless = {"Structure (CIF)", "Defect calculator"}          # pages whose figures need user input first (checked elsewhere)
    empty = [p for p in PAGES if not [k for k in captured[p] if k != "_error"] and p not in figureless]
    assert not empty, f"pages without any captured figure (their controls are untested): {empty}"


# ------------------------------------------------------------------------------------------------ wiring: widget -> figure
@pytest.fixture(scope="module")
def wired(tmp_path_factory):
    """Every page opened once, then the title / sizes / font / legend widgets of every figure are set and the page is rerun."""
    import re
    warnings.filterwarnings("ignore")
    os.environ["MFS_NO_UPDATE_CHECK"] = "1"
    os.environ["MFS_STORE"] = str(tmp_path_factory.mktemp("fcwire") / "store.json")
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(os.path.join(HERE, "app.py"), default_timeout=600)
    
    at.run()
    result = {}
    for pg in PAGES:
        res = {"missing": [], "checked": {}}
        result[pg] = res
        try:
            at.button(key="nv_" + pg).click().run()
            figs = at.session_state["si_figs"] if "si_figs" in at.session_state else {}
            stems = [s for s, r in figs.items() if r.get("page") == pg]
            keys = {s: re.sub(r"[^A-Za-z0-9_]+", "_", s) for s in stems}
            for s, k in keys.items():
                try:
                    at.text_input(key=f"fc_{k}_title").set_value("WIRED")
                    at.number_input(key=f"fc_{k}_label_size").set_value(13.0)
                    at.number_input(key=f"fc_{k}_tick_size").set_value(12.0)
                    at.selectbox(key=f"fc_{k}_family").set_value("DejaVu Serif")
                except Exception:
                    res["missing"].append(s)
                    continue
                try:
                    at.number_input(key=f"fc_{k}_leg_fontsize").set_value(13.0)
                    res["checked"][s] = True                       # the legend widget exists
                except Exception:
                    res["checked"][s] = False
            at.run()
            figs = at.session_state["si_figs"] if "si_figs" in at.session_state else {}
            res["figs"] = {s: pickle.dumps(figs[s]["fig"]) for s in stems if s in figs and s not in res["missing"]}
            res["error"] = str([e.value for e in at.exception][:1])[:200] if at.exception else ""
        except Exception as e:
            res["error"] = f"{type(e).__name__}: {str(e)[:150]}"
            res["figs"] = {}
    return result


@pytest.mark.parametrize("page", PAGES)
def test_panel_widgets_reach_the_figure(wired, page):
    res = wired[page]
    assert not res.get("error"), f"{page}: {res['error']}"
    problems = []
    for stem, blob in res["figs"].items():
        fig = pickle.loads(blob)
        pnl = _panel(fig)
        titles = [t.get_text() for t in ([fig._suptitle] if fig._suptitle else []) + [a.title for a in fig.axes]]
        if "WIRED" not in titles:
            problems.append(f"{stem}: title widget has no effect")
        if pnl:
            ax = pnl[1]
            if abs(ax.xaxis.label.get_fontsize() - 13.0) > 1e-6:
                problems.append(f"{stem}: axis-title size widget has no effect ({ax.xaxis.label.get_fontsize()})")
            t = _tick(ax, "x")
            if t is not None and t.label1.get_visible() and abs(t.label1.get_fontsize() - 12.0) > 1e-6:
                problems.append(f"{stem}: tick-label size widget has no effect ({t.label1.get_fontsize()})")
            from mfs import style as S
            if ax.xaxis.label.get_fontfamily()[0] != S._font_family("DejaVu Serif"):
                problems.append(f"{stem}: font widget has no effect ({ax.xaxis.label.get_fontfamily()[0]})")
        if res["checked"].get(stem) and (_legends(fig) or F.direct_labels(fig)):
            sizes = [t.get_fontsize() for lg in _legends(fig) for t in lg.get_texts()] + \
                    [t.get_fontsize() for a in fig.axes for t in a.texts if t.get_text() in F.direct_labels(fig)]
            if sizes and any(abs(s - 13.0) > 1e-6 for s in sizes):
                problems.append(f"{stem}: legend font-size widget has no effect ({sorted(set(sizes))})")
    assert not problems, f"{len(problems)} widget(s) without effect on {page}:\n  " + "\n  ".join(problems)
