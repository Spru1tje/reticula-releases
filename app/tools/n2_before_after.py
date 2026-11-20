"""Before/after for the N2 reader fix: v1 reader+plot vs v2 on the same MicroActive export.
usage: python tools/n2_before_after.py [file.xls] [out_prefix]"""
import os, sys, importlib.util
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import matplotlib.pyplot as plt
import mfs
from mfs import readers as R, analysis as A, plots as P, style as S, composer as C


def _legacy(name):
    spec = importlib.util.spec_from_file_location(f"mfs.{name}", os.path.join(ROOT, "tools", "legacy", f"{name}.py"))
    m = importlib.util.module_from_spec(spec); m.__package__ = "mfs"; spec.loader.exec_module(m); return m


def main(path, prefix):
    name = os.path.basename(path)
    old = _legacy("readers_v1").read_sorption(open(path, "rb"), name)
    P1 = _legacy("plots_v1")
    new = R.read_sorption(open(path, "rb"), name)
    nm_old = old["meta"].get("sample", "sample"); nm = new["meta"].get("sample", "sample")
    out = {}
    # ---- before (v1)
    S.apply_style()
    fig, ax = plt.subplots(figsize=(3.25, 2.6))
    P1.draw_isotherms(ax, [(nm_old, old["iso"])])
    a = old["iso"][old["iso"].branch == "ads"]
    b = A.bet_rouquerol(a.p_rel.values, a.q.values)
    ax.set_title(f"v1 reader: {len(old['iso'])} points merged from every table\nBET {b['area']:.0f} m²/g (file: {old['meta'].get('bet_reported', '–')})",
                 fontsize=6.5)
    out["before"] = S.save_all(fig, prefix + "_before_v1", ("png",))[0]; plt.close(fig)
    # ---- after (v2): linear + log isotherm and DFT PSD
    iso = new["iso"]; a = iso[iso.branch == "ads"]
    b = A.bet_rouquerol(a.p_rel.values, a.q.values)
    fns = [lambda ax: P.draw_isotherms(ax, [(nm, iso)]), lambda ax: P.draw_isotherms(ax, [(nm, iso)], logx=True, legend=False)]
    if new["psd"] is not None:
        fns.append(lambda ax: P.draw_psd(ax, [(nm, new["psd"])], xlim=(0.3, 3)))
    else:   # no PSD table in the file (e.g. TriStar full report): show the BET plot of the app's own fit instead
        fns.append(lambda ax: P.draw_bet(ax, b, nm))
    fig, axs = C.compose(fns, ncols=len(fns), width=7.0, panel_aspect=0.85, wspace=0.42)
    m = new["meta"]
    rep = m.get("bet_reported")
    txt = (f"$S_\\mathrm{{BET}}$ = {b['area']:.0f} m$^2$ g$^{{-1}}$ (Rouquerol, $p/p_0$ {b['p_lo']:.2g}–{b['p_hi']:.2g}, "
           f"$C$ = {b['C']:.0f}, $R^2$ = {b['R2']:.5f}, {b['n']} pts{'' if b['all_criteria'] else ', criteria not all met'})")
    if rep:
        rr = m.get("bet_reported_range")
        txt += f"\\nfile (MicroActive): {float(rep):.1f} m$^2$ g$^{{-1}}$" + (f" at $p/p_0$ {rr[0]:.2g}–{rr[1]:.2g}" if rr else "")
    axs[-1].text(0.97 if new["psd"] is not None else 0.03, 0.96, txt.replace("\\n", "\n").replace(" (Rouquerol, ", "\n(Rouquerol, ").replace(": 335", ":\n335").replace(" at p", "\nat p"), transform=axs[-1].transAxes, ha="right" if new["psd"] is not None else "left", va="top", fontsize=5.5)
    out["after"] = S.save_all(fig, prefix + "_after_v2", ("png",))[0]; plt.close(fig)
    out.update(n_old=len(old["iso"]), n_new=len(iso), n_ads=int((iso.branch == "ads").sum()), n_des=int((iso.branch == "des").sum()),
               bet_new=b["area"], bet_all=b["all_criteria"], bet_range=(b["p_lo"], b["p_hi"]), meta=m, tables=new["tables"])
    return out


if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "example_data", "Sorption", "SampleA_3Flex_full_report_EXAMPLE.xls")
    pre = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "example_figures", "N2_fullreport_EXAMPLE")
    r = main(f, pre)
    print({k: v for k, v in r.items() if k != "tables"})
