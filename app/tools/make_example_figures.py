"""Run every module headlessly on the EXAMPLE data and write sample figures (PNG 600 dpi + SVG + PDF + TIFF)."""
from __future__ import annotations
import os, sys, glob, json, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mfs import readers as R, analysis as A, plots as P, style as S, defects as D, composer as C, report
from mfs.chem import SPECIES, oxide_mass, AW

ED = os.path.join(ROOT, "example_data")
import argparse
_ap = argparse.ArgumentParser(description="Render every example figure headlessly")
_ap.add_argument("--out", default=os.path.join(ROOT, "example_figures"), help="output folder")
_ap.add_argument("--png-only", action="store_true", help="write PNG only (fast test mode)")
ARGS = _ap.parse_args()
OUT = ARGS.out
os.makedirs(OUT, exist_ok=True)
NAMES = {"SampleA": "Sample A", "SampleB": "Sample B", "ZrO2": "ZrO$_2$", "PM": "PM",
         "SampleC": "Sample C"}
ORDER = ["SampleA", "SampleB", "ZrO2", "PM", "SampleC"]
written = []


def save(fig, stem, formats=("png", "svg", "pdf", "tiff")):
    if ARGS.png_only:
        formats = ("png",)
    p = S.save_all(fig, os.path.join(OUT, stem), formats)
    written.append(p[0]); plt.close(fig)
    return p


def label_name(stem):
    return NAMES.get(stem, stem)


# ---------------------------------------------------------------- data loading
def load_pxrd():
    tr = []
    for k in ORDER:
        df = R.read_pxrd(os.path.join(ED, "PXRD", f"{k}_STOE.xy"))
        y = A.process(df.x, df.y, baseline="ALS", smooth_method="none", norm="max", lam=1e6, p=0.005)
        tr.append((label_name(k), df.x.values, y))
    sim, sticks, f = A.simulate_pxrd(os.path.join(ROOT, "references", "MOF-801_Zr_COD4121459.cif"), (3, 50), fwhm=0.12,
                                     substitute={"Zr": "Hf"})
    return tr, ("Simulated Sample A", sim.x.values, sim.y.values), sticks


def draw_pxrd(ax, **kw):
    tr, ref, st = load_pxrd()
    P.draw_stacked(ax, tr, "pxrd", offset=1.05, xlim=(4, 50), reference=ref, **kw)


def draw_ftir(ax, **kw):
    tr = []
    for k in ["Sulfosuccinic_acid", "ZrO2", "SampleA", "SampleB", "PM", "SampleC"]:
        df = R.read_ftir(os.path.join(ED, "FTIR", f"{k}_ATR.dpt"))
        y = A.process(df.x, df.y, norm="min-max")
        tr.append(("2-Sulfosuccinic acid" if k.startswith("Sulfos") else label_name(k), df.x.values, y))
    P.draw_stacked(ax, tr, "ftir_T", offset=1.0, xlim=(400, 2000), highlight=[(1700, 1740, "C=O"), (1190, 1230, "S=O"), (1030, 1055, "")], **kw)


def draw_nmr(ax, **kw):
    tr = []
    for k in ["SampleA", "SampleB", "SampleC"]:
        df = R.read_nmr(os.path.join(ED, "NMR", f"{k}_digest_totxt.txt"))
        m = (df.x > 1.5) & (df.x < 9.0) & ~((df.x > 4.55) & (df.x < 5.05))
        y = df.y.where(m, np.nan)
        y = y / np.nanmax(y)
        tr.append((label_name(k), df.x.values, y.values))
    P.draw_stacked(ax, tr, "nmr", offset=1.1, xlim=(9.0, 1.5),
                   annotations=[(8.44, "HCOO$^-$"), (6.52, "fumarate"), (3.62, "SS CH"), (2.68, "SS CH$_2$")], **kw)


def load_sorption():
    out = {}
    for k in ORDER:
        out[k] = R.read_sorption(os.path.join(ED, "Sorption", f"{k}_3Flex_MicroActive.xlsx"))
    return out


def draw_iso(ax, logx=False, **kw):
    s = load_sorption()
    P.draw_isotherms(ax, [(label_name(k), v["iso"]) for k, v in s.items()], logx=logx, **kw)


def draw_psd(ax, **kw):
    s = load_sorption()
    P.draw_psd(ax, [(label_name(k), v["psd"]) for k, v in s.items() if k != "ZrO2"], xlim=(0.3, 2.0), legend=False, **kw)


def load_tga():
    out = []
    for f, nm in [("SampleA_Netzsch.txt", "Sample A"), ("SampleB_TA_UniversalAnalysis.txt", "Sample B"),
                  ("SampleC_Mettler_STARe.txt", "Sample C")]:
        t = R.read_tga(os.path.join(ED, "TGA", f))
        out.append((nm, A.tga_analysis(t["T"], t["w"])))
    return out


def draw_tga(ax, **kw):
    P.draw_tga(ax, load_tga(), legend=True, **kw)


def load_kinetics():
    df = R.read_csv_table(os.path.join(ED, "Kinetics", "GlyGly_hydrolysis_60C_pD7.4.csv"), (r"sample", r"time"))
    sets, fits = [], {}
    for s, g in df.groupby("sample", sort=False):
        agg = g.groupby("time_h")["GlyGly_mM"].agg(["mean", "std"]).reset_index()
        fit = A.first_order_fit(g["time_h"], g["GlyGly_mM"])
        fits[s] = fit
        sets.append((s.replace("ZrO2", "ZrO$_2$"), agg["time_h"].values, agg["mean"].values, agg["std"].values, fit))
    return sets, fits


def draw_kin(ax, **kw):
    sets, _ = load_kinetics()
    P.draw_kinetics(ax, [s for s in sets if "Blank" not in s[0]] + [s for s in sets if "Blank" in s[0]], ylabel=r"[Gly-Gly]/[Gly-Gly]$_0$",
                    annotate_k=False, legend=kw.pop("legend", True), **kw)


def draw_kbars(ax, **kw):
    sets, fits = load_kinetics()
    names = [s[0] for s in sets]
    P.draw_k_bars(ax, names, [s[4]["k"] for s in sets], [s[4]["k_err"] for s in sets])


def main():
    S.apply_style()
    # 1 PXRD
    fig, ax = S.new_figure(S.SINGLE, 1.05); draw_pxrd(ax); save(fig, "01_PXRD_stacked")
    # 2 FTIR
    fig, ax = S.new_figure(S.SINGLE, 1.15); draw_ftir(ax); save(fig, "02_FTIR_stacked")
    # 2b Raman
    tr = []
    for k in ["SampleA", "SampleB", "SampleC"]:
        df = R.read_ftir(os.path.join(ED, "Raman", f"{k}_RAMII.dpt"))
        tr.append((label_name(k), df.x.values, A.process(df.x, df.y, baseline="ALS", norm="max", lam=1e6)))
    fig, ax = S.new_figure(S.SINGLE, 0.85); P.draw_stacked(ax, tr, "raman", xlim=(200, 1800), invert=False); save(fig, "02b_Raman_stacked")
    # 3 NMR
    fig, ax = S.new_figure(S.SINGLE, 0.9); draw_nmr(ax); save(fig, "03_NMR_digest_stacked")
    # 3b native TopSpin pdata
    df = R.read_nmr(os.path.join(ED, "NMR", "SampleA_digest_TopSpin_expno10.zip"))
    fig, ax = S.new_figure(S.SINGLE, 0.6)
    P.draw_stacked(ax, [("Sample A digest (TopSpin pdata)", df.x.values, (df.y / df.y.max()).values)], "nmr", xlim=(9.5, -0.5))
    save(fig, "03b_NMR_TopSpin_native")
    # 4 Sorption (a) isotherm (b) PSD
    fig, axs = C.compose([draw_iso, draw_psd], ncols=2, width=S.DOUBLE, panel_aspect=0.8)
    save(fig, "04_N2_isotherms_PSD")
    # 5 BET Rouquerol
    s = load_sorption()["SampleA"]["iso"]; a = s[s.branch == "ads"]
    bet = A.bet_rouquerol(a.p_rel.values, a.q.values)
    fig, axs = C.compose([lambda ax: P.draw_rouquerol(ax, bet), lambda ax: P.draw_bet(ax, bet, "Sample A")], ncols=2, width=S.DOUBLE * 0.8, panel_aspect=0.8, wspace=0.5)
    save(fig, "05_BET_Rouquerol")
    # 6 TGA + DTG
    fig = P.fig_tga(load_tga(), S.SINGLE, "panel"); save(fig, "06_TGA_DTG")
    # 7 kinetics
    fig, axs = C.compose([draw_kin, draw_kbars], ncols=2, width=S.DOUBLE, panel_aspect=0.78)
    save(fig, "07_GlyGly_kinetics")
    # 8 Langmuir
    ad = R.read_csv_table(os.path.join(ED, "Kinetics", "GlyGly_adsorption_RT.csv"), (r"sample", r"ce"))
    sets = []
    for nm, g in ad.groupby("sample", sort=False):
        sets.append((nm, g.Ce_mM.values, g.qe_umol_per_umol.values, A.langmuir_fit(g.Ce_mM, g.qe_umol_per_umol)))
    fig, ax = S.new_figure(S.SINGLE, 0.8); P.draw_langmuir(ax, sets); save(fig, "08_Langmuir")
    # 9 ICP bars (wt%)
    icp = R.read_icp(os.path.join(ED, "ICP", "ICP-OES_Syngistix_export.csv"))
    icp["wt"] = icp.conc * 0.050 / 10.0 * 100; icp["wt_sd"] = icp.sd * 0.050 / 10.0 * 100
    icp["sample"] = icp["sample"].str.replace("ZrO2", "ZrO$_2$")
    fig, axs = C.compose([lambda ax: P.draw_grouped_bars(ax, icp[icp.element == el], "wt", "wt_sd", ylabel=f"{el} (wt%)", legend=False) for el in ["Hf", "Zr", "S"]],
                         ncols=3, width=S.DOUBLE, panel_aspect=0.95, wspace=0.45)
    save(fig, "09_ICP_bars")
    # 10 SEM panel
    ims = []
    for f, t in [("SampleA_lowmag_SEM.tif", "Sample A"), ("SampleA_SEM.tif", "Sample A"), ("SampleB_SEM.tif", "Sample B"),
                 ("SampleC_SEM.tif", "Sample C")]:
        m = R.read_sem(os.path.join(ED, "SEM", f))
        ims.append({"image": m["image"], "pixel_nm": m["pixel_nm"], "crop_row": m["databar_row"], "title": t})
    fig = P.fig_sem_panel(ims, 2, S.DOUBLE * 0.75); save(fig, "10_SEM_panel")
    # 11 SDS-PAGE
    b = R.read_imagelab(os.path.join(ED, "SDS-PAGE", "ImageLab_LanesAndBands_SampleC_HEWL_60C.csv"))
    b["band_pct"] = 100 * b.band_pct / b.groupby("lane").band_pct.transform("sum")
    fig, ax = S.new_figure(S.SINGLE, 0.8); P.draw_bands(ax, b); save(fig, "11_SDS-PAGE_bands")
    # 12 Zeta + DLS
    z = R.read_zetasizer(os.path.join(ED, "Zeta-DLS", "Zetasizer_export_pD7.4_D2O.csv"))
    zs = z["summary"].groupby("sample", sort=False).agg(zeta=("zeta", "mean"), sd=("zeta", "std")).reset_index()
    zs["sample"] = zs["sample"].str.replace("ZrO2", "ZrO$_2$")
    def draw_zeta(ax):
        cols = S.colours_for(zs["sample"])
        xs = np.arange(len(zs))
        ax.bar(xs, zs.zeta, 0.65, yerr=zs.sd, color=cols, edgecolor="k", lw=0.6, error_kw=dict(lw=0.7, capsize=2.5))
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(xs, zs["sample"], rotation=30, ha="right"); ax.tick_params(axis="x", which="minor", bottom=False, top=False)
        ax.set_ylabel(r"$\zeta$-potential (mV)")
    def draw_dls(ax):
        dd = z["dist"].groupby(["sample", "size_nm"], sort=False).intensity.mean().reset_index()
        tr = [(s.replace("ZrO2", "ZrO$_2$"), g.size_nm.values, g.intensity.values) for s, g in dd.groupby("sample", sort=False)]
        P.draw_overlay(ax, tr, "Hydrodynamic diameter (nm)", "Intensity (%)", logx=True, legend=False)
        ax.set_xlim(5, 5000); ax.set_ylim(bottom=0)
    fig, axs = C.compose([draw_zeta, draw_dls], ncols=2, width=S.DOUBLE, panel_aspect=0.8)
    save(fig, "12_Zeta_DLS")
    # 13 XPS (Hf 4f stacked + Zr 3d composite)
    tr = []
    for k in ["SampleA", "SampleB", "SampleC"]:
        x = R.read_xps(os.path.join(ED, "XPS", f"{k}_Hf4f_CasaXPS.csv"))
        tr.append((label_name(k), x.x.values, x.y.values))
    zx = R.read_xps(os.path.join(ED, "XPS", "SampleC_Zr3d_CasaXPS.csv"))
    def d_hf(ax):
        P.draw_xps(ax, tr, region=(13, 23)); ax.text(0.03, 0.96, "Hf 4f", transform=ax.transAxes, va="top")
    def d_zr(ax):
        P.draw_xps(ax, [(label_name("SampleC"), zx.x.values, zx.y.values)], region=(178, 190)); ax.text(0.03, 0.96, "Zr 3d", transform=ax.transAxes, va="top")
    fig, axs = C.compose([d_hf, d_zr], ncols=2, width=S.DOUBLE * 0.85, panel_aspect=0.85)
    save(fig, "13_XPS_Hf4f_Zr3d")
    # 14 Defect calculator - composite + parent from the example files
    res_list = []
    for key in ["SampleA", "SampleB", "SampleC"]:
        inp = example_defect_input(key)
        if not ARGS.png_only:   # --png-only (tests) must not rewrite the tracked example inputs
            with open(os.path.join(ED, "Defects", f"defect_input_{key}.json"), "w", encoding="utf-8") as f:
                json.dump(report.input_to_dict(inp), f, indent=1)
        res = D.run_all(inp)
        res_list.append(res)
        fig = P.fig_defect_comparison(res); save(fig, f"14_defects_{key}")
        report.pdf_report(res, os.path.join(OUT, f"14_defect_report_{key}.pdf"))
        D.results_table(res).to_csv(os.path.join(OUT, f"14_defect_table_{key}.csv"), index=False)
    # 15 composed main-text figure
    fig, axs = C.compose([draw_pxrd, draw_iso, (draw_tga, {}), (draw_kin, {"legend": False})], ncols=2, width=S.DOUBLE, panel_aspect=0.82)
    save(fig, "15_composed_Figure1")
    v2_demos()
    gallery()
    print("\n".join(written))
    return res_list


def example_defect_input(key):
    """Build calculator inputs FROM THE EXAMPLE FILES (TGA trace, NMR digest, ICP), like a user would."""
    truth = json.load(open(os.path.join(ED, "TGA", "example_design_truth.json"), encoding="utf-8"))
    nm = {"SampleA": "Sample A", "SampleB": "Sample B", "SampleC": "Sample C"}[key]
    tfile = {"SampleA": "SampleA_Netzsch.txt", "SampleB": "SampleB_TA_UniversalAnalysis.txt", "SampleC": "SampleC_Mettler_STARe.txt"}[key]
    t = R.read_tga(os.path.join(ED, "TGA", tfile)); ta = A.tga_analysis(t["T"], t["w"])
    nmr = R.read_nmr(os.path.join(ED, "NMR", f"{key}_digest_totxt.txt"))
    I = lambda a, b: A.integrate(nmr.x.values, nmr.y.values, a, b)
    I_std = I(-0.15, 0.15)
    norm = 9 / I_std
    def mk(name, integ, at_plateau=False, main=False):
        f, z, role, nH, _ = SPECIES[name]
        return D.Species(name, f, z, role, nH, integ * norm, at_plateau, main)
    sp = [mk("Fumarate", I(6.4, 6.65), main=True), mk("Formate", I(8.35, 8.55))]
    ssI = I(2.4, 2.9)
    if ssI * norm > 0.05:
        s = mk("2-Sulfosuccinate (dicarboxylate, SO3H kept)", ssI); s.nH = 2; sp.append(s)
    icp = R.read_icp(os.path.join(ED, "ICP", "ICP-OES_Syngistix_export.csv"))
    g = icp[icp["sample"] == nm].set_index("element")
    wt = lambda el: float(g.loc[el, "conc"]) * 0.050 / 10.0 * 100
    tr = truth[nm]
    # Mexp for the NMR mass (the 'user' weighed 1.000 µmol-equivalent in 0.60 mL -> [MOF] = 1.0 mM in the design)
    from tools_helpers import design_mass
    inp = D.DefectInput(sample=nm, framework="MOF-801", metal="Hf", extraframework_Zr=(key == "SampleC"), species=sp,
                        tga_W_initial=ta["W_initial"], tga_W_200=float(np.interp(200, ta["T"], ta["w"])),
                        tga_W_plateau=ta["W_plateau"], tga_W_end=ta["W_end"],
                        nmr_std_conc_mM=2.0, nmr_std_integral=9.0, nmr_sample_mg=design_mass(tr) * 0.6e-3 * 1.0, nmr_volume_mL=0.6,
                        cl_per_node=tr["Cl"], icp_Zr_wt=wt("Zr"), icp_Hf_wt=wt("Hf"), icp_state="as-measured (RT)",
                        xrf_S_per_M=None, V_meas={"SampleA": 0.192, "SampleB": 0.185, "SampleC": 0.165}[key], V_ideal=0.178,
                        linker_umol_per_mg=(A.integrate(nmr.x.values, nmr.y.values, 6.4, 6.65) * norm / 2) / 9 * 9 * 2.0 * 0.6 / (design_mass(tr) * 0.6e-3),
                        tit_mass_mg=50.0 if key == "SampleA" else None, tit_c_NaOH=0.1,
                        tit_V_eq=find_eq(key), tit_HCl_mmol=0.050)
    return inp


def find_eq(key):
    if key != "SampleA":
        return None
    df = R.read_table(os.path.join(ED, "Titration", "SampleA_titration_NaOH.csv"))
    eq = A.equivalence_points(df.iloc[:, 0].values, df.iloc[:, 1].values, n=3)
    return [float(v) for v in eq]


def v2_demos():
    """v2: scale-bar options, custom ticks/legend/units, N2 full-report before/after."""
    from mfs import figctl as F
    from PIL import Image, ImageDraw, ImageFont
    # 16 SEM scale bar options
    d = R.read_sem(os.path.join(ED, "SEM", "SampleA_SEM.tif"))
    pix_mag = R.pixel_from_magnification(25000, d["image"].shape[1])
    specs = [(d["pixel_nm"], dict(mode="auto"), "auto: 500 nm"),
             (d["pixel_nm"], dict(mode="value", value=1, unit="µm", position="lower left", color="black", box=True, box_alpha=0.85), "value 1 µm, black, box"),
             (d["pixel_nm"], dict(mode="length", length_frac=0.35, unit="nm", position="lower left", thickness=0.03, font_size=9), "length 35 % → label"),
             (pix_mag, dict(mode="value", value=200, unit="nm", position="lower right", box=True, thickness=0.012), "mag. 25 000× override, 200 nm")]
    ims = [dict(image=d["image"], pixel_nm=px, crop_row=d["databar_row"], bar=b, title=t) for px, b, t in specs]
    fig = P.fig_sem_panel(ims, 2, S.DOUBLE)
    save(fig, "16_SEM_scalebar_options", ("png",))
    # 17 default vs customised (same data)
    s = load_sorption()
    isos = [(label_name(k), s[k]["iso"]) for k in ["SampleA", "SampleB", "SampleC"]]
    S.apply_style()
    f1, a1 = S.new_figure(S.SINGLE, 0.85); P.draw_isotherms(a1, isos)
    f2, a2 = S.new_figure(S.SINGLE, 0.85); P.draw_isotherms(a2, isos)
    cfg = F.default_settings()
    cfg["title"] = "N_2 sorption, 77 K"
    cfg["ticks"].update(preset="2 sides (bottom + left)", direction="out", frame="match ticks", minor=True, major_len=4)
    cfg["panels"]["0"] = {"y_family": "uptake", "y_unit": "mmol g⁻¹", "ymajor": 2, "yminor": 0.5, "ymax": 13,
                          "xmajor": 0.2, "xminor": 0.05, "xlabel": "p/p_0"}
    cfg["legend"].update(custom=True, loc="outside top", ncol=2, frame=True, markerscale=1.4, fontsize=6.5,
                         entries=[{"label": "Sample C", "text": "composite", "show": True},
                                  {"label": "Sample B", "text": "Sample B", "show": True},
                                  {"label": "Sample A", "text": "Sample A (parent)", "show": True},
                                  {"label": "adsorption", "text": "ads. (filled)", "show": True},
                                  {"label": "desorption", "text": "des. (open)", "show": True}])
    F.apply(f2, cfg)
    p1 = save(f1, "17a_default", ("png",))[0]; p2 = save(f2, "17b_custom", ("png",))[0]
    with open(os.path.join(OUT, "17_custom_settings.json"), "w", encoding="utf-8") as fh:
        fh.write(F.to_json(cfg))
    A_, B_ = Image.open(p1).convert("RGB"), Image.open(p2).convert("RGB")
    h = max(A_.height, B_.height) + 90
    G = Image.new("RGB", (A_.width + B_.width + 80, h), "white")
    G.paste(A_, (20, 90 + (h - 90 - A_.height) // 2)); G.paste(B_, (A_.width + 60, 90 + (h - 90 - B_.height) // 2))
    dr = ImageDraw.Draw(G)
    from tools.fonts import sans                  # Liberation on Linux, Arial on Windows, Helvetica on macOS
    fnt = sans(44, bold=True)
    dr.text((30, 20), "default (v2)", fill="black", font=fnt)
    dr.text((A_.width + 70, 20), "customised: title, mmol/g, ticks out (2 sides), increments, legend outside top",
            fill="black", font=sans(40))
    G.save(os.path.join(OUT, "17_custom_ticks_legend.png"), dpi=(600, 600))
    os.remove(p1); os.remove(p2); written.remove(p1); written.remove(p2)
    written.append(os.path.join(OUT, "17_custom_ticks_legend.png"))
    # 18 N2 MicroActive full-report export: v1 vs v2 reader
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import n2_before_after
    r = n2_before_after.main(os.path.join(ED, "Sorption", "SampleA_3Flex_full_report_EXAMPLE.xls"), os.path.join(OUT, "18_N2_fullreport"))
    written.extend([r["before"], r["after"]])


def gallery():
    from PIL import Image
    files = sorted(glob.glob(os.path.join(OUT, "[01]*.png")))
    files = [f for f in files if "gallery" not in f]
    thumbs = []
    W = 900
    for f in files:
        im = Image.open(f).convert("RGB")
        im.thumbnail((W, 700))
        thumbs.append(im)
    cols = 4
    rows = (len(thumbs) + cols - 1) // cols
    cellh = 700
    G = Image.new("RGB", (cols * W + (cols + 1) * 30, rows * cellh + (rows + 1) * 30), "white")
    for i, im in enumerate(thumbs):
        r, c = divmod(i, cols)
        x = 30 + c * (W + 30) + (W - im.width) // 2
        y = 30 + r * (cellh + 30) + (cellh - im.height) // 2
        G.paste(im, (x, y))
    G.save(os.path.join(OUT, "gallery.png"), dpi=(150, 150))
    written.append(os.path.join(OUT, "gallery.png"))


if __name__ == "__main__":
    main()
