"""Full heating/cooling cycle EXAMPLE built from the user's REAL Mettler STARe heating curve (SG_DSC_1.1, 1 K/min,
25 -> 70 °C) + a SIMULATED 5 min isothermal at 70 °C + a SIMULATED cooling 70 -> 25 °C at 1 K/min with a
crystallisation exotherm (onset 40.8 °C, ΔHc = +98 J/g). The cooling part is synthetic and labelled as EXAMPLE.

Writes example_data/DSC/full_cycle/:
  SG_DSC_full-cycle_Mettler_EXAMPLE.txt       one curve block, continuous time (case a: one file, all segments)
  SG_DSC_heating_Mettler_REAL.txt             the real heating export (case b, part 1)
  SG_DSC_cooling_Mettler_EXAMPLE.txt          the simulated cooling alone, t restarting at 0, later time stamp (case b, part 2)
  SG_DSC_3-curves_Mettler_EXAMPLE.txt         three 'Curve Name:' blocks heat / iso / cool, t restarting (case c)"""
import os, shutil, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mfs import dsc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL = os.path.join(ROOT, "tests", "data", "Mettler_STARe_SG_DSC_1.1_1K_REAL.txt")
OUT = os.path.join(ROOT, "example_data", "DSC", "full_cycle")
TRUTH = dict(T_on=40.8, k=3.0, theta=0.6, dHc=98.0, beta=1.0, iso_min=5.0)   # peak at T_on - (k-1)*theta = 39.6 °C


def build(seed=1):
    rng = np.random.default_rng(seed)
    d = dsc.read_dsc(REAL); c = d["curves"][0]
    t_h, T_h, Tr_h, q_h = c["t"] * 60, c["T"], c["Tr"], c["q"]             # s, °C, °C, W/g
    m = (T_h < 45) | (T_h > 52)
    pb = np.polyfit(T_h[m & (T_h > 27)], q_h[m & (T_h > 27)], 1)
    noise = float(np.std(q_h[m & (T_h > 27)] - np.polyval(pb, T_h[m & (T_h > 27)])))
    dt = float(np.median(np.diff(t_h)))
    # isothermal 70 °C: relaxes from the heating baseline to the (mirrored) cooling baseline
    n_iso = int(TRUTH["iso_min"] * 60 / dt)
    t_i = t_h[-1] + dt * (1 + np.arange(n_iso))
    qc0 = -np.polyval(pb, 70.0) + 0.004
    q_i = qc0 + (q_h[-1] - qc0) * np.exp(-(t_i - t_i[0]) / 30.0) + rng.normal(0, noise, n_iso)
    T_i = np.full(n_iso, 70.0) + rng.normal(0, 0.002, n_iso); Tr_i = np.full(n_iso, 70.0)
    # cooling 70 -> 25 °C at 1 K/min
    n_c = int(45 * 60 / dt) + 1
    t_c = t_i[-1] + dt * (1 + np.arange(n_c))
    Tr_c = 70.0 - (t_c - t_c[0]) / 60.0 * TRUTH["beta"]
    x = np.clip(TRUTH["T_on"] - Tr_c, 0, None) / TRUTH["theta"]
    shape = x ** (TRUTH["k"] - 1) * np.exp(-x)
    shape /= np.trapezoid(shape, t_c)                                    # 1/s
    exo = TRUTH["dHc"] * shape                                           # W/g, integral over s = J/g
    T_c = Tr_c + 0.03 + 0.02 * shape / shape.max()                       # small lag + recalescence
    q_c = -np.polyval(pb, Tr_c) + 0.004 + exo + rng.normal(0, noise, n_c)
    seg = {"heat": (t_h, T_h, Tr_h, q_h), "iso": (t_i, T_i, Tr_i, q_i), "cool": (t_c, T_c, Tr_c, q_c)}
    return seg, d["meta"]


def _fmt(v, nd):
    return f"{v:.{nd}f}".replace(".", ",")


def write_mettler(path, blocks, title_time, method, remark):
    """blocks: list of (curve name, t s, Ts, Tr, value W/g)."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{'Title':<64}{title_time}\n" + "_" * 80 + "\n\nTitle\nEvaluation Categories:\n\n")
        for name, t, T, Tr, q in blocks:
            fh.write(f"Curve Name:\n  !]2[{name}\nCurve Values:\n")
            fh.write(f"{'Index':>15}{'t':>15}{'Ts':>15}{'Tr':>15}{'Value':>15}\n{'':>15}{'[s]':>15}{'[°C]':>15}{'[°C]':>15}{'[Wg^-1]':>15}\n")
            for k in range(len(t)):
                fh.write(f"{k:>15}{_fmt(t[k], 0):>15}{_fmt(T[k], 3):>15}{_fmt(Tr[k], 3):>15}{_fmt(q[k], 8):>15}\n")
            fh.write("\n")
        fh.write(f"Sample:\n  SG_DSC_1.1_1K, 5,6400 mg\n\nSample Remarks:\n  {remark}\n\nMethod:\n  {method}\n\n"
                 "Sample Holder:\n  Aluminum Standard 40ul\n\n" + "_" * 80 + "\n" + f"{'Not signed':<66}STARe SW 16.20\n")


def write_all():
    os.makedirs(OUT, exist_ok=True)
    seg, meta = build()
    rm = "EXAMPLE: real heating (SG_DSC_1.1) + SIMULATED isothermal and cooling exotherm - not measured data"
    h, i, c = seg["heat"], seg["iso"], seg["cool"]
    full = [np.concatenate([h[k], i[k], c[k]]) for k in range(4)]
    write_mettler(os.path.join(OUT, "SG_DSC_full-cycle_Mettler_EXAMPLE.txt"),
                  [("SG_DSC_1.1_1K full cycle EXAMPLE (simulated cooling)", *full)], "03.04.2024 09:18",
                  "DSC +25°C+70°C 1°C/min 2 Steps (EXAMPLE: +5 min 70°C, -1°C/min to 25°C simulated)", rm)
    shutil.copy(REAL, os.path.join(OUT, "SG_DSC_heating_Mettler_REAL.txt"))
    write_mettler(os.path.join(OUT, "SG_DSC_cooling_Mettler_EXAMPLE.txt"),
                  [("SG_DSC_1.1_1K cooling EXAMPLE (simulated)", c[0] - c[0][0], c[1], c[2], c[3])], "03.04.2024 10:10",
                  "DSC +70°C+25°C -1°C/min (EXAMPLE, simulated)", rm)
    write_mettler(os.path.join(OUT, "SG_DSC_3-curves_Mettler_EXAMPLE.txt"),
                  [("SG_DSC_1.1_1K heating", h[0] - h[0][0], h[1], h[2], h[3]),
                   ("SG_DSC_1.1_1K 70°C isothermal EXAMPLE (simulated)", i[0] - i[0][0], i[1], i[2], i[3]),
                   ("SG_DSC_1.1_1K cooling EXAMPLE (simulated)", c[0] - c[0][0], c[1], c[2], c[3])], "03.04.2024 09:18",
                  "DSC +25°C+70°C 1°C/min 2 Steps (EXAMPLE: cooling simulated)", rm)
    with open(os.path.join(OUT, "README_EXAMPLE.txt"), "w", encoding="utf-8") as fh:
        fh.write(__doc__ + f"\nTruth of the simulated crystallisation: onset {TRUTH['T_on']} °C, peak "
                 f"{TRUTH['T_on'] - (TRUTH['k'] - 1) * TRUTH['theta']:.1f} °C (programme T), ΔHc = +{TRUTH['dHc']} J/g.\n")
    return sorted(os.listdir(OUT))


if __name__ == "__main__":
    print(write_all())
