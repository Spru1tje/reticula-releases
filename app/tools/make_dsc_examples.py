"""Synthetic DSC example files (EXAMPLE data) in every supported vendor format, with known truth:
heating 25->200 °C at 10 K/min: Tg 60 °C (Δq 0.05 W/g), melting peak 120 °C (ΔH = -80 J/g, endothermic);
cooling 200->25 °C: crystallisation 95 °C (+70 J/g); second heating again. Plus an STA (TGA-DSC) file."""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "example_data", "DSC")
TRUTH = dict(beta=10.0, Tg=60.0, dq=0.05, Tm=120.0, dHm=-80.0, Tc=95.0, dHc=70.0, mass=5.0)


def _peak(T, Tc, area, beta, sT=2.0):
    st = sT * 60 / beta                              # s
    t = (T - Tc) * 60 / beta
    return area * np.exp(-0.5 * (t / st) ** 2) / (st * np.sqrt(2 * np.pi))


def _sig(T, T0, h, w=2.5):
    return h / (1 + np.exp(-(T - T0) / w))


def segments(seed=0):
    rng = np.random.default_rng(seed); b = TRUTH["beta"]
    T1 = np.arange(25, 200.001, 0.1)
    h1 = -0.02 - 0.0004 * (T1 - 25) - _sig(T1, TRUTH["Tg"], TRUTH["dq"]) + _peak(T1, TRUTH["Tm"], TRUTH["dHm"], b)
    T2 = T1[::-1]
    h2 = 0.03 + 0.0003 * (T2 - 25) + _peak(T2, TRUTH["Tc"], TRUTH["dHc"], b)
    h3 = h1.copy()
    segs = []
    t0 = 0.0
    for i, (T, h) in enumerate([(T1, h1), (T2, h2), (T1, h3)]):
        t = t0 + np.arange(len(T)) * 0.1 / b         # min
        t0 = t[-1] + 0.1 / b
        segs.append(pd.DataFrame({"t": t, "T": T, "q": h + rng.normal(0, 2e-4, len(T)), "seg": i + 1}))
    return segs


def write_all():
    os.makedirs(OUT, exist_ok=True)
    s = segments(); m = TRUTH["mass"]
    # TA TRIOS (exo up, W/g normalised), one [step] per segment
    with open(os.path.join(OUT, "MOF-PEG_TA-TRIOS_EXAMPLE.txt"), "w", encoding="utf-8") as fh:
        fh.write("Filename\tMOF-PEG_TA-TRIOS_EXAMPLE.tri\nInstrument\tDSC2500\nSize\t5.0000\tmg\n\n")
        for i, d in enumerate(s):
            fh.write("[step]\n" + ["Ramp 10.00 °C/min to 200.00 °C", "Ramp 10.00 °C/min to 25.00 °C", "Ramp 10.00 °C/min to 200.00 °C"][i] + "\n")
            fh.write("Time\tTemperature\tHeat Flow (Normalized)\nmin\t°C\tW/g\n")
            for r in d.itertuples():
                fh.write(f"{r.t:.4f}\t{r.T:.3f}\t{r.q:.6f}\n")
            fh.write("\n")
    # TA Universal Analysis (one heating, mW)
    d = s[0]
    with open(os.path.join(OUT, "MOF-PEG_TA-UniversalAnalysis_EXAMPLE.txt"), "w", encoding="utf-8") as fh:
        fh.write("Instrument\tDSC Q2000\nSample\tMOF-PEG\nSize\t5.0000 mg\nSig1\tTime (min)\nSig2\tTemperature (°C)\nSig3\tHeat Flow (mW)\nStartOfData\n")
        for r in d.itertuples():
            fh.write(f"{r.t:.4f}\t{r.T:.3f}\t{r.q * m:.6f}\n")
    # Netzsch Proteus (exo down: endothermic positive, decimal commas, ';'), all segments with a Segment column
    with open(os.path.join(OUT, "MOF-PEG_Netzsch_EXAMPLE.csv"), "w", encoding="cp1252") as fh:
        fh.write("#EXPORTTYPE:;DATA ALL\n#INSTRUMENT:;NETZSCH DSC 214 Polyma\n#SAMPLE:;MOF-PEG\n#SAMPLE MASS /mg:;5,000\n#SEGMENTS:;3\n")
        fh.write("##Temp./°C;Time/min;DSC/(mW/mg);Segment\n")
        for d in s:
            for r in d.itertuples():
                fh.write(f"{r.T:.3f};{r.t:.4f};{-r.q:.6f};{r.seg}\n".replace(".", ","))
    # Mettler STARe (endo negative, mW, decimal commas, fixed width), heating 1 only
    d = s[0]
    with open(os.path.join(OUT, "MOF-PEG_Mettler-STARe_EXAMPLE.txt"), "w", encoding="utf-8") as fh:
        fh.write("Title                                                           12.11.2026 10:00\n" + "_" * 80 + "\n\nCurve Name:\n  !]1[MOF-PEG\nCurve Values:\n")
        fh.write(f"{'Index':>15}{'t':>15}{'Ts':>15}{'Tr':>15}{'Value':>15}\n{'':>15}{'[s]':>15}{'[°C]':>15}{'[°C]':>15}{'[mW]':>15}\n")
        for k, r in enumerate(d.itertuples()):
            fh.write(f"{k:>15}{r.t * 60:>15.1f}{r.T:>15.3f}{r.T:>15.3f}{r.q * m:>15.6f}\n".replace(".", ","))
        fh.write("\nSample:\n  MOF-PEG, 5,0000 mg\nMethod:\n  DSC 25-200 °C 10 K/min\n" + "_" * 80 + "\nNot signed" + " " * 50 + "STARe SW 16.20\n")
    # PerkinElmer Pyris (endo up, mW)
    with open(os.path.join(OUT, "MOF-PEG_PerkinElmer-Pyris_EXAMPLE.txt"), "w", encoding="utf-8") as fh:
        fh.write("Pyris Software\nSample ID: MOF-PEG\nSample Weight: 5.000 mg\n\n")
        fh.write("Time\tUnsubtracted Heat Flow\tProgram Temperature\tSample Temperature\n(min)\t(mW)\t(°C)\t(°C)\n")
        for r in d.itertuples():
            fh.write(f"{r.t:.4f}\t{-r.q * m:.6f}\t{r.T:.3f}\t{r.T:.3f}\n")
    # generic CSV with segment column (exo up, W/g)
    pd.concat(s).rename(columns={"t": "Time (min)", "T": "Temperature (°C)", "q": "Heat flow (W/g)", "seg": "Segment"}) \
        .to_csv(os.path.join(OUT, "MOF-PEG_generic_EXAMPLE.csv"), index=False)
    # STA (TGA-DSC), Netzsch STA 449: desolvation endotherm + decomposition exotherm with mass loss
    T = np.arange(30, 700.001, 0.5); t = (T - 30) / 10
    w = 100 - 12 / (1 + np.exp(-(T - 110) / 12)) - 35 / (1 + np.exp(-(T - 520) / 15))
    q = -0.05 + _peak(T, 110, -250.0, 10, 18) + _peak(T, 525, 1500.0, 10, 12) - 0.0002 * (T - 30)
    with open(os.path.join(OUT, "MOF-801_STA_Netzsch_EXAMPLE.csv"), "w", encoding="cp1252") as fh:
        fh.write("#INSTRUMENT:;NETZSCH STA 449 F3\n#SAMPLE:;MOF-801\n#SAMPLE MASS /mg:;10,000\n##Temp./°C;Time/min;DSC/(mW/mg);Mass/%\n")
        for a, b, c_, d_ in zip(T, t, q, w):
            fh.write(f"{a:.2f};{b:.3f};{-c_:.6f};{d_:.4f}\n".replace(".", ","))


if __name__ == "__main__":
    write_all(); print(sorted(os.listdir(OUT)))
