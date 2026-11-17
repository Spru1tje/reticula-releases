"""v2.14: the two SYNTHETIC photoluminescence example files (example_data/PL). Run:  python tools/make_pl_examples.py

1. PL_thermometry_Er_upconversion_EXAMPLE.csv - Er3+ visible upconversion bands 2H11/2 (~524 nm) and 4S3/2 (~545 nm) of an Yb,Er-doped
   nanocrystal, 293.15-323.15 K in 5 K steps (wide table: wavelength + one column per temperature; photon counts with Poisson noise).
   The band ratio follows D = alpha exp(-dE/kB T) with dE = 741 cm-1 and alpha = 10.7 (values of the order reported for NaYF4:Yb,Er).
2. PL_thermal_quenching_phosphor_EXAMPLE.csv - one broad band (Eu2+-like, ~450 nm) between 80 K and 500 K with single-channel thermal
   quenching I0 / (1 + A exp(-Ea / kB T)), Ea = 0.25 eV, A = 2000 (T50 ~ 380 K), a slight red shift and broadening with temperature."""
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "example_data", "PL")
KB = 8.617333262e-5
EV_CM = 8065.54429


def g(x, c, s):
    return np.exp(-0.5 * ((x - c) / s) ** 2)


def thermometry(seed=7):
    rng = np.random.default_rng(seed)
    wl = np.arange(500.0, 700.01, 0.5)
    temps = np.arange(293.15, 323.16, 5.0)
    cols = []
    for T in temps:
        D = 10.7 * np.exp(-(741.0 / EV_CM) / (KB * T))
        i2 = 52000.0 * (1 - 0.0016 * (T - 293.15))              # 4S3/2 (545 nm), counts at the maximum
        i1 = D * i2 * (3.0 / 2.6)                               # 2H11/2 (524 nm); the widths differ, so the band AREAS follow D
        y = i1 * g(wl, 524.0, 2.6) + i2 * g(wl, 545.0, 3.0) + 9000.0 * g(wl, 658.0, 6.0) + 120.0
        cols.append(rng.poisson(y).astype(float))
    return wl, [f"{T:.2f} K" for T in temps], cols


def quenching(seed=11):
    rng = np.random.default_rng(seed)
    wl = np.arange(380.0, 620.01, 1.0)
    temps = np.array([80, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500], float)
    cols = []
    for T in temps:
        I = 1.0 / (1.0 + 2000.0 * np.exp(-0.25 / (KB * T)))
        c = 450.0 + 0.012 * (T - 80.0)                           # red shift
        s = 22.0 * np.sqrt(1.0 + 0.6 * (T / 300.0) ** 2) / 1.3     # broadening
        y = 41000.0 * I * (22.0 / s) * g(wl, c, s) + rng.normal(0, 60.0, wl.size) + 150.0     # area ~ I(T): broadening does not change it
        cols.append(np.clip(y, 0, None))
    return wl, [f"{T:.0f} K" for T in temps], cols


def write(path, wl, names, cols, note):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(f"# {note}\n")
        fh.write("Wavelength (nm)," + ",".join(names) + "\n")
        for k in range(len(wl)):
            fh.write(f"{wl[k]:.1f}," + ",".join(f"{c[k]:.0f}" for c in cols) + "\n")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    write(os.path.join(OUT, "PL_thermometry_Er_upconversion_EXAMPLE.csv"), *thermometry(),
          "SYNTHETIC EXAMPLE - Er3+ upconversion emission (2H11/2 and 4S3/2) 293-323 K, photon counts; not a measurement")
    write(os.path.join(OUT, "PL_thermal_quenching_phosphor_EXAMPLE.csv"), *quenching(),
          "SYNTHETIC EXAMPLE - single-band phosphor emission 80-500 K with thermal quenching (Ea 0.25 eV); not a measurement")
    print("written to", OUT)
