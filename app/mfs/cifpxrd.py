"""Simulated PXRD from a CIF (pymatgen kinematic, pseudo-Voigt broadening) with hkl labels, and the overlay on a
measured pattern with a least-squares scale + linear background and the residual."""
from __future__ import annotations
import re
import numpy as np
from . import analysis as A


def hkl_tex(h):
    """'11-2' / '1, 1, -2' -> mathtext with overbars."""
    parts = re.findall(r"-?\d", str(h).replace(",", " "))
    return "(" + "".join(r"$\bar{%s}$" % p[1] if p.startswith("-") else p for p in parts) + ")"


def simulate(src, two_theta=(3, 50), wavelength="CuKa1", fwhm=0.10):
    p, sticks, _ = A.simulate_pxrd(src, two_theta, wavelength=wavelength, fwhm=fwhm)
    return p, sticks.sort_values("two_theta").reset_index(drop=True)


def fit_to(meas_x, meas_y, sim):
    """Scale s and linear background (b0 + b1·2θ) so that s·sim + bg ≈ measured (least squares on the measured grid)."""
    ys = np.interp(meas_x, sim.x, sim.y, left=0, right=0)
    Xd = np.vstack([ys, np.ones_like(meas_x), meas_x]).T
    coef, *_ = np.linalg.lstsq(Xd, meas_y, rcond=None)
    fit = Xd @ coef
    rwp = float(np.sqrt(np.sum((meas_y - fit) ** 2) / np.sum(meas_y ** 2)))
    return {"scale": coef[0], "bg": coef[1:], "fit": fit, "resid": meas_y - fit, "Rwp": rwp, "sim_on_grid": ys}


def draw(ax, sim, sticks, label_top=12, label_min=5.0, color="#0072B2", meas=None, fitres=None, rax=None, name="simulated"):
    """Simulated pattern with hkl labels for the strongest reflections (or overlay on a measured one if `meas`)."""
    if meas is not None and fitres is not None:
        mx, my = meas
        ax.plot(mx, my, color="k", lw=0.8, label="measured")
        ax.plot(mx, fitres["fit"], color=color, lw=0.8, label=f"{name} (scaled + bg)")
        top = max(my.max(), fitres["fit"].max())
        ax.vlines(sticks.two_theta, -0.06 * top, -0.02 * top, color=color, lw=0.5)
        ax.set_ylim(-0.08 * top, 1.15 * top)
        if rax is not None:
            rax.plot(mx, fitres["resid"], color="#7F7F7F", lw=0.6); rax.axhline(0, color="k", lw=0.4)
            rax.set_ylabel("Δ"); rax.set_xlabel(r"2$\theta$ (°)")
            rax.text(0.99, 0.9, f"Rwp = {100 * fitres['Rwp']:.1f} %", transform=rax.transAxes, ha="right", va="top", fontsize=6)
        ax.legend(frameon=False, fontsize=6, loc="upper right")
        ys = np.interp(sticks.two_theta, mx, fitres["fit"])
    else:
        ax.plot(sim.x, sim.y * 100, color=color, lw=0.8, label=name)
        ax.vlines(sticks.two_theta, -6, -2, color=color, lw=0.5); ax.set_ylim(-8, 118)
        ys = np.interp(sticks.two_theta, sim.x, sim.y * 100)
    lab = sticks[sticks.I >= label_min].nlargest(label_top, "I").sort_values("two_theta")
    groups = []                                     # reflections closer than ~1.2 % of the range share one label
    for _, r in lab.iterrows():
        if groups and r.two_theta - groups[-1][-1].two_theta < 0.012 * (sim.x.max() - sim.x.min()):
            groups[-1].append(r)
        else:
            groups.append([r])
    for g in groups:
        tt = float(np.mean([r.two_theta for r in g]))
        y = max(float(np.interp(r.two_theta, sticks.two_theta, ys)) for r in g)
        ax.annotate("/".join(hkl_tex(r.hkl) for r in g), (tt, y), xytext=(0, 4), textcoords="offset points", rotation=90, ha="center",
                    va="bottom", fontsize=5)
    ax.set_ylabel("Intensity (a.u.)"); ax.set_yticks([])
    if rax is None:
        ax.set_xlabel(r"2$\theta$ (°)")
    ax.set_xlim(sim.x.min(), sim.x.max())


def table(sticks, n=None):
    d = sticks.copy(); d["hkl"] = d.hkl.astype(str)
    d = d.rename(columns={"two_theta": "2θ (°)", "I": "I (rel.)", "d": "d (Å)"})
    return (d.nlargest(n, "I (rel.)").sort_values("2θ (°)") if n else d).round(4)
