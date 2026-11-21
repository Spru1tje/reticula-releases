"""Electrophoretic mobility -> zeta potential (v2.18.11).

- Smoluchowski (valid for kappa*a >> 1, i.e. large particles in electrolyte): mu = eps * zeta / eta  ->  zeta = eta * mu / eps.
- Hückel (valid for kappa*a << 1, small particles / low ionic strength): mu = 2 eps zeta / (3 eta)  ->  zeta = 3 eta mu / (2 eps).
Both are the standard textbook limits (Smoluchowski 1906; Hückel 1924). Real particles lie between the limits; the function says which
model applies when kappa*a is given, and warns when the chosen model is outside its range. General practice: report the model used."""
from __future__ import annotations

EPS0 = 8.8541878128e-12  # F/m
ETA_WATER_25C = 0.89e-3  # Pa s
EPS_R_WATER_25C = 78.5


def zeta_from_mobility(mu, eta=ETA_WATER_25C, eps_r=EPS_R_WATER_25C, model="smoluchowski", ka=None):
    """mu in m^2/(V s). Returns dict with zeta in V and mV, the model used, and a note on validity."""
    eps = eps_r * EPS0
    if model == "smoluchowski":
        zeta = eta * mu / eps
    elif model == "huckel":
        zeta = 3.0 * eta * mu / (2.0 * eps)
    else:
        raise ValueError("model must be 'smoluchowski' or 'huckel'")
    note = ""
    if ka is not None:
        if model == "smoluchowski" and ka < 100:
            note = f"kappa*a = {ka:g} is below ~100: Smoluchowski underestimates; consider Henry / Ohshima."
        if model == "huckel" and ka > 1:
            note = f"kappa*a = {ka:g} is above ~1: Hückel underestimates; consider Smoluchowski."
    return {"zeta_V": zeta, "zeta_mV": zeta * 1e3, "model": model, "note": note}
