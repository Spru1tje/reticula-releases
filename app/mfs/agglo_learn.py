"""v2.21: learn which objects are agglomerates from the user's own tags.

The user tags objects on a micrograph as 'agglomerate' (label 1) or 'rod / particle' (label 0). Each tag is stored with the
measured features of that object (not with the picture), so tags from several micrographs can be combined.

Model: L2-regularised logistic regression on six standardised features, fitted by iteratively re-weighted least squares
(Newton). Small, transparent, and easy to inspect: the coefficients say which feature drives the decision.

Rules of use
  * The learned model is only applied when there are at least MIN_PER_CLASS tags of each kind.
  * Its quality is reported as k-fold cross-validated accuracy and balanced accuracy, on the tags only.
  * Without enough tags the rule of mfs.rods.flag_agglomerates is used.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

FEATURES = ["log10_area_nm2", "log10_length_nm", "log10_width_nm", "log10_aspect", "solidity", "log10_rel_area"]
MIN_PER_CLASS = 20
MATCH_PX = 3.0          # a tag belongs to the object whose centre is this close (pixels)
LAMBDA = 1.0            # L2 penalty on the standardised coefficients (not on the intercept)


def default_store(sample: str | None = None) -> Path:
    """Tag file of one sample. Tags are kept per sample, so a new sample starts without the tags of the previous one."""
    base = Path(os.environ.get("RETICULA_DATA", Path.home() / ".reticula")) / "tags"
    base.mkdir(parents=True, exist_ok=True)
    name = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in (sample or "default")).strip("_") or "default"
    return base / f"{name}.json"


def features(obj: dict, nm_per_px: float, median_area_nm2: float) -> dict:
    """The feature vector of one object (the same numbers are used when tagging and when applying the model)."""
    area = max(obj["area_px"] * nm_per_px ** 2, 1e-6)
    return {
        "log10_area_nm2": math.log10(area),
        "log10_length_nm": math.log10(max(obj["length_nm"], 1e-3)),
        "log10_width_nm": math.log10(max(obj["width_nm"], 1e-3)),
        "log10_aspect": math.log10(max(obj["aspect"], 1e-3)),
        "solidity": float(obj["solidity"]),
        "log10_rel_area": math.log10(area / max(median_area_nm2, 1e-6)),
    }


def _vec(f: dict) -> np.ndarray:
    return np.array([f[k] for k in FEATURES], float)


# ----------------------------------------------------------------------------- tag store
def load_tags(path=None) -> list[dict]:
    p = Path(path) if path else default_store()
    if not p.exists():
        return []
    return json.loads(p.read_text(encoding="utf-8")).get("tags", [])


def save_tags(tags: list[dict], path=None) -> Path:
    p = Path(path) if path else default_store()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"version": 1, "features": FEATURES, "tags": tags}, indent=1), encoding="utf-8")
    return p


def add_tag(tags: list[dict], image: str, obj: dict, nm_per_px: float, median_area_nm2: float, label: int) -> list[dict]:
    """Add or replace the tag of one object (same image and centre within MATCH_PX counts as the same object)."""
    out = [t for t in tags if not (t["image"] == image and math.hypot(t["cx"] - obj["cx"], t["cy"] - obj["cy"]) < MATCH_PX)]
    out.append({"image": image, "cx": float(obj["cx"]), "cy": float(obj["cy"]), "label": int(label),
                "features": features(obj, nm_per_px, median_area_nm2)})
    return out


def remove_tag(tags: list[dict], image: str, obj: dict) -> list[dict]:
    return [t for t in tags if not (t["image"] == image and math.hypot(t["cx"] - obj["cx"], t["cy"] - obj["cy"]) < MATCH_PX)]


def counts(tags) -> dict:
    y = [t["label"] for t in tags]
    return {"agglomerate": int(sum(1 for v in y if v == 1)), "rod": int(sum(1 for v in y if v == 0))}


# ----------------------------------------------------------------------------- model
def _fit_logistic(X: np.ndarray, y: np.ndarray, lam: float = LAMBDA, iters: int = 50):
    """Newton / IRLS for L2-regularised logistic regression. X is standardised (n, d). Returns (w, b)."""
    n, d = X.shape
    A = np.hstack([X, np.ones((n, 1))])
    R = np.eye(d + 1) * lam
    R[-1, -1] = 0.0                                   # no penalty on the intercept
    beta = np.zeros(d + 1)
    for _ in range(iters):
        z = A @ beta
        p = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        g = A.T @ (p - y) + R @ beta
        W = p * (1 - p) + 1e-9
        H = A.T @ (A * W[:, None]) + R
        step = np.linalg.solve(H, g)
        beta -= step
        if np.max(np.abs(step)) < 1e-8:
            break
    return beta[:-1], float(beta[-1])


def fit(tags: list[dict]) -> dict | None:
    """Fit the model on all tags. Returns None when there are not enough tags of both kinds."""
    c = counts(tags)
    if c["agglomerate"] < MIN_PER_CLASS or c["rod"] < MIN_PER_CLASS:
        return None
    X = np.array([_vec(t["features"]) for t in tags])
    y = np.array([t["label"] for t in tags], float)
    mu, sd = X.mean(0), X.std(0) + 1e-9
    w, b = _fit_logistic((X - mu) / sd, y)
    return {"w": w.tolist(), "b": b, "mu": mu.tolist(), "sd": sd.tolist(), "features": FEATURES, "n": len(tags),
            "counts": c, "coef": dict(zip(FEATURES, (w / 1.0).round(3).tolist()))}


def predict_proba(model: dict, f: dict) -> float:
    x = (_vec(f) - np.array(model["mu"])) / np.array(model["sd"])
    z = float(x @ np.array(model["w"]) + model["b"])
    return 1.0 / (1.0 + math.exp(-max(min(z, 30), -30)))


def cross_validate(tags: list[dict], k: int = 5, seed: int = 0) -> dict | None:
    """k-fold cross-validation on the tags. Returns accuracy, balanced accuracy and the number of folds used."""
    if counts(tags)["agglomerate"] < MIN_PER_CLASS or counts(tags)["rod"] < MIN_PER_CLASS:
        return None
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(tags))
    folds = np.array_split(idx, k)
    pred = np.zeros(len(tags), int)
    truth = np.array([t["label"] for t in tags])
    X_all = np.array([_vec(t["features"]) for t in tags])
    for f_idx in folds:
        train = np.setdiff1d(idx, f_idx)
        ytr = truth[train]
        if ytr.sum() == 0 or ytr.sum() == len(ytr):
            continue
        mu, sd = X_all[train].mean(0), X_all[train].std(0) + 1e-9
        w, b = _fit_logistic((X_all[train] - mu) / sd, ytr.astype(float))
        z = ((X_all[f_idx] - mu) / sd) @ w + b
        pred[f_idx] = (z > 0).astype(int)
    acc = float(np.mean(pred == truth))
    tpr = float(np.mean(pred[truth == 1] == 1)); tnr = float(np.mean(pred[truth == 0] == 0))
    return {"accuracy": acc, "balanced_accuracy": 0.5 * (tpr + tnr), "sensitivity_agglomerate": tpr,
            "specificity_rod": tnr, "k": k, "n": len(tags)}


def apply(objs: list[dict], model: dict | None, nm_per_px: float, fallback) -> dict:
    """Set obj['agglomerate'] with the learned model (threshold p >= 0.5) or with the fallback rule. Returns a summary."""
    if model is None:
        fallback(objs, nm_per_px)
        return {"method": "rule", "n_agglomerate": sum(o.get("agglomerate", False) for o in objs)}
    areas = np.array([o["area_px"] * nm_per_px ** 2 for o in objs]) if objs else np.array([1.0])
    med = float(np.median(areas))
    for o in objs:
        p = predict_proba(model, features(o, nm_per_px, med))
        o["p_agglomerate"] = round(p, 3)
        o["agglomerate"] = bool(p >= 0.5)
    return {"method": "learned", "n_agglomerate": sum(o["agglomerate"] for o in objs), "n_tags": model["n"]}
