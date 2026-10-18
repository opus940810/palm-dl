"""Evaluation metrics and bootstrap confidence intervals."""
from __future__ import annotations

import numpy as np
from sklearn.metrics import (average_precision_score, brier_score_loss, confusion_matrix,
                             f1_score, matthews_corrcoef, roc_auc_score)


def find_threshold(y: np.ndarray, p: np.ndarray) -> float:
    """Threshold that maximises MCC on the given (validation) set.

    Vectorised: sort once and evaluate every cut-off from cumulative counts
    instead of calling sklearn once per candidate threshold (which is O(n^2)).
    """
    y = np.asarray(y).astype(np.int64)
    p = np.asarray(p, dtype=float)
    order = np.argsort(-p, kind="mergesort")
    p_sorted = p[order]
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    n_pos, n_neg = int(y_sorted.sum()), int(len(y_sorted) - y_sorted.sum())
    tn, fn = n_neg - fp, n_pos - tp
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = np.where(denom > 0, (tp * tn - fp * fn) / np.maximum(denom, 1e-12), 0.0)
    return float(p_sorted[int(np.argmax(mcc))])


def ece(y: np.ndarray, p: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    total = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (p >= lo) & (p < hi if hi < 1 else p <= hi)
        if m.sum():
            total += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(total)


def metrics_at(y: np.ndarray, p: np.ndarray, t: float) -> dict:
    pred = (p >= t).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    sn = tp / (tp + fn) if tp + fn else 0.0
    sp = tn / (tn + fp) if tn + fp else 0.0
    return {
        "n": int(len(y)), "n_pos": int(y.sum()),
        "threshold": round(float(t), 4),
        "acc": (tp + tn) / len(y),
        "sn": sn, "sp": sp,
        "mcc": matthews_corrcoef(y, pred),
        "f1": f1_score(y, pred, zero_division=0),
        "auc": roc_auc_score(y, p) if len(np.unique(y)) > 1 else float("nan"),
        "aupr": average_precision_score(y, p),
        "brier": brier_score_loss(y, p),
        "ece": ece(y, p),
        "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
    }


def bootstrap_ci(y: np.ndarray, p: np.ndarray, t: float, metric: str = "mcc",
                 n_boot: int = 1000, seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    vals = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy, pp = y[idx], p[idx]
        if yy.min() == yy.max():
            continue
        vals.append(metrics_at(yy, pp, t)[metric])
    return tuple(np.percentile(vals, [2.5, 97.5]))
