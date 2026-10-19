"""Fit a Platt calibrator for the forest view alone.

The proteome-wide scan cannot run the ESM view at this scale on the available
hardware (the CPU-only torch model stalls on long inference runs), so the scan is
scored with the forest view. This calibrator maps forest logits to probabilities
using the same out-of-fold validation rows and prevalence as the ensemble.
"""
from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score

from config import MODELS, RESULTS
from final_model import logit


def main() -> None:
    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    va = d[(d.model == "rf") & (d.set == "val")]
    z = logit(va.p.values).reshape(-1, 1)
    y = va.y.values
    clf = LogisticRegression(C=1.0, max_iter=1000).fit(z, y)
    p = clf.predict_proba(z)[:, 1]
    meta = {"n": int(len(y)), "n_pos": int(y.sum()),
            "auc": float(roc_auc_score(y, va.p.values)),
            "brier_before": float(brier_score_loss(y, va.p.values)),
            "brier_after": float(brier_score_loss(y, p)),
            "coef": float(clf.coef_[0, 0]), "intercept": float(clf.intercept_[0])}
    joblib.dump(clf, MODELS / "calibrator_rf.joblib")
    (MODELS / "calibrator_rf_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
