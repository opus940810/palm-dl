"""Merge per-model out-of-fold predictions into a single table."""
from __future__ import annotations

import pandas as pd

from config import DATA_PROC, RESULTS

COLS = ["model", "fold", "set", "acc", "pos", "y", "p", "source", "is_human"]


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    mp = sites.set_index(["acc", "pos"]).is_human

    def load(name: str) -> pd.DataFrame:
        path = RESULTS / name
        if not path.exists():
            print("missing", name)
            return pd.DataFrame(columns=COLS)
        d = pd.read_csv(path)
        if "is_human" not in d.columns:
            d["is_human"] = [bool(mp.get((a, p), False)) for a, p in zip(d.acc, d.pos)]
        for c in COLS:
            if c not in d.columns:
                d[c] = None
        return d[COLS]

    cv = pd.concat([load("preds_cv_base.csv"), load("preds_cv_cnn.csv"),
                    load("preds_cv_esm.csv"), load("preds_cv_ensemble.csv")], ignore_index=True)
    cv.to_csv(RESULTS / "preds_cv_all.csv", index=False)
    print("cv_all:", cv.groupby(["model", "set"]).size().to_dict())

    h = pd.concat([load("preds_cvh_base.csv"), load("preds_cvh_cnn.csv"), load("preds_cvh_esm.csv"), load("preds_cvh_ensemble.csv")], ignore_index=True)
    h.to_csv(RESULTS / "preds_cvh_all.csv", index=False)
    print("cvh_all:", h.groupby(["model", "set"]).size().to_dict())


if __name__ == "__main__":
    main()
