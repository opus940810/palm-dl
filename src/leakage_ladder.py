"""Benchmark audit: how much does the train/test split protocol inflate MCC?

The same random-forest model is trained and evaluated under four protocols that
are increasingly strict about homology:

  1. site_random    - rows split at random (identical windows appear in both sets)
  2. protein        - proteins split at random, but no protein appears twice
  3. cluster        - 4-mer cosine >= 0.30 clusters split (near-homologues disjoint)
  4. species_human  - trained without human data, tested on human sites

For each protocol we also report the fraction of test-positive windows that are
exactly present among training positives, and the AUC / MCC.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, train_test_split

from config import DATA_PROC, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from metrics import find_threshold, metrics_at
from models import make_rf, set_seed
from run_cv import balanced_subset
from run_species import half_of


def eval_protocol(name: str, df: pd.DataFrame, train_idx, test_idx, X, ctx) -> dict:
    tr = df.iloc[train_idx]
    rng = np.random.default_rng(SEED)
    val_clusters = set(rng.choice(tr.cluster.unique(), size=max(2, int(0.15 * tr.cluster.nunique())),
                                 replace=False))
    is_val = tr.cluster.isin(val_clusters).values
    val_idx = train_idx[is_val]
    inner = df.iloc[train_idx[~is_val]]
    bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
    n_pos = int(inner.label.sum())
    keep_bg = bg.sample(n=min(len(bg), n_pos), random_state=SEED)
    fit_idx = np.concatenate([inner[inner.source != "background_cys"].index.values, keep_bg.index.values])

    set_seed(SEED)
    clf = make_rf()
    clf.fit(np.hstack([X[fit_idx], ctx[fit_idx]]), df.label.values[fit_idx])
    p_val = clf.predict_proba(np.hstack([X[val_idx], ctx[val_idx]]))[:, 1]
    p_test = clf.predict_proba(np.hstack([X[test_idx], ctx[test_idx]]))[:, 1]
    y_val, y_test = df.label.values[val_idx], df.label.values[test_idx]
    t_nat = find_threshold(y_val, p_val)
    bi = balanced_subset(df.iloc[test_idx], p_test, SEED + 11)
    bval = balanced_subset(df.iloc[val_idx], p_val, SEED + 5)
    t_bal = find_threshold(y_val[bval], p_val[bval])
    m = metrics_at(y_test, p_test, t_nat)
    mb = metrics_at(y_test[bi], p_test[bi], t_bal)
    train_windows = set(df.window.values[fit_idx][df.label.values[fit_idx] == 1])
    test_pos_w = df.window.values[test_idx][y_test == 1]
    leak = np.mean([w in train_windows for w in test_pos_w]) if len(test_pos_w) else 0.0
    row = {"protocol": name, "n_train": len(fit_idx), "n_train_pos": int(n_pos),
           "n_test": len(test_idx), "n_test_pos": int(y_test.sum()),
           "leaked_pos_windows": round(float(leak), 4),
           "auc": m["auc"], "mcc_nat": m["mcc"], "aupr": m["aupr"],
           "mcc_bal": mb["mcc"], "auc_bal": mb["auc"], "acc": m["acc"]}
    print(pd.Series(row).round(3).to_string(), "\n")
    return row


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    cv = pd.read_csv(DATA_PROC / "splits_cv.csv")
    sp = pd.read_csv(DATA_PROC / "split_species.csv")[["acc", "pos", "split_species"]]
    df = cv.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                  on=["acc", "pos"], validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    df = df.merge(sp, on=["acc", "pos"], how="left")
    X = encode_flat(df.window.tolist())
    ctx = context_features(df)
    print(f"{len(df)} rows / {int(df.label.sum())} positives\n")

    rows = []
    idx_all = np.arange(len(df))
    tr, te = train_test_split(idx_all, test_size=0.2, random_state=SEED,
                              stratify=df.label.values)
    rows.append(eval_protocol("site_random", df, tr, te, X, ctx))

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    tr, te = next(gss.split(idx_all, groups=df.acc.values))
    rows.append(eval_protocol("protein", df, tr, te, X, ctx))

    tr = np.where(df.fold.values != 0)[0]
    te = np.where(df.fold.values == 0)[0]
    rows.append(eval_protocol("cluster", df, tr, te, X, ctx))

    mask_tr = df.split_species.fillna("") == "train"
    bg_human = (df.split_species.fillna("") == "test") & (df.source == "background_cys")
    half = df.acc.map(half_of)
    mask_tr = mask_tr | (bg_human & (half == 0))
    mask_te = (df.split_species.fillna("") == "test") & ~(bg_human & (half == 0))
    rows.append(eval_protocol("species_human", df, np.where(mask_tr)[0], np.where(mask_te)[0], X, ctx))

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "leakage_ladder.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
