"""Models: deep CNN, ESM embedding head, and classical baselines."""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef

from config import SEED


def set_seed(seed: int = SEED) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)


class PalmCNN(nn.Module):
    """Multi-kernel 1-D CNN over the 31-aa window plus protein context."""

    def __init__(self, in_ch: int = 23, ctx_dim: int = 11, width: int = 32, dropout: float = 0.3):
        super().__init__()
        kernels = (3, 5, 7, 11)
        self.convs = nn.ModuleList([nn.Conv1d(in_ch, width, k, padding=k // 2) for k in kernels])
        self.bn = nn.BatchNorm1d(width * len(kernels))
        self.fc1 = nn.Linear(width * len(kernels) + ctx_dim, 64)
        self.drop = nn.Dropout(dropout)
        self.fc2 = nn.Linear(64, 1)

    def forward(self, x: torch.Tensor, ctx: torch.Tensor) -> torch.Tensor:
        hs = [F.relu(c(x)) for c in self.convs]
        h = torch.cat([F.adaptive_max_pool1d(v, 1).squeeze(-1) for v in hs], dim=1)
        h = self.bn(h)
        h = self.drop(F.relu(self.fc1(torch.cat([h, ctx], dim=1))))
        return self.fc2(h).squeeze(-1)


class EmbeddingMLP(nn.Module):
    """Small head on top of frozen protein-language-model embeddings."""

    def __init__(self, in_dim: int, ctx_dim: int = 11, hidden: int = 256, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim + ctx_dim, hidden), nn.BatchNorm1d(hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden // 2), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden // 2, 1),
        )

    def forward(self, x: torch.Tensor, ctx: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([x, ctx], dim=1)).squeeze(-1)


def train_torch(model, Xtr, Ctr, ytr, Xva, Cva, yva, *, epochs=50, bs=128, lr=1e-3,
                weight_decay=1e-4, patience=8, verbose=False):
    """Train with early stopping on validation MCC; return model, best info, history."""
    from metrics import find_threshold

    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    n_pos = float(ytr.sum())
    pos_weight = torch.tensor([(len(ytr) - n_pos) / max(n_pos, 1.0)], dtype=torch.float32)
    lossf = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    Xtr_t, Ctr_t, ytr_t = map(torch.as_tensor, (Xtr, Ctr, ytr))
    Xva_t, Cva_t = map(torch.as_tensor, (Xva, Cva))
    ds = torch.utils.data.TensorDataset(Xtr_t, Ctr_t, ytr_t)
    dl = torch.utils.data.DataLoader(ds, batch_size=bs, shuffle=True, drop_last=len(ds) > bs)
    best = {"score": -2.0, "state": None, "epoch": -1, "thr": 0.5}
    hist, bad = [], 0
    for ep in range(epochs):
        model.train()
        tot = 0.0
        for xb, cb, yb in dl:
            opt.zero_grad()
            loss = lossf(model(xb, cb), yb)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            tot += loss.item() * len(yb)
        model.eval()
        with torch.no_grad():
            p = torch.sigmoid(model(Xva_t, Cva_t)).numpy()
        thr = find_threshold(yva, p)
        try:
            from sklearn.metrics import roc_auc_score
            m = float(roc_auc_score(yva, p))  # prevalence-independent model selection
        except ValueError:
            m = float(matthews_corrcoef(yva, p >= thr))
        hist.append({"epoch": ep, "loss": tot / len(ds), "val_auc": m, "val_thr": thr})
        if verbose:
            print(f"    ep{ep:02d} loss={tot / len(ds):.4f} val_auc={m:.4f}")
        if m > best["score"] + 1e-4:
            best = {"score": m, "state": {k2: v.clone() for k2, v in model.state_dict().items()},
                    "epoch": ep, "thr": thr}
            bad = 0
        else:
            bad += 1
            if bad >= patience:
                break
    model.load_state_dict(best["state"])
    return model, best, hist


def predict_torch(model, X, C, bs=1024):
    model.eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(X), bs):
            out.append(torch.sigmoid(model(torch.as_tensor(X[i:i + bs]),
                                           torch.as_tensor(C[i:i + bs]))).numpy())
    return np.concatenate(out)


def make_logreg(C: float = 1.0):
    # n_jobs=1: joblib's process backend uses named pipes, which the sandbox blocks
    return LogisticRegression(C=C, max_iter=3000, class_weight="balanced", n_jobs=1)


def make_rf(n_estimators: int = 400):
    return RandomForestClassifier(n_estimators=n_estimators, min_samples_leaf=3,
                                  class_weight="balanced_subsample", n_jobs=1, random_state=SEED)
