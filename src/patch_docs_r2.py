"""Update README, Chinese notes and server README for the second-round numbers."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

r = (ROOT / "README.md")
t = r.read_text(encoding="utf-8")
reps = {
    "| Cluster-held-out 5-fold CV, all species | 4,424 | **0.944** (0.940–0.948) | **0.611** | 0.766 |":
        "| Cluster-held-out 5-fold CV, all species | 4,424 | **0.938** (cluster CI 0.917–0.954) | **0.590** | 0.747 |",
    "| Human test folds (multi-species training) | 557 | **0.914** (0.900–0.927) | 0.314 | 0.689 |":
        "| Human test folds (multi-species training) | 557 | **0.910** (0.890–0.929) | 0.297 | 0.662 |",
    "| Human-only training, human test (ESM-2 view) | 557 | 0.892 | 0.286 | 0.609 |":
        "| Human-only training, human test (nested) | 557 | 0.901 | 0.280 | 0.661 |",
    "| Hard negatives only (positives + same-protein Cys) | 4,424 | 0.956 | 0.678 | 0.775 |":
        "| Hard negatives only (positives + same-protein Cys) | 4,424 | 0.952 | 0.664 | 0.762 |",
    "PalmDL reaches 0.914, a marginal advantage whose 95% CI touches 0.900 (P(AUC ≤ 0.900) = 0.025).":
        "PalmDL reaches 0.910; this is context, not a paired comparison (within-study P(AUC ≤ 0.900) = 0.08).",
    "Ensemble gains over either single view are +0.013 AUC (paired bootstrap, CI excludes zero).":
        "Ensemble gains are +0.008 over ESM-2 (cluster CI 0.002–0.014) and +0.007 over the forest (cluster CI −0.004–0.018).",
    "| within species (cluster-disjoint 20% holdout) | 0.925 |": "| within species (cluster-disjoint 20% holdout) | 0.900 |",
    "| mammal → mammal | 0.898 |": "| mammal → mammal | 0.896 (range 0.834–0.948) |",
    "| mammal → *A. thaliana* | 0.974 |": "| mammal → *A. thaliana* | 0.978 |",
    "| *A. thaliana* → mammal | 0.774 |": "| *A. thaliana* → mammal | 0.765 (range 0.739–0.781) |",
    "* FDR-controlled operating points, thresholds chosen inside each fold: 1% → precision 0.953 /\n  recall 20.2%; 5% → 0.903 / 34.5%; 10% → 0.855 / 40.4%\n  (`results/error_analysis/fdr_operating_points.csv`, Figure 8b).":
        "* Empirical validation-target operating points (thresholds chosen inside each fold), with the\n  observed held-out false-discovery proportion: 1% → precision 0.953 / recall 0.201 / FDP 4.7%;\n  5% → 0.900 / 0.331 / 10.0%; 10% → 0.849 / 0.398 / 15.1%; 20% → 0.759 / 0.490 / 24.1%\n  (`results/error_analysis/fdr_operating_points.csv`, Figure 8b). These are not a demonstrated\n  FDR-controlling procedure.",
}
for a, b in reps.items():
    if a in t:
        t = t.replace(a, b)
    else:
        print("README MISS:", a[:60])
r.write_text(t, encoding="utf-8")

c = ROOT / "中文说明.md"
u = c.read_text(encoding="utf-8")
creps = {
    "**标签修正后的真实结果（v2）**：聚类留出 CV AUC = **0.944**（0.940–0.948），人类测试位点 AUC = **0.914**（0.900–0.927）——比已发表最佳 0.900 略高，但置信区间正好触及 0.900，因此表述为\"持平、略占优\"，不写\"显著超过\"。":
        "**标签修正 + 嵌套选型后的真实结果（v2.1）**：聚类留出 CV AUC = **0.938**（簇 bootstrap 0.917–0.954），人类测试位点 AUC = **0.910**（0.890–0.929）。与已发表 0.900 只能描述性并列（不同测试集），不做优效性声明。",
    "| 聚类留出 5 折（4,424 正位点） | **0.944**（0.940–0.948） | 0.611 | 0.766 |":
        "| 聚类留出 5 折（4,424 正位点，嵌套选型） | **0.938**（0.917–0.954） | 0.590 | 0.747 |",
    "| 人类测试折（多物种训练） | **0.914**（0.900–0.927） | 0.314 | 0.689 |":
        "| 人类测试折（多物种训练） | **0.910**（0.890–0.929） | 0.297 | 0.662 |",
    "| 仅人类数据训练（验证集选中 ESM 单视角） | 0.892 | 0.286 | 0.609 |":
        "| 仅人类数据训练（嵌套） | 0.901 | 0.280 | 0.661 |",
    "| 硬负样本基准 | **0.956** | 0.678 | 0.775 |": "| 硬负样本基准 | **0.952** | 0.664 | 0.762 |",
    "集成相对任一单视角提升 ΔAUC +0.013（配对 bootstrap，CI 不含 0）。":
        "集成相对 ESM 提升 +0.008（CI 0.002–0.014），相对森林 +0.007（CI −0.004–0.018，未分辨）。",
    "**我们的 0.914 在零共享窗口、聚类留出条件下取得，但 P(AUC≤0.900)=0.025、CI 下界正好 0.900，措辞应为\"略高于、统计上边缘\"。**":
        "**我们的 0.910 在零共享窗口、聚类留出条件下取得；这只是描述性并列，不做优效性推断。**",
}
for a, b in creps.items():
    if a in u:
        u = u.replace(a, b)
    else:
        print("CN MISS:", a[:60])
c.write_text(u, encoding="utf-8")

s = ROOT / "server" / "README.md"
v = s.read_text(encoding="utf-8")
v = v.replace("## Notes\n", "## Notes\n\n- The cross-validated model used within-fold **rank** averaging, while this service averages view\n  **logits** (needed for single-sequence scoring). On the same folds the logit average scores AUC\n  0.940 versus 0.944 for the rank average under the global combination, so service probabilities are\n  slightly less sharp than the cross-validated ensemble.\n")
s.write_text(v, encoding="utf-8")
print("docs patched")
