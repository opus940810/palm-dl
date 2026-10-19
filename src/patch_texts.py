"""Patch README strings in the GitHub-bundle and review-package generators."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

gh = ROOT / "src" / "make_github_bundle.py"
t = gh.read_text(encoding="utf-8")
old = ("Headline numbers: cluster-held-out AUC **0.956** (95% CI 0.953–0.959); human test\n"
       "sites AUC **0.901** (0.885–0.915), statistically indistinguishable from the best\n"
       "published human value 0.900 (bootstrap P = 0.46); strict species holdout AUC 0.822.")
new = ("Headline numbers: cluster-held-out AUC **0.944** (95% CI 0.940–0.948); human test\n"
       "sites AUC **0.914** (0.900–0.927), marginally above the best published human value\n"
       "0.900 (bootstrap P = 0.025); strict species holdout AUC 0.885.")
assert old in t
gh.write_text(t.replace(old, new), encoding="utf-8")

rp = ROOT / "src" / "make_review_package.py"
t = rp.read_text(encoding="utf-8")
reps = {
    "| `04_data/processed/sites.csv`（核心位点表，71,594 行）":
        "| `04_data/processed/sites.csv`（核心位点表，67,266 行）",
    "- 聚类留出 5 折，全物种：AUC **0.956**（95% CI 0.953–0.959），MCC 0.753":
        "- 聚类留出 5 折，全物种：AUC **0.944**（95% CI 0.940–0.948），MCC 0.611",
    "- 人类测试位点：AUC **0.901**（0.885–0.915），与已发表最佳 0.900 **统计不可区分**（bootstrap P = 0.46）":
        "- 人类测试位点：AUC **0.914**（0.900–0.927），略高于已发表最佳 0.900（bootstrap P = 0.025）",
    "- 严格物种留出：AUC 0.822":
        "- 严格物种留出：AUC 0.885",
    "- 泄漏审计：随机位点级划分 43.7% 测试正样本窗口出现在训练集（AUC 0.990）；聚类划分 0%（AUC 0.950）":
        "- 泄漏审计：随机位点级划分 41.1% 测试正样本窗口出现在训练集（AUC 0.980）；聚类划分 0%（AUC 0.944）",
    "- 全蛋白组扫描：262,286 个 Cys，5% FDR 时 322 个高置信位点（227 个新位点），胞质拓扑富集 16.3 倍":
        "- 全蛋白组扫描（森林视角）：262,286 个 Cys，5% 档 301 个位点（53 个新位点），胞质富集 11.1 倍",
    "- 跨物种矩阵：哺乳动物间 AUC 0.94；细菌↔哺乳动物 0.55":
        "- 跨物种矩阵：哺乳动物间 AUC 0.90；植物→哺乳动物 0.77（细菌已移出）",
    "- 结构特征：配对 5 折 AUC 0.890 → 0.910":
        "- 结构特征：配对 5 折 AUC 0.890 → 0.911",
    "- 93% 的正样本是\"by similarity\"推断注释；实验证据子集 AUC 0.908":
        "- 90% 的正样本是\"by similarity\"推断注释；实验证据子集 AUC 0.907",
}
for a, b in reps.items():
    if a in t:
        t = t.replace(a, b)
    else:
        print("MISS:", a[:50])
rp.write_text(t, encoding="utf-8")
print("generators patched")
