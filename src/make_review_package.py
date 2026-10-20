"""Build the reviewer package (full and light versions).

Stage layout:
  00_README_审阅说明.md
  01_manuscript/      manuscript DOCX/PDF/markdown
  02_supplementary/   supplementary DOCX/PDF
  03_figures/         9 PNG figures + legends
  04_data/            raw UniProt snapshot + processed tables + data dictionary
  05_results/         summary tables, OOF predictions, scan scores, statistics
  06_manifest/        MANIFEST.csv + checksums.sha256

Usage:
  python src/make_review_package.py --stage     # build the tree (PDFs converted separately)
  python src/make_review_package.py --zip       # checksums + full/light zips
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / "review"
STAGE = REL / f"PalmDL_review_{date.today().isoformat()}"

LIGHT_DIRS = ("00_README_审阅说明.md", "01_manuscript", "02_supplementary",
              "03_figures", "04_data/processed", "04_data/DATA_DICTIONARY.md",
              "05_results/tables", "05_results/leakage_ladder.csv",
              "05_results/species_matrix_auc.csv", "05_results/headline_stats.json",
              "05_results/stats_vs_published.json", "05_results/scan/fdr_summary.csv",
              "05_results/structure/paired_cv.csv", "06_manifest")

README = """# PalmDL 审阅包 / Reviewer package

**日期 / Date:** {today}
**稿件 / Manuscript:** PalmDL: a leakage-aware deep learning framework and benchmark for protein S-palmitoylation site prediction
**代码 / Code:** https://github.com/opus940810/palm-dl （v2，最新 commit）

> 说明：Zenodo 上的 v1.0.0 存档是标签修正前的版本，修正版正在重新存档；本次审阅请以本包内容和 GitHub 上的 v2 为准。

> 作者信息、单位与基金仍在填写中（稿中为 `[Author 1]` 等占位符），不影响审阅内容。

## 这个包里有什么 / Contents

| 目录 | 内容 |
|---|---|
| `01_manuscript/` | 主稿：DOCX（可编辑）、PDF（便于阅读）、Markdown 源稿 |
| `02_supplementary/` | 补充材料：Tables S1–S7（DOCX + PDF） |
| `03_figures/` | 9 张正文图（300 dpi PNG）+ 图注 `FIGURE_LEGENDS.md` |
| `04_data/raw/` | UniProt 2026_03 原始导出（脂化注释表 + 无脂化人蛋白表），带下载元数据 |
| `04_data/processed/` | 位点表、官方折分配、物种划分、蛋白表；字段说明见 `DATA_DICTIONARY.md` |
| `05_results/` | 交叉验证汇总、泄漏审计、物种矩阵、统计检验、全蛋白组 FDR 清单、结构分析、文献校订表；含 OOF 预测与全蛋白组打分 |
| `06_manifest/` | 文件清单 `MANIFEST.csv` 与校验和 `checksums.sha256` |

## 建议的审阅顺序 / Suggested order

1. `01_manuscript/PalmDL_manuscript.pdf` —— 先通读正文（约 25 页）
2. `02_supplementary/PalmDL_supplementary.pdf` —— 补充表 S1–S7
3. `03_figures/` —— 对照图注看图
4. 需要核对具体数字时：`05_results/tables/` 与 `05_results/*.json`
5. 需要复现时：`04_data/processed/splits_cv.csv` 是官方折文件，代码见 GitHub

## 核心数字（便于快速核对）/ Headline numbers

- 聚类留出 5 折，全物种（嵌套选型）：AUC **0.938**（簇 CI 0.917–0.954），MCC 0.590
- 人类测试位点：AUC **0.910**（0.890–0.929）；与已发表 0.900 仅描述性并列
- 严格物种留出：AUC 0.885（随机森林）
- 泄漏审计：随机位点级划分 41.1% 测试正样本窗口出现在训练集（AUC 0.980）；聚类划分 0%（AUC 0.944）
- 全蛋白组扫描（森林视角）：262,286 个 Cys，5% 档 301 个位点（53 个新位点），胞质富集 11.1 倍
- 跨物种矩阵：哺乳动物间 0.883–0.948（猪→人 0.834 例外）；植物→哺乳动物 0.739–0.781
- 结构特征：配对 5 折 AUC 0.890 → 0.911（探索性）

## 数据来源 / Provenance

- UniProt release **2026_03**（2026-10-18 检索），CC BY 4.0
- AlphaFold DB 模型（CC BY 4.0）——未随包分发，可由 `src/af_features.py` 重新下载
- Europe PMC 开放获取全文 992 篇——**未随包分发**（版权原因），仅保留抽取与校订结果（`05_results/literature/`）

## 限制 / Caveats（审阅时请留意）

- 90% 的正样本是"by similarity"推断注释；实验证据子集 AUC 0.907
- 无同批 head-to-head 对比（GPS-Palm 服务器已下线、独立版为 GUI）
- 文献挖掘集与 RFCM-PALM 集均未通过基序富集检查，故不作为基准使用
- 结构特征仅在 479 个蛋白上有模型，未并入最终发布模型

## 反馈 / Feedback

请在 DOCX 中用批注或修订模式；如对某个数字有疑问，`05_results/` 下有对应的逐折/汇总文件。
"""

DATA_DICT = """# 数据字典 / Data dictionary

所有文件均为 UTF-8、逗号分隔（CSV）；原始数据为制表符分隔（TSV）。

## 04_data/processed/sites.csv（核心位点表，67,266 行）

| 列 | 含义 |
|---|---|
| `acc` | UniProt 登录号 |
| `pos` | 位点在规范序列中的位置（1-based） |
| `label` | 1 = S-棕榈酰化位点；0 = 负样本 |
| `residue` | 该位置在序列中的残基（99.0% 为 C） |
| `source` | `uniprot_exp`（实验证据）/ `uniprot_inferred`（推断）/ `same_protein_cys`（同蛋白硬负样本）/ `background_cys`（背景负样本） |
| `evidence` | UniProt 证据码（ECO），推断位点多为 `ECO:0000250` |
| `compartment` | UniProt 拓扑注释：cytoplasmic / extracellular / lumenal / signal_peptide / other / unknown |
| `window` | 以该位点为中心的 31 残基窗口（两端用 X 补齐） |
| `protein_length` | 蛋白长度 |
| `organism` | NCBI 分类号 |
| `gene` | 基因名 |
| `is_human` | 是否人源 |

## splits_cv.csv（官方 5 折，聚类留出）

| 列 | 含义 |
|---|---|
| `acc`, `pos` | 位点 |
| `cluster` | 4-mer 余弦 ≥0.30 聚出的同源簇（同簇不跨折） |
| `fold` | 0–4，交叉验证折号 |

## split_species.csv（跨物种独立测试）

| 列 | 含义 |
|---|---|
| `split_species` | `train`（非人物种，且已剔除人物种同源簇与人源传播注释）/ `test`（人源）/ `drop`（因泄漏剔除） |

## proteins.csv

`acc, organism, protein, gene, length, n_sites`：带正位点的蛋白清单。

## 04_data/raw/

- `uniprot_lipid_reviewed.tsv`：Swiss-Prot 中含脂化特征的 9,089 条entry
- `uniprot_negatives.tsv`：19,622 条无任何脂化注释的人源 reviewed entry
- `download_meta.json`：检索语句、行数、UniProt 版本

## 05_results/ 下的关键表

| 文件 | 内容 |
|---|---|
| `tables/table1_cv.csv` | 全物种 5 折各模型指标 |
| `tables/table2_cv_human.csv` | 人类测试折 |
| `tables/table3_cv_human_only.csv` | 仅人类训练 |
| `tables/table4_species.csv` | 严格物种留出 |
| `tables/table5_hard.csv` | 硬负样本基准 |
| `tables/table6_leakage.csv` | 四种划分的泄漏审计 |
| `leakage_ladder.csv` | 同上（原始输出） |
| `species_matrix_auc.csv` | 8×8 跨物种 AUC 矩阵 |
| `headline_stats.json` | bootstrap 置信区间与配对检验 |
| `stats_vs_published.json` | 与 0.900 / 0.855 的 bootstrap 比较 |
| `scan/fdr_summary.csv` | 全蛋白组 5 个 FDR 档位 |
| `scan/shortlist_fdr5.csv` | 5% 档的 301 个候选位点 |
| `scan/all_scores.csv` | 262,286 个 Cys 的打分（全量） |
| `preds_cv_all.csv` | 全部位点的折外预测 |
| `structure/hybrid_summary.json` | 结构感知模型的配对 CV 结果 |
| `literature/curated_sites.csv` | 人工校订后的 99 个文献位点 |

## 复现

```bash
pip install -r requirements.txt
python src/download_data.py && python src/build_dataset.py && python src/splits.py
# 其余步骤见 docs/REPRODUCE.md
```
"""

FIGURE_LEGENDS = """# 图注 / Figure legends

**Figure 1. 数据构成。** (a) 贡献位点最多的物种（人/小鼠高亮）；(b) 证据层级（477 实验 vs 6,328 推断，对数坐标）；(c) 正位点的拓扑注释。

**Figure 2. 划分方式决定表观性能。** (a) 四种协议下的 AUC 与平衡 MCC；(b) 测试正样本中精确 31-aa 窗口出现在训练集的比例（随机/蛋白级 40–44%，聚类/物种 0%）。

**Figure 3. 聚类留出交叉验证的判别能力。** (a) ROC；(b) 精确率–召回，含 PWM、逻辑回归、随机森林、CNN、ESM-2 与 PalmDL 集成。

**Figure 4. 人类位点在三种评估机制下的表现。** (a) AUC；(b) 平衡 MCC；多物种训练 / 仅人类训练 / 物种留出。

**Figure 5. 校准。** 自然流行率下的可靠性曲线（RF、CNN、ESM-2、校准后的集成）。

**Figure 6. 学到的基序。** 以半胱氨酸为中心 −7…+8 的氨基酸位置特异性 log-odds（正样本 vs 背景负样本）。

**Figure 7. 棕榈酰化"密码"是真核共有的。** 单物种模型的 8×8 迁移矩阵（行=训练物种，列=测试物种）；对角线来自同物种 20% 聚类留出。

**Figure 8. 模型依赖什么、预测能信到什么程度。** (a) 300 个人类正位点的饱和突变效应；(b) FDR 控制下的精确率–召回操作点。

**Figure 9. 全蛋白质组扫描与 FDR 清单。** (a) 262,286 个人类 Cys 与 20,000 个乱序诱饵的打分分布；(b) 估计 FDR 随清单大小的变化；(c) 各档位中已知/新位点的构成。
"""


def copytree(rel: str) -> int:
    src, dst = ROOT / rel, STAGE / rel
    if not src.exists():
        print("  missing:", rel)
        return 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
        return sum(1 for p in dst.rglob("*") if p.is_file())
    shutil.copy2(src, dst)
    return 1


def stage() -> None:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True, exist_ok=True)
    (STAGE / "00_README_审阅说明.md").write_text(README.format(today=date.today().isoformat()),
                                                 encoding="utf-8")
    files = [
        ("manuscript/PalmDL_manuscript.docx", "01_manuscript"),
        ("manuscript/manuscript.md", "01_manuscript"),
        ("manuscript/PalmDL_supplementary.docx", "02_supplementary"),
        ("data/raw/uniprot_lipid_reviewed.tsv", "04_data/raw"),
        ("data/raw/uniprot_negatives.tsv", "04_data/raw"),
        ("data/raw/download_meta.json", "04_data/raw"),
        ("results/scan/fdr_summary.csv", "05_results/scan"),
        ("results/scan/compartment_enrichment.csv", "05_results/scan"),
        ("results/scan/scan_summary.json", "05_results/scan"),
        ("results/scan/shortlist_fdr1.csv", "05_results/scan"),
        ("results/scan/shortlist_fdr5.csv", "05_results/scan"),
        ("results/scan/shortlist_fdr20.csv", "05_results/scan"),
        ("results/scan/all_scores.csv", "05_results/scan"),
        ("results/structure/hybrid_summary.json", "05_results/structure"),
        ("results/structure/hybrid_cv.csv", "05_results/structure"),
        ("results/structure/thirdview_summary.json", "05_results/structure"),
        ("results/structure/structure_summary.json", "05_results/structure"),
        ("results/literature/curated_sites.csv", "05_results/literature"),
        ("results/literature/filter_table.csv", "05_results/literature"),
        ("results/literature/metrics.json", "05_results/literature"),
        ("results/leakage_ladder.csv", "05_results"),
        ("results/species_matrix_auc.csv", "05_results"),
        ("results/headline_stats.json", "05_results"),
        ("results/stats_vs_published.json", "05_results"),
        ("results/preds_cv_all.csv", "05_results"),
    ]
    dirs = [
        ("figures", "03_figures"),
        ("data/processed", "04_data/processed"),
        ("results/tables", "05_results/tables"),
    ]
    total = 0
    for src_rel, dst_rel in files:
        src = ROOT / src_rel
        if not src.exists():
            print("  missing:", src_rel)
            continue
        (STAGE / dst_rel).mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, STAGE / dst_rel / src.name)
        total += 1
    for src_rel, dst_rel in dirs:
        src = ROOT / src_rel
        if not src.exists():
            print("  missing:", src_rel)
            continue
        shutil.copytree(src, STAGE / dst_rel, dirs_exist_ok=True)
        total += sum(1 for p in (STAGE / dst_rel).rglob("*") if p.is_file())
    (STAGE / "04_data" / "DATA_DICTIONARY.md").write_text(DATA_DICT, encoding="utf-8")
    (STAGE / "03_figures" / "FIGURE_LEGENDS.md").write_text(FIGURE_LEGENDS, encoding="utf-8")
    print(f"staged {total} files -> {STAGE}")


def zip_package(src: Path, out: Path) -> float:
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(src.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(src.parent).as_posix())
    return out.stat().st_size / 1e6


def finalize() -> None:
    manifest = STAGE / "06_manifest"
    manifest.mkdir(exist_ok=True)
    rows = []
    for p in sorted(STAGE.rglob("*")):
        if p.is_file() and "06_manifest" not in p.parts:
            rel = p.relative_to(STAGE).as_posix()
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            rows.append((rel, p.stat().st_size, h))
    (manifest / "MANIFEST.csv").write_text(
        "path,bytes,sha256\n" + "\n".join(f"{r},{s},{h}" for r, s, h in rows) + "\n", encoding="utf-8")
    (manifest / "checksums.sha256").write_text(
        "\n".join(f"{h}  {r}" for r, s, h in rows) + "\n", encoding="utf-8")

    full = REL / f"PalmDL_review_FULL_{date.today().isoformat()}.zip"
    size = zip_package(STAGE, full)
    print(f"full package: {full.name} ({size:.1f} MB, {len(rows)} files)")

    light = REL / f"PalmDL_review_LIGHT_{date.today().isoformat()}.zip"
    with zipfile.ZipFile(light, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(STAGE.rglob("*")):
            if not p.is_file():
                continue
            rel = p.relative_to(STAGE).as_posix()
            if any(rel == d or rel.startswith(d + "/") for d in LIGHT_DIRS):
                z.write(p, p.relative_to(STAGE.parent).as_posix())
    print(f"light package: {light.name} ({light.stat().st_size / 1e6:.1f} MB)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["stage", "zip"])
    args = ap.parse_args()
    if args.mode == "stage":
        stage()
    else:
        finalize()


if __name__ == "__main__":
    main()
