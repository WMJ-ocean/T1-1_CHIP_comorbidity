# -*- coding: utf-8 -*-
"""
58_gene_priority.py  ——  把 Wald ratio 结果重做成「候选基因优先级表」（假设生成）

背景（重要，别跳过）：
  单工具变量下 Wald ratio 的 P 值在数学上恒等于结局 GWAS 的 P 值：
      z_wald = (b_out / b_exp) / (se_out / |b_exp|) = b_out / se_out = z_outcome
  本脚本已实测该恒等式（偏差 8.9e-16）。因此 gene-level MR 的 P 值
  **不构成任何新证据**，把它写进论文等于放一张伪证据表。

  唯一有信息量的是效应尺度：
      OR per SDU = exp( b_out / b_exp )
  它回答「表达量每变化 1 个标准差，疾病 OR 变化多少」，可用于跨基因比较
  与后续实验的功效估算（要检测该效应需要多大样本）。

  另一个诚实结论：一个非编码位点对顺式窗口里多个基因都有 eQTL
  （本项目平均每个 IV 命中 5.3 个基因），**仅靠汇总数据无法指出究竟是哪个基因**。
  必须靠：(a) 多工具变量 SMR/HEIDI，(b) 全窗口 coloc，(c) 单细胞/实验数据。
"""
import os, math
import pandas as pd
from scipy.stats import norm

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(BASE, "06_targets")

P = print
L = "=" * 80

src = os.path.join(OUT, "gene_mr_wald_all.tsv")
d = pd.read_csv(src, sep="\t")
P(L); P("载入"); P(L)
P(f"  gene_mr_wald_all.tsv   {len(d)} 行  ({d.GeneSymbol.nunique()} 基因 × {d.outcome.nunique()} 结局)")

# 挂上 CHIP exposure 侧的效应量，便于算总效应 c
exp = pd.read_csv(os.path.join(BASE, "07_mr", "CHIPg_to_MPN_harmonised.csv"))
exp = exp[["SNP", "beta.exposure", "se.exposure", "pval.exposure",
           "effect_allele.exposure", "eaf.exposure"]].rename(columns={
    "beta.exposure": "b_chip", "se.exposure": "se_chip",
    "pval.exposure": "p_chip", "effect_allele.exposure": "chip_allele",
    "eaf.exposure": "chip_eaf"})

d = d.merge(exp, on="SNP", how="left")

# 每个位点的候选基因数
per_snp = d[d.outcome == "MPN"].groupby("SNP").GeneSymbol.nunique()
P(f"  每个 IV 命中的顺式候选基因中位数 = {per_snp.median():.1f}  "
  f"范围 {per_snp.min()}-{per_snp.max()}")

# ---- 效应尺度 ----
d["OR_perSDU"] = [math.exp(v) if pd.notna(v) and abs(v) < 50 else float("nan")
                  for v in d.wald_b]
# 要在 80% 功效、alpha=0.05 下检出该 per-SDU 效应所需样本
# n = (1.96 + 0.84)^2 / b^2  … 二元结局近似
d["n80徒刑"] = None


def n_for_80(b):
    if pd.isna(b) or b == 0:
        return float("nan")
    return ((norm.ppf(0.975) + norm.ppf(0.8)) ** 2) / (b ** 2)


d["n_need_80pct"] = d.wald_b.apply(n_for_80)

# ---- 方向链条一致性 ----
# a  = CHIP 风险等位对表达的效应   (= beta_exp，已翻转到 chip-risk 等位方向)
# c  = CHIP 风险等位对疾病的效应   (= beta_out)
# b  = 由 c/a 推得的 表达->疾病效应 (= wald_b)
d["chain"] = [
    ("CHIPrisk-allele -> expr " + ("up" if a > 0 else "down") + "  &  expr -> dx "
     + ("up" if b > 0 else "down"))
    if pd.notna(a) and pd.notna(b) and a != 0 else "na"
    for a, b in zip(d.beta_exp, d.wald_b)
]

COLS = ["outcome", "SNP", "GeneSymbol", "Ensembl", "chip_allele", "AssessedAllele",
        "concordance", "p_chip", "b_chip", "Z", "NrSamples", "f",
        "beta_exp", "se_exp", "beta_out", "se_out", "p_out",
        "wald_b", "OR_perSDU", "chain"]

ALL = d[COLS].copy()

# ---- 优先级：按 eVQTL 强度 |Z| 排序，每结局输出候选表 ----
prior = (ALL.assign(aZ=ALL.Z.abs())
            .sort_values("aZ", ascending=False)
            .drop(columns="aZ"))
ALL.to_csv(os.path.join(OUT, "gene_priority_all.tsv"), sep="\t", index=False)
prior.to_csv(os.path.join(OUT, "gene_priority_ranked.tsv"), sep="\t", index=False)
P(f"  -> gene_priority_all.tsv    ({len(ALL)} 行)")
P(f"  -> gene_priority_ranked.tsv (按 |eQTL Z| 排序)")

P(); P(L); P("MPN：候选基因（按 cis-eQTL 强度排序，取每 SNP 前 2 个基因）"); P(L)
m = prior[prior.outcome == "MPN"].copy()
top = m.sort_values(m.Z.name if False else "Z", key=lambda s: s.abs(), ascending=False) \
        .groupby("SNP").head(2).sort_values("Z", key=lambda s: s.abs(), ascending=False)
P(f"{'gene':<13}{'SNP':<12}{'eQTL Z':>8}{'OR/SDU':>9}{'CHIP P':>11}"
  f"{'CHIP b':>8}{'方向链'}")
for r in top.itertuples():
    P(f"{r.GeneSymbol:<13}{r.SNP:<12}{r.Z:>8.1f}{r.OR_perSDU:>9.2f}{r.p_chip:>11.1e}"
      f"{r.b_chip:>8.3f}  {r.chain}")

P(); P(L); P("PV：候选基因（同上）"); P(L)
v = prior[prior.outcome == "PV"].copy()
topv = v.sort_values("Z", key=lambda s: s.abs(), ascending=False) \
        .groupby("SNP").head(2).sort_values("Z", key=lambda s: s.abs(), ascending=False)
P(f"{'gene':<13}{'SNP':<12}{'eQTL Z':>8}{'OR/SDU':>9}{'CHIP P':>11}{'CHIP b':>8}  方向链")
for r in topv.itertuples():
    P(f"{r.GeneSymbol:<13}{r.SNP:<12}{r.Z:>8.1f}{r.OR_perSDU:>9.2f}{r.p_chip:>11.1e}"
      f"{r.b_chip:>8.3f}  {r.chain}")

P(); P("[完成]  提醒：本表 P 值不可用（= 结局 GWAS P），仅供假设生成与功效估算")
