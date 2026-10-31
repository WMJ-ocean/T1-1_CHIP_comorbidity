# -*- coding: utf-8 -*-
"""
56_gene_mr.py  ——  遗留一「SMR + eQTL」的落地执行（v2，修正列名坑 + 全向量化）

目的：
  eQTLGen 只提供 Zscore 与 NrSamples，不含 beta/se。用官方换算得到近似的
  表达效应量与标准误，再与 FinnGen R12 的 MPN/PV/MDS/ITP/AA 结局做 Wald ratio，
  得到「每 1 个 SDU 基因表达变化的 MPN 风险 OR」。

换算公式（eQTLGen README / FAQ）：
    denom = sqrt( 2 * f * (1-f) * (n + Z^2) )
    beta  = Z / denom
    se    = 1 / denom
  f = **AssessedAllele 的频率**（不是随便一个等位），n = NrSamples。

坑位备忘：
  1. eQTLGen MAF 文件给的是 AlleleB 的频率。AssessedAllele == AlleleB 时 f = freqB，
     == AlleleA 时 f = 1 - freqB。
  2. harmonised 结局文件列名含点号（beta.outcome 等），pandas 的 itertuples 会把它们
     改名成 _18/_20 …… 所以禁止用 itertuples 取值，先 rename 成下划线再算。
  3. 每个基因在本设计里只有 1 个工具变量（即该位点的 CHIP lead SNP），因此只能做
     Wald ratio，无法做 HEIDI 异质性检验（那需要 LD 参考面板，本机暂无）。
"""
import os, gzip, math
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW  = os.path.join(BASE, "02_raw")
OUT  = os.path.join(BASE, "06_targets")

P = print
L = "=" * 78

P(L); P("Step 1  装载输入"); P(L)

iv = pd.read_csv(os.path.join(BASE, "05_iv", "CHIP_germline_iv.tsv"), sep="\t")
P(f"  胚系 CHIP IV        {len(iv)} 个")

tgt = pd.read_csv(os.path.join(OUT, "iv_eqtlgen_targets.tsv"), sep="\t")
P(f"  IV-基因 eQTL 记录   {len(tgt)} 行   "
  f"涉及 {tgt.SNP.nunique()} 个 SNP / {tgt.GeneSymbol.nunique()} 个基因")

# 结局列名统一改成下划线（关键：避免 itertuples 重命名坑）
RMAP = {"effect_allele.outcome": "out_eff", "other_allele.outcome": "out_oth",
        "beta.outcome": "beta_out", "se.outcome": "se_out",
        "pval.outcome": "p_out", "eaf.outcome": "eaf_out"}

outcomes = {}
for oc in ["MPN", "PV", "MDS", "ITP", "AA"]:
    p = os.path.join(BASE, "07_mr", f"CHIPg_to_{oc}_harmonised.csv")
    if not os.path.exists(p):
        P(f"  结局 {oc:4s}           缺失，跳过"); continue
    df = pd.read_csv(p)[["SNP"] + list(RMAP)].rename(columns=RMAP)
    df["outcome"] = oc
    outcomes[oc] = df
    P(f"  结局 {oc:4s}           {len(df)} SNP")

P(); P(L); P("Step 2  抽取这 21 个 SNP 的等位频率"); P(L)

want = set(iv.SNP)
rows, nline = [], 0
with gzip.open(os.path.join(RAW, "eqtlgen", "eQTLGen_SNP_AF_AlleleB.txt.gz"), "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    iS, iA, iB, iF = h.index("SNP"), h.index("AlleleA"), h.index("AlleleB"), h.index("AlleleB_all")
    for line in f:
        nline += 1
        fl = line.split("\t")
        if fl[iS] in want:
            rows.append((fl[iS], fl[iA], fl[iB], float(fl[iF])))
maf = pd.DataFrame(rows, columns=["SNP", "AlleleA", "AlleleB", "freqB"])
P(f"  扫描 {nline:,} 行，命中 {len(maf)} 个 SNP")
miss = want - set(maf.SNP)
if miss:
    P(f"  [警告] MAF 文件缺失: {sorted(miss)}")

P(); P(L); P("Step 3  补取 AssessedAllele 并换算 beta / se"); P(L)

# eQTLGen 主表按 (SNP, Gene) 定位 AssessedAllele / OtherAllele
need = set(zip(tgt.SNP, tgt.GeneSymbol))
al = {}
with gzip.open(os.path.join(RAW, "eqtlgen", "eQTLGen_cis_eQTL_FDR005.txt.gz"), "rt") as f:
    h = f.readline().rstrip("\n").split("\t")
    c = {n: i for i, n in enumerate(h)}
    for line in f:
        fl = line.rstrip("\n").split("\t")
        s = fl[c["SNP"]]
        if s in want:
            k = (s, fl[c["GeneSymbol"]])
            if k in need:
                al[k] = (fl[c["AssessedAllele"]], fl[c["OtherAllele"]])
P(f"  补到 {len(al)} 条等位信息")

d = tgt.merge(maf, on="SNP", how="left")
d["AssessedAllele"] = [al.get((s, g), (None, None))[0] for s, g in zip(d.SNP, d.GeneSymbol)]
d["OtherAllele"]    = [al.get((s, g), (None, None))[1] for s, g in zip(d.SNP, d.GeneSymbol)]

# f = AssessedAllele 频率
same_B = d.AssessedAllele == d.AlleleB
same_A = d.AssessedAllele == d.AlleleA
d["f"] = None
d.loc[same_B, "f"] = d.loc[same_B, "freqB"]
d.loc[same_A, "f"] = 1.0 - d.loc[same_A, "freqB"]
d["f"] = pd.to_numeric(d.f)
P(f"  可定频率 {d.f.notna().sum()} / {len(d)} 行   "
  f"(匹配 AlleleB {int(same_B.sum())} 条，AlleleA {int(same_A.sum())} 条)")

n  = pd.to_numeric(d.NrSamples)
z  = pd.to_numeric(d.Z)
den = (2.0 * d.f * (1.0 - d.f) * (n + z ** 2)) ** 0.5
d["beta_exp_raw"] = z / den
d["se_exp"]       = 1.0 / den

# 内部一致性校验：按构造 |beta/se| 应严格等于 |Z|
chk = (d.beta_exp_raw / d.se_exp - z).abs().max()
P(f"  一致性自检 |beta/se - Z| 最大偏差 = {chk:.2e}  " + ("OK" if chk < 1e-6 else "!! 异常"))
P(f"  |beta_exp| 中位数 {d.beta_exp_raw.abs().median():.4f}   "
  f"范围 {d.beta_exp_raw.abs().min():.4f} – {d.beta_exp_raw.abs().max():.4f}")
P(f"  se_exp   中位数 {d.se_exp.median():.4f}")

P(); P(L); P("Step 4  Wald ratio：基因表达 -> 各血液病结局"); P(L)

allres = []
for oc, od in outcomes.items():
    m = d.merge(od, on="SNP", how="inner")
    if not len(m):
        continue
    # 等位方向对齐：AssessedAllele 应与 FinnGen effect allele 一致，否则翻转表达效应
    flip = (m.AssessedAllele != m.out_eff).astype(int) * 2 - 1   # 一致 +1 / 不一致 -1
    m["concordance"] = ["same" if v > 0 else "flipped" for v in flip]
    m["beta_exp"] = m.beta_exp_raw * flip
    m["outcome"] = oc

    bx, sx, by, sy = m.beta_exp, m.se_exp, m.beta_out, m.se_out
    m["wald_b"] = by / bx
    m["wald_se"] = (sy / bx).abs()          # Delta 法一阶近似
    m["wald_z"] = m.wald_b / m.wald_se
    m["wald_p"] = [math.erfc(abs(v) / math.sqrt(2)) if pd.notna(v) else float("nan")
                   for v in m.wald_z]
    m["OR_perSDU"] = [math.exp(v) if pd.notna(v) else float("nan") for v in m.wald_b]
    allres.append(m)

R = pd.concat(allres, ignore_index=True)
P(f"  生成 {len(R)} 条检验  "
  f"({R.GeneSymbol.nunique()} 基因 × {R.outcome.nunique()} 结局)")


def stars(p):
    if pd.isna(p): return "na"
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


R["sig"] = R.wald_p.apply(stars)

best = (R.assign(_z=R.Z.abs())
          .sort_values("_z", ascending=False)
          .groupby(["GeneSymbol", "outcome"], as_index=False).head(1)
          .drop(columns="_z")
          .sort_values(["outcome", "wald_p"]))

COLS = ["outcome", "GeneSymbol", "SNP", "Ensembl", "out_eff", "out_oth",
        "AssessedAllele", "OtherAllele", "concordance", "f",
        "Z", "NrSamples", "beta_exp", "se_exp",
        "beta_out", "se_out", "p_out", "wald_b", "wald_se", "wald_p", "OR_perSDU", "sig"]

allp  = os.path.join(OUT, "gene_mr_wald_all.tsv")
bestp = os.path.join(OUT, "gene_mr_wald_best.tsv")
R[COLS].to_csv(allp, sep="\t", index=False)
best[COLS].to_csv(bestp, sep="\t", index=False)
P(f"  全表  -> gene_mr_wald_all.tsv   ({len(R)} 行)")
P(f"  精简  -> gene_mr_wald_best.tsv  ({len(best)} 行，每基因 x 结局取最强 eQTL)")

P(); P(L); P("Step 5  结果预览：基因表达 -> MPN / PV"); P(L)
for oc in ["MPN", "PV"]:
    sub = best[(best.outcome == oc) & best.wald_p.notna()].sort_values("wald_p")
    if not len(sub):
        continue
    P(f"\n--- {oc}   ({len(sub)} 个基因) ---")
    P(f"{'gene':<11}{'SNP':<13}{'Z_exp':>8}{'b_exp':>8}{'b_out':>8}"
      f"{'OR/SDU':>8}{'P_Wald':>11}  sig")
    for r in sub.itertuples():
        if pd.isna(r.wald_b):
            continue
        P(f"{r.GeneSymbol:<11}{r.SNP:<13}{r.Z:>8.1f}{r.beta_exp:>8.3f}{r.beta_out:>8.3f}"
          f"{r.OR_perSDU:>8.3f}{r.wald_p:>11.2e}  {r.sig}")

P(); P("[完成]")
