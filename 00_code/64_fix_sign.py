# -*- coding: utf-8 -*-
"""
64_fix_sign.py

修正 56_gene_mr.py 中等位翻转方向的系统性错误。

背景：
  flip = (AssessedAllele != out_eff).astype(int) * 2 - 1
  这条写反了 —— 当 AssessedAllele == out_eff（一致、不该翻）时算出 -1，
  不一致（该翻）时算出 +1。交叉表显示 100% 的行都在反对角线上，
  也就是说每一行都恰好相差因子 -1，因此可以精确还原。

校正内容：
  beta_exp      -> -beta_exp
  wald_b        -> -wald_b          (wald_se 不变)
  OR_perSDU     -> 1/OR_perSDU
  concordance   -> 'same' <-> 'flipped'
  wald_p / sig  不变（P 值不受符号影响；且单 IV 时 wald_p 恒等于 p_out）

另外新增：
  beta.exposure  来自 harmonised 文件（该等位对 CHIP 风险的效应，可能为负）
  exp_dir_chip   在「CHIP 易感性升高」方向上，表达是被推高(+)还是压低(-)
"""
import os, io
import pandas as pd, numpy as np

BASE = r'C:\Users\Hero\WorkBuddy\2026-10-14-23-35-01\outputs\pv_hotspot\T1-1_CHIP_comorbidity'
P = print

src = os.path.join(BASE, '06_targets', 'gene_mr_wald_all.tsv')
d = pd.read_csv(src, sep='\t')
P(f'读入 {len(d)} 行 × {d.shape[1]} 列')

# ---------- 1. 取 exposure beta ----------
hmap = {}
for oc in d.outcome.unique():
    f = os.path.join(BASE, '07_mr', f'CHIPg_to_{oc}_harmonised.csv')
    if not os.path.exists(f):
        P(f'  (缺 {f})'); continue
    h = pd.read_csv(f)
    for snp, b in zip(h.SNP, h['beta.exposure']):
        hmap[snp] = b
P(f'取到 {len(hmap)} 个 SNP 的 beta.exposure')

if 'beta_exp' in d.columns and not d.attrs.get('fixed'):
    pass

d['beta_chip'] = d.SNP.map(hmap)

# ---------- 2. 统一乘 -1 还原 ----------
for c in ['beta_exp', 'wald_b']:
    if c in d.columns:
        d[c] = -d[c]
if 'OR_perSDU' in d.columns:
    d['OR_perSDU'] = 1.0 / d['OR_perSDU']
if 'concordance' in d.columns:
    d['concordance'] = d['concordance'].map({'same': 'flipped', 'flipped': 'same'})

# ---------- 3. 校验：|beta_exp/se_exp| 应等于 |Z| ----------
chk = (d.beta_exp / d.se_exp - d.Z).abs().max()
P(f'校验 |beta_exp/se_exp - Z| 最大值 = {chk:.3e}  (应接近 0)')

# ---------- 4. CHIP 易感性升高方向上的表达方向 ----------
d['exp_dir_chip'] = np.sign(d.beta_chip) * np.sign(d.beta_exp)
d['exp_dir_chip'] = d['exp_dir_chip'].map({1.0: 'up', -1.0: 'down', 0.0: 'na'})

out = os.path.join(BASE, '06_targets', 'gene_mr_wald_corrected.tsv')
d.to_csv(out, sep='\t', index=False)
P(f'写出 -> {out}')

P('\n=== 校正前后对照（前 6 行） ===')
P(d[['outcome', 'GeneSymbol', 'SNP', 'Z', 'beta_exp', 'OR_perSDU',
     'beta_chip', 'exp_dir_chip']].head(6).to_string(index=False))

P('\n=== CHIP 升高方向上表达被推高的基因占比 ===')
P(d.groupby('outcome').exp_dir_chip.value_counts().unstack(fill_value=0).to_string())
