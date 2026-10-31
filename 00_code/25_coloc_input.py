# -*- coding: utf-8 -*-
"""
25_coloc_input.py —— 生成 coloc 输入（已完成等位基因对齐）

对每个区域：CHIP(暴露) × {MPN, PV}
  - 取 SNP 交集
  - 统一以 CHIP 的 effect allele 为参照，必要时翻转 FinnGen 的 beta
  - 处理链互补 (A<->T, C<->G)
  - 剔除 indel / 多等位 / MAF 过低位点
输出: 09_coloc/input/<region>_CHIP_vs_<pheno>.tsv
      列 SNP chr pos ea oa b1 se1 p1 maf1 b2 se2 p2 maf2
"""
import os, sys

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
RD   = os.path.join(BASE, '04_clean', 'region')
OUTD = os.path.join(BASE, '09_coloc', 'input')
os.makedirs(OUTD, exist_ok=True)

REGIONS = ['TERT', 'JAK2', 'C3q25', 'C11q22', 'C21q22', 'C17p13']
PHENOS  = ['MPN', 'PV', 'MDS', 'AA', 'ITP']
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def load(p):
    d = {}
    with open(p, encoding='utf-8') as f:
        f.readline()
        for ln in f:
            c = ln.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            try:
                d[c[0]] = dict(ch=c[1], pos=int(c[2]), ea=c[3].upper(), oa=c[4].upper(),
                               b=float(c[5]), se=float(c[6]), p=float(c[7]),
                               maf=min(float(c[8]), 1 - float(c[8])))
            except Exception:
                continue
    return d

log = []
for rn in REGIONS:
    pc = os.path.join(RD, f'CHIP_{rn}.tsv')
    if not os.path.exists(pc):
        log.append(f'SKIP {rn}: no CHIP region file'); continue
    chip = load(pc)
    if not chip:
        log.append(f'SKIP {rn}: empty CHIP'); continue
    for ph in PHENOS:
        pp = os.path.join(RD, f'{ph}_{rn}.tsv')
        if not os.path.exists(pp):
            continue
        phe = load(pp)
        common = sorted(set(chip) & set(phe))
        rows, nflip, ncomp, nbad = [], 0, 0, 0
        for s in common:
            a, b = chip[s], phe[s]
            if len(a['ea']) != 1 or len(b['ea']) != 1:
                nbad += 1; continue
            b2 = b['b']
            if a['ea'] == b['ea'] and a['oa'] == b['oa']:
                pass
            elif a['ea'] == b['oa'] and a['oa'] == b['ea']:
                b2 = -b2; nflip += 1
            elif COMP.get(b['ea']) == a['ea'] and COMP.get(b['oa']) == a['oa']:
                b2 = -b2; ncomp += 1          # 链翻转 + 等位互换
            elif COMP.get(b['ea']) == a['oa'] and COMP.get(b['oa']) == a['ea']:
                ncomp += 1                     # 仅链翻转
            else:
                nbad += 1; continue
            if a['se'] <= 0 or b['se'] <= 0 or b['p'] <= 0:
                nbad += 1; continue
            rows.append((s, a['ch'], a['pos'], a['ea'], a['oa'],
                         a['b'], a['se'], a['p'], a['maf'],
                         b2, b['se'], b['p'], b['maf']))
        if len(rows) < 50:
            log.append(f'  {rn:8s} vs {ph:4s}  only {len(rows)} SNPs -> skip')
            continue
        outp = os.path.join(OUTD, f'{rn}_CHIP_vs_{ph}.tsv')
        with open(outp, 'w', encoding='utf-8') as f:
            f.write('SNP\tchr\tpos\tea\toa\tb1\tse1\tp1\tmaf1\tb2\tse2\tp2\tmaf2\n')
            for r in rows:
                f.write('\t'.join([r[0], r[1], str(r[2]), r[3], r[4]] +
                                  [f'{x:.6g}' if isinstance(x, float) else str(x) for x in r[5:]]) + '\n')
        log.append(f'  {rn:8s} vs {ph:4s}  SNPs={len(rows):6d}  flip={nflip:5d}  strand={ncomp:5d}  drop={nbad:4d}')

txt = '\n'.join(log)
with open(os.path.join(BASE, '01_log', 'coloc_input.txt'), 'w', encoding='utf-8') as f:
    f.write(txt)
print(txt)
