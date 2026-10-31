# -*- coding: utf-8 -*-
"""
从 04_clean/*_p5e6.tsv 提取独立 lead SNP（1Mb 窗口、按 p 最小保留）
输出：05_iv/<tag>_iv.tsv  +  汇总表
"""
import os, csv, sys
from collections import defaultdict

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
CLEAN = os.path.join(BASE, '04_clean')
OUTD = os.path.join(BASE, '05_iv')
os.makedirs(OUTD, exist_ok=True)
WINDOW = 1_000_000

summary = []
for fn in sorted(os.listdir(CLEAN)):
    if not fn.endswith('_p5e6.tsv'):
        continue
    tag = fn.replace('_p5e6.tsv', '')
    rows = []
    with open(os.path.join(CLEAN, fn), encoding='utf-8') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for r in rd:
            try:
                r['pval'] = float(r['pval'])
                r['beta'] = float(r['beta'])
                r['sebeta'] = float(r['sebeta'])
                r['af'] = float(r['af']) if r['af'] not in ('', 'NA', 'nan') else float('nan')
                r['pos'] = int(r['pos'])
            except Exception:
                continue
            rows.append(r)
    # 按染色体分组，窗口内取最小 p
    byc = defaultdict(list)
    for r in rows:
        byc[r['chrom']].append(r)
    iv = []
    for c, rs in byc.items():
        rs.sort(key=lambda x: x['pos'])
        cur = []
        for r in rs:
            if cur and r['pos'] - cur[-1]['pos'] > WINDOW:
                iv.append(min(cur, key=lambda x: x['pval']))
                cur = [r]
            else:
                cur.append(r)
        if cur:
            iv.append(min(cur, key=lambda x: x['pval']))
    iv.sort(key=lambda x: x['pval'])

    dst = os.path.join(OUTD, '%s_iv.tsv' % tag)
    with open(dst, 'w', encoding='utf-8') as g:
        g.write('SNP\tchr\tpos\teffect_allele\tother_allele\tbeta\tse\tpval\teaf\n')
        for r in iv:
            g.write('\t'.join([r['rsid'], r['chrom'], str(r['pos']), r['alt'], r['ref'],
                               '%.6g' % r['beta'], '%.6g' % r['sebeta'],
                               '%.6g' % r['pval'], '%.4g' % r['af']]) + '\n')
    ngws = sum(1 for r in iv if r['pval'] < 5e-8)
    summary.append((tag, len(rows), len(iv), ngws, min(r['pval'] for r in iv) if iv else 1))
    print('%s: p5e6=%d  IV=%d  GWS_IV=%d  minP=%.3g' % (tag, len(rows), len(iv), ngws,
                                                        min(r['pval'] for r in iv) if iv else 1))

print()
print('%-55s %-10s %-8s %-8s %s' % ('dataset', 'p5e6', 'IV', 'GWS_IV', 'minP'))
for s in summary:
    print('%-55s %-10d %-8d %-8d %.3g' % s)

with open(os.path.join(OUTD, 'IV_summary.tsv'), 'w', encoding='utf-8') as g:
    g.write('dataset\tn_p5e6\tn_IV\tn_GWS_IV\tminP\n')
    for s in summary:
        g.write('%s\t%d\t%d\t%d\t%.6g\n' % s)
print()
print('written ->', OUTD)
