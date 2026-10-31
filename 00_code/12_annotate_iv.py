# -*- coding: utf-8 -*-
"""给 IV 加基因注释（GRCh38），标出 CHIP/MPN 关键基因区域"""
import os, csv

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
IVD = os.path.join(BASE, '05_iv')

# GRCh38 关键基因窗口 (chrom, start, end, gene)
GENES = [
    ('9', 4_980_000, 5_220_000, 'JAK2'),
    ('4', 105_140_000, 106_300_000, 'TET2'),
    ('2', 25_200_000, 25_340_000, 'DNMT3A'),
    ('20', 32_350_000, 32_425_000, 'ASXL1'),
    ('1', 43_300_000, 43_420_000, 'MPL'),
    ('19', 12_930_000, 13_060_000, 'CALR'),
    ('2', 197_370_000, 197_420_000, 'SF3B1'),
    ('17', 76_725_000, 76_745_000, 'SRSF2'),
    ('17', 7_650_000, 7_700_000, 'TP53'),
    ('5', 1_250_000, 1_320_000, 'TERT'),
    ('3', 169_000_000, 169_550_000, 'MECOM'),
    ('6', 135_000_000, 135_450_000, 'HBS1L-MYB'),
    ('1', 155_000_000, 155_300_000, 'GBA/GBAP1'),
    ('3', 46_800_000, 47_000_000, 'PPM1D'),
    ('12', 56_000_000, 56_400_000, 'STAT2/IL23A'),
    ('14', 105_400_000, 105_600_000, 'TMEM121/CDC42BPB'),
    ('7', 92_000_000, 92_400_000, 'CDK6'),
    ('10', 6_000_000, 6_500_000, 'PFKFB3/IL2RA'),
]

def annot(c, p):
    hits = [g[3] for g in GENES if g[0] == c and g[1] <= p <= g[2]]
    return '/'.join(hits) if hits else ''

for fn in sorted(os.listdir(IVD)):
    if not fn.endswith('_iv.tsv'):
        continue
    tag = fn.replace('_iv.tsv', '')
    rows = list(csv.DictReader(open(os.path.join(IVD, fn), encoding='utf-8'), delimiter='\t'))
    print('=' * 105)
    print('%s   n=%d' % (tag, len(rows)))
    print('%-14s %-4s %-13s %-4s %-4s %-11s %-11s %-6s %s' % (
        'SNP', 'chr', 'pos', 'EA', 'OA', 'beta', 'se', 'p', 'gene'))
    nkey = 0
    for r in rows:
        gene = annot(r['chr'], int(r['pos']))
        if gene:
            nkey += 1
        p = float(r['pval'])
        star = '*' if p < 5e-8 else ' '
        print('%-14s %-4s %-13s %-4s %-4s %-11s %-11s %-9.2e%s %s' % (
            r['SNP'], r['chr'], r['pos'], r['effect_allele'], r['other_allele'],
            r['beta'][:9], r['se'][:9], p, star, gene))
    print('  -> 落在关键基因区的 IV: %d / %d' % (nkey, len(rows)))
