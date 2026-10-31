# -*- coding: utf-8 -*-
"""合并 MPN/PV/AA 的 GWS IV 成 SNP 列表，供跨病提取"""
import os, csv

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
IVD = os.path.join(BASE, '05_iv')
OUTD = os.path.join(BASE, '06_targets')
os.makedirs(OUTD, exist_ok=True)

want = {}
for fn, tag in [('finngen_R12_MYELOPROF_NONCML_MPN_iv.tsv', 'MPN'),
                ('finngen_R12_POLYCYTVERA_PV_iv.tsv', 'PV'),
                ('finngen_R12_D3_OTHERAPLASTICANAEMIA_AA_iv.tsv', 'AA')]:
    p = os.path.join(IVD, fn)
    for r in csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'):
        if float(r['pval']) < 5e-8 and r['SNP'] and not r['SNP'].startswith('rs,'):
            for s in r['SNP'].split(','):
                want.setdefault(s.strip(), set()).add(tag)

out = os.path.join(OUTD, 'gws_snps.txt')
with open(out, 'w', encoding='utf-8') as g:
    for s in sorted(want):
        g.write(s + '\n')
print('n unique GWS SNPs =', len(want))
from collections import Counter
c = Counter()
for s, tags in want.items():
    c['+'.join(sorted(tags))] += 1
for k, v in c.most_common():
    print('  %-20s %d' % (k, v))
print('written ->', out)
