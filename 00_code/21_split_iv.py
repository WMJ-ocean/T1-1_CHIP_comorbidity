# -*- coding: utf-8 -*-
"""把 CHIP IV 分为胚系易感组 / 体细胞突变读取组"""
import csv, os

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
IV = os.path.join(BASE, '05_iv', 'CHIP_iv_gws.tsv')
TD = os.path.join(BASE, '06_targets')

rows = list(csv.DictReader(open(IV, encoding='utf-8'), delimiter='\t'))
hdr = ['SNP', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'eaf']

germ = [r for r in rows if abs(float(r['beta'])) < 1.0]
som = [r for r in rows if abs(float(r['beta'])) >= 1.0]

for name, sub in [('CHIP_germline', germ), ('CHIP_somatic', som)]:
    dst = os.path.join(BASE, '05_iv', name + '_iv.tsv')
    with open(dst, 'w', encoding='utf-8') as g:
        g.write('\t'.join(hdr) + '\n')
        for r in sub:
            g.write('\t'.join([r['SNP'], r['chr'], r['pos'], r['effect_allele'],
                               r['other_allele'], r['beta'], r['se'], r['pval'], r['eaf']]) + '\n')
    snpl = os.path.join(TD, name + '_snps.txt')
    with open(snpl, 'w', encoding='utf-8') as g:
        for r in sub:
            g.write(r['SNP'] + '\n')
    print('%-16s n=%2d -> %s' % (name, len(sub), dst))
print()
print('germline SNPs:', ' '.join(r['SNP'] for r in germ))
print()
print('somatic  SNPs:', ' '.join(r['SNP'] for r in som))
