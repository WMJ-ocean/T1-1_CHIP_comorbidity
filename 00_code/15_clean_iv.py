# -*- coding: utf-8 -*-
"""清理 IV：只保留 p<5e-8 且 rsID 唯一的位点"""
import os, csv, re

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
IVD = os.path.join(BASE, '05_iv')
OUTD = os.path.join(BASE, '05_iv')

for fn in sorted(os.listdir(IVD)):
    if not fn.endswith('_iv.tsv'):
        continue
    tag = fn.replace('_iv.tsv', '')
    out = []
    seen = set()
    for r in csv.DictReader(open(os.path.join(IVD, fn), encoding='utf-8'), delimiter='\t'):
        if float(r['pval']) >= 5e-8:
            continue
        snp = r['SNP'].split(',')[0].strip()
        if not re.match(r'^rs\d+$', snp) or snp in seen:
            continue
        seen.add(snp)
        out.append([snp, r['chr'], r['pos'], r['effect_allele'], r['other_allele'],
                    r['beta'], r['se'], r['pval'], r['eaf']])
    dst = os.path.join(OUTD, '%s_iv_gws.tsv' % tag)
    with open(dst, 'w', encoding='utf-8') as g:
        g.write('SNP\tchr\tpos\teffect_allele\tother_allele\tbeta\tse\tpval\teaf\n')
        for row in out:
            g.write('\t'.join(row) + '\n')
    print('%-52s GWS IV = %d' % (tag[:52], len(out)))
