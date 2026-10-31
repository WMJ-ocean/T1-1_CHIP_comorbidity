# -*- coding: utf-8 -*-
"""
P1-2: 为 MDS / ITP（无全基因组显著位点）构造 p<5e-6 的反向 MR 输入
  1) 合并所有需要反向提取的 rsID -> 06_targets/gws_snps.txt
  2) 从 CHIP_p5e6.tsv 构造放宽阈值(p<5e-6)的胚系 IV，用于正向功效检验
输出:
  06_targets/gws_snps.txt                 (供 16_parse_chip.py 用)
  05_iv/sens/CHIP_germline_p5e6.tsv       放宽阈值胚系 IV
  06_targets/CHIPg_p5e6_snps.txt          其 rsID 列表
"""
import csv, os

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
IVD  = os.path.join(BASE, '05_iv')
TD   = os.path.join(BASE, '06_targets')
SD   = os.path.join(BASE, '05_iv', 'sens')
os.makedirs(SD, exist_ok=True)

def rd(p):
    return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))

# ---------- 1) 合并反向 MR 需要的 rsID ----------
allrs = []
seen  = set()
srcs  = [
    ('germline', os.path.join(IVD, 'CHIP_germline_iv.tsv')),
    ('somatic',  os.path.join(IVD, 'CHIP_somatic_iv.tsv')),
    ('MPN',      os.path.join(IVD, 'finngen_R12_MYELOPROF_NONCML_MPN_iv_gws.tsv')),
    ('PV',       os.path.join(IVD, 'finngen_R12_POLYCYTVERA_PV_iv_gws.tsv')),
    ('MDS',      os.path.join(IVD, 'finngen_R12_C3_MYELODYSP_SYNDR_EXALLC_MDS_iv.tsv')),
    ('ITP',      os.path.join(IVD, 'finngen_R12_D3_ITP_ITP_iv.tsv')),
    ('AA',       os.path.join(IVD, 'finngen_R12_D3_OTHERAPLASTICANAEMIA_AA_iv.tsv')),
]
counts = {}
for name, p in srcs:
    if not os.path.exists(p):
        continue
    rs = [r['SNP'] for r in rd(p) if r['SNP'].startswith('rs')]
    counts[name] = len(rs)
    for s in rs:
        if s not in seen:
            seen.add(s); allrs.append(s)

with open(os.path.join(TD, 'gws_snps.txt'), 'w', encoding='utf-8') as g:
    for s in allrs:
        g.write(s + '\n')
print('gws_snps.txt:', {k: v for k, v in counts.items()}, '-> union', len(allrs))

# ---------- 2) 放宽阈值的 CHIP 胚系 IV ----------
rows = rd(os.path.join(BASE, '04_clean', 'CHIP_p5e6.tsv'))
rows = [r for r in rows if abs(float(r['beta'])) < 1.0]          # 去掉体细胞突变读取
# CHIP_p5e6.tsv 混有"反向 MR 目标位点"(无条件保留、无 p 阈值)，正向 IV 必须过滤
rows = [r for r in rows if float(r['pval']) < 5e-6]
rows.sort(key=lambda r: float(r['pval']))

kept, bychr = [], {}
WIN = 1_000_000
for r in rows:
    c = r['chrom']; pos = int(r['pos'])
    ok = True
    for kp in bychr.get(c, []):
        if abs(pos - kp) < WIN:
            ok = False; break
    if ok:
        kept.append(r)
        bychr.setdefault(c, []).append(pos)

hdr = ['SNP', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'eaf']
dst = os.path.join(SD, 'CHIP_germline_p5e6.tsv')
with open(dst, 'w', encoding='utf-8') as g:
    g.write('\t'.join(hdr) + '\n')
    for r in kept:
        # CHIP_p5e6.tsv 的列是 ref/alt：alt 才是效应等位（与 FinnGen targets 一致）
        # 写错会把 beta 符号整体翻转 -> IVW 的 OR 取倒数（实测 2.93 -> 0.32）
        g.write('\t'.join([r.get('SNP', r.get('rsid')), r['chrom'], r['pos'],
                           r['alt'], r['ref'], r['beta'], r['sebeta'], r['pval'], r['af']]) + '\n')
with open(os.path.join(TD, 'CHIPg_p5e6_snps.txt'), 'w', encoding='utf-8') as g:
    for r in kept:
        g.write(r.get('SNP', r.get('rsid')) + '\n')
print('relaxed germline IV (p<5e-6):', len(kept), '->', dst)
print('  minP=%.2e  maxP=%.2e' % (min(float(r['pval']) for r in kept),
                                  max(float(r['pval']) for r in kept)))
