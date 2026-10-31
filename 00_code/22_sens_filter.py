# -*- coding: utf-8 -*-
"""
22_sens_filter.py
生成敏感性分析用的 IV 子集：
  S1  去 JAK2 区 (chr9:4.8-5.4 Mb)      —— 排除 JAK2 46/1 单倍型 / V617F 连锁
  S2  去 TERT 区 (chr5:1.2-1.4 Mb)      —— 排除端粒长度通路
  S3  去 JAK2+TERT 双区                   —— 最严格
  S4  去 MHC 区 (chr6:25-35 Mb)          —— 排除 HLA 多效性
"""
import os, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IV   = os.path.join(BASE, '05_iv')
OUT  = os.path.join(BASE, '05_iv', 'sens')
os.makedirs(OUT, exist_ok=True)

REGIONS = {
    'JAK2': (9, 4_800_000, 5_400_000),
    'TERT': (5, 1_200_000, 1_400_000),
    'MHC':  (6, 25_000_000, 35_000_000),
}

def load(p):
    with open(p, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    return hdr, rows

def write(p, hdr, rows):
    with open(p, 'w', encoding='utf-8') as f:
        f.write('\t'.join(hdr) + '\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')

def drop(hdr, rows, regions):
    ci, pi = hdr.index('chr'), hdr.index('pos')
    keep, dropd = [], []
    for r in rows:
        c, p = r[ci], int(float(r[pi]))
        hit = [k for k, (rc, lo, hi) in regions.items() if c.lstrip('chr') == str(rc) and lo <= p <= hi]
        (dropd if hit else keep).append((r, hit))
    return keep, dropd

TARGETS = {
    'CHIP_germline': 'CHIP_germline_iv.tsv',
    'CHIP_all':      'CHIP_iv_gws.tsv',
    'MPN':           'finngen_R12_MYELOPROF_NONCML_MPN_iv_gws.tsv',
    'PV':            'finngen_R12_POLYCYTVERA_PV_iv_gws.tsv',
    'MDS':           'finngen_R12_C3_MYELODYSP_SYNDR_EXALLC_MDS_iv_gws.tsv',
    'AA':            'finngen_R12_D3_OTHERAPLASTICANAEMIA_AA_iv_gws.tsv',
    'ITP':           'finngen_R12_D3_ITP_ITP_iv_gws.tsv',
}

SCHEMES = {
    'noJAK2':      {'JAK2': REGIONS['JAK2']},
    'noTERT':      {'TERT': REGIONS['TERT']},
    'noJAK2TERT':  {'JAK2': REGIONS['JAK2'], 'TERT': REGIONS['TERT']},
    'noMHC':       {'MHC':  REGIONS['MHC']},
}

log = []
for tag, fn in TARGETS.items():
    p = os.path.join(IV, fn)
    if not os.path.exists(p):
        log.append(f'MISSING {fn}')
        continue
    hdr, rows = load(p)
    log.append(f'\n=== {tag}  (n={len(rows)}) ===')
    for sname, regs in SCHEMES.items():
        keep, dropd = drop(hdr, rows, regs)
        outp = os.path.join(OUT, f'{tag}_{sname}.tsv')
        write(outp, hdr, [r for r, _ in keep])
        removed = ';'.join(f"{r[hdr.index('SNP')]}({'/'.join(h)})" for r, h in dropd)
        log.append(f'  {sname:11s} keep={len(keep):3d}  removed={len(dropd):2d}  {removed}')

txt = '\n'.join(log)
with open(os.path.join(BASE, '01_log', 'sens_filter.txt'), 'w', encoding='utf-8') as f:
    f.write(txt)
print(txt)
