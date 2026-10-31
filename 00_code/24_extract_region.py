# -*- coding: utf-8 -*-
"""
24_extract_region.py —— 为 coloc 提取区域级完整 sumstats（不限 p 值）

用法: python 24_extract_region.py <tag>
      tag ∈ {CHIP, MPN, PV, MDS, AA, ITP}

输出: 04_clean/region/<tag>_<region>.tsv
      列: SNP chr pos ea oa beta se pval eaf
"""
import gzip, os, re, sys, time

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
OUTD = os.path.join(BASE, '04_clean', 'region')
os.makedirs(OUTD, exist_ok=True)

# 区域定义 (GRCh38)。±200 kb 围绕共享 lead SNP
REGIONS = {
    'TERT':   ('5',  1_050_000, 1_500_000),   # rs7705526 / rs2853677
    'JAK2':   ('9',  4_700_000, 5_500_000),   # rs1327494(46/1) / rs77375493(V617F)
    'C3q25':  ('3', 160_200_000, 160_600_000),  # rs2305407
    'C11q22': ('11', 108_200_000, 108_600_000),  # rs10890839
    'C21q22': ('21', 34_800_000, 35_200_000),   # rs2834708 / rs56068023 (RUNX1 区)
    'C17p13': ('17', 7_600_000, 8_000_000),     # rs7406485 / rs78378222 (TP53 区)
}

SRC = {
    'CHIP': os.path.join(BASE, '02_raw', 'GCST90165267_CHIP_all.h.tsv.gz'),
    'MPN':  os.path.join(BASE, '02_raw', 'finngen_R12_MYELOPROF_NONCML_MPN.gz'),
    'PV':   os.path.join(BASE, '02_raw', 'finngen_R12_POLYCYTVERA_PV.gz'),
    'MDS':  os.path.join(BASE, '02_raw', 'finngen_R12_C3_MYELODYSP_SYNDR_EXALLC_MDS.gz'),
    'AA':   os.path.join(BASE, '02_raw', 'finngen_R12_D3_OTHERAPLASTICANAEMIA_AA.gz'),
    'ITP':  os.path.join(BASE, '02_raw', 'finngen_R12_D3_ITP_ITP.gz'),
}

tag = sys.argv[1] if len(sys.argv) > 1 else 'CHIP'
src = SRC[tag]
if not os.path.exists(src):
    print('MISSING', src); sys.exit(1)

re_se = re.compile(r'REGENIE_SE=([-\d.eE+]+)')
re_beta = re.compile(r'REGENIE_BETA=([-\d.eE+]+)')

outs = {}
for rn in REGIONS:
    fh = open(os.path.join(OUTD, f'{tag}_{rn}.tsv'), 'w', encoding='utf-8')
    fh.write('SNP\tchr\tpos\tea\toa\tbeta\tse\tpval\teaf\n')
    outs[rn] = fh

counts = {rn: 0 for rn in REGIONS}
n = 0
t0 = time.time()
logf = open(os.path.join(BASE, '01_log', f'region_{tag}.log'), 'w', encoding='utf-8')

def w(s):
    logf.write(s + '\n'); logf.flush(); print(s, flush=True)

with gzip.open(src, 'rt', encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    ix = {c: i for i, c in enumerate(header)}

    if tag == 'CHIP':
        i_ch, i_pos = ix['chromosome'], ix['base_pair_location']
        i_ea, i_oa = ix['effect_allele'], ix['other_allele']
        i_p, i_af, i_ai, i_rs = ix['p_value'], ix['effect_allele_frequency'], ix['additional_info'], ix['rsid']
        i_or = ix['odds_ratio']
        is_chip = True
    else:
        i_ch, i_pos = ix['#chrom'], ix['pos']
        i_ea, i_oa = ix['alt'], ix['ref']
        i_p, i_af, i_rs = ix['pval'], ix['af_alt'], ix['rsids']
        i_beta, i_se = ix['beta'], ix['sebeta']
        is_chip = False

    for line in f:
        n += 1
        if n % 5000000 == 0:
            w(f'  ... {n/1e6:.0f}M  {time.time()-t0:.0f}s  ' +
              ' '.join(f'{k}={v}' for k, v in counts.items()))
        # 快速预筛：靠 tab 位置切分代价高，但需要 ch/pos，先粗切
        parts = line.split('\t', 3)
        if len(parts) < 4:
            continue
        ch = parts[i_ch]
        try:
            pos = int(parts[i_pos])
        except Exception:
            continue
        hit = None
        for rn, (rc, lo, hi) in REGIONS.items():
            if ch == rc and lo <= pos <= hi:
                hit = rn; break
        if hit is None:
            continue
        full = line.rstrip('\n').split('\t')
        if len(full) <= max(i_p, i_af, i_rs):
            continue
        rs = full[i_rs]
        if not rs or rs == '.':
            continue
        rs = rs.split(',')[0]
        p = full[i_p]
        if is_chip:
            ai = full[i_ai]
            ms = re_se.search(ai); mb = re_beta.search(ai)
            if not ms:
                continue
            se = ms.group(1)
            if mb:
                beta = mb.group(1)
            else:
                try:
                    beta = '%.6g' % __import__('math').log(float(full[i_or]))
                except Exception:
                    continue
        else:
            beta, se = full[i_beta], full[i_se]
        outs[hit].write('\t'.join([rs, ch, full[i_pos], full[i_ea], full[i_oa],
                                   beta, se, p, full[i_af]]) + '\n')
        counts[hit] += 1

for fh in outs.values():
    fh.close()

w(f'\n=== {tag} done ===')
w(f'scanned {n} lines in {time.time()-t0:.0f}s')
for k, v in counts.items():
    w(f'  {k:8s} {v}')
logf.close()
