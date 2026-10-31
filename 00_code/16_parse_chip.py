# -*- coding: utf-8 -*-
"""
解析 CHIP (GCST90165267) harmonised sumstats
注意: standard_error 列为 NA; 真实 SE 在 additional_info 的 REGENIE_SE=;
      odds_ratio = exp(REGENIE_BETA)
列: chromosome base_pair_location effect_allele other_allele odds_ratio standard_error
    effect_allele_frequency p_value variant_id ... additional_info ... rsid
输出: 04_clean/CHIP_p5e6.tsv  (chrom pos ref alt rsid pval beta sebeta af)
"""
import gzip, os, re, sys, time

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
SRC = os.path.join(BASE, '02_raw', 'GCST90165267_CHIP_all.h.tsv.gz')
DST = os.path.join(BASE, '04_clean', 'CHIP_p5e6.tsv')
DST2 = os.path.join(BASE, '06_targets', 'CHIP_targets.tsv')
RSFILE = os.path.join(BASE, '06_targets', 'gws_snps.txt')   # 反向 MR 需要这些位点全部拿到
LOGP = os.path.join(BASE, '01_log', 'parse_chip.txt')
THRESH = 5e-6

# 反向 MR 的目标位点（不设 p 阈值，必须全部提取）
want = set()
if os.path.exists(RSFILE):
    for ln in open(RSFILE, encoding='utf-8'):
        ln = ln.strip()
        if ln and ln.startswith('rs'):
            want.add(ln)

log = open(LOGP, 'w', encoding='utf-8')
def w(s):
    log.write(s + '\n'); log.flush(); print(s, flush=True)

if not os.path.exists(SRC):
    w('SRC missing: %s' % SRC); sys.exit(1)
w('src size = %.0f MB' % (os.path.getsize(SRC) / 1e6))

n = 0; kept = 0; nbad = 0; minp = 1.0; c5e8 = 0; c1e6 = 0; nhit = 0
t0 = time.time()
out = open(DST, 'w', encoding='utf-8')
out.write('chrom\tpos\tref\talt\trsid\tpval\tbeta\tsebeta\taf\n')
out2 = open(DST2, 'w', encoding='utf-8')
out2.write('chrom\tpos\tref\talt\trsid\tpval\tbeta\tsebeta\taf\n')
found = set()

with gzip.open(SRC, 'rt', encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    ix = {c: i for i, c in enumerate(header)}
    need = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
            'odds_ratio', 'p_value', 'effect_allele_frequency', 'additional_info', 'rsid']
    miss = [c for c in need if c not in ix]
    if miss:
        w('MISSING COLUMNS: %s' % miss); sys.exit(1)
    i_ch, i_pos, i_ea, i_oa = ix['chromosome'], ix['base_pair_location'], ix['effect_allele'], ix['other_allele']
    i_or, i_p, i_af, i_ai, i_rs = ix['odds_ratio'], ix['p_value'], ix['effect_allele_frequency'], ix['additional_info'], ix['rsid']
    re_se = re.compile(r'REGENIE_SE=([-\d.eE+]+)')
    re_beta = re.compile(r'REGENIE_BETA=([-\d.eE+]+)')

    for line in f:
        n += 1
        parts = line.rstrip('\n').split('\t')
        if len(parts) <= i_ai:
            nbad += 1; continue
        try:
            p = float(parts[i_p])
        except Exception:
            nbad += 1; continue
        if p < minp:
            minp = p
        if p < 5e-8:
            c5e8 += 1
        if p < 1e-6:
            c1e6 += 1
        rs = parts[i_rs]
        in_want = rs in want
        if p >= THRESH and not in_want:
            continue
        ai = parts[i_ai]
        ms = re_se.search(ai)
        mb = re_beta.search(ai)
        if ms:
            se = float(ms.group(1))
            beta = float(mb.group(1)) if mb else None
        else:
            try:
                se = float(parts[ix['standard_error']])
                beta = None
            except Exception:
                nbad += 1; continue
        if beta is None:
            try:
                beta = float(__import__('math').log(float(parts[i_or])))
            except Exception:
                nbad += 1; continue
        if not (se > 0):
            nbad += 1; continue
        row = '\t'.join([parts[i_ch], parts[i_pos], parts[i_oa], parts[i_ea], rs,
                         parts[i_p], '%.6g' % beta, '%.6g' % se, parts[i_af]]) + '\n'
        out.write(row)
        kept += 1
        if in_want:
            out2.write(row)
            found.add(rs)
            nhit += 1
        if n % 5000000 == 0:
            w('  ... %d M lines, kept %d (%.0fs)' % (n / 1e6, kept, time.time() - t0))

out.close()
out2.close()
w('')
w('total lines   = %d' % n)
w('p<5e-8       = %d' % c5e8)
w('p<1e-6       = %d' % c1e6)
w('kept (p<5e-6)= %d' % kept)
w('target hits  = %d / %d  (missing: %s)' % (len(found), len(want), ','.join(sorted(want - found)) or 'none'))
w('min p        = %.3g' % minp)
w('bad/skipped  = %d' % nbad)
w('elapsed      = %.0fs' % (time.time() - t0))
w('written -> %s' % DST)
w('written -> %s' % DST2)
