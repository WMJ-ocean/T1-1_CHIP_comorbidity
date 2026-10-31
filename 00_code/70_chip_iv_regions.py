# -*- coding: utf-8 -*-
"""70_chip_iv_regions.py —— 为 eQTL coloc 抽取 CHIP GWAS 在 21 个胚系 IV ±1Mb 窗口内的全量 sumstats

为什么要新抽：原有 04_clean/region/CHIP_*.tsv 只覆盖 6 个经典 CHIP 位点，
而 eQTL coloc 需要逐个 IV 位点做（eQTL 窗口已按同样 ±1Mb 抽好）。

CHIP GWAS (GCST90165267, Kessler 2022) 的 standard_error 列全为 NA，
真实 SE 在 additional_info 的 REGENIE_SE=... 里，beta 同理。
输出: 04_clean/region/CHIPIV_<rs>.tsv  (SNP chr pos ea oa beta se pval eaf)
"""
import gzip, io, math, os, re, time

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, '02_raw', 'GCST90165267_CHIP_all.h.tsv.gz')
IVF = os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv')
OUTD = os.path.join(BASE, '04_clean', 'region')
os.makedirs(OUTD, exist_ok=True)
WIN = 1_000_000

ivs = []
for ln in io.open(IVF, encoding='utf-8').readlines()[1:]:
    c = ln.rstrip('\n').split('\t')
    if len(c) >= 3:
        ivs.append((c[0], c[1], int(c[2])))
print(f'[setup] IV {len(ivs)} 个', flush=True)

wins = []
for rs, ch, pos in ivs:
    wins.append((rs, ch, pos - WIN, pos + WIN))
by_chr = {}
for w in wins:
    by_chr.setdefault(w[1], []).append(w)

re_se = re.compile(r'REGENIE_SE=([-\d.eE+]+)')
re_beta = re.compile(r'REGENIE_BETA=([-\d.eE+]+)')

outs = {}
for rs, ch, lo, hi in wins:
    fh = io.open(os.path.join(OUTD, f'CHIPIV_{rs}.tsv'), 'w', encoding='utf-8')
    fh.write('SNP\tchr\tpos\tea\toa\tbeta\tse\tpval\teaf\n')
    outs[rs] = fh
counts = {rs: 0 for rs, _, _, _ in wins}

n = 0
t0 = time.time()
with gzip.open(SRC, 'rt', encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    ix = {c: i for i, c in enumerate(header)}
    i_ch, i_pos = ix['chromosome'], ix['base_pair_location']
    i_ea, i_oa = ix['effect_allele'], ix['other_allele']
    i_p, i_af = ix['p_value'], ix['effect_allele_frequency']
    i_ai, i_rs, i_or = ix['additional_info'], ix['rsid'], ix['odds_ratio']
    need = max(i_p, i_af, i_ai, i_rs, i_or)
    for line in f:
        n += 1
        if n % 5_000_000 == 0:
            print(f'  ...{n/1e6:.0f}M  {time.time()-t0:.0f}s  '
                  f'命中合计 {sum(counts.values()):,}', flush=True)
        parts = line.split('\t', 4)
        if len(parts) < 5:
            continue
        ch = parts[i_ch]
        ws = by_chr.get(ch)
        if not ws:
            continue
        try:
            pos = int(parts[i_pos])
        except ValueError:
            continue
        hit = None
        for rs, _, lo, hi in ws:
            if lo <= pos <= hi:
                hit = rs; break
        if hit is None:
            continue
        full = line.rstrip('\n').split('\t')
        if len(full) <= need:
            continue
        s = full[i_rs]
        if not s or s == '.':
            continue
        s = s.split(',')[0]
        ai = full[i_ai]
        ms = re_se.search(ai)
        if not ms:
            continue
        mb = re_beta.search(ai)
        if mb:
            beta = mb.group(1)
        else:
            try:
                beta = '%.6g' % math.log(float(full[i_or]))
            except Exception:
                continue
        outs[hit].write('\t'.join([s, ch, full[i_pos], full[i_ea], full[i_oa],
                                   beta, ms.group(1), full[i_p], full[i_af]]) + '\n')
        counts[hit] += 1

for fh in outs.values():
    fh.close()
print(f'\n=== 完成 === 扫描 {n:,} 行, 用时 {time.time()-t0:.0f}s', flush=True)
for rs, c in sorted(counts.items(), key=lambda x: -x[1]):
    print(f'  {rs:<12} {c:>7,}', flush=True)
