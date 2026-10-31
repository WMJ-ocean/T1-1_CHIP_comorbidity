# -*- coding: utf-8 -*-
"""比较 22 个 GWS 位点在五个病中的效应方向与显著性"""
import os, math

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
TD = os.path.join(BASE, '06_targets')
TAGS = ['MPN', 'PV', 'MDS', 'AA', 'ITP']

data = {}   # tag -> {rsid: (beta, se, p)}
for t in TAGS:
    p = os.path.join(TD, '%s_targets.tsv' % t)
    d = {}
    with open(p, encoding='utf-8') as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            c = line.rstrip('\n').split('\t')
            try:
                d[c[4]] = (float(c[6]), float(c[7]), float(c[5]))
            except Exception:
                pass
    data[t] = d

snps = sorted(set().union(*[set(d.keys()) for d in data.values()]))
print('n SNPs =', len(snps))
print()
hdr = '%-14s' % 'SNP'
for t in TAGS:
    hdr += ' %-22s' % t
print(hdr)
print('-' * 125)

direction = {t: [] for t in TAGS}
for s in snps:
    line = '%-14s' % s
    for t in TAGS:
        if s in data[t]:
            b, se, p = data[t][s]
            sig = '**' if p < 5e-8 else ('*' if p < 1e-5 else '  ')
            line += ' %+8.3f(p%.1e)%s' % (b, p, sig)
            direction[t].append((s, 1 if b > 0 else -1))
        else:
            line += ' %-22s' % 'NA'
    print(line)

print()
print('=== 各病在这些位点上的方向分布 ===')
for t in TAGS:
    pos = sum(1 for _, d in direction[t] if d > 0)
    neg = len(direction[t]) - pos
    print('%-6s n=%2d  beta>0: %2d  beta<0: %2d' % (t, len(direction[t]), pos, neg))

# 与 MPN 的方向一致性
print()
print('=== 与 MPN 效应方向的一致性（仅同向/反向计数）===')
m = dict(direction['MPN'])
for t in TAGS[1:]:
    same = opp = 0
    for s, d in direction[t]:
        if s in m:
            if m[s] == d:
                same += 1
            else:
                opp += 1
    print('%-6s 同向 %2d  反向 %2d' % (t, same, opp))
