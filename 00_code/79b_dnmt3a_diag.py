# -*- coding: utf-8 -*-
"""79b_dnmt3a_diag.py —— DNMT3A 亚型 MR 的稳健性诊断
逐位点 Wald 比、留一法、加权中位数、I^2，判断异质性来源。
"""
import io, os, math, gzip, sys
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, '02_raw')
IVD = os.path.join(BASE, '05_iv')
FINN = {'MPN': 'finngen_R12_MYELOPROF_NONCML_MPN.gz',
        'PV': 'finngen_R12_POLYCYTVERA_PV.gz'}


def read_tsv(p):
    with io.open(p, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for ln in f:
            if ln.strip():
                yield dict(zip(hdr, ln.rstrip('\n').split('\t')))


def prune(rows, win):
    rows = sorted(rows, key=lambda r: float(r['pval']))
    keep = []
    for r in rows:
        pos = int(r['pos'])
        if all(not (k['chr'] == r['chr'] and abs(int(k['pos']) - pos) < win) for k in keep):
            keep.append(r)
    return keep


def load_finn(fn, keys):
    d = {}
    path = os.path.join(RAW, fn)
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        f.readline()
        for ln in f:
            c = ln.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            k = (c[0], c[1])
            if k in keys:
                d[k] = dict(ref=c[2], alt=c[3], beta=float(c[8]), se=float(c[9]), p=float(c[6]),
                            rsid=c[4])
    return d


def _norm(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def pz(z):
    return 2 * (1 - _norm(abs(z)))


def ivw(bx, by, sy):
    bx, by, sy = map(np.asarray, (bx, by, sy))
    w = 1.0 / sy**2
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    se = math.sqrt(1.0 / np.sum(w * bx**2))
    k = len(bx)
    q = float(np.sum(w * (by - b * bx)**2))
    return b, se, q, k


def wmedian(bx, by, sy):
    r = np.asarray(by) / np.asarray(bx)
    w = (np.asarray(bx)**2) / (np.asarray(sy)**2)
    o = np.argsort(r); r = r[o]; w = w[o]
    c = np.cumsum(w) - 0.5 * w
    c /= np.sum(w)
    j = np.searchsorted(c, 0.5)
    j = min(max(j, 1), len(r) - 1)
    return float(r[j - 1] + (r[j] - r[j - 1]) * (0.5 - c[j - 1]) / (c[j] - c[j - 1]))


rows = list(read_tsv(os.path.join(IVD, 'dnmt3a_chip_gws.tsv')))
cand = [r for r in rows
        if len(r['effect_allele']) == 1 and len(r['other_allele']) == 1
        and abs(float(r['beta'])) <= 1.3 and float(r['se']) > 0]
keep = prune(cand, 10_000_000)
keys = {(r['chr'], r['pos']) for r in keep}

for phe, fn in FINN.items():
    d = load_finn(fn, keys)
    recs = []
    for r in keep:
        o = d.get((r['chr'], r['pos']))
        if not o:
            continue
        ea, oa = r['effect_allele'], r['other_allele']
        s = 1 if (o['alt'] == ea and o['ref'] == oa) else (-1 if (o['alt'] == oa and o['ref'] == ea) else 0)
        if not s:
            continue
        bx = float(r['beta']); by = o['beta'] * s; sy = o['se']
        recs.append(dict(snp=r['snp'], chr=r['chr'], pos=r['pos'], bx=bx, by=by, sy=sy,
                         ratio=by / bx, rsid=o['rsid'], p_out=o['p']))
    bx = [x['bx'] for x in recs]; by = [x['by'] for x in recs]; sy = [x['sy'] for x in recs]
    b, se, q, k = ivw(bx, by, sy)
    bfe, sefe, _, _ = ivw(bx, by, sy)
    se_fe = math.sqrt(1.0 / np.sum(1.0 / np.asarray(sy)**2))
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    print(f'\n===== DNMT3A -> {phe} =====')
    print(f'  IVW(随机)  logOR={b:+.4f} se={se:.4f} p={pz(b/se):.3g}  OR={math.exp(b):.2f}')
    print(f'  IVW(固定)  logOR={bfe:+.4f} se={se_fe:.4f} p={pz(bfe/se_fe):.3g}')
    print(f'  Q={q:.1f} df={k-1} I2={i2:.1f}%   k={k}')
    print(f'  加权中位数 logOR={wmedian(bx, by, sy):+.4f}')
    recs.sort(key=lambda x: x['p_out'])
    print(f'  {"snp":22s} {"chr:pos":18s} {"bx":>8s} {"by":>8s} {"ratio":>8s}')
    for x in recs:
        print(f'  {x["snp"]:22s} {x["chr"]+":"+x["pos"]:18s} {x["bx"]:+8.3f} {x["by"]:+8.3f} {x["ratio"]:+8.3f}')
    # 留一法
    loo = []
    for i in range(len(recs)):
        sub = [x for j, x in enumerate(recs) if j != i]
        bb, ss, qq, kk = ivw([x['bx'] for x in sub], [x['by'] for x in sub], [x['sy'] for x in sub])
        loo.append((pz(bb / ss), bb, recs[i]['snp']))
    loo.sort()
    print('  留一法最极端 5 个（p 最小 -> 最大 的端点）:')
    for p_, bb, snp in loo[:2] + loo[-2:]:
        print(f'    去掉 {snp:22s} -> logOR={bb:+.4f} p={p_:.3g}')
    # 去掉 TERT 后的估计
    sub = [x for x in recs if not x['snp'].startswith('5:12')]
    if len(sub) >= 2:
        bb, ss, qq, kk = ivw([x['bx'] for x in sub], [x['by'] for x in sub], [x['sy'] for x in sub])
        print(f'  剔除 TERT 位点后: k={kk} logOR={bb:+.4f} se={ss:.4f} p={pz(bb/ss):.3g} Q={qq:.1f}')
