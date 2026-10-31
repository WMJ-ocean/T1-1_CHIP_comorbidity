# -*- coding: utf-8 -*-
"""对 22 个共享位点做效应方向一致性符号检验（二项检验）+ 效应相关性"""
import os, math

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
TD = os.path.join(BASE, '06_targets')
TAGS = ['MPN', 'PV', 'MDS', 'AA', 'ITP']

data = {}
for t in TAGS:
    d = {}
    with open(os.path.join(TD, '%s_targets.tsv' % t), encoding='utf-8') as f:
        f.readline()
        for line in f:
            c = line.rstrip('\n').split('\t')
            try:
                d[c[4]] = (float(c[6]), float(c[7]), float(c[5]))
            except Exception:
                pass
    data[t] = d

snps = [s for s in data['MPN'] if all(s in data[t] for t in TAGS)]

def binom_two_sided(k, n, p=0.5):
    # 精确二项检验（双边），用 log 组合数避免溢出
    def pmf(i):
        return math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                        + i * math.log(p) + (n - i) * math.log(1 - p))
    obs = pmf(k)
    tot = sum(pmf(i) for i in range(n + 1) if pmf(i) <= obs + 1e-15)
    return min(1.0, tot)

def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n; my = sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return float('nan'), float('nan')
    r = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)
    # t 检验
    if abs(r) >= 1:
        return r, 0.0
    t = r * math.sqrt((n - 2) / (1 - r * r))
    # 近似 p（t 分布，df=n-2），用正态近似即可
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))
    return r, p

print('=== 与 MPN 的效应方向一致性（n=%d 位点）===\n' % len(snps))
print('%-6s %-10s %-10s %-12s %-12s %s' % ('病种', '同向', '反向', '同向比例', '符号检验P', 'Pearson r (P)'))
print('-' * 82)
for t in TAGS:
    same = opp = 0
    xs, ys = [], []
    for s in snps:
        b1 = data['MPN'][s][0]; b2 = data[t][s][0]
        xs.append(b1); ys.append(b2)
        if b1 * b2 > 0:
            same += 1
        else:
            opp += 1
    p_sign = binom_two_sided(min(same, opp), len(snps))
    r, pr = pearson(xs, ys)
    print('%-6s %-10d %-10d %-12.1f%% %-12.4f %.3f (P=%.3g)' % (
        t, same, opp, 100.0 * same / len(snps), p_sign, r, pr))

print()
print('注：符号检验为零假设「方向随机(50/50)」的双边精确二项检验；')
print('    Pearson r 基于各病 beta 值（原始尺度，未加权）。')
print()
print('=== 逐位点：MPN 显著位点在 AA 中的表现 ===')
print('%-14s %-24s %-24s %s' % ('SNP', 'MPN beta (p)', 'AA beta (p)', '方向'))
for s in sorted(snps, key=lambda x: data['MPN'][x][2]):
    b1, se1, p1 = data['MPN'][s]
    b2, se2, p2 = data['AA'][s]
    if p1 < 5e-8:
        dirn = '相反' if b1 * b2 < 0 else '相同'
        print('%-14s %+8.3f (%.1e)        %+8.3f (%.2f)        %s' % (s, b1, p1, b2, p2, dirn))
