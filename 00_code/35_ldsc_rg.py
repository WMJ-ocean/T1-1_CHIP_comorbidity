# -*- coding: utf-8 -*-
"""
35_ldsc_rg.py —— 自包含的交叉性状 LD score 回归（仅依赖 numpy）

为什么自己写：官方 ldsc.py 依赖 Python2 时代的 pandas API，
在 numpy 2.x / pandas 3.x 上无法直接运行；这里用纯 numpy 实现
Bulik-Sullivan et al. (2015) 的 h2 / rg 估计，并用合成数据自检。

模型
    单性状 : E[z_j^2]  = 1 + (N h^2 / M) * l_j        -> slope = N h^2 / M
    交叉   : E[z1 z2]  = rho * sqrt(N1 N2) / M * l_j   -> slope12 = rho sqrt(N1 N2)/M
    rg = rho / sqrt(h1^2 h2^2) = slope12 / sqrt(slope1 * slope2)

用法
    python 35_ldsc_rg.py --ld-dir <eur_w_ld_chr_dir> --ld-prefix eur_w_ld_chr \
        --sumstats CHIP=.../CHIP.sumstats.tsv MPN=.../MPN.sumstats.tsv ... \
        --out ../09_ldsc/rg_matrix.csv
    python 35_ldsc_rg.py --selftest        # 用已知真值的合成数据验证估计器
"""
import argparse, gzip, os, sys
import math
import numpy as np

# numpy >= 2.0 移除了 np.math.* 别名，正态尾概率改用标准库 math.erf
def norm_p2(z):
    """双尾正态 P 值"""
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))


# ----------------------------------------------------------------- IO
def read_ldscores(ld_dir, prefix, chroms=range(1, 23)):
    """读 <prefix><chr>.l2.ldscore(.gz) -> dict SNP -> L2"""
    ld = {}
    for c in chroms:
        for ext in ['.l2.ldscore.gz', '.l2.ldscore']:
            p = os.path.join(ld_dir, f'{prefix}{c}{ext}')
            if os.path.exists(p):
                break
        else:
            continue
        op = gzip.open if p.endswith('.gz') else open
        with op(p, 'rt') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            i_snp = hdr.index('SNP'); i_l2 = hdr.index('L2')
            for ln in f:
                parts = ln.rstrip('\n').split('\t')
                ld[parts[i_snp]] = float(parts[i_l2])
    return ld


def read_sumstats(path):
    """TSV: SNP  A1  A2  Z  N   -> dict SNP -> (A1, A2, Z, N)"""
    d = {}
    with open(path, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        ix = {k: hdr.index(k) for k in ['SNP', 'A1', 'A2', 'Z', 'N']}
        for ln in f:
            p = ln.rstrip('\n').split('\t')
            if len(p) < len(hdr):
                continue
            try:
                z = float(p[ix['Z']])
            except ValueError:
                continue
            if not np.isfinite(z):
                continue
            d[p[ix['SNP']]] = (p[ix['A1']].upper(), p[ix['A2']].upper(), z, float(p[ix['N']]))
    return d


# ------------------------------------------------------- core estimator
def _wls(x, y, w):
    """加权最小二乘 y ~ 1 + x，返回 (intercept, slope)"""
    X = np.column_stack([np.ones_like(x), x])
    W = w / np.mean(w)
    XtW = X.T * W
    beta = np.linalg.solve(XtW @ X, XtW @ y)
    return beta[0], beta[1]


def _weights_uni(l, c, intercept=1.0):
    """官方 LDSC 的异方差权重: Var(z_j^2) ~ 2 (intercept + c*l_j)^2

    注意：不要额外乘 1/l。官方的 1/w_ld 来自另一套 "LD weights" 文件，
    不是 l 本身；而且 l 可能为负（本面板 min = -0.26），乘 1/l 会让
    接近 0 或负的 l 拿到天文数字权重，一个点就毁掉整条回归。
    """
    return 1.0 / (2.0 * np.square(intercept + c * l))


def _weights_cross(l, c1, c2, intercept=1.0):
    return 1.0 / ((intercept + c1 * l) * (intercept + c2 * l))


def _mask_good(l):
    """官方 ldsc 同样只保留 l > 1 的 SNP"""
    return np.asarray(l, float) > 1.0


def fit_uni(l, z2, N, M, nblocks=200):
    l = np.asarray(l, float); z2 = np.asarray(z2, float)
    m = _mask_good(l) & np.isfinite(z2)
    l, y = l[m], z2[m] - 1.0
    h2 = 0.1
    for _ in range(3):                                   # 迭代求权
        c = N * h2 / M
        b = _wls(l, y, _weights_uni(l, c))
        h2 = b[1] * M / N
    slope = b[1]
    # jackknife（对 slope），再换算到 h2 尺度
    idx = np.array_split(np.arange(len(l)), nblocks)
    jk = np.empty(nblocks)
    for k, ids in enumerate(idx):
        mm = np.ones(len(l), bool); mm[ids] = False
        jk[k] = _wls(l[mm], y[mm], _weights_uni(l[mm], N * h2 / M))[1]
    se_slope = np.sqrt((nblocks - 1) * np.mean((jk - jk.mean()) ** 2))
    return dict(h2=h2, se=se_slope * M / N, slope=slope, se_slope=se_slope,
                intercept=1.0 + b[0])       # b[0] 是 y=Z^2-1 的截距，加回 1


def fit_pair(l, z1, z2, N1, N2, M, h1, h2, nblocks=200):
    l = np.asarray(l, float)
    z1 = np.asarray(z1, float); z2 = np.asarray(z2, float)
    m = _mask_good(l) & np.isfinite(z1) & np.isfinite(z2)
    l, z1, z2 = l[m], z1[m], z2[m]
    y = z1 * z2
    c1, c2 = N1 * h1 / M, N2 * h2 / M
    b = _wls(l, y, _weights_cross(l, c1, c2))
    idx = np.array_split(np.arange(len(l)), nblocks)
    jk = np.empty(nblocks)
    for k, ids in enumerate(idx):
        mm = np.ones(len(l), bool); mm[ids] = False
        jk[k] = _wls(l[mm], y[mm], _weights_cross(l[mm], c1, c2))[1]
    se = np.sqrt((nblocks - 1) * np.mean((jk - jk.mean()) ** 2))
    return dict(slope=b[1], se=se, intercept=b[0])


def delta_method_rg(s12, se12, s1, se1, s2, se2):
    """rg = s12/sqrt(s1 s2) 的一阶 SE（三个斜率由同一批 SNP 估得，保守起见忽略协方差）"""
    rg = s12 / np.sqrt(s1 * s2)
    d12 = 1.0 / np.sqrt(s1 * s2)
    d1 = -0.5 * s12 / (s1 * np.sqrt(s1 * s2))
    d2 = -0.5 * s12 / (s2 * np.sqrt(s1 * s2))
    var = (d12 * se12) ** 2 + (d1 * se1) ** 2 + (d2 * se2) ** 2
    return rg, np.sqrt(max(var, 0.0))


# ------------------------------------------------------------- selftest
def selftest(n=120000, M=1_200_000):
    rng = np.random.default_rng(20261029)
    l = np.abs(rng.normal(20, 12, n)) + 1.0
    N1, N2 = 400_000, 380_000
    h1, h2, rg_true = 0.05, 0.02, 0.45
    c1, c2 = N1 * h1 / M, N2 * h2 / M
    rho = rg_true * np.sqrt(h1 * h2)
    cov12 = rho * np.sqrt(N1 * N2) / M * l
    cov = np.empty((n, 2, 2))
    cov[:, 0, 0] = c1 * l + 1.0
    cov[:, 1, 1] = c2 * l + 1.0
    cov[:, 0, 1] = cov[:, 1, 0] = cov12
    z = np.array([rng.multivariate_normal([0, 0], cov[i]) for i in range(n)])
    u1 = fit_uni(l, z[:, 0] ** 2, N1, M)
    u2 = fit_uni(l, z[:, 1] ** 2, N2, M)
    pr = fit_pair(l, z[:, 0], z[:, 1], N1, N2, M, u1['h2'], u2['h2'])
    rg, rg_se = delta_method_rg(pr['slope'], pr['se'],
                                u1['h2'] * N1 / M, u1['se'] * N1 / M,
                                u2['h2'] * N2 / M, u2['se'] * N2 / M)
    print('--- selftest (n=%d SNPs) ---' % n)
    print('  h1^2  true=%.3f  est=%.4f  (SE %.4f)' % (h1, u1['h2'], u1['se']))
    print('  h2^2  true=%.3f  est=%.4f  (SE %.4f)' % (h2, u2['h2'], u2['se']))
    print('  rg    true=%.3f  est=%.4f  (SE %.4f)' % (rg_true, rg, rg_se))
    zsc = (rg - rg_true) / rg_se if rg_se > 0 else np.nan
    ok = abs(zsc) < 4
    print('  z = %.2f  -> %s' % (zsc, 'PASS' if ok else 'FAIL'))
    return 0 if ok else 1


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ld-dir')
    ap.add_argument('--ld-prefix', default='eur_w_ld_chr')
    ap.add_argument('--sumstats', nargs='*', default=[],
                    help='TAG=path 形式')
    ap.add_argument('--M', type=float, default=None,
                    help='参考面板 SNP 总数（默认 = 读到的 LD score 条数）')
    ap.add_argument('--out')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()

    if a.selftest:
        return selftest()

    if not a.ld_dir or not a.sumstats:
        ap.error('需要 --ld-dir 与 --sumstats（或用 --selftest）')

    print('reading LD scores ...', flush=True)
    ld = read_ldscores(a.ld_dir, a.ld_prefix)
    M = a.M if a.M else float(len(ld))
    print('  %d SNPs, M = %.0f' % (len(ld), M), flush=True)

    ss = {}
    for spec in a.sumstats:
        tag, path = spec.split('=', 1)
        ss[tag] = read_sumstats(path)
        print('  %-6s %d SNPs' % (tag, len(ss[tag])), flush=True)

    # 以第一个表型为参照，统一等位基因方向后取交集
    tags = list(ss)
    ref = tags[0]
    snps = sorted(set(ld) & set(ss[ref]))
    print('intersection with LD panel: %d' % len(snps), flush=True)

    Z, Ns, Nav = {}, {}, {}
    l = np.array([ld[s] for s in snps])
    Z[ref] = np.array([ss[ref][s][2] for s in snps])
    Nav[ref] = np.median([ss[ref][s][3] for s in snps])
    ref_alleles = {s: (ss[ref][s][0], ss[ref][s][1]) for s in snps}

    for t in tags[1:]:
        z = np.empty(len(snps)); n = np.empty(len(snps)); miss = 0
        for i, s in enumerate(snps):
            if s not in ss[t]:
                z[i] = np.nan; n[i] = np.nan; miss += 1; continue
            a1, a2, zi, ni = ss[t][s]
            r1, r2 = ref_alleles[s]
            if (a1, a2) == (r2, r1):        # 等位基因对调 -> Z 取反
                zi = -zi
            z[i] = zi; n[i] = ni
        keep = np.isfinite(z)
        print('  %-6s matched %d / %d' % (t, keep.sum(), len(snps)), flush=True)
        Z[t] = z[keep]; Nav[t] = np.nanmedian(n[keep])
        for tt in Z:
            if tt != t and len(Z[tt]) != keep.sum():
                Z[tt] = Z[tt][keep]
        l = l[keep]
        snps = [s for s, k in zip(snps, keep) if k]

    print('\n--- univariate h2 ---')
    uni = {}
    for t in tags:
        u = fit_uni(l, Z[t] ** 2, Nav[t], M)
        uni[t] = dict(h2=u['h2'], se=u['se'], slope=u['h2'] * Nav[t] / M,
                      slope_se=u['se'] * Nav[t] / M, N=Nav[t])
        ok = u['h2'] > 0 and u['h2'] / u['se'] >= 2.0     # h2 至少 2 倍 SE 才算可估
        uni[t]['ok'] = ok
        print('  %-6s N=%7.0f  h2=%.4f (SE %.4f)  z=%.1f  intercept=%.3f  %s'
              % (t, Nav[t], u['h2'], u['se'], u['h2'] / u['se'], u['intercept'],
                 '' if ok else '<-- underpowered, rg not estimable'))

    print('\n--- genetic correlation (rg) ---')
    rows = []
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            t1, t2 = tags[i], tags[j]
            # h2 <= 0 表示该性状在此样本量下 LDSC 无功效，rg 不可估（分母为 0/负）
            if not (uni[t1]['ok'] and uni[t2]['ok']):
                bad = [x for x in (t1, t2) if not uni[x]['ok']]
                print('  %-6s vs %-6s  rg = NA  (underpowered: %s)'
                      % (t1, t2, ', '.join(bad)))
                rows.append(dict(trait1=t1, trait2=t2, rg=np.nan, se=np.nan,
                                 z=np.nan, pval=np.nan,
                                 h2_1=uni[t1]['h2'], h2_2=uni[t2]['h2'],
                                 n_snps=len(snps), intercept=np.nan))
                continue
            p = fit_pair(l, Z[t1], Z[t2], Nav[t1], Nav[t2], M,
                         uni[t1]['h2'], uni[t2]['h2'])
            rg, se = delta_method_rg(p['slope'], p['se'],
                                     uni[t1]['slope'], uni[t1]['slope_se'],
                                     uni[t2]['slope'], uni[t2]['slope_se'])
            zsc = rg / se if se > 0 else np.nan
            pval = norm_p2(zsc) if se > 0 else np.nan
            print('  %-6s vs %-6s  rg=%+.3f (SE %.3f)  P=%.3g  intercept=%+.4f'
                  % (t1, t2, rg, se, pval, p['intercept']))
            rows.append(dict(trait1=t1, trait2=t2, rg=rg, se=se, z=zsc, pval=pval,
                             h2_1=uni[t1]['h2'], h2_2=uni[t2]['h2'],
                             n_snps=len(snps), intercept=p['intercept']))
    if a.out and rows:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        import csv
        with open(a.out, 'w', newline='', encoding='utf-8') as g:
            wcsv = csv.DictWriter(g, fieldnames=list(rows[0]))
            wcsv.writeheader(); wcsv.writerows(rows)
        print('\nwritten ->', a.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
