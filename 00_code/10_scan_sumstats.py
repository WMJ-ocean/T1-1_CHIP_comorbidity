# -*- coding: utf-8 -*-
"""
扫描所有 sumstats，输出：
  - 总 SNP 数、最小 p、各阈值下显著位点数、基因组膨胀因子 lambda
  - 导出 p < 5e-6 的子集（供 MR / clump 用）
FinnGen 列: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
"""
import gzip, os, sys, json, math, time

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
RAW = os.path.join(BASE, '02_raw')
OUTD = os.path.join(BASE, '04_clean')
os.makedirs(OUTD, exist_ok=True)
LOGP = os.path.join(BASE, '01_log', 'scan_sumstats.txt')

THRESH = 5e-6

def scan_finngen(path, tag, out):
    n = 0
    minp = 1.0
    c = {5e-8: 0, 1e-6: 0, 1e-5: 0, 1e-4: 0}
    chis = []
    hits = []
    t0 = time.time()
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').lstrip('#').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            parts = line.rstrip('\n').split('\t')
            n += 1
            try:
                p = float(parts[idx['pval']])
            except Exception:
                continue
            if p < minp:
                minp = p
            if p < 5e-8:
                c[5e-8] += 1
            if p < 1e-6:
                c[1e-6] += 1
            if p < 1e-5:
                c[1e-5] += 1
            if p < 1e-4:
                c[1e-4] += 1
            if n % 7 == 0:  # 抽样 1/7 估算 lambda
                try:
                    b = float(parts[idx['beta']]); se = float(parts[idx['sebeta']])
                    if se > 0:
                        chis.append((b / se) ** 2)
                except Exception:
                    pass
            if p < THRESH:
                hits.append('\t'.join([parts[idx['#chrom']] if '#chrom' in idx else parts[idx['chrom']],
                                       parts[idx['pos']], parts[idx['ref']], parts[idx['alt']],
                                       parts[idx['rsids']], parts[idx['pval']],
                                       parts[idx['beta']], parts[idx['sebeta']], parts[idx['af_alt']]]))
            if n % 5000000 == 0:
                out.write('    ... %d M lines (%.0fs)\n' % (n / 1e6, time.time() - t0))
                out.flush()
    lam = (sum(chis) / len(chis)) if chis else float('nan')
    dst = os.path.join(OUTD, '%s_p5e6.tsv' % tag)
    with open(dst, 'w', encoding='utf-8') as g:
        g.write('chrom\tpos\tref\talt\trsid\tpval\tbeta\tsebeta\taf\n')
        g.write('\n'.join(hits) + '\n')
    return dict(tag=tag, n=n, minp=minp, c5e8=c[5e-8], c1e6=c[1e-6], c1e5=c[1e-5],
                c1e4=c[1e-4], lam=lam, nhits=len(hits), secs=time.time() - t0)

def main():
    with open(LOGP, 'w', encoding='utf-8') as out:
        res = []
        for fn in sorted(os.listdir(RAW)):
            if not fn.endswith('.gz') or 'manifest' in fn:
                continue
            tag = fn.replace('.gz', '')
            p = os.path.join(RAW, fn)
            out.write('=== %s (%.0f MB) ===\n' % (fn, os.path.getsize(p) / 1e6))
            out.flush()
            r = scan_finngen(p, tag, out)
            res.append(r)
            out.write(json.dumps(r, ensure_ascii=False) + '\n\n')
            out.flush()
        out.write('\n=== SUMMARY ===\n')
        out.write('%-58s %-12s %-12s %-8s %-8s %-8s %s\n' % ('dataset', 'nSNP', 'minP', 'p<5e-8', 'p<1e-6', 'lambda', 'kept'))
        for r in res:
            out.write('%-58s %-12d %-12.3g %-8d %-8d %-8.3f %d\n' % (
                r['tag'][:58], r['n'], r['minp'], r['c5e8'], r['c1e6'], r['lam'], r['nhits']))
    print('done ->', LOGP)

if __name__ == '__main__':
    main()
