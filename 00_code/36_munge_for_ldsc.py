# -*- coding: utf-8 -*-
"""
36_munge_for_ldsc.py —— 把本项目 6 份 sumstats 转成 35_ldsc_rg.py 需要的格式
用法:
    python 36_munge_for_ldsc.py <snplist.txt> <TAG> [<TAG2> ...]
输出:
    09_ldsc/<TAG>.sumstats.tsv    列: SNP  A1  A2  Z  N

说明
    * snplist 用 eur_w_ld_chr 里的 SNP 列即可(即 HapMap3 集合), 不必另下 w_hm3:
        zcat eur_w_ld_chr/eur_w_ld_chr1.l2.ldscore.gz | cut -f2 > ...
    * CHIP: standard_error 列全 NA, 真实 SE 在 additional_info 的 REGENIE_SE=
    * FinnGen: N 自动从 03_meta/finngen_R12_manifest.tsv 读取
    * 效应等位 = alt (与 FinnGen / GWAS Catalog harmonised 口径一致)
"""
import gzip, os, re, sys, time

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
RAW  = os.path.join(BASE, '02_raw')
OUTD = os.path.join(BASE, '09_ldsc')
os.makedirs(OUTD, exist_ok=True)

FILES = {
    'MPN': 'finngen_R12_MYELOPROF_NONCML_MPN.gz',
    'PV':  'finngen_R12_POLYCYTVERA_PV.gz',
    'MDS': 'finngen_R12_C3_MYELODYSP_SYNDR_EXALLC_MDS.gz',
    'ITP': 'finngen_R12_D3_ITP_ITP.gz',
    'AA':  'finngen_R12_D3_OTHERAPLASTICANAEMIA_AA.gz',
    'CHIP': 'GCST90165267_CHIP_all.h.tsv.gz',
}
N_CHIP = 418401

re_se   = re.compile(r'REGENIE_SE=([-\d.eE+]+)')
re_beta = re.compile(r'REGENIE_BETA=([-\d.eE+]+)')


def finngen_N(tag):
    """从 03_meta/finngen_R12_manifest.tsv 读 num_cases + num_controls"""
    mf = os.path.join(BASE, '03_meta', 'finngen_R12_manifest.tsv')
    if not os.path.exists(mf):
        return None
    key = {'MPN': 'MYELOPROF_NONCML', 'PV': 'POLYCYTVERA',
           'MDS': 'MYELODYSP', 'ITP': 'D3_ITP', 'AA': 'OTHERAPLASTIC'}[tag]
    with open(mf, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        i_pc = hdr.index('phenocode')
        i_ca = hdr.index('num_cases'); i_co = hdr.index('num_controls')
        for ln in f:
            p = ln.rstrip('\n').split('\t')
            if len(p) > i_co and key in p[i_pc]:      # MDS/AA 的表型码带 C3_/D3_ 前缀
                try:
                    return float(p[i_ca]) + float(p[i_co])
                except ValueError:
                    return None
    return None


def munge(tag, want, log):
    path = os.path.join(RAW, FILES[tag])
    N = N_CHIP if tag == 'CHIP' else (finngen_N(tag) or 300000)
    out = os.path.join(OUTD, tag + '.sumstats.tsv')
    n_hit = 0
    t0 = time.time()
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f, \
         open(out, 'w', encoding='utf-8') as g:
        g.write('SNP\tA1\tA2\tZ\tN\n')
        hdr = f.readline().rstrip('\n').split('\t')
        if tag == 'CHIP':
            ix = {c: hdr.index(c) for c in
                  ['rsid', 'other_allele', 'effect_allele', 'p_value',
                   'additional_info', 'standard_error', 'odds_ratio']}
            i_ai = ix['additional_info']
        else:
            ix = {c: hdr.index(c) for c in ['rsids', 'ref', 'alt', 'beta', 'sebeta']}
        for ln in f:
            p = ln.rstrip('\n').split('\t')
            if tag == 'CHIP':
                rs = p[ix['rsid']]
                if rs not in want:
                    continue
                ai = p[i_ai]
                ms = re_se.search(ai); mb = re_beta.search(ai)
                if not ms:
                    continue
                se = float(ms.group(1))
                if mb:
                    b = float(mb.group(1))
                else:
                    try:
                        import math
                        b = math.log(float(p[ix['odds_ratio']]))
                    except (ValueError, IndexError):
                        continue
                if se <= 0:
                    continue
                a1, a2 = p[ix['effect_allele']], p[ix['other_allele']]
            else:
                rs = p[ix['rsids']]
                if rs not in want:
                    continue
                try:
                    b = float(p[ix['beta']]); se = float(p[ix['sebeta']])
                except ValueError:
                    continue
                if se <= 0:
                    continue
                a1, a2 = p[ix['alt']], p[ix['ref']]
            g.write('%s\t%s\t%s\t%.6f\t%.0f\n' % (rs, a1, a2, b / se, N))
            n_hit += 1
    msg = '%-5s N=%-7d hits=%-7d %.0fs -> %s' % (tag, N, n_hit, time.time() - t0, out)
    print(msg, flush=True); log.write(msg + '\n')


def main():
    if len(sys.argv) < 3:
        print(__doc__); return 1
    snpfile, tags = sys.argv[1], sys.argv[2:]
    want = set()
    for ln in open(snpfile, encoding='utf-8'):
        ln = ln.strip()
        if ln and ln != 'SNP':
            want.add(ln)
    print('snplist = %d' % len(want))
    log = open(os.path.join(BASE, '01_log', 'munge_ldsc.txt'), 'w', encoding='utf-8')
    for t in tags:
        munge(t, want, log)
    log.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
