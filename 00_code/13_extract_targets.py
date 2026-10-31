# -*- coding: utf-8 -*-
"""
按 rsID 列表从 sumstats (gz) 中提取目标位点
用法: python 13_extract_targets.py <snp_list.txt> <tag1,tag2,...>
输出: 06_targets/<tag>_targets.tsv  (chrom pos ref alt rsid pval beta sebeta af)
"""
import gzip, os, sys, time

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
RAW = os.path.join(BASE, '02_raw')
OUTD = os.path.join(BASE, '06_targets')
os.makedirs(OUTD, exist_ok=True)

# 文件 tag -> 路径
FILES = {
    'MDS': 'finngen_R12_C3_MYELODYSP_SYNDR_EXALLC_MDS.gz',
    'ITP': 'finngen_R12_D3_ITP_ITP.gz',
    'AA':  'finngen_R12_D3_OTHERAPLASTICANAEMIA_AA.gz',
    'MPN': 'finngen_R12_MYELOPROF_NONCML_MPN.gz',
    'PV':  'finngen_R12_POLYCYTVERA_PV.gz',
}

def main():
    snpfile = sys.argv[1]
    tags = sys.argv[2].split(',') if len(sys.argv) > 2 else list(FILES)
    want = set()
    for ln in open(snpfile, encoding='utf-8'):
        ln = ln.strip()
        if not ln or ln.startswith('SNP') or ln.startswith('#'):
            continue
        want.add(ln.split()[0])
    print('target SNPs =', len(want), flush=True)

    for tag in tags:
        fn = FILES.get(tag)
        if not fn:
            print('unknown tag', tag); continue
        p = os.path.join(RAW, fn)
        t0 = time.time()
        found = {}
        n = 0
        with gzip.open(p, 'rt') as f:
            header = f.readline().rstrip('\n').lstrip('#').split('\t')
            ix = {c: i for i, c in enumerate(header)}
            ic = ix.get('chrom') if 'chrom' in ix else ix.get('#chrom')
            for line in f:
                n += 1
                # 快速预筛：rsID 列内容至少一个字符
                parts = line.rstrip('\n').split('\t')
                rs = parts[ix['rsids']]
                if not rs or rs == '.':
                    continue
                hits = [s for s in rs.split(',') if s in want]
                if not hits:
                    continue
                for s in hits:
                    if s in found:
                        continue
                    found[s] = [parts[ic], parts[ix['pos']], parts[ix['ref']], parts[ix['alt']],
                                s, parts[ix['pval']], parts[ix['beta']], parts[ix['sebeta']],
                                parts[ix['af_alt']]]
        dst = os.path.join(OUTD, '%s_targets.tsv' % tag)
        with open(dst, 'w', encoding='utf-8') as g:
            g.write('chrom\tpos\tref\talt\trsid\tpval\tbeta\tsebeta\taf\n')
            for s in sorted(found):
                g.write('\t'.join(found[s]) + '\n')
        print('%s: scanned %d lines in %.0fs, found %d / %d' % (
            tag, n, time.time() - t0, len(found), len(want)), flush=True)

if __name__ == '__main__':
    main()
