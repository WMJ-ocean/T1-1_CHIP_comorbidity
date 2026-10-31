#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
37_prep_ldsc_panel.py
解压 eur_w_ld_chr 参考面板，并从 LD score 文件本身反推 HapMap3 snplist。

背景：官方 w_hm3.snplist.bz2 与 eur_w_ld_chr.tar.bz2 同在 Broad requester-pays 桶
（本网络不可得）。但 eur_w_ld_chr 本身就是 HapMap3 子集，其 .l2.ldscore 文件的
第一列即 rsID —— 可直接作为 snplist 使用，从而绕开 w_hm3.snplist。

用法:
    python 37_prep_ldsc_panel.py <tar.gz 路径> <输出目录>
"""
import gzip
import io
import os
import sys
import tarfile

sys.stdout.reconfigure(encoding='utf-8', errors='replace')


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    tgz, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    print('extracting %s -> %s' % (tgz, outdir), flush=True)
    with tarfile.open(tgz, 'r:gz') as tf:
        # 防路径穿越：只接受普通文件，剥离 ../ 前缀
        members = []
        for m in tf.getmembers():
            if not m.isfile():
                continue
            m.name = os.path.basename(m.name.replace('\\', '/'))
            if m.name:
                members.append(m)
        tf.extractall(outdir, members=members)
    print('  extracted %d members' % len(members), flush=True)

    # 找 1..22 的 l2.ldscore
    found = []
    for c in range(1, 23):
        for cand in ['%d.l2.ldscore.gz' % c, '%d.l2.ldscore' % c]:
            p = os.path.join(outdir, cand)
            if os.path.exists(p):
                found.append((c, p))
                break
    print('  found %d chromosome ldscore files' % len(found), flush=True)
    if not found:
        print('ERROR: no *.l2.ldscore files found; contents:')
        for f in sorted(os.listdir(outdir))[:40]:
            print('   ', f)
        sys.exit(2)

    # 检查是否带 M / ldscore 多列
    with gzip.open(found[0][1], 'rt') if found[0][1].endswith('.gz') \
            else io.open(found[0][1], encoding='utf-8') as f:
        head = [next(f).rstrip('\n').split('\t') for _ in range(3)]
    print('  header/rows sample: %s' % (head[0][:6],), flush=True)
    print('  ncols = %d' % len(head[0]), flush=True)

    snp_out = os.path.join(outdir, 'hm3_snplist.txt')
    n = 0
    per_chr = {}
    with open(snp_out, 'w', encoding='utf-8') as g:
        for c, p in found:
            cnt = 0
            op = gzip.open(p, 'rt', encoding='utf-8', errors='replace') \
                if p.endswith('.gz') else io.open(p, encoding='utf-8', errors='replace')
            with op as f:
                for ln in f:
                    if not ln.strip():
                        continue
                    p = ln.rstrip('\n').split('\t')
                    if len(p) < 2:
                        continue
                    if p[0].strip().upper() == 'CHR':      # header 行
                        continue
                    # 列顺序: CHR  SNP  BP  CM  MAF  L2  -> rsID 在第 2 列
                    rs = p[1].strip()
                    if not rs or rs.upper() in ('SNP', 'RSID', 'RSID1', 'NA', '.'):
                        continue
                    g.write(rs + '\n')
                    cnt += 1
            per_chr[c] = cnt
            n += cnt
            print('    chr%-3d %8d SNPs' % (c, cnt), flush=True)
    print('TOTAL %d SNPs -> %s' % (n, snp_out), flush=True)


if __name__ == '__main__':
    main()
