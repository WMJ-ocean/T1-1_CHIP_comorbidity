# -*- coding: utf-8 -*-
"""79_dnmt3a_mr2.py —— 补 ④（修正版）：DNMT3A 亚型 CHIP -> MPN / PV 的 MR

与 78 的差别（78 是照 GWAS Catalog 通用 schema 写的，实测列不符）：
  GCST90165271 实际列 = name(chr:pos:ref:alt，无 rsID) / chromosome / base_pair_location /
  other_allele / effect_allele / trait / cohort / model / odds_ratio / ci_lower / ci_upper /
  p_value / effect_allele_frequency / num_cases ... / additional_info
  → 无 beta / standard_error 列；真值在 additional_info 的 REGENIE_BETA= / REGENIE_SE=
  → 无法按 rsID 与总体 CHIP 的 21 个 IV 匹配，故改用 chr:pos 匹配 FinnGen

分阶段：scan（流式全扫，抽 p<5e-8）-> mr（修剪 + 匹配 + IVW）
"""
import io, os, sys, time, math, zlib
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, '02_raw')
IVD = os.path.join(BASE, '05_iv')
DEST = r'D:\pqtl_tmp\GCST90165271_buildGRCh38.tsv.gz'
FINN = {'MPN': 'finngen_R12_MYELOPROF_NONCML_MPN.gz',
        'PV': 'finngen_R12_POLYCYTVERA_PV.gz'}
FINN_DIR = os.path.join(RAW, 'sumstats')
PV_SIG = 5e-8
SOMATIC_BETA = 1.3
GWS_F = os.path.join(IVD, 'dnmt3a_chip_gws.tsv')
OUTF = os.path.join(BASE, '07_mr', 'DNMT3A_subtype_to_outcomes.tsv')


def lines(path):
    do = zlib.decompressobj(31); tail = b''
    with open(path, 'rb') as fi:
        while True:
            c = fi.read(1 << 22)
            if not c:
                break
            buf = do.decompress(c)
            while do.eof:
                nd = zlib.decompressobj(31)
                rest = do.unused_data
                buf += nd.decompress(rest); do = nd
            tail += buf
            if b'\n' in tail:
                ps = tail.split(b'\n'); tail = ps.pop()
                for L in ps:
                    yield L
    tail += do.flush()
    if tail:
        for L in tail.split(b'\n'):
            yield L


def parse_info(s):
    """additional_info -> (REGENIE_BETA, REGENIE_SE, MAC, INFO)"""
    b = se = mac = info = None
    for kv in s.split(';'):
        if '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        try:
            if k == 'REGENIE_BETA':
                b = float(v)
            elif k == 'REGENIE_SE':
                se = float(v)
            elif k == 'MAC':
                mac = float(v)
            elif k == 'INFO':
                info = float(v)
        except ValueError:
            pass
    return b, se, mac, info


def scan():
    ip = 11  # p_value
    n = ng = 0
    t0 = time.time()
    loci = []
    with io.open(GWS_F, 'w', encoding='utf-8') as fg:
        fg.write('snp\tchr\tpos\tother_allele\teffect_allele\tbeta\tse\tpval\teaf\tMAC\tINFO\n')
        for L in lines(DEST):
            n += 1
            if n == 1:
                continue  # header
            if (n % 5_000_000) == 0:
                print(f'   {n/1e6:.0f}M 行  GWS={ng}  {time.time()-t0:.0f}s', flush=True)
            p = L.split(b'\t')
            if len(p) < 22:
                continue
            try:
                pv = float(p[ip])
            except ValueError:
                continue
            if pv >= PV_SIG:
                continue
            name = p[0].decode()
            b, se, mac, info = parse_info(p[21].decode())
            if b is None or se is None:
                continue
            fg.write(f'{name}\t{p[1].decode()}\t{p[2].decode()}\t{p[3].decode()}\t'
                     f'{p[4].decode()}\t{b:.6g}\t{se:.6g}\t{pv:.4g}\t'
                     f'{p[12].decode()}\t{mac if mac is not None else ""}\t'
                     f'{info if info is not None else ""}\n')
            ng += 1
            loci.append((pv, name, b))
    print(f'[scan] 扫描 {n:,} 行；p<5e-8 且含 BETA/SE 共 {ng:,} 行；{time.time()-t0:.0f}s', flush=True)
    loci.sort()
    print('[scan] 最显著 15 个：', flush=True)
    for pv, nm, b in loci[:15]:
        print(f'   {nm:24s} p={pv:.3g} beta={b:+.3f}', flush=True)
    print('   ->', GWS_F, flush=True)


def read_tsv(p):
    with io.open(p, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for ln in f:
            if ln.strip():
                yield dict(zip(hdr, ln.rstrip('\n').split('\t')))


def _norm(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def ivw(bx, by, sy):
    bx, by, sy = map(np.asarray, (bx, by, sy))
    w = 1.0 / sy**2
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    se = math.sqrt(1.0 / np.sum(w * bx**2))
    k = len(bx)
    q = float(np.sum(w * (by - b * bx)**2))
    if k > 1 and q > (k - 1):
        se *= math.sqrt(q / (k - 1))
    z = b / se
    return dict(k=k, b=b, se=se, z=z, p=2 * (1 - _norm(abs(z))), q=q)


def prune(rows, win):
    rows = sorted(rows, key=lambda r: float(r['pval']))
    keep = []
    for r in rows:
        pos = int(r['pos'])
        if all(not (k['chr'] == r['chr'] and abs(int(k['pos']) - pos) < win) for k in keep):
            keep.append(r)
    return keep


def load_finn(fn, keys):
    """按 (chr,pos) 索引 FinnGen；返回 {(chr,pos): row}"""
    d = {}
    path = os.path.join(FINN_DIR, fn)
    if not os.path.exists(path):
        path = os.path.join(RAW, fn)
    import gzip
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        f.readline()
        for ln in f:
            c = ln.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            k = (c[0], c[1])
            if k in keys:
                d[k] = dict(ref=c[2], alt=c[3], beta=float(c[8]), se=float(c[9]), p=c[6])
    return d


def mr():
    rows = list(read_tsv(GWS_F))
    print(f'[mr] DNMT3A GWS 行 {len(rows):,}', flush=True)
    cand = []
    for r in rows:
        try:
            b = float(r['beta']); se = float(r['se'])
        except (KeyError, ValueError):
            continue
        if len(r['effect_allele']) != 1 or len(r['other_allele']) != 1:
            continue
        if abs(b) > SOMATIC_BETA or se <= 0:
            continue
        cand.append(r)
    print(f'[mr] 胚系候选（剔除 |beta|>{SOMATIC_BETA} 与 indel）{len(cand):,}', flush=True)
    if not cand:
        print('[mr] 无候选'); return

    results = {}
    for win in (10_000_000, 1_000_000):
        keep = prune(cand, win)
        print(f'[mr] 窗口 {win//1_000_000} Mb -> 工具 {len(keep)} 个: '
              f'{[r["snp"] for r in keep][:20]}', flush=True)
        keys = {(r['chr'], r['pos']) for r in keep}
        for phe, fn in FINN.items():
            d = load_finn(fn, keys)
            bx, by, sy, used = [], [], [], []
            for r in keep:
                o = d.get((r['chr'], r['pos']))
                if not o:
                    continue
                ea, oa = r['effect_allele'], r['other_allele']
                if o['alt'] == ea and o['ref'] == oa:
                    s = 1
                elif o['alt'] == oa and o['ref'] == ea:
                    s = -1
                else:
                    continue
                bx.append(float(r['beta'])); by.append(o['beta'] * s); sy.append(o['se'])
                used.append(r['snp'])
            tag = f'{phe}|{win//1_000_000}Mb'
            if len(bx) < 2:
                print(f'[mr] {tag}: 可用 {len(bx)}，跳过'); continue
            res = ivw(bx, by, sy)
            results[tag] = (res, used)
            print(f'[mr] DNMT3A -> {tag}: k={res["k"]} logOR={res["b"]:+.4f} '
                  f'({res["se"]:.4f}) p={res["p"]:.3g} OR={math.exp(res["b"]):.2f} Q={res["q"]:.2f}',
                  flush=True)

    with io.open(OUTF, 'w', encoding='utf-8') as f:
        f.write('exposure\toutcome\twindow\tk\tbeta\tse\tz\tp\tOR\tQ\tsnps\n')
        for tag, (r, used) in results.items():
            phe, win = tag.split('|')
            f.write(f'DNMT3A-CHIP\t{phe}\t{win}\t{r["k"]}\t{r["b"]:.6g}\t{r["se"]:.6g}\t'
                    f'{r["z"]:.4f}\t{r["p"]:.4g}\t{math.exp(r["b"]):.4f}\t{r["q"]:.4g}\t'
                    f'{",".join(used)}\n')
    print('[mr] ->', OUTF, flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if cmd in ('scan', 'all'):
        scan()
    if cmd in ('mr', 'all'):
        mr()
