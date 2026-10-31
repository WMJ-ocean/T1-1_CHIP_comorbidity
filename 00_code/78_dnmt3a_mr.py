# -*- coding: utf-8 -*-
"""78_dnmt3a_mr.py —— 补 ④：用 DNMT3A 亚型 CHIP 工具变量复现 CHIP -> MPN / PV

背景：Kessler 2022 (PMID 36450978) 公开了分亚型 CHIP GWAS（同一批 43,187,107 SNP）：
  GCST90165271 = CHIP (DNMT3A mutation)，16,219 例 / 342,869 对照  ← 唯一有功效的亚型
  GCST90165281 = CHIP (TET2 mutation)，3,918 例（次要，可选）
  GCST90165273 = CHIP (JAK2 mutation)，186 例（功效为零，不可用于 MR）
若 DNMT3A 亚型工具变量也能复现 CHIP->MPN，则「独立于 JAK2」从
「我们的工具变量里没有 46/1」升级为「效应跨不同克隆起源一致」。

分阶段：
  dl    多线程下载 DNMT3A sumstats 到 D:（约 3.52 GiB）
  scan  流式扫描：抽出 (a) p<5e-8 的全部行、(b) 21 个总体 CHIP IV 所在行
  mr    选工具（10 Mb 距离修剪、剔除 |beta|>1.3 的体细胞读出）-> 从本地 FinnGen 取结局 beta -> IVW
"""
import io, os, sys, time, threading, math
import numpy as np
import requests, urllib3
urllib3.disable_warnings()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, '02_raw')
IVD = os.path.join(BASE, '05_iv')
OUT = os.path.join(BASE, '07_mr_protein')
TMP = r'D:\pqtl_tmp'
os.makedirs(TMP, exist_ok=True)

ACC = 'GCST90165271'
URL = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
       'GCST90165001-GCST90166000/GCST90165271/GCST90165271_buildGRCh38.tsv.gz')
DEST = os.path.join(TMP, f'{ACC}_buildGRCh38.tsv.gz')
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research'}
FINN = {'MPN': 'finngen_R12_MYELOPROF_NONCML_MPN.gz',
        'PV': 'finngen_R12_POLYCYTVERA_PV.gz'}
PV_SIG = 5e-8
PRUNE = 10_000_000
SOMATIC_BETA = 1.3


def download(nth=16, chunk=8 << 20):
    size = None
    for _ in range(3):
        try:
            r = requests.head(URL, headers=UA, timeout=60, verify=False)
            if r.headers.get('Content-Length'):
                size = int(r.headers['Content-Length']); break
        except Exception:
            pass
        time.sleep(2)
    if not size:
        raise RuntimeError('no content-length')
    print(f'[dl] remote {size/2**30:.2f} GiB', flush=True)
    ok = DEST + '.ok'
    if os.path.exists(ok) and os.path.getsize(DEST) == size:
        print('[dl] cached', flush=True); return size
    blocks = [(i, min(i + chunk, size)) for i in range(0, size, chunk)]
    lock = threading.Lock(); done = [0]; fails = []
    with open(DEST, 'wb') as f:
        f.truncate(size)

    def work(blk):
        s, e = blk
        for att in range(6):
            try:
                h = dict(UA); h['Range'] = f'bytes={s}-{e-1}'
                got = 0
                with requests.get(URL, headers=h, stream=True, timeout=180,
                                  verify=False) as rr:
                    if rr.status_code not in (200, 206):
                        raise RuntimeError(f'HTTP {rr.status_code}')
                    with open(DEST, 'r+b') as fo:
                        fo.seek(s)
                        for c in rr.iter_content(1 << 20):
                            fo.write(c); got += len(c)
                            if got >= e - s:
                                break
                if got != e - s:
                    raise RuntimeError(f'short {got}/{e-s}')
                with lock:
                    done[0] += got
                return
            except Exception as ex:
                with lock:
                    fails.append(str(ex)[:60])
                time.sleep(3 + att * 5)
        raise RuntimeError('block fail ' + (fails[-1] if fails else ''))

    idx = [0]; t0 = time.time()

    def runner():
        while True:
            with lock:
                if idx[0] >= len(blocks):
                    return
                b = blocks[idx[0]]; idx[0] += 1
            work(b)

    ths = [threading.Thread(target=runner) for _ in range(nth)]
    for t in ths: t.start()
    for t in ths: t.join()
    if done[0] != size:
        raise RuntimeError(f'byte mismatch {done[0]}/{size}')
    with open(ok, 'w') as fh:
        fh.write(str(size))
    print(f'[dl] done {size/2**30:.2f} GiB in {time.time()-t0:.0f}s', flush=True)
    return size


def lines(path):
    import zlib
    do = zlib.decompressobj(31); tail = b''
    with open(path, 'rb') as fi:
        while True:
            c = fi.read(1 << 22)
            if not c:
                break
            buf = do.decompress(c)
            while do.eof:
                nd = zlib.decompressobj(31)
                buf += nd.decompress(do.unused_data); do = nd
            tail += buf
            if b'\n' in tail:
                ps = tail.split(b'\n'); tail = ps.pop()
                for L in ps:
                    yield L
    tail += do.flush()
    if tail:
        for L in tail.split(b'\n'):
            yield L


def scan():
    ivs = set()
    for ln in io.open(os.path.join(IVD, 'CHIP_germline_iv.tsv'), encoding='utf-8').readlines()[1:]:
        c = ln.split('\t')
        if c and c[0].startswith('rs'):
            ivs.add(c[0])
    gws_f = os.path.join(IVD, 'dnmt3a_chip_gws.tsv')
    iv_f = os.path.join(IVD, 'dnmt3a_chip_at_iv.tsv')
    it = lines(DEST)
    hdr = next(it).decode('utf-8', 'replace').rstrip('\r\n').split('\t')
    gi = {n: i for i, n in enumerate(hdr)}
    print('[scan] header:', hdr, flush=True)
    ip, ir = gi['p_value'], gi['rsid']
    n = 0; ng = 0; niv = 0
    t0 = time.time()
    with io.open(gws_f, 'w', encoding='utf-8') as fg, \
            io.open(iv_f, 'w', encoding='utf-8') as fi2:
        fg.write('\t'.join(hdr) + '\n'); fi2.write('\t'.join(hdr) + '\n')
        for L in it:
            n += 1
            if (n % 5_000_000) == 0:
                print(f'   {n/1e6:.0f}M 行  GWS={ng}  IV={niv}  {time.time()-t0:.0f}s', flush=True)
            s = L.decode('utf-8', 'replace').rstrip('\r\n')
            p = s.split('\t')
            if len(p) <= max(ip, ir):
                continue
            rs = p[ir]
            if rs in ivs:
                fi2.write(s + '\n'); niv += 1
            try:
                pv = float(p[ip])
            except ValueError:
                continue
            if pv < PV_SIG:
                fg.write(s + '\n'); ng += 1
    print(f'[scan] 扫描 {n:,} 行；GWS {ng:,}；21 IV 命中 {niv}；{time.time()-t0:.0f}s', flush=True)
    print('   ->', gws_f, '/', iv_f, flush=True)


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


def mr():
    gws = list(read_tsv(os.path.join(IVD, 'dnmt3a_chip_gws.tsv')))
    print(f'[mr] DNMT3A GWS 行 {len(gws):,}', flush=True)
    # 1) 剔除体细胞读出（|beta|>1.3）与 indel
    inst = []
    for r in gws:
        try:
            b = float(r['beta']); se = float(r['standard_error'])
        except (KeyError, ValueError):
            continue
        if any(len(a) != 1 for a in (r.get('effect_allele', 'X'), r.get('other_allele', 'X'))):
            continue
        if abs(b) > SOMATIC_BETA or se <= 0:
            continue
        inst.append(r)
    print(f'[mr] 胚系候选 {len(inst):,}', flush=True)
    # 2) 10 Mb 距离修剪
    inst.sort(key=lambda r: float(r['p_value']))
    keep = []
    for r in inst:
        ch = r['chromosome']; pos = int(r['base_pair_location'])
        if all(not (k['chromosome'] == ch and abs(int(k['base_pair_location']) - pos) < PRUNE) for k in keep):
            keep.append(r)
    print(f'[mr] 距离修剪后工具 {len(keep)}: {[r["rsid"] for r in keep]}', flush=True)
    if not keep:
        print('[mr] 无可用工具，终止'); return
    want = {r['rsid'] for r in keep}
    out = {}
    for phe, fn in FINN.items():
        d = {}
        with io.open(os.path.join(RAW, fn), encoding='utf-8') as f:
            f.readline()
            for ln in f:
                c = ln.rstrip('\n').split('\t')
                if len(c) < 10:
                    continue
                ids = c[4]
                hit = [x for x in ids.split(',') if x in want] if ',' in ids else ([ids] if ids in want else [])
                if hit:
                    d[hit[0]] = dict(chrom=c[0], pos=c[1], ref=c[2], alt=c[3],
                                     beta=float(c[8]), se=float(c[9]), p=c[6])
        out[phe] = d
        print(f'   {phe}: 命中 {len(d)}/{len(want)}', flush=True)

    rows = []
    for phe in FINN:
        bx, by, sy, used = [], [], [], []
        for r in keep:
            o = out[phe].get(r['rsid'])
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
            used.append(r['rsid'])
        if len(bx) < 2:
            print(f'[mr] {phe}: 可用 {len(bx)}，跳过'); continue
        res = ivw(bx, by, sy)
        rows.append((phe, res, used))
        print(f'[mr] DNMT3A -> {phe}: k={res["k"]} logOR={res["b"]:+.4f} '
              f'({res["se"]:.4f}) p={res["p"]:.3g}  OR={math.exp(res["b"]):.2f}', flush=True)
    with io.open(os.path.join(BASE, '07_mr', 'DNMT3A_subtype_to_outcomes.tsv'),
                 'w', encoding='utf-8') as f:
        f.write('exposure\toutcome\tk\tbeta\tse\tz\tp\tOR\tsnps\n')
        for phe, r, used in rows:
            f.write(f'DNMT3A-CHIP\t{phe}\t{r["k"]}\t{r["b"]:.6g}\t{r["se"]:.6g}\t'
                    f'{r["z"]:.4f}\t{r["p"]:.4g}\t{math.exp(r["b"]):.4f}\t{",".join(used)}\n')
    print('[mr] -> 07_mr/DNMT3A_subtype_to_outcomes.tsv', flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if cmd in ('dl', 'all'):
        download()
    if cmd in ('scan', 'all'):
        scan()
    if cmd in ('mr', 'all'):
        mr()
