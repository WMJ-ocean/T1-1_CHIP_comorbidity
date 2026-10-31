# -*- coding: utf-8 -*-
"""61f_dl_pqtl_mt.py —— pQTL 抽取 v2：先多线程落盘（校验字节数），再本地解析

为什么改：
  旧版 61_pqtl_extract.py 用 requests 单流 + iter_content 边下边解。即使修好了
  多成员 gzip，**HTTP 长连接仍会被中途静默掐断**（实测只读到文件前几 MB），
  于是 20 个蛋白有 18 个抽到 0 行、CXCL8 只抽到 2 个显著位点。
  改为：Range 8 线程分块落盘 -> 校验 size == Content-Length -> 本地多成员解压。

用法：
  python 61f_dl_pqtl_mt.py            # 全部 20 个蛋白
  python 61f_dl_pqtl_mt.py IL6 CXCL8  # 只跑指定蛋白
"""
import io, os, sys, time, zlib, threading
import requests, urllib3
urllib3.disable_warnings()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, '02_raw', 'pqtl')
META = os.path.join(BASE, '03_meta')
TMP = r'D:\pqtl_tmp'                       # 大文件落 D 盘
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

FTP = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
       'GCST90274001-GCST90275000/GCST{acc}/GCST{acc}.tsv.gz')
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research'}
NTH, CHUNK = 8, 8 << 20
PV_SIG = 5e-8

PANEL = [
    ('GCST90274815', 'IL6', 'Interleukin-6'),
    ('GCST90274839', 'TNF', 'Tumor necrosis factor'),
    ('GCST90274805', 'IL1A', 'Interleukin-1-alpha'),
    ('GCST90274817', 'CXCL8', 'Interleukin-8'),
    ('GCST90274821', 'CCL2', 'MCP-1'),
    ('GCST90274825', 'CCL3', 'MIP-1alpha'),
    ('GCST90274818', 'TGFB1_LAP', 'Latency-associated peptide TGF-beta 1'),
    ('GCST90274838', 'TGFA', 'TGF-alpha'),
    ('GCST90274826', 'MMP1', 'Matrix metalloproteinase-1'),
    ('GCST90274827', 'MMP10', 'Matrix metalloproteinase-10'),
    ('GCST90274848', 'VEGFA', 'VEGF-A'),
    ('GCST90274830', 'TNFRSF11B', 'Osteoprotegerin'),
    ('GCST90274794', 'IFNG', 'Interferon gamma'),
    ('GCST90274795', 'IL10', 'Interleukin-10'),
    ('GCST90274840', 'LTA', 'TNF-beta'),
    ('GCST90274801', 'IL17A', 'Interleukin-17A'),
    ('GCST90274803', 'IL18', 'Interleukin-18'),
    ('GCST90274831', 'OSM', 'Oncostatin-M'),
    ('GCST90274819', 'LIF', 'Leukemia inhibitory factor'),
    ('GCST90274833', 'KITLG', 'Stem cell factor'),
]

IVS = set()
with io.open(os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv'), encoding='utf-8') as f:
    for ln in f:
        p = ln.split('\t')[0].strip()
        if p.startswith('rs') and len(p) > 2:
            IVS.add(p)
print(f'[setup] CHIP 胚系 IV = {len(IVS)} 个   PV_SIG = {PV_SIG}', flush=True)


def download(acc, dest):
    """8 线程 Range 分块；区间半开 [s, e) 严格不重叠。返回 size 或 None"""
    acc = acc.replace('GCST', '')          # 模板里已含 GCST 前缀（勿写成 GCSTGCST）
    url = FTP.format(acc=acc)
    size = None
    for _ in range(3):
        try:
            r = requests.head(url, headers=UA, timeout=60, verify=False,
                              allow_redirects=True)
            if r.headers.get('Content-Length'):
                size = int(r.headers['Content-Length']); break
            print(f'    [head] status={r.status_code} keys={list(r.headers)[:6]}', flush=True)
        except Exception as e:
            print(f'    [head-exc] {type(e).__name__}: {str(e)[:90]}', flush=True)
        time.sleep(2)
    if size is None:                       # 有的镜像 HEAD 不给长度 -> 用 Range 探
        h = dict(UA); h['Range'] = 'bytes=0-0'
        r = requests.get(url, headers=h, stream=True, timeout=60, verify=False)
        cr = r.headers.get('Content-Range', '')
        r.close()
        if '/' in cr:
            size = int(cr.split('/')[-1])
    if not size:
        raise RuntimeError('no content-length')
    ok = dest + '.ok'
    if os.path.exists(ok) and os.path.exists(dest) and os.path.getsize(dest) == size:
        return size
    blocks = [(i, min(i + CHUNK, size)) for i in range(0, size, CHUNK)]
    lock = threading.Lock()
    done = [0]
    fails = []

    def work(blk):
        s, e = blk
        for att in range(6):
            try:
                h = dict(UA); h['Range'] = f'bytes={s}-{e - 1}'
                got = 0
                with requests.get(url, headers=h, stream=True, timeout=180,
                                  verify=False) as rr:
                    if rr.status_code not in (200, 206):
                        raise RuntimeError(f'HTTP {rr.status_code}')
                    with open(dest, 'r+b') as fo:
                        fo.seek(s)
                        for c in rr.iter_content(1 << 20):
                            fo.write(c); got += len(c)
                            if got >= e - s:
                                break
                if got != e - s:
                    raise RuntimeError(f'短读 {got}/{e-s}')
                with lock:
                    done[0] += got
                return
            except Exception as ex:
                with lock:
                    fails.append(str(ex)[:70])
                time.sleep(3 + att * 5)
        raise RuntimeError('block fail: ' + (fails[-1] if fails else '?'))

    with open(dest, 'wb') as fo:
        fo.truncate(size)
    t0 = time.time()
    q = list(blocks); idx = [0]
    def runner():
        while True:
            with lock:
                if idx[0] >= len(q):
                    return
                blk = q[idx[0]]; idx[0] += 1
            work(blk)
    ths = [threading.Thread(target=runner) for _ in range(NTH)]
    for t in ths: t.start()
    for t in ths: t.join()
    if done[0] != size or os.path.getsize(dest) != size:
        raise RuntimeError(f'byte mismatch {done[0]}/{size}')
    with open(ok, 'w') as fh:
        fh.write(str(size))
    print(f'    [dl] {size/1e6:.0f} MB  {time.time()-t0:.0f}s '
          f'({size/1e6/max(time.time()-t0,1):.1f} MB/s)', flush=True)
    return size


def lines_file(path):
    """多成员 gzip 本地读"""
    do = zlib.decompressobj(31); tail = b''
    with open(path, 'rb') as fi:
        while True:
            chunk = fi.read(1 << 22)
            if not chunk:
                break
            buf = do.decompress(chunk)
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


def one(acc, sym):
    f_iv = os.path.join(OUT, f'{sym}_step1_rows.tsv')
    f_ins = os.path.join(OUT, f'{sym}_instruments.tsv')
    dest = os.path.join(TMP, f'GCST{acc}.tsv.gz')
    t0 = time.time()
    try:
        download(acc, dest)
    except Exception as e:
        return sym, 'DLFAIL:' + str(e)[:40], 0, 0, 0, 1.0

    try:
        return _scan(dest, sym, f_iv, f_ins, t0)
    except Exception as e:
        try:
            os.replace(dest, dest + '.bad')
            if os.path.exists(dest + '.ok'):
                os.replace(dest + '.ok', dest + '.ok.bad')
        except Exception:
            pass
        return sym, 'SCANFAIL:' + type(e).__name__, 0, 0, 0, 1.0


def _scan(dest, sym, f_iv, f_ins, t0):
    it = lines_file(dest)
    hdr = next(it).decode('utf-8', 'replace').rstrip('\r\n').split('\t')
    gi = {n: i for i, n in enumerate(hdr)}
    for need in ('p_value', 'rsid'):
        if need not in gi:
            return sym, 'NOHDR:' + need, 0, 0, 0, 1.0
    ip, ir = gi['p_value'], gi['rsid']
    hit, ins, n = [], [], 0
    minp = (1.0, '')
    for L in it:
        n += 1
        s = L.decode('utf-8', 'replace').rstrip('\r\n')
        if (n & 0xFFFFF) == 0:
            print(f'      {sym} 已扫 {n/1e6:.1f}M 行  {time.time()-t0:.0f}s', flush=True)
        p = s.split('\t')
        if len(p) <= max(ip, ir):
            continue
        try:
            pv = float(p[ip])
        except ValueError:
            continue
        if pv < minp[0]:
            minp = (pv, p[ir])
        if p[ir] in IVS:
            hit.append(s)
        if pv < PV_SIG:
            ins.append(s)

    with io.open(f_iv, 'w', encoding='utf-8') as fo:
        fo.write('\t'.join(hdr) + '\n' + ('\n'.join(hit) + '\n' if hit else ''))
    with io.open(f_ins, 'w', encoding='utf-8') as fo:
        fo.write('\t'.join(hdr) + '\n' + ('\n'.join(ins) + '\n' if ins else ''))
    print(f'    [scan] {sym} 行={n:,}  IV命中={len(hit)}  显著={len(ins):,}  '
          f'最小p={minp[0]:.3g}@{minp[1]}  {time.time()-t0:.0f}s', flush=True)
    return sym, 'OK', len(hit), len(ins), n, minp[0]


if __name__ == '__main__':
    want = [a.upper() for a in sys.argv[1:]]
    todo = [p for p in PANEL if not want or p[1].upper() in want]
    if want:                                  # 按命令行给出的顺序执行（可指定优先级）
        todo.sort(key=lambda p: want.index(p[1].upper()))
    print(f'[setup] 待处理 {len(todo)} 个蛋白: {[t[1] for t in todo]}', flush=True)
    rows = []
    for acc, sym, _ in todo:
        print(f'  --- {sym} ({acc}) ---', flush=True)
        out = one(acc, sym)
        rows.append(out)
        print(f'  {sym:<10} {out[1]:<20} IV={out[2]:<3} sig={out[3]:<6} '
              f'rows={out[4]:,} minp={out[5]:.3g}', flush=True)
    with io.open(os.path.join(META, 'pqtl_extract_summary_v2.tsv'), 'w',
                 encoding='utf-8') as fo:
        fo.write('symbol\tstatus\tn_iv_rows\tn_sig\tn_rows\tmin_p\n')
        for r in rows:
            fo.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.3g}\n')
    print('\n完成', flush=True)
