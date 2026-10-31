# -*- coding: utf-8 -*-
"""61e_verify_il6.py —— 落盘 IL-6 pQTL 并用两种方式核对行数/信号，定位「0 显著」根因"""
import os, zlib, time, requests, urllib3
urllib3.disable_warnings()

D = r'D:\pqtl_tmp'
os.makedirs(D, exist_ok=True)
f = os.path.join(D, 'GCST90274815.tsv.gz')
URL = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST90274815/GCST90274815.tsv.gz'
UA = {'User-Agent': 'Mozilla/5.0 research'}

if not os.path.exists(f) or os.path.getsize(f) < 300_000_000:
    t0 = time.time(); n = 0
    with requests.get(URL, headers=UA, stream=True, timeout=300, verify=False) as r, open(f, 'wb') as fo:
        for c in r.iter_content(1 << 20):
            fo.write(c); n += len(c)
    print(f'[dl] {n/1e6:.1f} MB  {time.time()-t0:.0f}s', flush=True)
else:
    print('[dl] cached', os.path.getsize(f), flush=True)


def lines_file(path):
    """从本地文件多成员解压"""
    do = zlib.decompressobj(31); tail = b''
    with open(path, 'rb') as fi:
        while True:
            chunk = fi.read(1 << 20)
            if not chunk:
                break
            buf = do.decompress(chunk)
            while do.eof:
                nd = zlib.decompressobj(31)
                buf += nd.decompress(do.unused_data); do = nd
            tail += buf
            if b'\n' in tail:
                parts = tail.split(b'\n'); tail = parts.pop()
                for L in parts:
                    yield L


n = 0; mins = (1.0, None); cis = (1.0, None); nchr7 = 0
t0 = time.time()
for L in lines_file(f):
    if n == 0:
        print('[hdr]', L[:200].decode('utf-8', 'replace'), flush=True)
        n += 1; continue
    p = L.decode('utf-8', 'replace').rstrip('\r').split('\t')
    n += 1
    if len(p) < 11:
        continue
    try:
        pv = float(p[7])
    except ValueError:
        continue
    if pv < mins[0]:
        mins = (pv, p)
    if p[0] == '7':
        nchr7 += 1
        try:
            pos = int(p[1])
        except ValueError:
            continue
        if 22_600_000 <= pos <= 22_800_000 and pv < cis[0]:
            cis = (pv, p)
print(f'\n[reader] 行数 = {n:,}  (含表头)   用时 {time.time()-t0:.0f}s', flush=True)
print(f'[reader] 全文件最小 p = {mins[0]:.3g}  行 = {mins[1]}', flush=True)
print(f'[reader] chr7 行数 = {nchr7:,}', flush=True)
print(f'[reader] IL6 cis 区最小 p = {cis[0]:.3g}  行 = {cis[1]}', flush=True)
print('[SIZE]', os.path.getsize(f), flush=True)
