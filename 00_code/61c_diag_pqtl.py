# -*- coding: utf-8 -*-
"""61c_diag_pqtl.py —— 诊断 pQTL 抽取为 0 的根因（只读前 30 MB，不落盘）"""
import io, zlib, requests, urllib3
urllib3.disable_warnings()

FTP = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST{acc}/GCST{acc}.tsv.gz'
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research'}
TARGETS = [('90274826', 'MMP1'), ('90274815', 'IL6'), ('90274817', 'CXCL8')]


def lines_from(resp, max_lines=None):
    do = zlib.decompressobj(31); tail = ''
    n = 0
    for chunk in resp.iter_content(1 << 18):
        buf = do.decompress(chunk)
        while do.eof:
            nd = zlib.decompressobj(31)
            buf += nd.decompress(do.unused_data); do = nd
        tail += buf.decode('utf-8', 'replace')
        if '\n' in tail:
            parts = tail.split('\n'); tail = parts.pop()
            for L in parts:
                yield L; n += 1
                if max_lines and n >= max_lines:
                    return
    tail += do.flush().decode('utf-8', 'replace')
    if tail:
        yield tail


for acc, sym in TARGETS:
    url = FTP.format(acc=acc)
    r = requests.get(url, headers=UA, stream=True, timeout=180, verify=False)
    hdr = None; rows = 0; sig = []; colnames = []
    for i, L in enumerate(lines_from(r)):
        if i == 0:
            hdr = L.split('\t'); continue
        p = L.split('\t'); rows += 1
        if len(p) < 6: continue
        # 按列名定位
        if not colnames:
            colnames = {c: k for k, c in enumerate(hdr)}
            print(f'\n=== {sym} ({acc}) ===')
            print('  header:', hdr)
            print('  rsid idx =', colnames.get('rsid'), ' p idx =', colnames.get('p_value'))
        try:
            pv = float(p[colnames['p_value']])
        except Exception:
            continue
        if pv < 5e-8:
            sig.append((round(pv, 12), p[colnames.get('rsid', 0)]))
        if rows >= 300000:
            break
    r.close()
    sig.sort()
    print(f'  读到 {rows} 行（前 30 MB 内）; p<5e-8 = {len(sig)}; 最小 p 前3 = {sig[:3]}')
