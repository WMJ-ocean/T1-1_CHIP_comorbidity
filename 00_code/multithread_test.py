# -*- coding: utf-8 -*-
"""测试 EBI 多线程 Range 下载能否提速"""
import requests, time, threading
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
U = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165267/harmonised/GCST90165267.h.tsv.gz'

# 1) Range 支持?
r = requests.get(U, headers=dict(UA, Range='bytes=0-1048575'), timeout=60, verify=False, stream=True)
n = sum(len(c) for c in r.iter_content(65536))
print('Range test: code=%s  got=%d bytes  accept-ranges=%s' % (
    r.status_code, n, r.headers.get('Accept-Ranges', '?')), flush=True)
print('Content-Range:', r.headers.get('Content-Range', '?'), flush=True)

def dl_range(idx, start, size, results):
    end = start + size - 1
    h = dict(UA, Range='bytes=%d-%d' % (start, end))
    t0 = time.time()
    try:
        rr = requests.get(U, headers=h, timeout=120, verify=False, stream=True)
        b = 0
        for c in rr.iter_content(1 << 20):
            b += len(c)
        results[idx] = (b, time.time() - t0)
    except Exception as e:
        results[idx] = (-1, str(e)[:50])

for nthreads in [1, 4, 8]:
    results = {}
    ths = []
    each = 8 * 1024 * 1024
    t0 = time.time()
    for i in range(nthreads):
        # 不同偏移，避免缓存
        st = 100 * 1024 * 1024 + i * 20 * 1024 * 1024
        t = threading.Thread(target=dl_range, args=(i, st, each, results))
        t.start()
        ths.append(t)
    for t in ths:
        t.join()
    el = time.time() - t0
    tot = sum(v[0] for v in results.values() if isinstance(v[0], int))
    print('threads=%d  total=%.1f MB in %.1fs -> %.2f MB/s' % (nthreads, tot / 1e6, el, tot / el / 1e6), flush=True)
print('DONE')
