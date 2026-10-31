# -*- coding: utf-8 -*-
"""61g_bench_ebi.py —— 测 EBI pQTL 下载带宽（协议 / 并发）"""
import time, threading, requests, urllib3
urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research'}
P = ('ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
     'GCST90274001-GCST90275000/GCST90274839/GCST90274839.tsv.gz')


def one(url, s, n, out):
    h = dict(UA); h['Range'] = f'bytes={s}-{s+n-1}'
    r = requests.get(url, headers=h, stream=True, timeout=60, verify=False)
    got = 0
    for c in r.iter_content(1 << 16):
        got += len(c)
    r.close()
    out[0] += got


def bench(scheme, nth, mb=2):
    url = f'{scheme}://{P}'
    out = [0]; t0 = time.time()
    ts = [threading.Thread(target=one, args=(url, i * mb * 2**20, mb * 2**20, out))
          for i in range(nth)]
    for t in ts: t.start()
    for t in ts: t.join()
    el = time.time() - t0
    print(f'{scheme} x{nth:<2}: {out[0]/2**20:.1f} MB / {el:.1f}s = '
          f'{out[0]/2**20/max(el,.01):.2f} MB/s', flush=True)


if __name__ == '__main__':
    for sch, n in [('https', 8), ('http', 8), ('https', 16), ('http', 24)]:
        try:
            bench(sch, n)
        except Exception as e:
            print(f'{sch} x{n} ERR {type(e).__name__}: {str(e)[:80]}', flush=True)
