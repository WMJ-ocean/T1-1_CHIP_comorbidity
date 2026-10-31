# -*- coding: utf-8 -*-
"""实测 FinnGen 下载速度（只下前 12MB 用于测速，不落盘）"""
import requests, time
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

u = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_MYELODYSP_SYNDR_EXALLC.gz'
t0 = time.time()
n = 0
try:
    r = requests.get(u, headers=UA, timeout=90, verify=False, stream=True)
    print('code', r.status_code, 'total', r.headers.get('Content-Length'))
    for chunk in r.iter_content(1 << 20):
        n += len(chunk)
        el = time.time() - t0
        if n >= 12 * 1024 * 1024:
            break
        if el > 45:
            break
    el = time.time() - t0
    sp = n / el / 1024 / 1024
    print('downloaded %.1f MB in %.1f s -> %.2f MB/s' % (n / 1024 / 1024, el, sp))
    print('estimated time for 800MB: %.1f min' % (800 / sp / 60) if sp > 0 else 'NA')
except Exception as e:
    print('ERR', type(e).__name__, str(e)[:120])
print('DONE')
