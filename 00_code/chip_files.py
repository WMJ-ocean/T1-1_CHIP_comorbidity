# -*- coding: utf-8 -*-
"""列 CHIP (GCST90165267 等) FTP 目录文件名 + 测 EBI 下载速度"""
import requests, re, time, sys
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
BASE = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/%s/'

for acc in ['GCST90165267', 'GCST90165261', 'GCST90165259', 'GCST90165270', 'GCST90165276', 'GCST90165282']:
    d = BASE % acc
    try:
        r = requests.get(d, headers=UA, timeout=25, verify=False)
        if r.status_code != 200:
            print('%s dir code=%s' % (acc, r.status_code), flush=True)
            continue
        files = [x for x in re.findall(r'href="([^"?][^"]*)"', r.text) if '?' not in x and not x.startswith('/')]
        print('%s files=%s' % (acc, files), flush=True)
        for fn in files:
            if fn.endswith('.gz') or fn.endswith('.tsv'):
                fu = d + fn
                h = requests.head(fu, headers=UA, timeout=25, verify=False)
                sz = h.headers.get('Content-Length', '?')
                print('    %-46s %s bytes' % (fn, sz), flush=True)
    except Exception as e:
        print('%s ERR %s' % (acc, str(e)[:70]), flush=True)

print()
print('=== EBI 下载测速 (12MB) ===', flush=True)
d = BASE % 'GCST90165267'
r0 = requests.get(d, headers=UA, timeout=25, verify=False)
files = [x for x in re.findall(r'href="([^"?][^"]*)"', r0.text) if x.endswith('.gz')]
if files:
    u = d + files[0]
    t0 = time.time()
    r = requests.get(u, headers=UA, timeout=90, verify=False, stream=True)
    n = 0
    for chunk in r.iter_content(1 << 20):
        n += len(chunk)
        if n >= 12 * 1024 * 1024:
            break
    el = time.time() - t0
    print('%s -> %.2f MB/s (%.1f MB in %.1fs)' % (files[0], n / el / 1e6, n / 1e6, el), flush=True)
print('DONE')
