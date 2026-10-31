# -*- coding: utf-8 -*-
"""探测 CHIP harmonised 目录与文件大小"""
import requests, re
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
B = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/%s/'

for acc in ['GCST90165267', 'GCST90165261', 'GCST90165259', 'GCST90165270']:
    d = (B % acc) + 'harmonised/'
    try:
        r = requests.get(d, headers=UA, timeout=25, verify=False)
        print('%s harmonised code=%s' % (acc, r.status_code), flush=True)
        if r.status_code == 200:
            files = [x for x in re.findall(r'href="([^"?][^"]*)"', r.text) if '?' not in x and not x.startswith('/')]
            for fn in files:
                fu = d + fn
                h = requests.head(fu, headers=UA, timeout=20, verify=False)
                sz = h.headers.get('Content-Length', '?')
                mb = '%.0f MB' % (int(sz) / 1e6) if sz.isdigit() else sz
                print('    %-52s %s' % (fn, mb), flush=True)
    except Exception as e:
        print('%s ERR %s' % (acc, str(e)[:70]), flush=True)

print()
print('=== README (GCST90165267) ===', flush=True)
r = requests.get((B % 'GCST90165267') + 'README.txt', headers=UA, timeout=25, verify=False)
print(r.text[:1500], flush=True)
print()
print('=== meta yaml ===', flush=True)
r = requests.get((B % 'GCST90165267') + 'GCST90165267_buildGRCh38.tsv.gz-meta.yaml', headers=UA, timeout=25, verify=False)
print(r.text[:2000], flush=True)
print('DONE')
