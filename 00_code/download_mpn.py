# -*- coding: utf-8 -*-
"""下载 FinnGen R12 MPN 相关 endpoint"""
import requests, os, time
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
OUT = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw'
os.makedirs(OUT, exist_ok=True)

TARGETS = [
    ('MYELOPROF_NONCML', 'MPN'),
    ('POLYCYTVERA', 'PV'),
]
for pheno, tag in TARGETS:
    url = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz' % pheno
    dst = os.path.join(OUT, 'finngen_R12_%s_%s.gz' % (pheno, tag))
    if os.path.exists(dst) and os.path.getsize(dst) > 100 * 1024 * 1024:
        print('SKIP %s' % os.path.basename(dst), flush=True)
        continue
    t0 = time.time()
    r = requests.get(url, headers=UA, timeout=120, verify=False, stream=True)
    total = int(r.headers.get('Content-Length', 0))
    print('DOWN %s %.0f MB' % (tag, total / 1e6), flush=True)
    n = 0
    with open(dst, 'wb') as f:
        for chunk in r.iter_content(1 << 22):
            if not chunk:
                continue
            f.write(chunk)
            n += len(chunk)
    el = time.time() - t0
    print('OK %s %.0f MB in %.0fs (%.1f MB/s)' % (tag, n / 1e6, el, n / el / 1e6), flush=True)
print('DONE')
