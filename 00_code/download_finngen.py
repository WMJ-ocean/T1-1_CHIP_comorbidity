# -*- coding: utf-8 -*-
"""下载 FinnGen R12 sumstats (MDS / AA / ITP / PV / ET)"""
import requests, os, sys, time
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

OUT = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw'
os.makedirs(OUT, exist_ok=True)

TARGETS = [
    ('C3_MYELODYSP_SYNDR_EXALLC', 'MDS'),
    ('D3_OTHERAPLASTICANAEMIA',   'AA'),
    ('D3_ITP',                    'ITP'),
]

for pheno, tag in TARGETS:
    url = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz' % pheno
    dst = os.path.join(OUT, 'finngen_R12_%s_%s.gz' % (pheno, tag))
    if os.path.exists(dst) and os.path.getsize(dst) > 100 * 1024 * 1024:
        print('SKIP exists: %s (%.0f MB)' % (os.path.basename(dst), os.path.getsize(dst) / 1e6))
        continue
    t0 = time.time()
    try:
        r = requests.get(url, headers=UA, timeout=120, verify=False, stream=True)
        total = int(r.headers.get('Content-Length', 0))
        print('DOWN %s  total=%.0f MB' % (tag, total / 1e6), flush=True)
        n = 0
        with open(dst, 'wb') as f:
            for chunk in r.iter_content(1 << 22):
                if not chunk:
                    continue
                f.write(chunk)
                n += len(chunk)
                if n % (200 * 1024 * 1024) < (1 << 22):
                    el = time.time() - t0
                    print('   %.0f MB / %.0f MB  (%.1f MB/s, %.0fs)' % (
                        n / 1e6, total / 1e6, n / el / 1e6, el), flush=True)
        el = time.time() - t0
        print('OK %s  %.0f MB in %.0fs  (%.2f MB/s)' % (tag, n / 1e6, el, n / el / 1e6), flush=True)
    except Exception as e:
        print('ERR %s: %s' % (tag, str(e)[:120]), flush=True)
print('ALL DONE')
