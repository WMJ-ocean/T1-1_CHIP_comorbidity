# -*- coding: utf-8 -*-
"""62_dl_msigdb.py  下载 MSigDB v2024.1.Hs symbols GMT（后台）。"""
import requests, urllib3, time, os
urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
u = 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/msigdb.v2024.1.Hs.symbols.gmt'
out = r'C:\Users\Hero\WorkBuddy\2026-10-14-23-35-01\outputs\pv_hotspot\analysis\00_data\msigdb_v2024.1_Hs.gmt'
try:
    r = requests.get(u, headers=UA, stream=True, timeout=120, verify=False)
    n = 0
    with open(out, 'wb') as fo:
        for c in r.iter_content(1 << 18):
            if not c:
                continue
            fo.write(c); n += len(c)
    print(f'OK {n/1e6:.1f} MB -> {out}')
except Exception as e:
    print('ERR', type(e).__name__, e)
