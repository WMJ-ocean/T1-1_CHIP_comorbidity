# -*- coding: utf-8 -*-
"""
60_probe_pqtl.py
摸清 GWAS Catalog 上 Zheng 2020 Nat Genet 的 91 个炎症蛋白 pQTL:
  1) accession -> 蛋白名 (REST v2)
  2) 各文件的真实体积 / 是否支持 Range
"""
import requests, re, io, sys, time
sys.stdout.reconfigure(encoding='utf-8')
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research'}
FTP = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
REST = 'https://www.ebi.ac.uk/gwas/rest/api/studies/'

ACC = list(range(90274758, 90274849))          # 91 个
P = lambda m: print(m, flush=True)

# ---------- 1) 体积探测（抽查 3 个）+ 是否支持 Range ----------
P('=== 1. 文件体积 / Range 支持 ===')
size_map = {}
for acc in ACC[:3] + [ACC[-1]]:
    for sub in ['GCST90274001-GCST90275000']:
        u = f'{FTP}{sub}/GCST{acc}/GCST{acc}.tsv.gz'
        try:
            r = requests.get(u, headers=UA, stream=True, timeout=40)
            if r.status_code == 200:
                cl = r.headers.get('Content-Length')
                P(f'  {acc}: size={int(cl)/1e6 if cl else "?"} MB  range={r.headers.get("Accept-Ranges")}')
                if cl:
                    size_map[acc] = int(cl)
                r.close()
                break
            r.close()
        except Exception as e:
            P(f'  {acc}: ERR {type(e).__name__}')
print()

# ---------- 2) accession -> 蛋白名（REST v2 批量） ----------
P('=== 2. 获取蛋白名 ===')
names = {}
ok = fail = 0
for acc in ACC:
    u = f'{REST}GCST{acc}'
    try:
        r = requests.get(u, headers=UA, timeout=30)
        if r.status_code == 200:
            j = r.json()
            trait = j.get('diseaseTrait', {}).get('trait', '')
            title = j.get('initialSampleDescription', '') or ''
            names[acc] = trait
            ok += 1
        else:
            fail += 1
    except Exception:
        fail += 1
P(f'  成功 {ok} / 失败 {fail}')

import json
out = r'C:\Users\Hero\WorkBuddy\2026-10-14-23-35-01\outputs\pv_hotspot\T1-1_CHIP_comorbidity\03_meta\pqtl_map.json'
io.open(out, 'w', encoding='utf-8').write(json.dumps(names, ensure_ascii=False, indent=1))
P(f'  写入 {out}')
P('  前 25 个：')
for k, v in list(names.items())[:25]:
    P(f'    {k}  {v}')
