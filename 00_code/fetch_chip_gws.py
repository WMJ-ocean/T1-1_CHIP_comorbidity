# -*- coding: utf-8 -*-
"""从 GWAS Catalog REST API v2 直接取 CHIP (GCST90165267) 的全基因组显著关联"""
import requests, json, sys, os
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
OUT = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/04_clean'
os.makedirs(OUT, exist_ok=True)

ACC = 'GCST90165267'

# 探测 associations 端点
base = 'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/%s/associations' % ACC
for q in ['?page=0&size=5', '?size=5', '']:
    r = requests.get(base + q, headers=UA, timeout=40, verify=False)
    print('--- %s -> code=%s' % (q or '(none)', r.status_code), flush=True)
    if r.status_code == 200:
        j = r.json()
        print('   KEYS:', sorted(j.keys()), flush=True)
        emb = j.get('_embedded', {})
        print('   EMB KEYS:', list(emb.keys()), flush=True)
        asso = emb.get('associations', [])
        print('   n=', len(asso), flush=True)
        if asso:
            print('   sample assoc:', json.dumps(asso[0], ensure_ascii=False)[:1200], flush=True)
        print('   page:', j.get('page'), flush=True)
        break
    else:
        print('   ', r.text[:200], flush=True)
print('DONE')
