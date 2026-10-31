# -*- coding: utf-8 -*-
"""76_chip_subtypes.py —— 核查 Kessler 2022 是否公开了 CHIP 分亚型 sumstats"""
import requests, urllib3
urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0 research'}
API = 'https://www.ebi.ac.uk/gwas/rest/api/v2/studies'

r = requests.get(f'{API}/GCST90165267', headers=UA, timeout=60, verify=False).json()
pm = r.get('pubmed_id')
print('overall CHIP GWAS: GCST90165267 | pubmed:', pm,
      '| N:', r.get('initial_sample_size'))

r2 = requests.get(API, params={'pubmedId': pm, 'size': 200}, headers=UA,
                  timeout=90, verify=False).json()
items = r2.get('_embedded', {}).get('studies', [])
print('same-PMID studies:', len(items))

KEY = ('DNMT3A', 'TET2', 'ASXL1', 'JAK2', 'SF3B1', 'SRSF2', 'TP53', 'PPM1D',
       'spliceosome', 'mosaic')
hits = []
for s in items:
    t = s.get('disease_trait') or ''
    if any(k.lower() in t.lower() for k in KEY):
        hits.append((s.get('accession_id'), t, s.get('initial_sample_size')))
print(f'\nsubtype-specific CHIP studies: {len(hits)}')
for a, t, n in hits:
    print(f'  {a} | {t[:80]} | {(n or "")[:55]}')
