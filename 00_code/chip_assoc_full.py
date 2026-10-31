# -*- coding: utf-8 -*-
"""拉取 GCST90165267 在 GWAS Catalog 中的全部 curated associations"""
import requests, json, os
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
OUT = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/04_clean'
os.makedirs(OUT, exist_ok=True)

ACC = 'GCST90165267'
allrows = []
page = 0
while True:
    u = 'https://www.ebi.ac.uk/gwas/rest/api/v2/associations?studyAccession=%s&size=1000&page=%d' % (ACC, page)
    r = requests.get(u, headers=UA, timeout=60, verify=False)
    if r.status_code != 200:
        print('page %d code=%s  %s' % (page, r.status_code, r.text[:150]), flush=True)
        break
    j = r.json()
    emb = j.get('_embedded', {}).get('associations', [])
    pg = j.get('page', {})
    allrows.extend(emb)
    print('page %d: got %d (total so far %d); pageInfo=%s' % (
        page, len(emb), len(allrows), {k: pg.get(k) for k in ['size', 'totalElements', 'totalPages', 'number']}), flush=True)
    if len(emb) < 1000 or pg.get('number', 0) >= pg.get('totalPages', 1) - 1:
        break
    page += 1
    if page > 30:
        break

print()
print('TOTAL associations =', len(allrows), flush=True)
if allrows:
    print('FIELDS:', sorted(allrows[0].keys()), flush=True)
    print()
    print(json.dumps(allrows[0], ensure_ascii=False, indent=1)[:2200], flush=True)
    ps = [a.get('p_value') for a in allrows if a.get('p_value') is not None]
    if ps:
        print()
        print('min p =', min(ps), ' n with p<5e-8 =', sum(1 for p in ps if p < 5e-8))
    json.dump(allrows, open(os.path.join(OUT, 'CHIP_assoc_raw.json'), 'w', encoding='utf-8'), ensure_ascii=False)
    print('saved ->', os.path.join(OUT, 'CHIP_assoc_raw.json'), flush=True)
print('DONE')
