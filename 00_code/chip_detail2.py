# -*- coding: utf-8 -*-
import requests, json, re
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

print('=== raw keys of GCST90165257 ===')
r = requests.get('https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90165257', headers=UA, timeout=30, verify=False)
j = r.json()
print('KEYS:', sorted(j.keys()))
print()
print(json.dumps(j, ensure_ascii=False, indent=1)[:2500])
