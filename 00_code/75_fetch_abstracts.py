# -*- coding: utf-8 -*-
"""75_fetch_abstracts.py —— 取关键文献摘要，核对「CHIP 与 MPN」的准确表述与数字"""
import requests, urllib3, re, textwrap
urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0 research'}
EPMC = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'

PMIDS = ['25931582', '32736382', '37221411', '32731609', '24321397', '24325356']
for pm in PMIDS:
    r = requests.get(EPMC, params={'query': f'EXT_ID:{pm} AND SRC:MED', 'format': 'json',
                                   'resultType': 'core'}, headers=UA, timeout=60)
    res = r.json().get('resultList', {}).get('result', [])
    if not res:
        print(f'--- PMID {pm}: NOT FOUND'); continue
    d = res[0]
    j = d.get('journalInfo', {}) or {}
    print('=' * 100)
    print(f"PMID {pm} | {d.get('authorString','')[:60]}")
    print(f"  {(d.get('title') or '')[:150]}")
    print(f"  {j.get('journal',{}).get('medlineAbbreviation') or j.get('journal',{}).get('title')} "
          f"{j.get('yearOfPublication')};{j.get('volume')}({j.get('issue')}):{d.get('pageInfo')}")
    ab = re.sub(r'<[^>]+>', '', d.get('abstractText') or '')
    for para in textwrap.wrap(ab[:1100], 110):
        print('   ', para)
