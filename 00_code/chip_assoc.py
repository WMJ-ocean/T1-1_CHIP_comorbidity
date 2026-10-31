# -*- coding: utf-8 -*-
"""用 GWAS Catalog solr association 检索 CHIP (GCST90165267) 显著位点"""
import requests, json, re
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

def show(name, u, n=900):
    try:
        r = requests.get(u, headers=UA, timeout=45, verify=False)
        print('--- %s code=%s' % (name, r.status_code), flush=True)
        if r.status_code != 200:
            print('    ', r.text[:180], flush=True)
            return None
        j = r.json()
        resp = j.get('response', {})
        docs = resp.get('docs', [])
        print('    numFound=%d  returned=%d' % (resp.get('numFound', 0), len(docs)), flush=True)
        if docs:
            print('    FIELDS:', sorted(docs[0].keys()), flush=True)
            print('    sample:', json.dumps(docs[0], ensure_ascii=False)[:900], flush=True)
        return j
    except Exception as e:
        print('--- %s ERR %s' % (name, str(e)[:100]), flush=True)
        return None

# 1) 探测 association 资源
show('assoc-all', 'https://www.ebi.ac.uk/gwas/api/search?q=*:*&rows=2&wt=json&resourcename=association')
# 2) 按 study 过滤（字段名需确认）
show('assoc-study', 'https://www.ebi.ac.uk/gwas/api/search?q=GCST90165267&rows=3&wt=json&resourcename=association')
# 3) V2 API 资源列表
try:
    r = requests.get('https://www.ebi.ac.uk/gwas/rest/api/v2/', headers=UA, timeout=30, verify=False)
    print('--- v2 root code=%s' % r.status_code, flush=True)
    j = r.json()
    print('    links:', list((j.get('_links') or {}).keys()), flush=True)
except Exception as e:
    print('v2 root ERR', str(e)[:80], flush=True)

# 4) 试探常见 association 端点
for p in ['https://www.ebi.ac.uk/gwas/rest/api/v2/associations?studyAccession=GCST90165267&size=3',
          'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90165267/snps?size=3',
          'https://www.ebi.ac.uk/gwas/rest/api/v2/singleNucleotidePolymorphisms?size=2']:
    try:
        r = requests.get(p, headers=UA, timeout=30, verify=False)
        print('    %-88s code=%s' % (p[-70:], r.status_code), flush=True)
        if r.status_code == 200:
            print('      ', r.text[:300], flush=True)
    except Exception as e:
        print('    ERR', p[-60:], str(e)[:60], flush=True)
print('DONE')
