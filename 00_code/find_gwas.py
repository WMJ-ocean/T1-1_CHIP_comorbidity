# -*- coding: utf-8 -*-
"""用 GWAS Catalog solr 搜索 API 定位目标 GWAS 的 accession"""
import requests, json, time
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

def solr(q, rows=40, extra=''):
    u = 'https://www.ebi.ac.uk/gwas/api/search?q=%s&rows=%d&wt=json&resourcename=study%s' % (
        requests.utils.quote(q), rows, extra)
    try:
        r = requests.get(u, headers=UA, timeout=40, verify=False)
        j = r.json()
        return j['response']['numFound'], j['response']['docs']
    except Exception as e:
        print('   ERR', str(e)[:80])
        return 0, []

QUERIES = [
    'clonal hematopoiesis',
    'clonal haematopoiesis',
    'CHIP clonal hematopoiesis of indeterminate potential',
    'myeloproliferative neoplasm',
    'aplastic anemia',
    'myelodysplastic syndrome',
    'immune thrombocytopenia',
]

for q in QUERIES:
    n, docs = solr(q, 30)
    print('=' * 100)
    print('QUERY: %-70s numFound=%d' % (q, n))
    for d in docs[:14]:
        acc = d.get('accessionId', '?')
        title = (d.get('title') or '')[:95]
        pubmed = d.get('pubmedId', '?')
        fp = d.get('fullPvalueSet', '?')
        nassoc = d.get('associationCount', '?')
        date = (d.get('publicationDate') or '')[:10]
        trait = d.get('traitName_s') or d.get('reportedTrait_s') or ''
        if isinstance(trait, list):
            trait = '; '.join(trait)[:60]
        print('  %s | fpv=%-5s | n=%-6s | %s | PMID %s' % (acc, fp, nassoc, date, pubmed))
        print('        T: %s' % title)
        print('        trait: %s' % trait)
    time.sleep(0.4)

print()
print('=== 定向验证已知 accession ===')
for acc in ['GCST90000032', 'GCST90018794', 'GCST90208187', 'GCST90208188']:
    u = 'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/%s' % acc
    try:
        r = requests.get(u, headers=UA, timeout=30, verify=False)
        if r.status_code == 200:
            j = r.json()
            print('%s OK  title=%s' % (acc, (j.get('title') or '')[:110]))
            print('     author=%s  date=%s  fpv=%s  nAssoc=%s' % (
                j.get('author'), j.get('publicationDate'), j.get('fullPvalueSet'), j.get('associationCount')))
        else:
            print('%s code=%s' % (acc, r.status_code))
    except Exception as e:
        print('%s ERR %s' % (acc, str(e)[:60]))
print('DONE')
