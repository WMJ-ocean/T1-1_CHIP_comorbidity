# -*- coding: utf-8 -*-
"""74_verify_chip_mpn.py —— 核实「CHIP 与 MPN 共病」这一说法的文献依据"""
import requests, urllib3
urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0 research'}
EPMC = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'


def q(label, query, n=6):
    r = requests.get(EPMC, params={'query': query, 'format': 'json',
                                   'resultType': 'lite', 'pageSize': n,
                                   'sort': 'CITED desc'}, headers=UA, timeout=60)
    print('##', label)
    for d in r.json().get('resultList', {}).get('result', []):
        print(f"   PMID {d.get('pmid'):>9} {d.get('pubYear')} "
              f"{(d.get('journalTitle') or '')[:20]:22s} | {(d.get('title') or '')[:100]}")


q('A. CHIP 定义（无显性血液病）',
  'TITLE:"clonal hematopoiesis of indeterminate potential" AND AUTH:"Steensma"')
q('B. CHIP 与血液肿瘤风险（前瞻队列）',
  'TITLE:"clonal hematopoiesis" AND TITLE:"hematologic" AND TITLE:"risk"', 8)
q('C. CHIP 与 MPN 同时出现的直接研究',
  'TITLE:"clonal hematopoiesis" AND (TITLE:"myeloproliferative" OR TITLE:"myeloproliferative neoplasm")', 8)
q('D. MPN 患者中 CHIP 型突变（DNMT3A/TET2/ASXL1）共存',
  '(TITLE:"myeloproliferative") AND (ABSTRACT:"DNMT3A" AND ABSTRACT:"TET2" AND ABSTRACT:"ASXL1") AND TITLE:"clonal"', 8)
q('E. JAK2 V617F 在一般人群中的前驱性',
  'TITLE:"JAK2 V617F" AND (TITLE:"population" OR TITLE:"general population" OR TITLE:"prevalence")', 8)
