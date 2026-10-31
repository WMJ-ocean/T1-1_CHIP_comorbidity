# -*- coding: utf-8 -*-
"""查 CHIP (Kessler 2022) 各 accession 的具体表型 + ftp 上是否有文件"""
import requests
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

ACCS = ['GCST90165%03d' % i for i in range(250, 275)]

print('=== CHIP study 详情 (GCST90165250-274) ===')
rows = []
for acc in ACCS:
    u = 'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/%s' % acc
    try:
        r = requests.get(u, headers=UA, timeout=30, verify=False)
        if r.status_code != 200:
            continue
        j = r.json()
        traits = j.get('reportedTraits') or []
        tn = '; '.join([t.get('trait', '') for t in traits]) if traits else (j.get('diseaseTrait', {}) or {}).get('trait', '?')
        ss = j.get('initialSampleSize', '')
        fpv = j.get('fullPvalueSet')
        na = j.get('associationCount')
        print('%s | fpv=%-5s | nAssoc=%-5s' % (acc, fpv, na))
        print('    trait : %s' % tn[:150])
        print('    title : %s' % (j.get('title') or '')[:110])
        print('    sample: %s' % ss[:170])
        rows.append((acc, tn, fpv, na))
    except Exception as e:
        pass

print()
print('=== ftp 文件探测 (GCST90165001-GCST90166000) ===')
for acc, tn, fpv, na in rows:
    if na and na >= 1:
        base = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/%s/' % acc
        try:
            r = requests.get(base, headers=UA, timeout=25, verify=False)
            if r.status_code == 200:
                import re
                files = re.findall(r'href="([^"?/][^"]*)"', r.text)
                files = [f for f in files if not f.startswith('/') and '?' not in f]
                print('%s (%s) -> %s' % (acc, tn[:40], files))
            else:
                print('%s code=%s' % (acc, r.status_code))
        except Exception as e:
            print('%s ERR %s' % (acc, str(e)[:50]))
print('DONE')
