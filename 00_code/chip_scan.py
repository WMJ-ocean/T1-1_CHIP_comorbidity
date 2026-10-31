# -*- coding: utf-8 -*-
"""扫描 Kessler 2022 (PMID 36450978) 全部子表型 accession"""
import requests, json, re, time, sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

found = []
for i in range(252, 340):
    acc = 'GCST90165%03d' % i
    try:
        r = requests.get('https://www.ebi.ac.uk/gwas/rest/api/v2/studies/%s' % acc, headers=UA, timeout=20, verify=False)
        if r.status_code != 200:
            continue
        j = r.json()
        if j.get('pubmed_id') != 36450978:
            continue
        dt = j.get('disease_trait', '')
        ss = (j.get('initial_sample_size') or '')[:60]
        full = j.get('full_summary_stats_available')
        nsnp = j.get('snp_count')
        na = j.get('snp_count')
        print('%s | full=%-5s | snp=%-10s | PMID=%s' % (acc, full, nsnp, j.get('pubmed_id')))
        print('     trait : %s' % dt)
        print('     sample: %s' % ss)
        print('     ftp   : %s' % (j.get('full_summary_stats') or '-'))
        found.append((acc, dt, full, nsnp, j.get('full_summary_stats')))
    except Exception:
        pass
    time.sleep(0.05)

print()
print('n studies =', len(found))
print()
print('=== ftp 文件名探测（对有 full stats 的）===')
for acc, dt, full, nsnp, ftp in found:
    if not full or not ftp:
        continue
    d = ftp.replace('http://', 'https://').replace('ftp.ebi.ac.uk', 'ftp.ebi.ac.uk') + '/'
    for attempt in [d, d.replace('/pub/', '/pub/')]:
        try:
            r = requests.get(attempt, headers=UA, timeout=25, verify=False)
            if r.status_code == 200:
                files = [f for f in re.findall(r'href="([^"?][^"]*)"', r.text) if '?' not in f and not f.startswith('/')]
                print('%s  %-55s -> %s' % (acc, dt[:55], files[:6]))
                break
            else:
                print('%s  dir code=%s' % (acc, r.status_code))
                break
        except Exception as e:
            print('%s ERR %s' % (acc, str(e)[:50]))
            break
print('DONE')
