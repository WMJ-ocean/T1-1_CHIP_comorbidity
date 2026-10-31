# -*- coding: utf-8 -*-
"""扫描 Kessler 2022 CHIP 子表型 + 列 FTP 文件大小，结果写文件"""
import requests, re, sys, os
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

OUT = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/01_log/chip_scan.txt'
f = open(OUT, 'w', encoding='utf-8')

def w(s=''):
    f.write(s + '\n')

studies = []
for i in range(252, 340):
    acc = 'GCST90165%03d' % i
    try:
        r = requests.get('https://www.ebi.ac.uk/gwas/rest/api/v2/studies/%s' % acc, headers=UA, timeout=15, verify=False)
        if r.status_code != 200:
            continue
        j = r.json()
        if j.get('pubmed_id') != 36450978:
            continue
        dt = j.get('disease_trait', '')
        ss = j.get('initial_sample_size') or ''
        full = j.get('full_summary_stats_available')
        nsnp = j.get('snp_count')
        ftp = j.get('full_summary_stats') or ''
        anc = ''
        m = re.search(r'(\d[\d,]*) (European|African|East Asian|South Asian)', ss)
        if m:
            anc = m.group(2)
        mc = re.search(r'([\d,]+) .*?cases', ss)
        ncase = mc.group(1) if mc else '?'
        studies.append(dict(acc=acc, trait=dt, full=full, nsnp=nsnp, ftp=ftp, anc=anc, ncase=ncase, ss=ss))
    except Exception:
        pass

w('=== Kessler 2022 (PMID 36450978) CHIP 子表型清单 ===')
w('%-15s %-6s %-8s %-11s %s' % ('acc', 'anc', 'ncase', 'nSNP', 'trait'))
for s in studies:
    w('%-15s %-6s %-8s %-11s %s' % (s['acc'], s['anc'], s['ncase'], s['nsnp'], s['trait'][:95]))

w()
w('=== FTP 文件名与大小（仅欧洲人群、非 gene burden）===')
sel = [s for s in studies if s['anc'] == 'European' and 'gene burden' not in s['trait'] and s['full']]
for s in sel:
    d = s['ftp'].replace('http://', 'https://') + '/'
    try:
        r = requests.get(d, headers=UA, timeout=25, verify=False)
        if r.status_code == 200:
            files = [x for x in re.findall(r'href="([^"?][^"]*)"', r.text) if '?' not in x and not x.startswith('/')]
            w('%-15s %-70s' % (s['acc'], s['trait'][:70]))
            for fn in files:
                fu = d + fn
                try:
                    h = requests.head(fu, headers=UA, timeout=20, verify=False)
                    sz = h.headers.get('Content-Length', '?')
                    mb = '%.0f MB' % (int(sz) / 1e6) if sz.isdigit() else sz
                    w('      %-42s %s' % (fn, mb))
                except Exception as e:
                    w('      %-42s ERR %s' % (fn, str(e)[:40]))
        else:
            w('%-15s dir code=%s' % (s['acc'], r.status_code))
    except Exception as e:
        w('%-15s ERR %s' % (s['acc'], str(e)[:60]))

w()
w('n studies total = %d ; selected EUR non-burden = %d' % (len(studies), len(sel)))
f.close()
print('written', OUT, os.path.getsize(OUT))
