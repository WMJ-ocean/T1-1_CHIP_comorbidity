# -*- coding: utf-8 -*-
"""下载 CHIP harmonised (HapMap3) sumstats"""
import requests, os, time, re, sys
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}
RAW = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw'
B = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/%s/'
LOG = open(r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/01_log/dl_chip.log', 'w', encoding='utf-8')
def w(s=''):
    LOG.write(s + '\n'); LOG.flush()

# 候选文件名（GWAS Catalog harmonised 命名规则）
CAND = ['harmonised/%s.h.tsv.gz']

TARGETS = [
    ('GCST90165267', 'CHIP_all'),       # CHIP (with/without mCA), EUR, 25657 cases
    ('GCST90165261', 'CHIP_nomCA'),     # CHIP (without mCA), EUR, 20606 cases
]

for acc, tag in TARGETS:
    d = B % acc
    # 先列 harmonised 目录
    try:
        r = requests.get(d + 'harmonised/', headers=UA, timeout=40, verify=False)
        files = [x for x in re.findall(r'href="([^"?][^"]*)"', r.text) if '?' not in x and not x.startswith('/')]
        w('%s harmonised files=%s' % (acc, files))
    except Exception as e:
        w('%s list ERR %s' % (acc, str(e)[:80]))
        files = []

    got = False
    for fn in files:
        if not fn.endswith('.gz'):
            continue
        fu = d + 'harmonised/' + fn
        try:
            h = requests.head(fu, headers=UA, timeout=30, verify=False)
            sz = int(h.headers.get('Content-Length', 0))
        except Exception as e:
            w('  HEAD ERR %s %s' % (fn, str(e)[:50]))
            continue
        w('  -> %s  %.0f MB' % (fn, sz / 1e6))
        if 'h.tsv' not in fn:
            continue
        dst = os.path.join(RAW, '%s_%s.h.tsv.gz' % (acc, tag))
        t0 = time.time()
        try:
            rr = requests.get(fu, headers=UA, timeout=180, verify=False, stream=True)
            n = 0
            with open(dst, 'wb') as f:
                for ch in rr.iter_content(1 << 22):
                    if not ch:
                        continue
                    f.write(ch); n += len(ch)
            el = time.time() - t0
            w('  OK %s %.0f MB in %.0fs (%.1f MB/s)' % (fn, n / 1e6, el, n / el / 1e6))
            got = True
        except Exception as e:
            w('  DL ERR %s' % str(e)[:80])
    if not got:
        w('%s 未拿到 harmonised 文件' % acc)
w('DONE')
