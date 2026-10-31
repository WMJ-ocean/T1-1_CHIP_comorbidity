# -*- coding: utf-8 -*-
"""1) 查看 FinnGen 文件头  2) 拉 manifest 找 MPN/PV/ET/PMF endpoint"""
import gzip, os, requests, sys, re, io
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

RAW = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw'
LOG = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/01_log'
out = open(os.path.join(LOG, 'inspect.txt'), 'w', encoding='utf-8')
def w(s=''):
    out.write(s + '\n')

print('--- inspect_finngen start ---', flush=True)
for fn in sorted(os.listdir(RAW)):
    p = os.path.join(RAW, fn)
    w('=== %s  (%.0f MB) ===' % (fn, os.path.getsize(p) / 1e6))
    with gzip.open(p, 'rt') as f:
        for i, line in enumerate(f):
            w(line.rstrip('\n')[:600])
            if i >= 2:
                break
    w()

print('--- manifest ---', flush=True)
u = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv'
r = requests.get(u, headers=UA, timeout=120, verify=False)
txt = r.text
mp = os.path.join(RAW, '..', '03_meta', 'finngen_R12_manifest.tsv')
os.makedirs(os.path.dirname(mp), exist_ok=True)
open(mp, 'w', encoding='utf-8').write(txt)
w('manifest rows=%d' % len(txt.strip().split('\n')))

KEYS = ['MYELOPROLIF', 'POLYCYTHEMIA', 'POLYCYTHAEMIA', 'THROMBOCYTHEMIA', 'MYELOFIBROSIS',
        'MYELODYSP', 'APLASTIC', 'THROMBOCYTOPEN', 'MYELOID', 'LEUKEMIA', 'LEUKAEMIA',
        'HAEMATOPOIET', 'HEMATOPOIET', 'BLOOD']
w()
w('=== FinnGen R12 manifest 血液系统相关 endpoint ===')
lines = txt.strip().split('\n')
hdr = lines[0].split('\t')
w('header: %s' % hdr)
for ln in lines[1:]:
    parts = ln.split('\t')
    if len(parts) < 5:
        continue
    blob = (parts[0] + ' ' + parts[1]).upper()
    if any(k in blob for k in KEYS):
        w('%-42s %-58s cases=%-8s ctrl=%s' % (parts[0], parts[1][:58], parts[3], parts[4]))
out.close()
print('done', flush=True)
