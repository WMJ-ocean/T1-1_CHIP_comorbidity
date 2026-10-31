# -*- coding: utf-8 -*-
"""61d_diag_pqtl2.py —— 看 p 值分布 + 元数据，判断文件是否为「仅显著」子集"""
import zlib, requests, urllib3, json
urllib3.disable_warnings()
FTP = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST{acc}/GCST{acc}.tsv.gz'
UA = {'User-Agent': 'Mozilla/5.0 research'}
acc = '90274817'


def lines_from(resp):
    do = zlib.decompressobj(31); tail = ''
    for chunk in resp.iter_content(1 << 18):
        buf = do.decompress(chunk)
        while do.eof:
            nd = zlib.decompressobj(31); buf += nd.decompress(do.unused_data); do = nd
        tail += buf.decode('utf-8', 'replace')
        if '\n' in tail:
            parts = tail.split('\n'); tail = parts.pop()
            for L in parts:
                yield L


r = requests.get(FTP.format(acc=acc), headers=UA, stream=True, timeout=180, verify=False)
ps = []; sample = []
for i, L in enumerate(lines_from(r)):
    if i == 0:
        print('HDR:', L)
        continue
    p = L.split('\t')
    if len(p) < 11: continue
    if len(sample) < 3: sample.append(L)
    try: ps.append(float(p[7]))
    except Exception: pass
    if len(ps) >= 400000: break
r.close()
print('\n样本行:')
for s in sample: print('  ', s)
ps.sort()
print(f'\n共 {len(ps)} 行')
print('最小 p  :', ps[0], ' 第10小:', ps[9])
print('中位 p  :', ps[len(ps)//2])
print('p<0.05  :', sum(1 for x in ps if x < 0.05))
print('p<1e-5  :', sum(1 for x in ps if x < 1e-5))
print('p<5e-8  :', sum(1 for x in ps if x < 5e-8))

m = requests.get(f'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST{acc}',
                 headers=UA, timeout=60, verify=False).json()
print('\n=== 元数据 ===')
for k in ('accessionId', 'diseaseTrait', 'initialSampleSize', 'fullPvalueSet', 'snpCount', 'pubmedId'):
    v = m.get(k)
    if isinstance(v, dict):
        v = {kk: v.get(kk) for kk in ('trait', 'totalSampleSize') if kk in v}
    print(f'  {k}: {v}')
