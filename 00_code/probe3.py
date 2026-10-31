# -*- coding: utf-8 -*-
import requests, re, time
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

def get(u, n=200000, t=40):
    r = requests.get(u, headers=UA, timeout=t, verify=False, stream=True)
    b = b''
    try:
        for c in r.iter_content(8192):
            b += c
            if len(b) > n:
                break
    except Exception:
        pass
    return r.status_code, b.decode('utf-8', 'ignore')

print('=== A. ftp.ebi summary_statistics 目录（找 GCST9020xxxx 范围）===')
st, txt = get('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/', 400000)
dirs = sorted(set(re.findall(r'href="(GCST\d+-GCST\d+)/"', txt)))
print('n ranges =', len(dirs))
print('last 12:', dirs[-12:])
hits = [d for d in dirs if d.startswith('GCST9020') or d.startswith('GCST9021')]
print('9020/9021 ranges:', hits)

print()
print('=== B. FinnGen R12 (release 子目录) ===')
for pheno in ['C3_MYELODYSP_SYNDR_EXALLC', 'D3_OTHERAPLASTICANAEMIA', 'D3_ITP', 'C3_MYELOPROLIFERATIVE_NEOPLASM']:
    u = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz' % pheno
    try:
        r = requests.head(u, headers=UA, timeout=30, verify=False)
        print('%-34s code=%s len=%s' % (pheno, r.status_code, r.headers.get('Content-Length', '?')))
    except Exception as e:
        print('%-34s ERR %s' % (pheno, str(e)[:60]))

print()
print('=== C. FinnGen R12 manifest (前几行) ===')
st, txt = get('https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv', 300000, 60)
print('code=', st, 'chars=', len(txt))
for line in txt.split('\n')[:3]:
    print('   ', line[:300])
if 'MYELODYS' in txt or 'APLASTIC' in txt:
    for line in txt.split('\n'):
        if 'MYELODYS' in line or 'APLASTIC' in line or 'MYELOPROLIF' in line:
            print('   HIT:', line[:250])
else:
    print('   (manifest 未包含这些表型，可能内容不同)')

print()
print('=== D. GWAS Catalog V2 API ===')
for u in ['https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90000032',
          'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/search/findByEfoTrait?efoTrait=clonal%20hematopoiesis',
          'https://www.ebi.ac.uk/gwas/api/search?q=clonal+hematopoiesis']:
    st, txt = get(u, 3000, 30)
    print('%-95s code=%s' % (u[-70:], st))
    print('   ', txt[:220].replace('\n', ' '))

print()
print('DONE')
