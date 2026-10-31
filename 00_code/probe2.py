# -*- coding: utf-8 -*-
import requests
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

urls = [
 ('ebi-home',        'https://www.ebi.ac.uk/'),
 ('gwas-catalog-home','https://www.ebi.ac.uk/gwas/'),
 ('gwas-rest-search', 'https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByEfoTrait?efoTrait=myeloproliferative%20neoplasm'),
 ('gwas-rest-acc',    'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90208187'),
 ('ftp-ebi-root',     'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'),
 ('ftp-ebi-range1',   'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/'),
 ('ftp-ebi-range2',   'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90208001-GCST90209000/'),
 ('gcs-fg-r12-root',  'https://storage.googleapis.com/storage/v1/b/finngen-public-data-r12/o?prefix=summary_stats/finngen_R12_C3_MYELODYS&maxResults=20'),
 ('gcs-fg-r11-root',  'https://storage.googleapis.com/storage/v1/b/finngen-public-data-r11/o?prefix=summary_stats/finngen_R11_C3_MYELODYS&maxResults=20'),
 ('ieu-home',         'https://gwas-api.mrcieu.ac.uk/'),
 ('ieu-docs',         'https://gwas-api.mrcieu.ac.uk/docs'),
 ('google',           'https://www.google.com/'),
]
for name, u in urls:
    try:
        r = requests.get(u, headers=UA, timeout=25, verify=False, stream=True)
        body = b''
        try:
            for c in r.iter_content(4096):
                body += c
                if len(body) > 1500:
                    break
        except Exception:
            pass
        txt = body.decode('utf-8', 'ignore')
        print('%-20s code=%s len=%s' % (name, r.status_code, r.headers.get('Content-Length', '?')))
        print('    ', txt[:400].replace('\n', ' '))
    except Exception as e:
        print('%-20s ERR %s %s' % (name, type(e).__name__, str(e)[:80]))
