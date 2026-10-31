# -*- coding: utf-8 -*-
"""测试数据下载通道是否可用（T1-1 第一步）"""
import requests, sys, time
requests.packages.urllib3.disable_warnings()

UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'}

def probe(name, url, method='HEAD', timeout=30):
    try:
        if method == 'HEAD':
            r = requests.head(url, headers=UA, timeout=timeout, allow_redirects=True, verify=False)
        else:
            r = requests.get(url, headers=UA, timeout=timeout, stream=True, allow_redirects=True, verify=False)
        cl = r.headers.get('Content-Length', '?')
        ct = r.headers.get('Content-Type', '?')
        print('%-28s %s  code=%s  len=%s  type=%s' % (name, 'OK ' if r.status_code < 400 else 'BAD', r.status_code, cl, ct))
        if method != 'HEAD' and r.status_code < 400:
            n = 0
            t0 = time.time()
            try:
                for chunk in r.iter_content(65536):
                    n += len(chunk)
                    if n > 300000:
                        break
            except Exception as e:
                print('    read err', e)
            dt = time.time() - t0
            sp = n / dt / 1024 if dt > 0 else 0
            print('    speed: %.1f KB/s (read %d bytes in %.1fs)' % (sp, n, dt))
        return r.status_code
    except Exception as e:
        print('%-28s ERR  %s' % (name, type(e).__name__), str(e)[:90])
        return None

print('=== 1. GWAS Catalog REST API ===')
probe('catalog-study-GCST90000032', 'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90000032')
probe('catalog-sumstats-list', 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative')

print()
print('=== 2. GWAS Catalog FTP (HTTPS mirror) ===')
# Kessler 2022 CHIP: GCST90208187 ~ GCST90208194 range; 先探测目录
for acc in ['GCST90208187', 'GCST90208188', 'GCST90208189', 'GCST90208190',
            'GCST90208191', 'GCST90208192', 'GCST90208193', 'GCST90208194']:
    u = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/%s/%s/' % (acc[0:11], acc)
    probe(acc + '-dir', u, 'GET')

print()
print('=== 3. FinnGen R12 ===')
probe('finngen-R12-MDS', 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_C3_MYELODYSP_SYNDR_EXALLC.gz', 'GET')
probe('finngen-R12-AA', 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_D3_OTHERAPLASTICANAEMIA.gz', 'GET')

print()
print('=== 4. IEU OpenGWAS API (需 JWT) ===')
probe('ieu-gwas-info', 'https://gwas-api.mrcieu.ac.uk/gwasinfo/ieu-b-4965')

print()
print('=== 5. LDSC reference (ld score) ===')
probe('ldsc-eur-chr1', 'https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2', 'HEAD', 40)

print()
print('DONE')
