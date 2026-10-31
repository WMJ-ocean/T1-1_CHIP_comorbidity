# -*- coding: utf-8 -*-
"""51_probe_eastasian.py — 东亚 GWAS 数据源可达性实测"""
import requests, json
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})

def go(tag, url, **kw):
    try:
        r = S.get(url, timeout=30, **kw)
        cl = r.headers.get("Content-Length")
        print(f"\n### {tag}\n  URL {url}\n  -> {r.status_code}  len={cl}")
        txt = r.text[:600].replace("\n", " ")
        print(f"  body: {txt}")
        return r
    except Exception as e:
        print(f"\n### {tag}\n  URL {url}\n  FAIL {type(e).__name__}: {e}")
        return None

# 1. BBJ PheWeb variant API（东亚最大 biobank）
go("BBJ PheWeb variant API", "https://pheweb.jp/api/variant/rs7705526")
# 2. BBJ PheWeb manifest
go("BBJ PheWeb phenotypes manifest", "https://pheweb.jp/api/phenotypes/0", timeout=60)
# 3. Pan-UKBB (含 EAS) 注册表
go("Pan-UKBB manifest", "https://pan.ukbb.broadinstitute.org/api/v1/statistics/manifest")
# 4. GWAS Catalog: Zheng 2020 炎症蛋白 pQTL GCST90274758
go("GWAS Catalog study GCST90274758 (炎症蛋白 pQTL)",
   "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90274758")
# 5. GWAS Catalog download 目录
go("GWAS Catalog summary stats dir",
   "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/")

print("\nDONE")
