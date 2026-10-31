# -*- coding: utf-8 -*-
"""49b_eqtlgen_sizes.py  — 查各文件实际体积 + SMR_results 目录"""
import requests
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
BASE = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"

files = [
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
    "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
    "cis-eQTL_significant_20181017.txt.gz",
    "cis-eQTLs_full_20180905.txt.gz",
    "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz",
    "SMR_formatted/README_cis_full_SMR",
    "SMR_formatted/cis-eQTL-SMR_20191212.tar.gz",
]
print("########## 文件体积 ##########")
for f in files:
    try:
        r = S.head(BASE + f, timeout=25, allow_redirects=True)
        cl = r.headers.get("Content-Length")
        ar = r.headers.get("Accept-Ranges")
        tag = f"{int(cl)/1e6:.0f} MB" if cl else "?"
        print(f"  [{r.status_code}] {tag:>9}  ranges={ar}  {f}")
    except Exception as e:
        print(f"  [ERR] {f}  {type(e).__name__}: {e}")

print("\n########## SMR_results 目录 ##########")
try:
    r = S.get(BASE + "SMR_results/", timeout=25)
    print("status", r.status_code)
    import re
    links = re.findall(r'href="([^"]+)"', r.text)
    for l in links[:40]:
        print("   ", l)
except Exception as e:
    print("FAIL", e)

print("\n########## README 内容 ##########")
try:
    r = S.get(BASE + "SMR_formatted/README_cis_full_SMR", timeout=25)
    print(r.text[:1500])
except Exception as e:
    print("FAIL", e)

print("\nDONE")
