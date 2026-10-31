# -*- coding: utf-8 -*-
"""52_probe_pqtl.py — 查 GWAS Catalog FTP 上的炎症蛋白 pQTL（Zheng 2020, 91 proteins）"""
import requests, re
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})

DIR = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/"
print("### 列出目录并找 GCST90274xxx (Zheng 2020 pQTL 区间)")
r = S.get(DIR, timeout=60)
print("status", r.status_code, "len", len(r.text))
rows = re.findall(r'<a href="([^"]+)">([^<]+)</a>.*?(\d{4}-\d{2}-\d{2} \d{2}:\d{2})', r.text[:200000])
hits = [x for x in rows if "GCST9027" in x[1] and x[1].endswith((".gz", ".tsv", ".txt"))]
print(f"命中 {len(hits)} 个文件")
for fn, _, date in hits[:12]:
    print("   ", fn)
if hits:
    # 查第一个的大小
    fn = hits[0][0]
    rr = S.head(DIR + fn, timeout=30)
    print(f"样例 HEAD {fn} -> {rr.headers.get('Content-Length')} bytes")

print("\n### GWAS Catalog v2 REST API")
for u in [
    "https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/GCST90274758",
    "https://gwas-catalog-api.ebi.ac.uk/gwas/summary-statistics/api/studies/GCST90274758",
]:
    try:
        rr = S.get(u, timeout=30)
        print(f"\n{u}\n  -> {rr.status_code}")
        print("  ", rr.text[:400])
    except Exception as e:
        print(f"\n{u}\n  FAIL {type(e).__name__}: {e}")

print("\nDONE")
