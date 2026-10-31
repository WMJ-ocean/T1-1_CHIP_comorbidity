# -*- coding: utf-8 -*-
"""52b_pqtl_dir.py — 正确解析 GWAS Catalog FTP 列表，定位 pQTL 文件"""
import requests, re
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})

DIR = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/"
r = S.get(DIR, timeout=60)
hrefs = re.findall(r'href="([^"]+)"', r.text)
print(f"href 总数 {len(hrefs)}")
gz = [h for h in hrefs if h.endswith(".gz") or h.endswith(".tsv")]
print(f"文件类 {len(gz)}")
print("前 8 个示例：")
for h in gz[:8]:
    print("   ", h)

target = [h for h in gz if "90274758" in h or "9027475" in h or "9027484" in h]
print(f"\n含 90274758 的文件：{target}")

if target:
    for t in target[:3]:
        rr = S.head(DIR + t, timeout=30)
        print(f"  HEAD {t} -> {rr.headers.get('Content-Length')} bytes, ranges={rr.headers.get('Accept-Ranges')}")
else:
    # 试 harmonised 子目录
    HD = DIR + "harmonised/"
    rr = S.get(HD, timeout=60)
    print(f"\nharmonised/ -> {rr.status_code}")
    if rr.status_code == 200:
        hh = [h for h in re.findall(r'href="([^"]+)"', rr.text) if h.endswith(".gz")]
        ht = [h for h in hh if "90274758" in h]
        print(f"  harmonised 下 {len(hh)} 个文件；命中 {len(ht)}")
        if ht:
            h = S.head(HD + ht[0], timeout=30)
            print(f"  HEAD {ht[0]} -> {h.headers.get('Content-Length')} bytes ranges={h.headers.get('Accept-Ranges')}")
        else:
            print("  样例:", hh[:5])

print("\nDONE")
