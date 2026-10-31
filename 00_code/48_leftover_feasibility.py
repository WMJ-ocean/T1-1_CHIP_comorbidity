# -*- coding: utf-8 -*-
"""48_leftover_feasibility.py
实测三项遗留的可行性：eQTLGen 下载速度 / GTEx API / IEU OpenGWAS eQTL
"""
import sys, time, json
import requests

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-script"}
S = requests.Session()
S.headers.update(UA)

def probe(name, url, method="HEAD", **kw):
    print(f"\n=== {name} ===")
    print(f"URL: {url}")
    t0 = time.time()
    try:
        r = S.request(method, url, timeout=25, allow_redirects=True, **kw)
        dt = time.time() - t0
        print(f"status={r.status_code}  {dt:.1f}s")
        cl = r.headers.get("Content-Length")
        ct = r.headers.get("Content-Type")
        ar = r.headers.get("Accept-Ranges")
        print(f"content-length={cl}  ({int(cl)/1e6:.1f} MB)" if cl else f"content-length=None")
        print(f"content-type={ct}  accept-ranges={ar}")
        return r
    except Exception as e:
        print(f"FAIL {type(e).__name__}: {e}")
        return None

# ---------- 1. eQTLGen SMR 格式（比全量 txt 小） ----------
probe("eQTLGen SMR-formatted",
      "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz")

# ---------- 2. 实测下载速度（Range 取 8 MB，允许 Bacon 再掐） ----------
url = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz"
print("\n=== 速度实测：Range 取前 8 MB ===")
t0 = time.time()
n = 0
rate = 0
try:
    with S.get(url, stream=True, timeout=30,
               headers={"Range": "bytes=0-8388607"}) as r:
        print(f"status={r.status_code} (206=支持断点)")
        for chunk in r.iter_content(1024 * 256):
            n += len(chunk)
            if time.time() - t0 > 20 or n >= 8 * 1024 * 1024:
                break
    dt = time.time() - t0
    rate = n / 1e6 / dt
    print(f"下载 {n/1e6:.2f} MB / {dt:.1f}s = {rate:.2f} MB/s")
    if rate > 0.01:
        print(f">>> 按此速度，1 GB 需要 {1000/rate/60:.0f} 分钟")
except Exception as e:
    print(f"FAIL {type(e).__name__}: {e}")

# ---------- 3. GTEx Portal API v2：按 rsID 查单组织 eQTL ----------
probe("GTEx API single-tissue eQTL for variant",
      "https://gtexportal.org/api/v2/association/singleTissueEqtl",
      method="GET",
      params={"variantId": "chr17_7673386_G_A_b38", "gencodeId": "ENSG00000141510.17"})

# ---------- 4. IEU OpenGWAS (ieugwasr 后端) ----------
probe("IEU OpenGWAS API root", "https://gwas-api.mrcieu.ac.uk/status")

print("\n=== DONE ===")
