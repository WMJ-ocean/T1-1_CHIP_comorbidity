# 50_eqtlgen_probe.py — 远程窥探 eQTLGen 文件：真实大小 + 表头 + 列名
import zlib, requests, io

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
B = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"

FILES = [
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
    "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
    "cis-eQTL_significant_20181017.txt.gz",
    "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz",
]

for f in FILES:
    url = B + f
    h = requests.get(url, headers=UA, stream=True, timeout=30)
    total = h.headers.get("Content-Length")
    print("=" * 78)
    print(f, "\n  size:", f"{int(total)/1e6:.1f} MB" if total else "?")
    # 取前 1 MB 并解压（raw deflate, wbits=-15）
    buf = b""
    for chunk in h.iter_content(32768):
        buf += chunk
        if len(buf) >= 1 << 20:
            break
    h.close()
    d = zlib.decompressobj(-15)
    txt = d.decompress(buf).decode("utf-8", "replace")
    lines = txt.split("\n")
    print("  columns:", lines[0].strip()[:300])
    if len(lines) > 1:
        print("  row1   :", lines[1].strip()[:300])
        if len(lines) > 2:
            print("  row2   :", lines[2].strip()[:300])
print("=" * 78)
print("DONE")
