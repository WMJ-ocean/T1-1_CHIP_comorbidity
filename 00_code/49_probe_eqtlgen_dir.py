# -*- coding: utf-8 -*-
"""49_probe_eqtlgen_dir.py
查 molgenis26 上有没有按染色体分卷的 eQTLGen 文件（比全量 1.37 GB 省得多）
"""
import requests, re, time
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})

BASE = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"

cands = [
    "",                                                    # 目录索引
    "RAW/",
    "SMR_formatted/",
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
]

for c in cands:
    u = BASE + c
    try:
        r = S.get(u, timeout=25)
        print(f"\n===== {c or '(root)'}  -> {r.status_code} =====")
        if r.status_code == 200:
            ht = r.text
            if "<html" in ht.lower()[:1000] or "<!doctype" in ht.lower()[:1000]:
                links = re.findall(r'href="([^"]+)"', ht)
                hits = [l for l in links if any(k in l.lower() for k in (".gz", ".zip", "chr", "smr", ".txt", ".besd", "/"))]
                print("  links(%d):" % len(hits))
                for l in hits[:60]:
                    print("   ", l)
            else:
                print("  (非 HTML, len=%d) 前 200 字节:" % len(ht))
                print("  ", ht[:200].replace("\n", " | ")[:300])
        elif r.status_code == 206:
            print("  206 partial")
    except Exception as e:
        print(f"\n===== {c} FAIL {type(e).__name__}: {e}")

# 尝试猜测 per-chromosome 文件名（常见命名）
print("\n\n########## 试猜 per-chr 文件 ##########")
guesses = []
for ch in [str(i) for i in range(1, 23)]:
    guesses += [
        f"{BASE}RAW/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded-chr{ch}.txt.gz",
        f"{BASE}SMR_formatted/cis-eQTL-SMR-chr{ch}.tar.gz",
        f"{BASE}cis-eQTL-chr{ch}.txt.gz",
    ]
for g in guesses[:14]:
    try:
        r = S.head(g, timeout=15, allow_redirects=True)
        cl = r.headers.get("Content-Length")
        tag = f"{int(cl)/1e6:.1f} MB" if cl else "-"
        print(f"  [{r.status_code}] {tag:>10}  {g.split('/')[-1]}")
    except Exception as e:
        print(f"  [ERR] {g.split('/')[-1]}  {type(e).__name__}")

print("\nDONE")
