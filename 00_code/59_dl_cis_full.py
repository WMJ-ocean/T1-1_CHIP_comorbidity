# -*- coding: utf-8 -*-
"""
59_dl_cis_full.py  —— 下载 eQTLGen 全量 cis-eQTL（3480 万+ SNP-基因对）
这是做真正的 SMR / coloc 所必需的：FDR<0.05 精简版只保留显著结果，
窗口内的零效应 SNP 被丢掉会让 coloc 的 PP.H4 严重虚高，必须拿全量。
"""
import os, threading, time, shutil, requests

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
BASE = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
FN = "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
DEST = r"C:\Users\Hero\WorkBuddy\2026-10-14-23-35-01\outputs\pv_hotspot\T1-1_CHIP_comorbidity\02_raw\eqtlgen"
os.makedirs(DEST, exist_ok=True)
PARTS = os.path.join(DEST, "_parts_full")
os.makedirs(PARTS, exist_ok=True)
outp = os.path.join(DEST, "eQTLGen_cis_eQTL_full.txt.gz")

URL = BASE + FN
r = requests.head(URL, headers=UA, timeout=60)
size = int(r.headers["Content-Length"])
print(f"远端 {FN}  size={size/1e6:.0f} MB  range={r.headers.get('Accept-Ranges')}")

N = 16                      # 提高并发到 16 线程
edges = [int(size * i / N) for i in range(N)] + [size - 1]
paths = [os.path.join(PARTS, f"full.p{i:03d}") for i in range(N)]
prog = [0] * N

def worker(i):
    a, b = edges[i], edges[i + 1]
    if os.path.exists(paths[i]) and os.path.getsize(paths[i]) == b - a + 1:
        prog[i] = b - a + 1
        return
    for attempt in range(5):
        try:
            h = dict(UA, Range=f"bytes={a}-{b}")
            rr = requests.get(URL, headers=h, stream=True, timeout=180)
            rr.raise_for_status()
            tmp = paths[i] + f".tmp{attempt}"
            got = 0
            with open(tmp, "wb") as fo:
                for chunk in rr.iter_content(1 << 20):
                    fo.write(chunk); got += len(chunk); prog[i] = min(got, b - a + 1)
            os.replace(tmp, paths[i])
            return
        except Exception as e:
            print(f"  [retry {i} attempt {attempt}] {type(e).__name__}: {str(e)[:60]}")
            prog[i] = 0
            time.sleep(3)

ths = [threading.Thread(target=worker, args=(i,)) for i in range(N)]
t0 = time.time()
for t in ths: t.start()
while any(t.is_alive() for t in ths):
    el = time.time() - t0
    tot = sum(prog)
    sp = tot / 1e6 / max(el, 0.1)
    eta = (size - tot) / 1e6 / max(sp, 1e-6) / 60
    print(f"\r  {tot/1e6:9.1f}/{size/1e6:.0f} MB  {sp:6.2f} MB/s  "
          f"{100*tot/size:5.1f}%  ETA {eta:5.1f} min   ", end="", flush=True)
    time.sleep(5)
for t in ths: t.join()

print()
print("合并 …")
with open(outp, "wb") as fo:
    for pth in paths:
        with open(pth, "rb") as fi:
            shutil.copyfileobj(fi, fo, 1 << 22)
print(f"  完成 {os.path.getsize(outp)/1e6:.1f} MB  (期望 {size/1e6:.1f})")
print(f"-> {outp}")
