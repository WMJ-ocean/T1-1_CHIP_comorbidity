# -*- coding: utf-8 -*-
"""
57_dl_besd.py  —— 下载 eQTLGen SMR 格式化数据（cis-eQTL-SMR_20191212，约 1.37 GB）

为什么要这份：
  FDR<0.05 那份只有显著位点，做 coloc / SMR 会严重高估共定位概率
  （窗口里的零效应 SNP 全被丢掉了）。SMR 格式化版本按 probe 组织，
  保留 cis 窗口内的全部 SNP 及其 beta/se，才能做真正的 approx-BF coloc。

注意：沙箱禁止 os.remove()，合并后分片用 os.replace() 移进 _parts 目录。
"""
import os, io, gzip, time, shutil
import requests

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
BASE = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/SMR_formatted/"
FN = "cis-eQTL-SMR_20191212.tar.gz"

DEST = r"C:\Users\Hero\WorkBuddy\2026-10-14-23-35-01\outputs\pv_hotspot\T1-1_CHIP_comorbidity\02_raw\eqtlgen"
os.makedirs(DEST, exist_ok=True)
PARTS = os.path.join(DEST, "_parts_besd")
os.makedirs(PARTS, exist_ok=True)

URL = BASE + FN
r = requests.get(URL, headers=UA, stream=True, timeout=60)
if r.status_code != 200:
    print("HTTP", r.status_code, URL); raise SystemExit(1)
size = int(r.headers["Content-Length"])
supports_range = r.headers.get("Accept-Ranges") == "bytes"
r.close()
print(f"远端 {FN}")
print(f"  size = {size/1e6:.1f} MB   Range 支持 = {supports_range}")

outp = os.path.join(DEST, "eQTLGen_SMR_20191212.tar.gz")

if not supports_range:
    t0 = time.time()
    with requests.get(URL, headers=UA, stream=True, timeout=120) as rr, open(outp, "wb") as fo:
        done = 0
        for chunk in rr.iter_content(1 << 20):
            fo.write(chunk); done += len(chunk)
            el = time.time() - t0
            print(f"\r  {done/1e6:8.1f}/{size/1e6:.1f} MB  {done/1e6/el:6.2f} MB/s  "
                  f"{100*done/size:5.1f}%", end="", flush=True)
else:
    N = 8
    edges = [int(size * i / N) for i in range(N)] + [size - 1]
    paths = [os.path.join(PARTS, f"besd.p{i:02d}") for i in range(N)]

    import threading
    prog = [0] * N

    def worker(i):
        a, b = edges[i], edges[i + 1]
        if os.path.exists(paths[i]) and os.path.getsize(paths[i]) == b - a + 1:
            prog[i] = b - a + 1
            return
        h = dict(UA, Range=f"bytes={a}-{b}")
        rr = requests.get(URL, headers=h, stream=True, timeout=120)
        rr.raise_for_status()
        tmp = paths[i] + ".tmp"
        with open(tmp, "wb") as fo:
            for chunk in rr.iter_content(1 << 20):
                fo.write(chunk); prog[i] += len(chunk)
        os.replace(tmp, paths[i])

    ths = [threading.Thread(target=worker, args=(i,)) for i in range(N)]
    t0 = time.time()
    for t in ths: t.start()
    while any(t.is_alive() for t in ths):
        el = time.time() - t0
        tot = sum(prog)
        print(f"\r  {tot/1e6:8.1f}/{size/1e6:.1f} MB  "
              f"{tot/1e6/max(el,0.1):6.2f} MB/s  {100*tot/size:5.1f}%  "
              f"ETA {max((size-tot)/1e6/max(tot/1e6/max(el,0.1),1e-6),0)/60:5.1f} min",
              end="", flush=True)
        time.sleep(2)
    for t in ths: t.join()
    print()
    print("  合并 …")
    with open(outp, "wb") as fo:
        for pth in paths:
            with open(pth, "rb") as fi:
                shutil.copyfileobj(fi, fo, 1 << 22)
    ok = os.path.getsize(outp) == size
    print(f"  合并完成 {os.path.getsize(outp)/1e6:.1f} MB   大小对齐: {ok}")

print(f"\n成品: {outp}")
