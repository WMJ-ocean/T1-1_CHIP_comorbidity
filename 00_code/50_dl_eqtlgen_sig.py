# -*- coding: utf-8 -*-
"""50_dl_eqtlgen_sig.py
多线程 Range 断点续传下载 eQTLGen 显著 cis-eQTL（260 MB）
用法: python 50_dl_eqtlgen_sig.py
支持 Ctrl-C 中断后重跑续传（依赖 .part 文件大小）
"""
import os, sys, time, threading
import requests

URL = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
OUT = os.path.join(os.path.dirname(__file__), "..", "02_raw", "eqtlgen")
os.makedirs(OUT, exist_ok=True)
DEST = os.path.join(OUT, "eQTLGen_cis_eQTL_FDR005.txt.gz")

NTHREAD = 8
CHUNK = 8 * 1024 * 1024          # 每块 8 MB
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

def total_size():
    for _ in range(3):
        try:
            r = requests.head(URL, timeout=30, headers=UA, allow_redirects=True)
            if r.status_code == 200 and r.headers.get("Content-Length"):
                return int(r.headers["Content-Length"])
        except Exception:
            time.sleep(2)
    raise RuntimeError("无法获取远端文件大小")

def fetch_range(idx, start, end, done_bytes):
    """下载 [start,end] 到 .part{idx}; 已有时从断点继续"""
    part = DEST + f".part{idx}"
    have = 0
    if os.path.exists(part):
        have = os.path.getsize(part)
    if have >= end - start + 1:
        with glob_lock:
            done_bytes[0] += (end - start + 1)
        return True
    pos = start + have
    attempt = 0
    while attempt < 5:
        try:
            hd = dict(UA)
            hd["Range"] = f"bytes={pos}-{end}"
            with requests.get(URL, headers=hd, stream=True, timeout=60) as r:
                if r.status_code not in (200, 206):
                    raise RuntimeError(f"status {r.status_code}")
                with open(part, "ab") as f:
                    for c in r.iter_content(1024 * 512):
                        if not c:
                            break
                        f.write(c)
                        n = len(c)
                        pos += n
                        with glob_lock:
                            done_bytes[0] += n
            if pos >= end:
                return True
            have = pos - start
        except Exception as e:
            attempt += 1
            print(f"  [part{idx}] retry {attempt}: {type(e).__name__} {e}", flush=True)
            time.sleep(3 * attempt)
    return False

glob_lock = threading.Lock()

def main():
    total = total_size()
    print(f"远端文件 {total/1e6:.1f} MB, {NTHREAD} 线程, 断点续传已支持", flush=True)
    t0 = time.time()
    done = [0]
    # 已有分片计入
    for i in range(0, total, CHUNK):
        p = DEST + f".part{i//CHUNK}"
        if os.path.exists(p):
            done[0] += os.path.getsize(p)

    blocks = [(i, min(i + CHUNK - 1, total - 1)) for i in range(0, total, CHUNK)]
    print(f"分片数 {len(blocks)}", flush=True)
    ths = []
    for idx, (s, e) in enumerate(blocks):
        t = threading.Thread(target=fetch_range, args=(idx, s, e, done), daemon=True)
        t.start()
        ths.append(t)
        if len(ths) >= NTHREAD:
            ths = [x for x in ths if x.is_alive()]

    while any(t.is_alive() for t in ths):
        time.sleep(5)
        el = time.time() - t0
        sp = done[0] / 1e6 / max(el, 1)
        print(f"\r  {done[0]/1e6:8.1f}/{total/1e6:.1f} MB  {sp:6.2f} MB/s  "
              f"{done[0]/total*100:5.1f}%  eta {max((total-done[0])/1e6/max(sp,1e-3),0)/60:5.1f} min   ",
              end="", flush=True)
    print(flush=True)

    parts = sorted([(i // CHUNK, DEST + f".part{i//CHUNK}") for i, _ in blocks])
    exp = sum(os.path.getsize(p) if os.path.exists(p) else 0 for _, p in parts)
    if exp != total:
        print(f"!! 分片合计 {exp} != 远端 {total}，请重跑本脚本续传")
        sys.exit(2)
    print("合并分片...")
    with open(DEST, "wb") as out:
        for _, p in parts:
            with open(p, "rb") as f:
                while True:
                    b = f.read(32 * 1024 * 1024)
                    if not b:
                        break
                    out.write(b)
            os.remove(p)
    print(f"完成: {DEST}  {os.path.getsize(DEST)/1e6:.1f} MB  {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
