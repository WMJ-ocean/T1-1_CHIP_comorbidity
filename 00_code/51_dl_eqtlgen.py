# 51_dl_eqtlgen.py — eQTLGen 多线程 Range 断点续传下载器
import os, sys, time, threading, queue
import requests

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
B = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"

FILES = [
    "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
    "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz",
]

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "02_raw", "eqtlgen")
os.makedirs(OUT, exist_ok=True)

NTHREAD = 8
CHUNK = 16 * 1024 * 1024   # 16 MB per task

lock = threading.Lock()
state = {}


def head_size(url):
    for _ in range(3):
        try:
            r = requests.head(url, headers=UA, timeout=30, allow_redirects=True)
            if r.status_code == 200:
                return int(r.headers["Content-Length"]), bool(r.headers.get("Accept-Ranges") == "bytes")
        except Exception:
            time.sleep(2)
    return None, False


def worker(url, dst):
    while True:
        try:
            item = q.get_nowait()
        except queue.Empty:
            return
        start, end, idx = item
        part = f"{dst}.part{idx:04d}"
        # 已经下完的跳过（断点续传）
        if os.path.exists(part):
            have = os.path.getsize(part)
            if have >= end - start + 1:
                with lock:
                    state["done"] += end - start + 1
                continue
            start = start + have
            mode = "ab"
        else:
            mode = "wb"
        for attempt in range(6):
            try:
                h = dict(UA)
                h["Range"] = f"bytes={start}-{end}"
                r = requests.get(url, headers=h, stream=True, timeout=60)
                if r.status_code not in (200, 206):
                    raise RuntimeError(f"HTTP {r.status_code}")
                pos = start
                with open(part, mode) as f:
                    for c in r.iter_content(256 * 1024):
                        f.write(c)
                        pos += len(c)
                        with lock:
                            state["done"] += len(c)
                r.close()
                if pos >= end:
                    break
                start = pos
                mode = "ab"
            except Exception as e:
                print(f"  [warn] part{idx} {type(e).__name__}: {e} retry {attempt+1}", flush=True)
                time.sleep(3 * (attempt + 1))
        else:
            print(f"  [FAIL] part{idx} 放弃", flush=True)


for fn in FILES:
    url = B + fn
    dst = os.path.join(OUT, fn)
    total, ranged = head_size(url)
    print(f"\n### {fn}\n  size={total/1e6:.1f} MB  accept-ranges={ranged}", flush=True)
    if total is None:
        print("  无法获取大小，跳过"); continue
    if os.path.exists(dst) and os.path.getsize(dst) == total:
        print("  已存在且完整，跳过"); continue

    q = queue.Queue()
    if ranged:
        i = 0
        for s in range(0, total, CHUNK):
            q.put((s, min(s + CHUNK - 1, total - 1), i)); i += 1
        state = {"done": 0, "total": total}
        t0 = time.time()
        ths = [threading.Thread(target=worker, args=(url, dst), daemon=True) for _ in range(NTHREAD)]
        [t.start() for t in ths]
        while any(t.is_alive() for t in ths):
            time.sleep(5)
            el = time.time() - t0
            d = state["done"]
            rate = d / el / 1e6
            eta = (total - d) / max(rate * 1e6, 1)
            print(f"\r  {d/1e6:7.1f}/{total/1e6:.1f} MB  {rate:5.2f} MB/s  "
                  f"{d/total*100:5.1f}%  ETA {eta/60:4.1f} min   ", end="", flush=True)
        [t.join() for t in ths]
        print(flush=True)
        # 合并分片
        with open(dst, "wb") as out:
            for idx in range(i):
                p = f"{dst}.part{idx:04d}"
                if not os.path.exists(p):
                    print(f"  [FAIL] 缺少分片 {p}"); break
                with open(p, "rb") as f:
                    while True:
                        b = f.read(8 * 1024 * 1024)
                        if not b: break
                        out.write(b)
                os.remove(p)
        print(f"  合并完成: {os.path.getsize(dst)/1e6:.1f} MB", flush=True)
    else:
        # 单线程兜底
        t0 = time.time()
        with requests.get(url, headers=UA, stream=True, timeout=60) as r, open(dst, "wb") as f:
            for c in r.iter_content(1024 * 1024):
                f.write(c)
        print(f"  单线程完成 {os.path.getsize(dst)/1e6:.1f} MB in {(time.time()-t0)/60:.1f} min", flush=True)

print("\nALL DONE")
