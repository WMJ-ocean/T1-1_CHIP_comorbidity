# 55_dl_maf.py — 只下载 eQTLGen 等位频率文件（240 MB），用于 Z->beta/se 换算
# 设计要点：全程不使用 os.remove（沙箱禁止删除），分片合并后用 os.replace 移入 _parts 子目录
import os, sys, time, threading, queue, gzip, shutil
import requests

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
URL = ("https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
       "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "02_raw", "eqtlgen")
os.makedirs(OUT, exist_ok=True)
DST = os.path.join(OUT, "eQTLGen_SNP_AF_AlleleB.txt.gz")
TRASH = os.path.join(OUT, "_parts")
os.makedirs(TRASH, exist_ok=True)

NTHREAD = 8
CHUNK = 16 * 1024 * 1024

lock = threading.Lock()
STATE = {}


def head_size():
    for _ in range(3):
        try:
            r = requests.head(URL, headers=UA, timeout=30, allow_redirects=True)
            if r.status_code == 200:
                return int(r.headers["Content-Length"])
        except Exception:
            time.sleep(2)
    return None


def worker(q):
    while True:
        try:
            start, end, idx = q.get_nowait()
        except queue.Empty:
            return
        part = os.path.join(TRASH, f"maf.part{idx:04d}")
        mode = "wb"
        if os.path.exists(part):
            have = os.path.getsize(part)
            if have >= end - start + 1:
                with lock:
                    STATE["done"] += end - start + 1
                continue
            start += have
            mode = "ab"
        for attempt in range(6):
            try:
                h = dict(UA)
                h["Range"] = f"bytes={start}-{end}"
                r = requests.get(URL, headers=h, stream=True, timeout=60)
                if r.status_code not in (200, 206):
                    raise RuntimeError(f"HTTP {r.status_code}")
                pos = start
                with open(part, mode) as f:
                    for c in r.iter_content(256 * 1024):
                        f.write(c)
                        pos += len(c)
                        with lock:
                            STATE["done"] += len(c)
                r.close()
                if pos >= end:
                    break
                start, mode = pos, "ab"
            except Exception as e:
                print(f"  [warn] part{idx} {type(e).__name__}: {str(e)[:80]} retry {attempt+1}", flush=True)
                time.sleep(3 * (attempt + 1))
        else:
            print(f"  [FAIL] part{idx} 放弃", flush=True)


total = head_size()
print(f"远端 {URL.split('/')[-1]}\n  size={total/1e6:.1f} MB", flush=True)

if os.path.exists(DST) and os.path.getsize(DST) == total:
    print("  已存在且完整，无需下载")
else:
    if os.path.exists(DST):
        os.replace(DST, os.path.join(TRASH, "maf.superseded"))
    q = queue.Queue()
    i = 0
    for s in range(0, total, CHUNK):
        q.put((s, min(s + CHUNK - 1, total - 1), i))
        i += 1
    STATE = {"done": 0}
    t0 = time.time()
    ths = [threading.Thread(target=worker, args=(q,), daemon=True) for _ in range(NTHREAD)]
    [t.start() for t in ths]
    while any(t.is_alive() for t in ths):
        time.sleep(5)
        el = time.time() - t0
        d = STATE["done"]
        rate = d / el / 1e6
        print(f"\r  {d/1e6:7.1f}/{total/1e6:.1f} MB  {rate:5.2f} MB/s  {d/total*100:5.1f}%  ETA {(total-d)/max(rate*1e6,1)/60:4.1f} min   ",
              end="", flush=True)
    [t.join() for t in ths]
    print(flush=True)

    with open(DST, "wb") as out:
        for idx in range(i):
            p = os.path.join(TRASH, f"maf.part{idx:04d}")
            if not os.path.exists(p):
                print(f"  [FAIL] 缺分片 {p}")
                sys.exit(2)
            with open(p, "rb") as f:
                shutil.copyfileobj(f, out, 32 * 1024 * 1024)
            # 不删除，改移入归档目录（沙箱禁止 os.remove）
            os.replace(p, os.path.join(TRASH, f"done.maf.part{idx:04d}"))
    print(f"  合并完成 {os.path.getsize(DST)/1e6:.1f} MB", flush=True)

print("  gzip CRC 校验…", flush=True)
ok, n = 0, 0
with gzip.open(DST, "rb") as g:
    while True:
        b = g.read(64 * 1024 * 1024)
        if not b:
            break
        ok += len(b)
        n += b.count(b"\n")
print(f"  OK  解压 {ok/1e6:.1f} MB / {n:,} 行")
print(f"\n成品: {DST}")
