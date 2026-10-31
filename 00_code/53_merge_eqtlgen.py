# -*- coding: utf-8 -*-
"""53_merge_eqtlgen.py — 合并分片 + gzip 完整性校验（不删除任何分片）"""
import os, gzip, time, hashlib

D = os.path.dirname(os.path.abspath(__file__)) + "/../02_raw/eqtlgen"
DEST = os.path.join(D, "eQTLGen_cis_eQTL_FDR005.txt.gz")
REMOTE = 322_751_000   # 实际以后面获取为准；此处用于提示

parts = sorted([f for f in os.listdir(D) if ".part" in f],
               key=lambda x: int(x.split(".part")[-1]))
print(f"分片 {len(parts)} 个")
if os.path.exists(DEST):
    os.replace(DEST, DEST + ".incomplete")

t0 = time.time()
n_in = 0
with open(DEST, "wb") as out:
    for p in parts:
        pp = os.path.join(D, p)
        with open(pp, "rb") as f:
            while True:
                b = f.read(32 * 1024 * 1024)
                if not b:
                    break
                out.write(b)
                n_in += len(b)
print(f"合并完成 {n_in/1e6:.1f} MB  {time.time()-t0:.0f}s")

# gzip CRC 完整性 + 原始字节数统计
print("开始流式校验（读取整个 gzip 流）…")
t1 = time.time()
raw = 0
nline = 0
h = hashlib.md5()
with gzip.open(DEST, "rb") as g:
    while True:
        b = g.read(64 * 1024 * 1024)
        if not b:
            break
        raw += len(b)
        nline += b.count(b"\n")
        h.update(b[:10_000_000])  # 只算开头做指纹，省时间
print(f"OK  解压 {raw/1e6:.1f} MB / {nline:,} 行  {time.time()-t1:.0f}s")

if os.path.exists(DEST + ".incomplete"):
    os.remove(DEST + ".incomplete")

print(f"\n成品: {DEST}  {os.path.getsize(DEST)/1e6:.1f} MB")
print(f"远端应约 {REMOTE/1e6:.1f} MB")
print("差异超过 1 MB 说明分片不齐，需重跑 50 号脚本续传")
