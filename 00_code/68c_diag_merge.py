# -*- coding: utf-8 -*-
"""68c_diag_merge.py —— 定位 3.88 GB 合并文件的解压失败字节偏移，并检查分片接缝"""
import os, zlib

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, '02_raw', 'eqtlgen', '_parts_full')
DST = os.path.join(BASE, '02_raw', 'eqtlgen', 'eQTLGen_cis_eQTL_full.txt.gz')

parts = sorted([f for f in os.listdir(SRC) if f.startswith('full.p')],
               key=lambda s: int(s.split('.p')[-1]))
sizes = [os.path.getsize(os.path.join(SRC, p)) for p in parts]
print('分片大小:', sizes, flush=True)
cum = [0]
for s in sizes:
    cum.append(cum[-1] + s)
print('累计边界:', cum, flush=True)

# 接缝检查：后一片首字节 vs 前一片末字节
print('\n=== 接缝字节 ===', flush=True)
for i in range(len(parts) - 1):
    with open(os.path.join(SRC, parts[i]), 'rb') as f:
        f.seek(-1, os.SEEK_END); last = f.read(1)
    with open(os.path.join(SRC, parts[i + 1]), 'rb') as f:
        first = f.read(1)
    print(f'  p{i:03d}|p{i+1:03d}  last={last.hex()} first={first.hex()} '
          f'{"**相同(疑似重复字节)**" if last == first else ""}', flush=True)

# 流式解压，记录失败位置
print('\n=== 解压定位 ===', flush=True)
do = zlib.decompressobj(31)
off = 0
nline = 0
err_off = None
with open(DST, 'rb') as fi:
    while True:
        chunk = fi.read(1 << 22)
        if not chunk:
            break
        try:
            buf = do.decompress(chunk)
            while do.eof:
                nd = zlib.decompressobj(31)
                buf += nd.decompress(do.unused_data); do = nd
        except zlib.error as e:
            err_off = off
            print(f'  [ERR] {e} @ 约 {off/1e6:.1f} MB  行数={nline:,}', flush=True)
            break
        off += len(chunk)
        nline += buf.count(b'\n')
print(f'  读到 {off/1e6:.1f} MB / {os.path.getsize(DST)/1e6:.1f} MB, 行={nline:,}', flush=True)
if err_off is not None:
    print(f'  失败偏移 {err_off:,} B', flush=True)
    for i, c in enumerate(cum):
        if abs(c - err_off) < 1_000_000:
            print(f'    距边界 cum[{i}]={c:,} 差 {err_off-c:+,} B', flush=True)
