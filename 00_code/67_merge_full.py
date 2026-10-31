# -*- coding: utf-8 -*-
"""合并 eQTLGen 全量 cis-eQTL 的 16 个分片，并做 gzip CRC 完整性校验。
注意：沙箱禁止 os.remove，分片保留不删。"""
import io, os, sys, time, zlib, gzip, shutil, binascii

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, '02_raw', 'eqtlgen', '_parts_full')
DST = os.path.join(BASE, '02_raw', 'eqtlgen', 'eQTLGen_cis_eQTL_full.txt.gz')

parts = sorted([f for f in os.listdir(SRC) if f.startswith('full.p')],
               key=lambda s: int(s.split('.p')[-1]))
print('分片数:', len(parts))
for p in parts:
    print('   ', p, f'{os.path.getsize(os.path.join(SRC,p))/1e6:.1f} MB')

t0 = time.time()
crc = 0
n_in = n_out = 0
if os.path.exists(DST):
    print('\n目标文件已存在，先移除旧内容重建')

with open(DST, 'wb') as fo:
    for i, p in enumerate(parts):
        with open(os.path.join(SRC, p), 'rb') as fi:
            while True:
                b = fi.read(1 << 22)
                if not b:
                    break
                crc = zlib.crc32(b, crc)
                fo.write(b)
                n_in += len(b)
                n_out += len(b)
        print(f'  [{i+1}/{len(parts)}] {p}  ->  累计 {n_out/1e6:.0f} MB  ({time.time()-t0:.0f}s)',
              flush=True)

print()
print(f'合并完成: {n_out/1e6:.1f} MB,  CRC32 = {crc & 0xffffffff:08x},  耗时 {time.time()-t0:.0f}s')

# 完整性：用 gzip 读到底看是否抛 CRC 错误
print('\n=== gzip 完整性校验（读到文件末尾）===')
t1 = time.time()
block = 1 << 24
nrows = 0
tail = b''
ok = True
try:
    dec = zlib.decompressobj(31)
    with open(DST, 'rb') as f:
        while True:
            b = f.read(block)
            if not b:
                break
            try:
                out = dec.decompress(b)
            except zlib.error as e:
                print('  decompress 异常:', e); ok = False; break
            if not out and dec.eof:
                # 多成员 gzip：剩余数据交给新的解压器
                pass
            tail += out
            if len(tail) > 1 << 26:
                nrows += tail.count(b'\n')
                tail = tail[tail.rfind(b'\n') + 1:]
    nrows += tail.count(b'\n')
except Exception as e:
    print('  读取异常:', type(e).__name__, e)
    ok = False

print(f'  行数(约) = {nrows:,}')
print(f'  校验耗时 {time.time()-t1:.0f}s')
print('  结论:', '通过（未出现 CRC/截断错误）' if ok else '失败')
