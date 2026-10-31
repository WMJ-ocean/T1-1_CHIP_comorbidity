# -*- coding: utf-8 -*-
"""68d_extract_windows_fixed.py —— 修正分片接缝后重抽 coloc 窗口

根因（2026-10-29）：59_dl_cis_full.py 里
    edges = [int(size*i/N) ...] + [size-1]
    a, b = edges[i], edges[i+1]        # b 是「下一个分片的首字节」
    Range: bytes={a}-{b}               # 闭区间 -> 每个接缝多下 1 个字节
=> 15 个分片接缝各重复 1 字节（实测 15/15 接缝 last==first，概率上不可能是巧合）。
   合并文件 3.88 GB 在 1937 MB 处炸「invalid bit length repeat」，且接缝之后的
   解压输出自 242 MB 起就已污染 —— 68b 之前抽出的 236k 行全部作废。

本脚本**不重建 3.88 GB 文件**，直接按分片流式读取，跳过每个接缝的重复首字节。
"""
import io, os, gzip, json, time, zlib

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, '02_raw', 'eqtlgen', '_parts_full')
F_IV = os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv')
F_COORD = os.path.join(BASE, '03_meta', 'iv_coords.json')
F_OUT = os.path.join(BASE, '06_targets', 'eqtlgen_windows_full.tsv.gz')
F_LOG = os.path.join(BASE, '01_log', '68d_extract_windows.log')
WIN = 1_000_000

parts = sorted([f for f in os.listdir(SRC) if f.startswith('full.p')],
               key=lambda s: int(s.split('.p')[-1]))
print(f'[setup] 分片 {len(parts)} 个', flush=True)


def raw_bytes():
    """按分片产出字节流，跳过接缝重复字节"""
    prev_last = None
    nskip = 0
    for i, p in enumerate(parts):
        with open(os.path.join(SRC, p), 'rb') as fi:
            first = True
            while True:
                b = fi.read(1 << 22)
                if not b:
                    break
                if first:
                    first = False
                    if prev_last is not None:
                        if b[0] == prev_last:
                            b = b[1:]; nskip += 1          # 丢掉重复字节
                        else:
                            print(f'  [warn] p{i:03d} 接缝不等 '
                                  f'({prev_last:#04x}!={b[0]:#04x})', flush=True)
                prev_last = b[-1]
                yield b
    print(f'[setup] 跳过接缝重复字节 {nskip} 个', flush=True)


def lines():
    do = zlib.decompressobj(31); tail = b''
    for chunk in raw_bytes():
        buf = do.decompress(chunk)
        while do.eof:
            nd = zlib.decompressobj(31)
            buf += nd.decompress(do.unused_data); do = nd
        tail += buf
        if b'\n' in tail:
            ps = tail.split(b'\n'); tail = ps.pop()
            for L in ps:
                yield L.decode('utf-8', 'replace')
    tail += do.flush()
    if tail:
        for L in tail.split(b'\n'):
            yield L.decode('utf-8', 'replace')


def read_ivs():
    out = []
    for ln in io.open(F_IV, encoding='utf-8'):
        rs = ln.split()[0].strip()
        if rs.startswith('rs'):
            out.append(rs)
    return out


def main():
    ivs = read_ivs()
    coords = json.load(io.open(F_COORD, encoding='utf-8'))
    print(f'[setup] IV = {len(ivs)} 个，坐标 {len(coords)} 个', flush=True)

    winmap = {}
    for rs, c in coords.items():
        winmap.setdefault(str(c['chr']), []).append((c['pos'] - WIN, c['pos'] + WIN, rs))
    for ch in winmap:
        winmap[ch].sort()
    print('[setup] 目标染色体:', sorted(winmap.keys(), key=lambda s: (len(s), s)), flush=True)

    it = lines()
    hdr = next(it).split('\t')
    print('[hdr]', hdr, flush=True)
    gi = {n: i for i, n in enumerate(hdr)}
    c_chr, c_pos = gi['SNPChr'], gi['SNPPos']

    nscan = nkeep = 0
    t0 = time.time()
    with gzip.open(F_OUT, 'wt', compresslevel=1) as fo:
        fo.write('\t'.join(hdr) + '\n')
        for L in it:
            nscan += 1
            if (nscan % 2_000_000) == 0:
                el = time.time() - t0
                print(f'  扫 {nscan/1e6:>6.1f}M 行  抽出 {nkeep:>9,}  {el:>6.0f}s', flush=True)
            p = L.split('\t')
            if len(p) <= max(c_chr, c_pos):
                continue
            w = winmap.get(p[c_chr].strip())
            if not w:
                continue
            try:
                ps = int(p[c_pos])
            except ValueError:
                continue
            lo, hi = 0, len(w) - 1
            while lo <= hi:
                mid = (lo + hi) // 2
                a, b, _ = w[mid]
                if ps < a:
                    hi = mid - 1
                elif ps > b:
                    lo = mid + 1
                else:
                    fo.write(L.rstrip('\r\n') + '\n'); nkeep += 1
                    break
    el = time.time() - t0
    print(f'=== 完成 === 扫描 {nscan:,} 行，抽出 {nkeep:,} 行，{el/60:.1f} min', flush=True)
    with io.open(F_LOG, 'a', encoding='utf-8') as f:
        f.write(f'rows_scanned={nscan}\nrows_extracted={nkeep}\nelapsed_sec={el:.0f}\n')


if __name__ == '__main__':
    main()
