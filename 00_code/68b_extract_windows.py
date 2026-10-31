# -*- coding: utf-8 -*-
"""从 eQTLGen 全量 cis-eQTL (3.88 GB) 抽取 21 个胚系 IV 周围 ±1 Mb 窗口内的
全部 SNP-基因对（含零效应者）——标准 coloc.abf 的必需输入。

设计要点：
  1. 该文件是**多成员 gzip**，必须用循环新建 decompressobj 消费 unused_data
     （同 61_pqtl_extract.py 的修复；原实现只解出第一个 member）。
  2. IV 的 chr/pos 先用 Ensembl REST 拿到，从而**只需扫描全表一遍**。
"""
import io, os, sys, time, gzip, json
import requests, urllib3
urllib3.disable_warnings()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_FULL = os.path.join(BASE, '02_raw', 'eqtlgen', 'eQTLGen_cis_eQTL_full.txt.gz')
F_IV = os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv')
F_COORD = os.path.join(BASE, '03_meta', 'iv_coords.json')
F_OUT = os.path.join(BASE, '06_targets', 'eqtlgen_windows_full.tsv.gz')
F_LOG = os.path.join(BASE, '01_log', '68_extract_windows.txt')
WIN = 1_000_000
UA = {'User-Agent': 'Mozilla/5.0', 'Content-Type': 'application/json'}


def read_ivs():
    out = []
    for ln in io.open(F_IV, encoding='utf-8'):
        rs = ln.split()[0].strip()
        if rs.startswith('rs'):
            out.append(rs)
    return out


def get_coords(ivs):
    if os.path.exists(F_COORD):
        return json.load(io.open(F_COORD, encoding='utf-8'))
    coords = {}
    for i in range(0, len(ivs), 200):
        batch = ivs[i:i + 200]
        try:
            r = requests.post('https://rest.ensembl.org/variation/human',
                              headers=UA, timeout=60, verify=False,
                              data=json.dumps({'ids': batch}))
            if r.status_code != 200:
                continue
            for k, v in r.json().items():
                m = v.get('mappings', [])
                if m:
                    coords[k] = {'chr': str(m[0].get('seq_region_name')),
                                 'pos': int(m[0].get('start'))}
        except Exception:
            continue
        time.sleep(0.4)
    json.dump(coords, io.open(F_COORD, 'w', encoding='utf-8'), indent=1)
    return coords


def lines():
    """多成员 gzip 流式行生成器。"""
    import zlib
    do = zlib.decompressobj(31)
    tail = b''
    with open(F_FULL, 'rb') as f:
        while True:
            chunk = f.read(1 << 22)
            if not chunk:
                break
            buf = do.decompress(chunk)
            while do.eof:
                nd = zlib.decompressobj(31)
                buf = buf + nd.decompress(do.unused_data)
                do = nd
            tail += buf
            if b'\n' in tail:
                ps = tail.split(b'\n')
                tail = ps.pop()
                for L in ps:
                    yield L.decode('utf-8', 'replace')
    tail += do.flush()
    if tail:
        for L in tail.split(b'\n'):
            yield L.decode('utf-8', 'replace')


def main():
    ivs = read_ivs()
    print(f'[setup] 胚系 IV = {len(ivs)} 个', flush=True)
    coords = get_coords(ivs)
    print(f'[setup] 拿到坐标 {len(coords)}/{len(ivs)}', flush=True)
    if not coords:
        print('无坐标，退出'); return

    # chr -> [(lo,hi,snp)]
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
    need = ['SNP', 'GeneSymbol', 'SNPChr', 'SNPPos', 'Pvalue', 'Zscore',
            'AssessedAllele', 'OtherAllele', 'NrSamples', 'GeneChr', 'GenePos']
    miss = [n for n in need if n not in gi]
    if miss:
        print('!! 缺少列:', miss, '→ 请用现有列兜底'); print('可用列:', hdr); return

    c_chr, c_pos = gi['SNPChr'], gi['SNPPos']
    nscan = nkeep = 0
    t0 = time.time()
    with gzip.open(F_OUT, 'wt', compresslevel=1) as fo:
        fo.write('\t'.join(hdr) + '\n')
        for i, L in enumerate(it):
            nscan += 1
            p = L.split('\t')
            if len(p) <= max(c_chr, c_pos):
                continue
            ch = p[c_chr].strip()
            w = winmap.get(ch)
            if not w:
                continue
            try:
                ps = int(p[c_pos])
            except ValueError:
                continue
            # 二分查找落在哪个窗口
            lo_j, hi_j = 0, len(w) - 1
            hit = False
            while lo_j <= hi_j:
                mid = (lo_j + hi_j) // 2
                a, b, _ = w[mid]
                if ps < a:
                    hi_j = mid - 1
                elif ps > b:
                    lo_j = mid + 1
                else:
                    hit = True
                    break
            if hit:
                fo.write(L.rstrip('\r\n') + '\n')
                nkeep += 1
            if nscan % 5_000_000 == 0:
                el = time.time() - t0
                print(f'  扫 {nscan/1e6:>7.1f}M 行  抽出 {nkeep:>9,}  '
                      f'{el:>6.0f}s  预计剩余 {(3880/(nscan/max(nscan,1)))*0+el/nscan*(35e6-nscan)/60:.0f}min',
                      flush=True)
    el = time.time() - t0
    print()
    print(f'=== 完成 ===  扫描 {nscan:,} 行，抽出 {nkeep:,} 行，耗时 {el/60:.1f} min')
    print('输出:', F_OUT)
    with io.open(F_LOG, 'w', encoding='utf-8') as f:
        f.write(f'rows_scanned={nscan}\nrows_extracted={nkeep}\nelapsed_sec={el:.0f}\n')


if __name__ == '__main__':
    main()
