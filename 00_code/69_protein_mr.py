# -*- coding: utf-8 -*-
"""69_protein_mr.py —— P0 两步蛋白孟德尔随机化

Stage1  CHIP(胚系 IV) -> 20 个炎症蛋白      : 用 pQTL 中 21 个 IV 所在行，Wald + IVW
Stage2  蛋白 -> MPN / PV                    : cis 工具(p<5e-8, ±1Mb, 250kb 修剪) + IVW
Med     中介比例 = a*b / c                  : c 为 CHIP -> 结局的总效应(已有 07_mr/*_MR.csv)

用法:
  python 69_protein_mr.py step1
  python 69_protein_mr.py cis
  python 69_protein_mr.py subset
  python 69_protein_mr.py step2
  python 69_protein_mr.py med
  python 69_protein_mr.py all
"""
import io, os, sys, gzip, json, math, time
import numpy as np
import requests, urllib3
urllib3.disable_warnings()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PQ = os.path.join(BASE, '02_raw', 'pqtl')
RAW = os.path.join(BASE, '02_raw')
META = os.path.join(BASE, '03_meta')
OUT = os.path.join(BASE, '07_mr_protein')
os.makedirs(OUT, exist_ok=True)

SKIP = set()          # TNF 的 pQTL 首次抽取走的是被掐断的流，正在重抽；重抽完删掉这一行
SYMS = ['IL6', 'TNF', 'IL1A', 'CXCL8', 'CCL2', 'CCL3', 'TGFB1_LAP', 'TGFA',
        'MMP1', 'MMP10', 'VEGFA', 'TNFRSF11B', 'IFNG', 'IL10', 'LTA', 'IL17A',
        'IL18', 'OSM', 'LIF', 'KITLG']
GENE_OF = {'TGFB1_LAP': 'TGFB1'}
FINN = {'MPN': 'finngen_R12_MYELOPROF_NONCML_MPN.gz',
        'PV':  'finngen_R12_POLYCYTVERA_PV.gz'}
CIS = 1_000_000
PRUNE = 250_000
PV_SIG = 5e-8
UA = {'User-Agent': 'Mozilla/5.0'}


def read_tsv(p):
    with io.open(p, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for ln in f:
            if ln.strip():
                yield dict(zip(hdr, ln.rstrip('\n').split('\t')))


def chip_iv():
    d = {}
    for r in read_tsv(os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv')):
        d[r['SNP']] = {'chr': r['chr'], 'pos': int(r['pos']),
                       'ea': r['effect_allele'], 'oa': r['other_allele'],
                       'b': float(r['beta']), 'se': float(r['se']),
                       'p': float(r['pval'])}
    return d


def ivw(bx, by, sy):
    """固定效应 IVW + 乘性随机效应校正"""
    bx, by, sy = map(np.asarray, (bx, by, sy))
    w = 1.0 / sy ** 2
    num = np.sum(w * bx * by)
    den = np.sum(w * bx ** 2)
    b = num / den
    se = math.sqrt(1.0 / den)
    k = len(bx)
    q = float(np.sum(w * (by - b * bx) ** 2))
    q_p = 1.0 - _chi2cdf(q, max(k - 1, 1))
    if k > 1 and q > (k - 1):
        se *= math.sqrt(q / (k - 1))
    z = b / se if se > 0 else 0.0
    p = 2 * (1 - _normcdf(abs(z)))
    return dict(k=k, b=b, se=se, z=z, p=p, q=q, q_p=q_p)


def _normcdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def _chi2cdf(x, k):
    """正则化下不完全 gamma P(k/2, x/2) — 级数法，k 小时够用"""
    if x <= 0:
        return 0.0
    a, xx = k / 2.0, x / 2.0
    if xx < a + 1:
        term = 1.0 / a; s = term; n = 0
        while abs(term) > 1e-12 * abs(s) and n < 500:
            n += 1
            term *= xx / (a + n); s += term
        return s * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
    # 连分式
    tiny = 1e-300
    b = xx + 1 - a; c = 1 / tiny; d = 1 / b; h = d
    for i in range(1, 300):
        an = -i * (i - a)
        b += 2
        d = an * d + b; d = tiny if abs(d) < tiny else d
        c = b + an / c; c = tiny if abs(c) < tiny else c
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return 1 - h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))


# ---------------- Stage 1 ----------------
def stage1():
    iv = chip_iv()
    rows = []
    detail = []
    for sym in SYMS:
        if sym in SKIP:
            print(f'  [skip] {sym} 暂缓（数据待重抽）'); continue
        f = os.path.join(PQ, f'{sym}_step1_rows.tsv')
        if not os.path.exists(f):
            print(f'  [skip] {sym} 无 step1 文件'); continue
        if os.path.getsize(f) < 300:
            print(f'  [skip] {sym} step1 文件只有表头（抽取失败）'); continue
        bx, by, sy, used = [], [], [], []
        for r in read_tsv(f):
            rs = r.get('rsid', '')
            if rs not in iv:
                continue
            c = iv[rs]
            ea, oa = r['effect_allele'], r['other_allele']
            if ea == c['ea'] and oa == c['oa']:
                s = 1
            elif ea == c['oa'] and oa == c['ea']:
                s = -1
            else:
                continue
            bp = float(r['beta']) * s
            sp = float(r['standard_error'])
            wal = bp / c['b']
            sew = sp / abs(c['b'])
            bx.append(c['b']); by.append(bp); sy.append(sp)
            used.append(rs)
            detail.append(dict(symbol=sym, SNP=rs, chip_b=c['b'], chip_se=c['se'],
                               prot_b=bp, prot_se=sp, wald=wal, wald_se=sew))
        if len(bx) < 2:
            print(f'  [skip] {sym} 可用 IV 仅 {len(bx)}'); continue
        res = ivw(bx, by, sy)
        res.update(symbol=sym, n_iv=len(bx), snps=','.join(used))
        rows.append(res)
        print(f'  {sym:<10} k={res["k"]:<3} beta={res["b"]:+.4f} se={res["se"]:.4f} '
              f'p={res["p"]:.3g}  Qp={res["q_p"]:.2g}', flush=True)
    with io.open(os.path.join(OUT, 'step1_chip_to_protein.tsv'), 'w', encoding='utf-8') as f:
        f.write('symbol\tn_iv\tbeta_SD_per_logOR\tse\tz\tp\tQ\tQ_p\tsnps\n')
        for r in rows:
            f.write(f'{r["symbol"]}\t{r["n_iv"]}\t{r["b"]:.6g}\t{r["se"]:.6g}\t'
                    f'{r["z"]:.4f}\t{r["p"]:.4g}\t{r["q"]:.3f}\t{r["q_p"]:.4g}\t{r["snps"]}\n')
    with io.open(os.path.join(OUT, 'step1_detail.tsv'), 'w', encoding='utf-8') as f:
        f.write('symbol\tSNP\tchip_b\tchip_se\tprot_b\tprot_se\twald\twald_se\n')
        for d in detail:
            f.write(f'{d["symbol"]}\t{d["SNP"]}\t{d["chip_b"]:.6g}\t{d["chip_se"]:.6g}\t'
                    f'{d["prot_b"]:.6g}\t{d["prot_se"]:.6g}\t{d["wald"]:.6g}\t{d["wald_se"]:.6g}\n')
    print(f'[step1] 完成 {len(rows)} 个蛋白 -> {OUT}')


# ---------------- cis 工具 ----------------
def gene_coords():
    f = os.path.join(META, 'protein_gene_coords.json')
    if os.path.exists(f):
        return json.load(io.open(f, encoding='utf-8'))
    out = {}
    for sym in SYMS:
        g = GENE_OF.get(sym, sym)
        try:
            r = requests.get(f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}',
                             headers=UA, timeout=60, verify=False,
                             params={'expand': 0, 'content-type': 'application/json'})
            if r.status_code == 200:
                j = r.json()
                out[sym] = {'gene': g, 'chr': str(j['seq_region_name']),
                            'start': int(j['start']), 'end': int(j['end'])}
                print(f'  {sym}: chr{j["seq_region_name"]}:{j["start"]}-{j["end"]}')
        except Exception as e:
            print(f'  [coord-fail] {sym} {type(e).__name__}')
        time.sleep(0.2)
    json.dump(out, io.open(f, 'w', encoding='utf-8'), indent=1)
    return out


def stage_cis():
    co = gene_coords()
    allrows, meta = [], []
    for sym in SYMS:
        f = os.path.join(PQ, f'{sym}_instruments.tsv')
        if not (os.path.exists(f) and sym in co):
            continue
        c = co[sym]
        cand = []
        for r in read_tsv(f):
            if r['chromosome'] != c['chr']:
                continue
            try:
                pos = int(r['base_pair_location']); pv = float(r['p_value'])
            except ValueError:
                continue
            if abs(pos - (c['start'] + c['end']) // 2) <= CIS:
                cand.append((pv, pos, r))
        cand.sort(key=lambda x: x[0])
        keep, kept_pos = [], []
        for pv, pos, r in cand:
            if all(abs(pos - q) > PRUNE for q in kept_pos):
                keep.append(r); kept_pos.append(pos)
        meta.append((sym, len(cand), len(keep)))
        for r in keep:
            allrows.append((sym, r))
        print(f'  {sym:<10} cis 候选 {len(cand):<4} 修剪后 {len(keep)}', flush=True)
    with io.open(os.path.join(OUT, 'protein_instruments.tsv'), 'w', encoding='utf-8') as f:
        f.write('symbol\tSNP\tchr\tpos\tea\toa\tbeta\tse\tp\teaf\tn\n')
        for sym, r in allrows:
            f.write(f'{sym}\t{r["rsid"]}\t{r["chromosome"]}\t{r["base_pair_location"]}\t'
                    f'{r["effect_allele"]}\t{r["other_allele"]}\t{r["beta"]}\t'
                    f'{r["standard_error"]}\t{r["p_value"]}\t{r["effect_allele_frequency"]}\t'
                    f'{r["n"]}\n')
    print(f'[cis] 共 {len(allrows)} 个工具 SNP（{len(meta)} 个蛋白）')


# ---------------- FinnGen 子集 ----------------
def stage_subset():
    need = set()
    for r in read_tsv(os.path.join(OUT, 'protein_instruments.tsv')):
        need.add(r['SNP'])
    print(f'[subset] 需要 {len(need)} 个 rsID')
    for phe, fn in FINN.items():
        p = os.path.join(RAW, fn)
        outp = os.path.join(RAW, f'finngen_subset_{phe}.tsv')
        if os.path.exists(outp):
            print(f'  {phe} 已存在，跳过'); continue
        t0 = time.time(); n = 0; hit = 0
        with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as fi, \
                io.open(outp, 'w', encoding='utf-8') as fo:
            fo.write('#chrom\tpos\tref\talt\trsids\tbeta\tsebeta\tpval\taf_alt\n')
            hdr = fi.readline()
            for ln in fi:
                n += 1
                f = ln.rstrip('\n').split('\t')
                if len(f) < 13:
                    continue
                ids = f[4]
                if ',' in ids:
                    hitid = [x for x in ids.split(',') if x in need]
                    id0 = hitid[0] if hitid else None
                else:
                    id0 = ids if ids in need else None
                if id0:
                    fo.write(f'{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}\t{id0}\t'
                             f'{f[8]}\t{f[9]}\t{f[6]}\t{f[10]}\n')
                    hit += 1
        print(f'  {phe}: 扫 {n:,} 行, 命中 {hit:,}, {time.time()-t0:.0f}s', flush=True)


# ---------------- Stage 2 ----------------
def stage2():
    ins = list(read_tsv(os.path.join(OUT, 'protein_instruments.tsv')))
    out = {}
    for phe in FINN:
        out[phe] = {r['rsids']: r for r in read_tsv(os.path.join(RAW, f'finngen_subset_{phe}.tsv'))}
    rows = []
    for sym in SYMS:
        sub = [r for r in ins if r['symbol'] == sym]
        if not sub:
            continue
        for phe in FINN:
            bx, by, sy = [], [], []
            for r in sub:
                o = out[phe].get(r['SNP'])
                if not o:
                    continue
                # 统一到 pQTL 的 effect allele
                if o['alt'] == r['ea'] and o['ref'] == r['oa']:
                    s = 1
                elif o['alt'] == r['oa'] and o['ref'] == r['ea']:
                    s = -1
                else:
                    continue
                bo = float(o['beta']) * s
                so = float(o['sebeta'])
                bp = float(r['beta'])
                if abs(bp) < 1e-12:
                    continue
                bx.append(bp); by.append(bo); sy.append(so)
            if len(bx) < 1:
                continue
            res = ivw(bx, by, sy) if len(bx) > 1 else dict(
                k=1, b=by[0] / bx[0], se=sy[0] / abs(bx[0]), z=0, p=0, q=0, q_p=1)
            if len(bx) == 1:
                z = res['b'] / res['se']; res['z'] = z
                res['p'] = 2 * (1 - _normcdf(abs(z)))
            res.update(symbol=sym, outcome=phe)
            rows.append(res)
            print(f'  {sym:<10} -> {phe:<4} k={res["k"]:<3} beta={res["b"]:+.4f} '
                  f'p={res["p"]:.3g}', flush=True)
    with io.open(os.path.join(OUT, 'step2_protein_to_outcome.tsv'), 'w', encoding='utf-8') as f:
        f.write('symbol\toutcome\tk\tbeta_logOR_per_SD\tse\tz\tp\tQ_p\n')
        for r in rows:
            f.write(f'{r["symbol"]}\t{r["outcome"]}\t{r["k"]}\t{r["b"]:.6g}\t'
                    f'{r["se"]:.6g}\t{r["z"]:.4f}\t{r["p"]:.4g}\t{r["q_p"]:.4g}\n')
    print(f'[step2] 完成 {len(rows)} 条')


# ---------------- 中介 ----------------
def stage_med():
    s1 = {r['symbol']: r for r in read_tsv(os.path.join(OUT, 'step1_chip_to_protein.tsv'))}
    s2 = list(read_tsv(os.path.join(OUT, 'step2_protein_to_outcome.tsv')))
    # 总效应 c：CHIP -> MPN/PV（读已有 IVW 结果）
    c = {}
    for phe in FINN:
        f = os.path.join(BASE, '07_mr', f'CHIPg_to_{phe}_MR.csv')
        if os.path.exists(f):
            import csv
            with io.open(f, encoding='utf-8') as fh:      # 注意：这是逗号分隔的 CSV
                for r in csv.DictReader(fh):
                    if r.get('method', '').startswith('Inverse variance'):
                        c[phe] = {'b': float(r['b']), 'se': float(r['se']),
                                  'p': float(r['pval'])}
                        break
    fout = os.path.join(OUT, 'mediation.tsv')
    with io.open(fout, 'w', encoding='utf-8') as f:
        f.write('symbol\toutcome\ta_chip_to_protein\tse_a\tb_protein_to_outcome\tse_b\t'
                'indirect\tsd_indirect\tprop_mediated\tc_total\tp_a\tp_b\n')
        for r in s2:
            sym = r['symbol']
            if sym not in s1 or r['outcome'] not in c:
                continue
            a, sa = float(s1[sym]['beta_SD_per_logOR']), float(s1[sym]['se'])
            b, sb = float(r['beta_logOR_per_SD']), float(r['se'])
            ind = a * b
            sd = math.sqrt((a * sb) ** 2 + (b * sa) ** 2)
            ct = c[r['outcome']]['b']
            prop = ind / ct if abs(ct) > 1e-9 else float('nan')
            f.write(f'{sym}\t{r["outcome"]}\t{a:.6g}\t{sa:.6g}\t{b:.6g}\t{sb:.6g}\t'
                    f'{ind:.6g}\t{sd:.6g}\t{prop:.4f}\t{ct:.6g}\t'
                    f'{s1[sym]["p"]}\t{r["p"]}\n')
    print('[med] ->', fout)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if cmd in ('step1', 'all'): stage1()
    if cmd in ('cis', 'all'): stage_cis()
    if cmd in ('subset', 'all'): stage_subset()
    if cmd in ('step2', 'all'): stage2()
    if cmd in ('med', 'all'): stage_med()
