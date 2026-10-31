# -*- coding: utf-8 -*-
"""
61_pqtl_extract.py

流式下载 Zheng 2020 Nat Genet 的 91 个炎症蛋白 pQTL，边解压边抽：
  A) 21 个 CHIP 胚系 IV 所在行      -> Step 1 (CHIP -> protein) 用
  B) p < 5e-8 的全部行              -> Step 2 (protein -> MPN) 的工具变量候选

不落原始 320MB/蛋白；不调用 os.remove（沙箱禁止）。
"""
import os, sys, io, json, time, gzip, zlib
import concurrent.futures as cf
import requests, urllib3

urllib3.disable_warnings()

ROOT = r'C:\Users\Hero\WorkBuddy\2026-10-14-23-35-01\outputs\pv_hotspot\T1-1_CHIP_comorbidity'
OUT  = os.path.join(ROOT, '02_raw', 'pqtl')
META = os.path.join(ROOT, '03_meta')
os.makedirs(OUT, exist_ok=True)

FTP = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST{acc}/GCST{acc}.tsv.gz'
UA  = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research'}

# ---------------- CHIP 胚系 IV ----------------
IVS = set()
with io.open(os.path.join(ROOT, '05_iv', 'CHIP_germline_iv.tsv'), encoding='utf-8') as f:
    for line in f:
        p = line.split('\t')
        if len(p) > 1 and p[0].startswith('rs'):
            IVS.add(p[0].strip())
# 去掉表头残留
IVS = {s for s in IVS if s.startswith('rs') and len(s) > 2}
print(f'[setup] CHIP 胚系 IV = {len(IVS)} 个', flush=True)

# ---------------- 预设蛋白 panel ----------------
# 依据 MPN -> 炎症微环境 -> 骨髓纤维化 主线预先设定，非事后挑选
PANEL = [
    # Tier 1 : MPN 炎症轴 / 纤维化轴核心
    ('GCST90274815', 'IL6',          'Interleukin-6'),
    ('GCST90274839', 'TNF',          'Tumor necrosis factor'),
    ('GCST90274805', 'IL1A',         'Interleukin-1-alpha'),
    ('GCST90274817', 'CXCL8',        'Interleukin-8'),
    ('GCST90274821', 'CCL2',         'MCP-1'),
    ('GCST90274825', 'CCL3',         'MIP-1alpha'),
    ('GCST90274818', 'TGFB1_LAP',    'Latency-associated peptide TGF-beta 1'),
    ('GCST90274838', 'TGFA',         'TGF-alpha'),
    ('GCST90274826', 'MMP1',         'Matrix metalloproteinase-1'),
    ('GCST90274827', 'MMP10',        'Matrix metalloproteinase-10'),
    ('GCST90274848', 'VEGFA',        'VEGF-A'),
    ('GCST90274830', 'TNFRSF11B',    'Osteoprotegerin'),
    # Tier 2 : 强生物学依据
    ('GCST90274794', 'IFNG',         'Interferon gamma'),
    ('GCST90274795', 'IL10',         'Interleukin-10'),
    ('GCST90274840', 'LTA',          'TNF-beta'),
    ('GCST90274801', 'IL17A',        'Interleukin-17A'),
    ('GCST90274803', 'IL18',         'Interleukin-18'),
    ('GCST90274831', 'OSM',          'Oncostatin-M'),
    ('GCST90274819', 'LIF',          'Leukemia inhibitory factor'),
    ('GCST90274833', 'KITLG',        'Stem cell factor'),
]
# ---------------- 按带宽现实收缩的执行子集 ----------------
# 修复后实测：单个 320 MB 文件约 4.7 分钟即可扫完 ~650 万行（而非原先的空转 21 分钟）。
# 20 个蛋白 ≈ 1.6 小时，可接受 -> 恢复全集。
CORE6 = None          # None = 全部执行
if CORE6:
    PANEL = [p for p in PANEL if p[1] in CORE6]
# Tier2/Tier3 完整清单保留在 00_code/_pqtl_panel_full.py 备查。

HDR = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
       'beta', 'standard_error', 'effect_allele_frequency', 'p_value',
       'variant_id', 'rsid', 'n']


def lines_from(resp):
    """把 HTTP gzip 流解压成文本行生成器。

    ★ 关键修正（2026-10-29）：GWAS Catalog 的 .tsv.gz 是**多成员 gzip**
      ——由大量约 21 KB 的小 member 串接而成（实测 25 MB 内含 1,194 个成员）。
      单一 decompressobj 解完第一个 member（仅 0.07 MB / 907 行）后
      `eof=True` 且不再产出，这正是原先「抽取 0 行」的根因。
      必须循环用 unused_data 开新解压器。
    """
    do = zlib.decompressobj(31)
    tail = ''
    for chunk in resp.iter_content(1 << 18):
        if not chunk:
            continue
        buf = do.decompress(chunk)
        while do.eof:                       # 当前 member 已结束 -> 接下一个
            nd = zlib.decompressobj(31)
            buf = buf + nd.decompress(do.unused_data)
            do = nd
        tail += buf.decode('utf-8', 'replace')
        if '\n' in tail:
            parts = tail.split('\n')
            tail = parts.pop()
            for L in parts:
                yield L
    tail += do.flush().decode('utf-8', 'replace')
    if tail:
        for L in tail.split('\n'):
            yield L


def one(item, attempts=3):
    acc, sym, name = item
    acc = acc.replace('GCST', '')          # 路径模板里已含 GCST 前缀
    f_iv  = os.path.join(OUT, f'{sym}_step1_rows.tsv')
    f_ins = os.path.join(OUT, f'{sym}_instruments.tsv')
    if os.path.exists(f_iv) and os.path.exists(f_ins):
        # 上一轮 bug 版产生的只有表头的空文件（<200B）必须剔除，否则会被误判为已完成。
        # 沙箱禁止 os.remove，改为 os.replace 归档。
        sz = os.path.getsize(f_iv) + os.path.getsize(f_ins)
        if sz < 400:
            arc = os.path.join(OUT, '_bad')
            os.makedirs(arc, exist_ok=True)
            for f in (f_iv, f_ins):
                try:
                    os.replace(f, os.path.join(arc, os.path.basename(f)))
                except Exception as e:
                    print('  [warn] 归档失败', f, e, flush=True)
            print(f'  [clean] {sym} 空结果已归档，本轮重抽', flush=True)
        else:
            return sym, 'CACHED', 0, 0

    for att in range(attempts):
        try:
            r = requests.get(FTP.format(acc=acc), headers=UA, stream=True,
                             timeout=120, verify=False)
            if r.status_code not in (200, 206):
                r.close(); time.sleep(3); continue

            hit, ins, nrow = [], [], 0
            for i, L in enumerate(lines_from(r)):
                if i == 0:
                    continue                      # 表头
                nrow += 1
                parts = L.split('\t')
                if len(parts) < 11:
                    continue
                try:
                    pv = float(parts[7])
                except ValueError:
                    continue
                if parts[9] in IVS:
                    hit.append(L.rstrip('\r'))
                if pv < 5e-8:
                    ins.append(L.rstrip('\r'))
            r.close()

            with io.open(f_iv, 'w', encoding='utf-8') as fo:
                fo.write('\t'.join(HDR) + '\n')
                fo.write('\n'.join(hit) + ('\n' if hit else ''))
            with io.open(f_ins, 'w', encoding='utf-8') as fo:
                fo.write('\t'.join(HDR) + '\n')
                fo.write('\n'.join(ins) + ('\n' if ins else ''))
            return sym, 'OK', len(hit), len(ins), nrow
        except Exception as e:
            print(f'    [{sym}] attempt{att} {type(e).__name__}: {e}', flush=True)
            time.sleep(4)
    return sym, 'FAIL', 0, 0, 0


if __name__ == '__main__':
    t0 = time.time()
    res = []
    with cf.ProcessPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(one, it): it for it in PANEL}
        for fu in cf.as_completed(futs):
            it = futs[fu]
            try:
                out = fu.result()
            except Exception as e:
                out = (it[1], 'ERR:' + type(e).__name__, 0, 0, 0)
            res.append(out)
            print(f'  {out[0]:<12} {out[1]:<12} step1_rows={out[2]:<4} '
                  f'instruments={out[3]:<6} ({time.time()-t0:.0f}s)', flush=True)

    # 保留原来的顺序写汇总
    order = {s: i for i, (_, s, _) in enumerate(PANEL)}
    res.sort(key=lambda x: order.get(x[0], 999))
    with io.open(os.path.join(META, 'pqtl_extract_summary.tsv'), 'w', encoding='utf-8') as fo:
        fo.write('symbol\tstatus\tn_step1_rows\tn_instruments\n')
        for r in res:
            fo.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n')
    print(f'\n完成，用时 {time.time()-t0:.0f}s')
