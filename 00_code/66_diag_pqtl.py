# -*- coding: utf-8 -*-
"""pQTL 抽取诊断：限定行数，逐环节统计丢弃数，定位 0 行的根因。"""
import io, os, sys, time, zlib
import requests, urllib3
urllib3.disable_warnings()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = {'User-Agent': 'Mozilla/5.0'}
FTP = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
       'GCST90274001-GCST90275000/GCST{acc}/GCST{acc}.tsv.gz')

# 胚系 IV
iv = set()
for ln in io.open(os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv'), encoding='utf-8'):
    rs = ln.split()[0].strip()
    if rs.startswith('rs'):
        iv.add(rs)
print('[setup] IV =', len(iv), '个', flush=True)

MAXLINES = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
ACC = sys.argv[2] if len(sys.argv) > 2 else '90274839'   # TNF

t0 = time.time()
url = FTP.format(acc=ACC)
r = requests.get(url, headers=UA, stream=True, timeout=120, verify=False)
print('[http] status', r.status_code, 'size', r.headers.get('Content-Length'), flush=True)

do = zlib.decompressobj(31)
tail = ''
nb_total = nb_short = nb_floatfail = nb_ivhit = nb_sig = 0
minp = 1.0
seen_rs = 0
first_rows = []


def handle(L):
    global nb_total, nb_short, nb_floatfail, nb_ivhit, nb_sig, minp, seen_rs
    nb_total += 1
    p = L.split('\t')
    if len(p) < 11:
        nb_short += 1
        return
    try:
        pv = float(p[7])
    except ValueError:
        nb_floatfail += 1
        return
    if pv < minp:
        minp = pv
    if pv < 5e-8:
        nb_sig += 1
    if p[9].startswith('rs'):
        seen_rs += 1
    if p[9] in iv:
        nb_ivhit += 1


for ci, chunk in enumerate(r.iter_content(1 << 18)):
    tail += do.decompress(chunk).decode('utf-8', 'replace')
    if '\n' in tail:
        ps = tail.split('\n')
        tail = ps.pop()
        for L in ps:
            if len(first_rows) < 3 and nb_total > 0:
                first_rows.append(L)
            handle(L)
    if nb_total >= MAXLINES:
        break
    if ci % 200 == 0:
        print(f'  chunk{ci:>5}  lines={nb_total:>8}  sig={nb_sig:>4}  '
              f'ivhit={nb_ivhit}  floatfail={nb_floatfail}  short={nb_short}  '
              f'elapsed={time.time()-t0:.0f}s', flush=True)
r.close()

print()
print('=========== 诊断结果 ===========')
print(f'  已处理行数        : {nb_total}')
print(f'  列数<11 被丢弃    : {nb_short}')
print(f'  float(p) 失败     : {nb_floatfail}')
print(f'  rsid 以 rs 开头   : {seen_rs}')
print(f'  IV 命中行数       : {nb_ivhit}')
print(f'  p<5e-8 行数       : {nb_sig}')
print(f'  已见最小 p 值     : {minp:.3e}')
print(f'  耗时              : {time.time()-t0:.0f} s')
print()
if nb_total == 0:
    print('  ==> 结论：解压/迭代环节完全没产出，lines_from 或 iter_content 坏了')
elif nb_sig == 0 and minp > 1e-5:
    print('  ==> 结论：能读到行，但扫到的区间里没有强信号（可能只是扫描不够远，')
    print('      因为文件按位置排序，chr1 开头通常无显著位点）')
elif nb_sig > 0:
    print('  ==> 结论：解析正常。若完整跑仍为 0，则是流式读取被中途中断。')
