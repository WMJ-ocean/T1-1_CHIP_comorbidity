# -*- coding: utf-8 -*-
"""77_bench_ftp.py —— 测 EBI FTP 协议带宽（HTTP 被限速时试另一条路）"""
import ftplib, time

FN = 'GCST90165271_buildGRCh38.tsv.gz'
DIR = '/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165271/'


class Stop(Exception):
    pass


ftp = ftplib.FTP('ftp.ebi.ac.uk', timeout=60)
ftp.login()
ftp.cwd(DIR)
size = ftp.size(FN)
print(f'FTP ok, size = {size/2**30:.2f} GiB')

got = [0]
t0 = time.time()
LIMIT = 24 * 2**20


def cb(b):
    got[0] += len(b)
    if got[0] >= LIMIT:
        raise Stop


try:
    ftp.retrbinary(f'RETR {FN}', cb, blocksize=1 << 16)
except Stop:
    pass
el = time.time() - t0
rate = got[0] / 2**20 / max(el, .01)
print(f'FTP 单连接: {got[0]/2**20:.1f} MB / {el:.1f}s = {rate:.2f} MB/s')
print(f'→ 全文件预计 {size/2**20/rate/60:.0f} min')
ftp.close()
