# -*- coding: utf-8 -*-
"""61h_rescan_local.py —— 用已落盘的 pQTL 文件重跑抽取（不重新下载）

场景：61f 的 TNF 任务下载完整（313 MB），扫描也读完了 12.6M 行，
但写结果文件时抛 PermissionError（文件被占用/权限瞬时丢失），随后把数据
文件改名成 .bad。此脚本把 .bad 改回来，直接调用 61f 的 _scan 重新抽取。

用法: python 61h_rescan_local.py TNF GCST90274839
"""
import os, sys, importlib.util

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, '00_code'))

spec = importlib.util.spec_from_file_location(
    'pqtl', os.path.join(BASE, '00_code', '61f_dl_pqtl_mt.py'))
pqtl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pqtl)

sym, acc = sys.argv[1], sys.argv[2]
dest = os.path.join(pqtl.TMP, f'GCST{acc}.tsv.gz')
bad = dest + '.bad'
if not os.path.exists(dest) and os.path.exists(bad):
    os.replace(bad, dest)
    print(f'[restore] {os.path.basename(bad)} -> {os.path.basename(dest)} '
          f'({os.path.getsize(dest)/1e6:.0f} MB)', flush=True)
print('[rescan]', sym, dest, flush=True)
print(pqtl._scan(dest, sym,
                 os.path.join(pqtl.OUT, f'{sym}_step1_rows.tsv'),
                 os.path.join(pqtl.OUT, f'{sym}_instruments.tsv'),
                 __import__('time').time()), flush=True)
