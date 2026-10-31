# -*- coding: utf-8 -*-
"""16 线程分块下载 CHIP harmonised sumstats（支持断点续传）"""
import requests, os, time, threading, sys
requests.packages.urllib3.disable_warnings()
UA = {'User-Agent': 'Mozilla/5.0'}

U = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165267/harmonised/GCST90165267.h.tsv.gz'
BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/02_raw'
TMP = os.path.join(BASE, '_chip_parts')
os.makedirs(TMP, exist_ok=True)
FINAL = os.path.join(BASE, 'GCST90165267_CHIP_all.h.tsv.gz')
LOG = open(os.path.join(BASE, '..', '01_log', 'dl_chip_mt.log'), 'a', encoding='utf-8')

NTH = 16
TOTAL = 4167276418  # 已知 Content-Length
CHUNK = (TOTAL + NTH - 1) // NTH

lock = threading.Lock()
progress = {'done': 0}

def log(s):
    with lock:
        LOG.write(s + '\n'); LOG.flush()

def dl(i):
    st = i * CHUNK
    en = min(TOTAL - 1, st + CHUNK - 1)
    part = os.path.join(TMP, 'part%02d' % i)
    want = en - st + 1
    if os.path.exists(part):
        cur = os.path.getsize(part)
        if cur >= want:
            with lock:
                progress['done'] += want
            log('part%02d SKIP (%.0f MB)' % (i, cur / 1e6))
            return
        st = st + cur  # 续传
    t0 = time.time()
    got = os.path.getsize(part) if os.path.exists(part) else 0
    try:
        r = requests.get(U, headers=dict(UA, Range='bytes=%d-%d' % (st, en)),
                         timeout=600, verify=False, stream=True)
        with open(part, 'ab') as f:
            for c in r.iter_content(1 << 20):
                if not c:
                    continue
                f.write(c)
                got += len(c)
                with lock:
                    progress['done'] += len(c)
        el = time.time() - t0
        log('part%02d OK %.0f MB in %.0fs (%.2f MB/s)' % (i, got / 1e6, el, got / el / 1e6))
    except Exception as e:
        log('part%02d ERR %s (got %.0f MB)' % (i, str(e)[:90], got / 1e6))

def monitor(ths):
    while any(t.is_alive() for t in ths):
        time.sleep(30)
        with lock:
            d = progress['done']
        log('PROGRESS %.0f / %.0f MB (%.1f%%)' % (d / 1e6, TOTAL / 1e6, 100.0 * d / TOTAL))

log('=== start 16-thread download %s ===' % time.strftime('%H:%M:%S'))
ths = [threading.Thread(target=dl, args=(i,)) for i in range(NTH)]
for t in ths:
    t.start()
monitor(ths)
for t in ths:
    t.join()
log('all threads done')

# 校验并合并
sizes = []
for i in range(NTH):
    p = os.path.join(TMP, 'part%02d' % i)
    sizes.append(os.path.getsize(p) if os.path.exists(p) else 0)
log('part sizes: %s' % sizes)
if sum(sizes) == TOTAL:
    log('size OK, merging...')
    with open(FINAL, 'wb') as out:
        for i in range(NTH):
            p = os.path.join(TMP, 'part%02d' % i)
            with open(p, 'rb') as f:
                while True:
                    b = f.read(1 << 23)
                    if not b:
                        break
                    out.write(b)
    log('merged -> %s (%.0f bytes)' % (FINAL, os.path.getsize(FINAL)))
else:
    log('SIZE MISMATCH sum=%d want=%d' % (sum(sizes), TOTAL))
log('DONE')
