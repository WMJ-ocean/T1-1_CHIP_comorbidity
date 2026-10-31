# -*- coding: utf-8 -*-
"""
31_annotate_ensembl.py —— 用 Ensembl REST 精确注释 IV 位点最近基因
输出 05_iv/CHIP_germline_iv_annotated.tsv（含 gene, distance, F_stat）
"""
import os, sys, time, json
import requests

BASE = r'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
SRC = os.path.join(BASE, '05_iv', 'CHIP_germline_iv.tsv')
DST = os.path.join(BASE, '05_iv', 'CHIP_germline_iv_annotated.tsv')
SESS = requests.Session()
SESS.headers.update({'Content-Type': 'application/json'})

def load(p):
    with open(p, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    return hdr, rows

def nearest(chrom, pos, pad=300000):
    """在 ±300kb 窗口内取所有基因，返回距离最近的那个"""
    u = (f'https://rest.ensembl.org/overlap/region/human/{chrom}:'
         f'{max(1,pos-pad)}-{pos+pad}?feature=gene;content-type=application/json')
    for attempt in range(3):
        try:
            r = SESS.get(u, timeout=45)
            if r.status_code != 200:
                return None, None, None
            gs = r.json()
            if not gs:
                return None, None, None
            best = min(gs, key=lambda g: min(abs(pos - g['start']), abs(pos - g['end']),
                                            0 if g['start'] <= pos <= g['end'] else 1e12))
            dist = 0 if best['start'] <= pos <= best['end'] else min(
                abs(pos - best['start']), abs(pos - best['end']))
            return best.get('external_name') or best.get('id'), dist, (best['start'], best['end'])
        except Exception as e:
            if attempt == 2:
                return None, None, None
            time.sleep(1.5)

hdr, rows = load(SRC)
ci, pi, bi, si = hdr.index('chr'), hdr.index('pos'), hdr.index('beta'), hdr.index('se')
out = []
log = []
for r in rows:
    snp, ch, pos = r[0], r[ci], int(float(r[pi]))
    beta, se = float(r[bi]), float(r[si])
    F = (beta / se) ** 2
    gene, dist, span = nearest(ch, pos)
    out.append(r + [gene or '—', str(dist) if dist is not None else '—', f'{F:.0f}'])
    log.append(f'{snp}\t{ch}:{pos}\t{gene}\tdist={dist}\tF={F:.0f}')
    print(log[-1], flush=True)

with open(DST, 'w', encoding='utf-8') as f:
    f.write('\t'.join(hdr + ['nearest_gene', 'distance_bp', 'F_stat']) + '\n')
    for o in out:
        f.write('\t'.join(o) + '\n')
print('\nwritten ->', DST)
