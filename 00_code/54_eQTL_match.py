# -*- coding: utf-8 -*-
"""54_eQTL_match.py — 从 eQTLGen 中抽出 21 个胚系 IV 的靶基因（流式，不全量入内存）"""
import os, gzip, csv, time, sys

BASE = os.path.dirname(os.path.abspath(__file__)) + "/.."
IV   = os.path.join(BASE, "05_iv", "CHIP_germline_iv.tsv")
EQTL = os.path.join(BASE, "02_raw", "eqtlgen", "eQTLGen_cis_eQTL_FDR005.txt.gz")
OUT  = os.path.join(BASE, "06_targets", "iv_eqtlgen_targets.tsv")

# 读 IV
ivs = {}
with open(IV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        ivs[row["SNP"]] = row.get("nearest_gene", "")
print(f"读入 {len(ivs)} 个胚系 IV", flush=True)

hit = []
t0 = time.time()
n = 0
with gzip.open(EQTL, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    i_snp, i_gene = idx["SNP"], idx["GeneSymbol"]
    i_p, i_z, i_ens = idx["Pvalue"], idx["Zscore"], idx["Gene"]
    i_nc, i_ns, i_fdr = idx["NrCohorts"], idx["NrSamples"], idx["FDR"]
    for line in f:
        n += 1
        if n % 2_000_000 == 0:
            print(f"  扫描 {n/1e6:.0f}M 行…", flush=True)
        p = line.split("\t")
        s = p[i_snp]
        if s in ivs:
            hit.append({
                "SNP": s, "GeneSymbol": p[i_gene], "Ensembl": p[i_ens],
                "P": p[i_p], "Z": p[i_z], "FDR": p[i_fdr],
                "NrCohorts": p[i_nc], "NrSamples": p[i_ns],
            })
print(f"扫描 {n:,} 行, 命中 {len(hit)} 条 SNP-基因对, {time.time()-t0:.0f}s", flush=True)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(hit[0].keys()), delimiter="\t")
    w.writeheader()
    for r in sorted(hit, key=lambda x: (x["SNP"], float(x["P"]))):
        w.writerow(r)
print(f"写入 {OUT}", flush=True)

# 每个 IV 的最佳靶基因 + 汇总
print("\n=== 每个 IV 的最强 eQTL 靶基因 ===")
print(f"{'SNP':<14}{'top gene':<16}{'P':<14}{'Z':<12}{'#genes'}")
from collections import defaultdict
g = defaultdict(list)
for r in hit:
    g[r["SNP"]].append(r)
for s in ivs:
    if s not in g:
        print(f"{s:<14}{'(eQTLGen 无显著关联)':<16}")
        continue
    rs = sorted(g[s], key=lambda x: abs(float(x["Z"])), reverse=True)
    top = rs[0]
    print(f"{s:<14}{top['GeneSymbol']:<16}{top['P']:<14}{top['Z']:<12}{len(rs)}")

nohit = [s for s in ivs if s not in g]
print(f"\n有 eQTL 的 IV: {len(g)}/{len(ivs)}；无显著 eQTL 的: {len(nohit)}")
if nohit:
    print("无匹配:", ", ".join(nohit))
