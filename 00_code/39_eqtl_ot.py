#!/usr/bin/env python
# 39_eqtl_ot.py —— 用 Open Targets Genetics GraphQL 为 CHIP 胚系 IV 做
#   (1) 位点->基因映射(V2G)  (2) 全血/血液组织 eQTL 方向与显著性
# 输出: 06_targets/eqtl_ot.tsv
import json, sys, time, os
import requests

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IV   = os.path.join(BASE, "05_iv", "CHIP_germline_iv_annotated.tsv")
OUT  = os.path.join(BASE, "06_targets", "eqtl_ot.tsv")
URL  = "https://api.genetics.opentargets.org/graphql"

Q = """
query v2g($id: String!) {
  genesForVariant(variantId: $id) {
    gene { id symbol }
    overallScore
    distance
    qtls {
      type
      source
      tissue { name }
      beta
      pval
    }
  }
}
"""

def rows(path):
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            yield dict(zip(hdr, p))

def main():
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0 (research; T1-1 CHIP MR)"})
    out = []
    ivs = list(rows(IV))
    print(f"[i] {len(ivs)} variants to query", flush=True)
    for i, r in enumerate(ivs, 1):
        chrom, pos = r["chr"], r["pos"]
        ea, oa = r["effect_allele"], r["other_allele"]
        vid = f"{chrom}_{pos}_{oa}_{ea}"
        payload = {"query": Q, "variables": {"id": vid}}
        ok = False
        for attempt in range(3):
            try:
                resp = session.post(URL, json=payload, timeout=60)
                if resp.status_code == 200:
                    js = resp.json()
                    if "errors" in js:
                        print(f"[{i}] {r['SNP']} GRAPHQL ERR: {str(js['errors'])[:120]}")
                        break
                    gv = (js.get("data") or {}).get("genesForVariant") or []
                    ok = True
                    break
            except Exception as e:
                time.sleep(2 + attempt * 3)
                gv = []
        if not ok:
            gv = []
        # 选 overallScore 最高的基因
        gv_sorted = sorted(gv, key=lambda g: (g.get("overallScore") or 0), reverse=True)
        top = gv_sorted[0] if gv_sorted else None
        gene = top["gene"]["symbol"] if top else ""
        score = top.get("overallScore") if top else None
        dist = top.get("distance") if top else None
        # 汇总 eQTL：优先全血/血液类组织
        best = None
        blood_hits = []
        for g in gv_sorted[:3]:
            for q in (g.get("qtls") or []):
                if q.get("type") != "eqtl":
                    continue
                tis = ((q.get("tissue") or {}).get("name") or "").lower()
                src = (q.get("source") or "").upper()
                if ("blood" in tis) or ("EQTLGEN" in src) or ("NATHANSON" in src):
                    blood_hits.append({
                        "gene": g["gene"]["symbol"],
                        "source": q.get("source"),
                        "tissue": (q.get("tissue") or {}).get("name"),
                        "beta": q.get("beta"),
                        "pval": q.get("pval"),
                    })
        if blood_hits:
            best = min(blood_hits, key=lambda x: x["pval"] if x["pval"] is not None else 1)
        rec = {
            "SNP": r["SNP"],
            "nearest_gene_ensembl": r.get("nearest_gene", ""),
            "ot_gene": gene,
            "ot_v2g_score": score,
            "ot_distance": dist,
            "eqtl_gene": best["gene"] if best else "",
            "eqtl_source": best["source"] if best else "",
            "eqtl_tissue": best["tissue"] if best else "",
            "eqtl_beta": best["beta"] if best else "",
            "eqtl_pval": best["pval"] if best else "",
            "n_blood_eqtl": len(blood_hits),
        }
        out.append(rec)
        print(f"[{i:2d}] {r['SNP']:12s} V2G={gene:12s} score={score} "
              f"eQTL={rec['eqtl_gene']:10s} beta={rec['eqtl_beta']} p={rec['eqtl_pval']}",
              flush=True)

    hdr = list(out[0].keys())
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for o in out:
            f.write("\t".join(str(o.get(k, "")) for k in hdr) + "\n")
    print(f"\n[OK] written -> {OUT}  ({len(out)} rows)")

if __name__ == "__main__":
    main()
