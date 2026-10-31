#!/usr/bin/env python
# 42_vep_annot.py —— Ensembl VEP REST：为 CHIP 胚系 IV 做功能后果注释
#   同时用 Open Targets Platform v4 GraphQL 取 V2G 基因与全血 eQTL
# 输出: 06_targets/iv_functional_annot.tsv
import json, os, time
import requests

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IV   = os.path.join(BASE, "05_iv", "CHIP_germline_iv_annotated.tsv")
OUT  = os.path.join(BASE, "06_targets", "iv_functional_annot.tsv")
OT   = "https://api.platform.opentargets.org/api/v4/graphql"
VEP  = "https://rest.ensembl.org/vep/human/id"

def rows(path):
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= len(hdr):
                yield dict(zip(hdr, p))

def vep_annot(snps):
    sess = requests.Session()
    sess.headers.update({"User-Agent": "Mozilla/5.0 research",
                         "Content-Type": "application/json",
                         "Accept": "application/json"})
    out = {}
    for i in range(0, len(snps), 200):
        chunk = snps[i:i+200]
        for attempt in range(3):
            try:
                r = sess.post(VEP, data=json.dumps({"ids": chunk}), timeout=120)
                if r.status_code == 200:
                    for it in r.json():
                        cs = []
                        genes = set()
                        for t in (it.get("transcript_consequences") or []):
                            if t.get("gene_symbol"):
                                genes.add(t["gene_symbol"])
                            for c in (t.get("consequence_terms") or []):
                                cs.append(c)
                        # 优先级最高的后果
                        order = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant",
                                 "stop_gained","frameshift_variant","stop_lost","start_lost",
                                 "missense_variant","protein_altering_variant","inframe_insertion",
                                 "inframe_deletion","splice_region_variant","5_prime_UTR_variant",
                                 "3_prime_UTR_variant","coding_sequence_variant","synonymous_variant",
                                 "regulatory_region_variant","TF_binding_site_variant","intron_variant",
                                 "non_coding_transcript_exon_variant","upstream_gene_variant",
                                 "downstream_gene_variant","intergenic_variant"]
                        top = next((o for o in order if o in cs), (cs[0] if cs else ""))
                        out[it.get("id", "")] = {
                            "consequence": top,
                            "vep_genes": ";".join(sorted(genes))[:80],
                            "impact": it.get("most_severe_consequence", ""),
                            "regulatory": "yes" if any("regulatory" in c or "TF_binding" in c for c in cs) else "no",
                        }
                    break
            except Exception as e:
                print("VEP retry:", str(e)[:80]); time.sleep(3)
        time.sleep(0.5)
    return out

def ot_annot(ivs):
    sess = requests.Session()
    sess.headers.update({"User-Agent": "Mozilla/5.0 research"})
    Q = """
    query v($id: String!) {
      variant(variantId: $id) {
        id
        rsIds
        variantEffect {
          method
          score
          assessment
        }
        genome { chromosome position }
      }
    }
    """
    out = {}
    for r in ivs:
        vid = f"{r['chr']}_{r['pos']}_{r['other_allele']}_{r['effect_allele']}"
        try:
            resp = sess.post(OT, json={"query": Q, "variables": {"id": vid}}, timeout=60)
            if resp.status_code == 200:
                js = resp.json()
                d = (js.get("data") or {}).get("variant")
                if d:
                    ve = d.get("variantEffect") or []
                    sc = max([v.get("score") or 0 for v in ve], default=0)
                    out[r["SNP"]] = {"ot_variant": d.get("id", ""), "ot_max_score": sc}
        except Exception as e:
            print("OT err:", str(e)[:60])
        time.sleep(0.2)
    return out

def main():
    ivs = list(rows(IV))
    snps = [r["SNP"] for r in ivs]
    print("[i] VEP annotation for", len(snps), "variants ...", flush=True)
    vep = vep_annot(snps)
    print("[i] Open Targets V2G ...", flush=True)
    ot = ot_annot(ivs)
    hdr = ["SNP", "chr", "pos", "nearest_gene_ensembl", "vep_consequence", "vep_most_severe",
           "vep_genes", "regulatory", "ot_variant", "ot_v2g_score"]
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for r in ivs:
            v = vep.get(r["SNP"], {})
            o = ot.get(r["SNP"], {})
            f.write("\t".join([r["SNP"], r["chr"], r["pos"], r.get("nearest_gene", ""),
                               v.get("consequence", ""), v.get("impact", ""),
                               v.get("vep_genes", ""), v.get("regulatory", ""),
                               o.get("ot_variant", ""), str(o.get("ot_v2g_score", ""))]) + "\n")
    print(f"[OK] -> {OUT}")
    for r in ivs:
        v = vep.get(r["SNP"], {})
        print(f"  {r['SNP']:12s} {r.get('nearest_gene',''):16s} {v.get('consequence',''):35s} {v.get('vep_genes','')[:40]}")

if __name__ == "__main__":
    main()
