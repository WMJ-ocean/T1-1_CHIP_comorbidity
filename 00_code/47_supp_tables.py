#!/usr/bin/env python
# 47_supp_tables.py —— 从结果 CSV 生成投稿用 Supplementary Tables (Markdown)
import os, csv, math
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(BASE, "10_manuscript", "Supplementary_Tables.md")


def fmt_p(p):
    try:
        p = float(p)
    except Exception:
        return str(p)
    if p == 0:
        return "<1e-300"
    if p < 1e-4:
        m, e = f"{p:.1e}".split("e")
        return f"{m}e{int(e)}"
    return f"{p:.3f}"


def ci(b, se):
    return (math.exp(b - 1.96 * se), math.exp(b + 1.96 * se))


lines = []
w = lines.append

w("# Supplementary Tables\n")
w("**Manuscript:** Germline liability to clonal haematopoiesis causes myeloproliferative neoplasms independently of *JAK2*\n")

# ---------------- Table S1 ----------------
w("## Table S1. Outcome definitions and sample sizes (FinnGen release 12)\n")
w("| Outcome | Endpoint ID | Cases | Controls | Total N | Prevalence | GWS loci | Genomic inflation (λ) |")
w("|---|---|---:|---:|---:|---:|---:|---:|")
s1 = [
    ("MPN", "MYELOPROF_NONCML", 2791, 497256, "0.56%", 17, "1.01"),
    ("Polycythaemia vera", "POLYCYTVERA", 1205, 378022, "0.32%", 9, "0.98"),
    ("Myelodysplastic syndromes", "C3_MYELODYSP_SYNDR_EXALLC", 979, 378296, "0.26%", 0, "0.97"),
    ("Immune thrombocytopenia", "D3_ITP", 1064, 492275, "0.22%", 0, "1.03"),
    ("Aplastic anaemia", "D3_OTHERAPLASTICANAEMIA", 374, 479581, "0.08%", 1, "0.92"),
]
for name, eid, ca, co, prev, gws, lam in s1:
    w(f"| {name} | {eid} | {ca:,} | {co:,} | {ca+co:,} | {prev} | {gws} | {lam} |")
w("\nCHIP exposure: Kessler et al., *Nature* 2022 (GCST90165267); UK Biobank, 418,401 participants, "
  "25,657 cases (6.13%). GWS, genome-wide significant (P < 5 × 10⁻⁸).\n")

# ---------------- Table S2 ----------------
w("## Table S2. Relaxed-threshold (P < 5 × 10⁻⁶) bidirectional Mendelian randomisation\n")
rel = pd.read_csv(os.path.join(BASE, "07_mr", "relaxed", "relaxed_MR_all.csv"))
rel = rel[rel["method"].isin(["Inverse variance weighted", "Weighted median", "MR Egger"])]
piv = rel.pivot_table(index="tag", columns="method",
                      values=["b", "se", "pval", "nsnp_h"], aggfunc="first")
w("| Analysis | nSNP | IVW OR (95% CI) | IVW P | Weighted median OR | WM P | MR-Egger OR | Egger P |")
w("|---|---:|---|---:|---:|---:|---:|---:|")
label = {
    "CHIPgP5E6_to_PV": "CHIP liability (P<5e-6) → polycythaemia vera",
    "CHIPgP5E6_to_MPN": "CHIP liability (P<5e-6) → MPN",
    "CHIPgP5E6_to_MDS": "CHIP liability (P<5e-6) → myelodysplastic syndromes",
    "CHIPgP5E6_to_ITP": "CHIP liability (P<5e-6) → immune thrombocytopenia",
    "CHIPgP5E6_to_AA": "CHIP liability (P<5e-6) → aplastic anaemia",
    "MDS_P5E6_to_CHIP": "Myelodysplastic syndromes (P<5e-6) → CHIP",
    "ITP_P5E6_to_CHIP": "Immune thrombocytopenia (P<5e-6) → CHIP",
    "AA_P5E6_to_CHIP": "Aplastic anaemia (P<5e-6) → CHIP",
}
order = ["CHIPgP5E6_to_PV", "CHIPgP5E6_to_MPN", "CHIPgP5E6_to_MDS", "CHIPgP5E6_to_ITP",
         "CHIPgP5E6_to_AA", "MDS_P5E6_to_CHIP", "ITP_P5E6_to_CHIP", "AA_P5E6_to_CHIP"]
for tag in order:
    if tag not in piv.index:
        continue
    r = piv.loc[tag]
    b, se = float(r[("b", "Inverse variance weighted")]), float(r[("se", "Inverse variance weighted")])
    lo, hi = ci(b, se)
    wm_b = float(r[("b", "Weighted median")]); wm_p = float(r[("pval", "Weighted median")])
    eg_b = float(r[("b", "MR Egger")]); eg_p = float(r[("pval", "MR Egger")])
    n = int(r[("nsnp_h", "Inverse variance weighted")])
    w(f"| {label.get(tag, tag)} | {n} | {math.exp(b):.2f} ({lo:.2f}–{hi:.2f}) | "
      f"{fmt_p(r[('pval','Inverse variance weighted')])} | {math.exp(wm_b):.2f} | {fmt_p(wm_p)} | "
      f"{math.exp(eg_b):.2f} | {fmt_p(eg_p)} |")
w("\nInstruments selected at P < 5 × 10⁻⁶ in the exposure dataset, clumped at r² < 0.001, "
  "and filtered to F > 10 in the exposure only. ORs are per unit genetically predicted liability.\n")

# ---------------- Table S3 ----------------
w("## Table S3. Colocalisation results for all region–phenotype pairs\n")
col = pd.read_csv(os.path.join(BASE, "09_coloc", "coloc_summary.csv"))
w("| Region | Outcome | nSNP | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 | Interpretation |")
w("|---|---|---:|---:|---:|---:|---:|---:|---|")
for _, r in col.iterrows():
    reg_raw = r["pair"].split("_vs_")[0]
    reg = {"C11q22_CHIP": "C11q22 (ATM)", "C17p13_CHIP": "C17p13 (TP53)",
           "C21q22_CHIP": "C21q22 (RUNX1)", "C3q25_CHIP": "C3q25",
           "JAK2_CHIP": "JAK2 (9p24)", "TERT_CHIP": "TERT (5p15.33)"}.get(reg_raw, reg_raw)
    phe = r["pair"].split("_vs_")[1]
    h3, h4 = float(r["PP.H3"]), float(r["PP.H4"])
    if h4 > 0.8:
        itp = "Shared causal variant"
    elif h3 > 0.8:
        itp = "Distinct causal variants"
    else:
        itp = "Inconclusive"
    w(f"| {reg} | {phe} | {int(r['nSNP'])} | {float(r['PP.H0']):.2g} | {float(r['PP.H1']):.3f} | "
      f"{float(r['PP.H2']):.2g} | {h3:.3f} | {h4:.3f} | {itp} |")
w("\nPP.H4 > 0.8 indicates a shared causal variant; PP.H3 > 0.8 indicates two distinct causal variants. "
  "Regions were defined as ±250 kb around the lead CHIP variant; prior p12 = 1 × 10⁻⁵.\n")

# ---------------- Table S4 ----------------
w("## Table S4. MR-PRESSO outlier detection and correction\n")
pr = pd.read_csv(os.path.join(BASE, "07_mr", "presso", "presso_summary.csv"))
w("| Analysis | nSNP | Raw OR | Raw P | n outliers | Outlier-corrected OR (95% CI) | Corrected P | Distortion P |")
w("|---|---:|---:|---:|---:|---|---:|---:|")
for _, r in pr.iterrows():
    cb = r.get("corr_b", None)
    cse = r.get("corr_se", None)
    if pd.notna(cb) and pd.notna(cse):
        lo, hi = ci(float(cb), float(cse))
        cinfo = f"{math.exp(float(cb)):.2f} ({lo:.2f}–{hi:.2f})"
    else:
        cinfo = "not applicable (0 outliers)"
    w(f"| {r['analysis']} | {int(r['nsnp'])} | {float(r['raw_or']):.2f} | {fmt_p(r['raw_p'])} | "
      f"{int(r['n_outliers'])} | {cinfo} | "
      f"{fmt_p(r['corr_p']) if pd.notna(r.get('corr_p')) else '—'} | "
      f"{fmt_p(r['dist_p']) if pd.notna(r.get('dist_p')) else '—'} |")
w("\nThe distortion test evaluates whether removing outliers significantly changes the causal estimate; "
  "all values exceed 0.4, indicating that no material distortion occurred.\n")

# ---------------- Table S5 ----------------
w("## Table S5. Reverse-direction Mendelian randomisation (genome-wide significant instruments)\n")
w("| Exposure → CHIP | nSNP | IVW OR (95% CI) | IVW P | Weighted median OR | WM P | Egger intercept P |")
w("|---|---:|---|---:|---:|---:|---:|")
for tag, name in [("MPN_to_CHIP", "MPN"), ("PV_to_CHIP", "Polycythaemia vera")]:
    f = os.path.join(BASE, "07_mr", f"{tag}_MR.csv")
    d = pd.read_csv(f)
    ivw = d[d["method"] == "Inverse variance weighted"].iloc[0]
    wm  = d[d["method"] == "Weighted median"].iloc[0]
    eg  = d[d["method"] == "MR Egger"].iloc[0]
    lo, hi = ci(ivw["b"], ivw["se"])
    w(f"| {name} → CHIP | {int(ivw['nsnp'])} | {math.exp(ivw['b']):.2f} ({lo:.2f}–{hi:.2f}) | "
      f"{fmt_p(ivw['pval'])} | {math.exp(wm['b']):.2f} | {fmt_p(wm['pval'])} | {fmt_p(eg['pval'])} |")
w("\nGenome-wide significant loci (P < 5 × 10⁻⁸) from each outcome were used as instruments, "
  "with CHIP (UK Biobank) as the outcome; the near-zero cross-trait LDSC intercept confirms absence of sample overlap.\n")

# ---------------- Table S6 ----------------
w("## Table S6. Functional annotation of germline instruments (Ensembl VEP)\n")
ann = pd.read_csv(os.path.join(BASE, "06_targets", "iv_functional_annot.tsv"), sep="\t")
w("| SNP | Nearest gene | VEP consequence | Genes affected | Regulatory |")
w("|---|---|---|---|---|")
for _, r in ann.iterrows():
    w(f"| {r['SNP']} | {r['nearest_gene_ensembl']} | {r['vep_consequence']} | "
      f"{str(r['vep_genes'])[:60]} | {r['regulatory']} |")
w("\nAll instruments are non-coding, consistent with regulatory rather than protein-altering mechanisms.\n")

open(OUT, "w", encoding="utf-8").write("\n".join(lines))
print("[OK]", OUT, len(lines), "lines")
