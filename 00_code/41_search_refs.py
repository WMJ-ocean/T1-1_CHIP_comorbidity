#!/usr/bin/env python
# 41_search_refs.py —— 按标题精确检索 PubMed，取回正确 PMID 与题录
import json, os, time, urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "03_meta", "refs_final.tsv")

QUERIES = [
 "Inherited causes of clonal haematopoiesis in 97,691 whole genomes",
 "Inherited myeloproliferative neoplasm risk affects haematopoietic stem cells",
 "A common JAK2 haplotype confers susceptibility to the myeloproliferative neoplasms",
 "Chromosomal alterations among age-related haematopoietic clones in Japan",
 "The MR-Base platform supports systematic causal inference across the human phenome",
 "Classification and personalized prognosis in myeloproliferative neoplasms",
 "The evolutionary dynamics and fitness landscape of clonal hematopoiesis",
 "The 5th edition of the World Health Organization Classification of Haematolymphoid Tumours: Myeloid and Histiocytic/Dendritic Neoplasms",
 "Clonal Hematopoiesis and Risk of Atherosclerotic Cardiovascular Disease",
 "Prediction of acute myeloid leukaemia risk in healthy elderly individuals",
 "Leukemia-Associated Somatic Mutations Drive Distinct Patterns of Age-Related Clonal Hemopoiesis",
 "Germ line variants predispose to both JAK2 V617F clonal hematopoiesis and myeloproliferative neoplasms",
 "Distinction of lymphoid and myeloid clonal hematopoiesis",
 "Clonal hematopoiesis, with and without candidate driver mutations, is common in the elderly",
 "Clonal Hematopoiesis and Evolution to Hematopoietic Malignancies",
 "Genome-wide analyses of 200,453 individuals yield new insights into the causes and consequences of clonal mosaicism",
 "Computationally efficient whole-genome regression for quantitative and binary traits",
 "Orienting the causal relationship between imprecisely measured traits using GWAS summary data",
 "Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate gene expression",
 "Mendelian randomization",
 "Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization",
 "Interpreting findings from Mendelian randomization using the MR-Egger method",
 "The UK Biobank resource with deep phenotyping and genomic data",
 "Estimating the proportion of variation in susceptibility to schizophrenia explained by common SNPs",
 "Clinical implications of clonal hematopoiesis",
 "Genome-wide association study identifies susceptibility loci for myeloproliferative neoplasms",
 "Aplastic anemia",
 "Sensitive detection of somatic point mutations and improved germline variant calling",
]

def esearch(q):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=3&term="
           + urllib.parse.quote('"%s"[Title]' % q))
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["esearchresult"].get("idlist", [])

def esummary(ids):
    if not ids:
        return []
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
           + ",".join(ids))
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["result"]

rows = []
for q in QUERIES:
    ids = esearch(q)
    time.sleep(0.4)
    if not ids:
        print("[MISS]", q[:60])
        continue
    res = esummary(ids[:1])
    d = res.get(ids[0], {})
    au = d.get("authors") or []
    first = au[0]["name"] if au else ""
    rows.append({
        "pmid": ids[0], "first_author": first, "title": d.get("title", ""),
        "journal": d.get("fulljournalname", d.get("source", "")),
        "date": d.get("pubdate", ""), "vol": d.get("volume", ""),
        "pages": d.get("pages", ""),
        "doi": next((x["value"] for x in d.get("articleids", []) if x["idtype"] == "doi"), ""),
    })
    print(f"{ids[0]:>9}  {first[:16]:16s} {d.get('source','')[:32]:32s} {d.get('pubdate','')[:12]:12s} {d.get('title','')[:60]}")
    time.sleep(0.4)

with open(OUT, "w", encoding="utf-8") as f:
    f.write("pmid\tfirst_author\ttitle\tjournal\tdate\tvol\tpages\tdoi\n")
    for r in rows:
        f.write("\t".join([r["pmid"], r["first_author"], r["title"], r["journal"],
                           r["date"], r["vol"], r["pages"], r["doi"]]) + "\n")
print(f"\n[OK] {len(rows)} -> {OUT}")
