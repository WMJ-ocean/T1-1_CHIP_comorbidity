#!/usr/bin/env python
# 41b_search_refs2.py —— 作者+关键词检索 PubMed，取回题录（宽松匹配，人工可选）
import json, os, time, urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "03_meta", "refs_final.tsv")

# (标签, PubMed query)
Q = [
 ("Bick2020_inheritedCHIP",      'Bick AG[Author] AND clonal haematopoiesis AND whole genomes'),
 ("Bao2020_inheritedMPN",        'Bao EL[Author] AND myeloproliferative neoplasm risk'),
 ("Olcaydu2009_JAK2hap",         'Olcaydu D[Author] AND JAK2 haplotype'),
 ("Terao2020_JapanCHIP",         'Terao C[Author] AND chromosomal alterations AND haematopoietic clones'),
 ("Hemani2018_MRBase",           'Hemani G[Author] AND MR-Base'),
 ("Grinfeld2018_MPNclass",       'Grinfeld J[Author] AND myeloproliferative neoplasms AND classification'),
 ("Watson2023_CHIPdynamics",     'Watson CJ[Author] AND clonal hematopoiesis AND fitness'),
 ("Khoury2022_WHO5",             'Khoury JD[Author] AND World Health Organization Classification AND Haematolymphoid'),
 ("Jaiswal2017_CHD",             'Jaiswal S[Author] AND Clonal Hematopoiesis AND Atherosclerotic'),
 ("Abelson2018_AMLrisk",         'Abelson S[Author] AND acute myeloid leukaemia risk'),
 ("McKerrell2015_CHIPdrivers",   'McKerrell T[Author] AND clonal hemopoiesis'),
 ("Hinds2016_JAK2germline",      'Hinds DA[Author] AND JAK2 V617F AND clonal hematopoiesis'),
 ("Niroula2021_lymphoidmyeloid", 'Niroula A[Author] AND clonal hematopoiesis'),
 ("Zink2017_CHIPelderly",        'Zink F[Author] AND clonal hematopoiesis'),
 ("Bowman2018_CHIPreview",       'Bowman RL[Author] AND Clonal Hematopoiesis AND Evolution'),
 ("Kar2022_mosaicism",           'Kar S[Author] AND clonal mosaicism AND genome-wide'),
 ("Mbatchou2021_REGENIE",        'Mbatchou J[Author] AND REGENIE'),
 ("Hemani2017_Steiger",          'Hemani G[Author] AND Orienting the causal relationship'),
 ("Vosa2021_eQTLGen",            'Vosa U[Author] AND eQTL AND polygenic scores'),
 ("Sanderson2022_MRprimer",      'Sanderson E[Author] AND Mendelian randomization AND Nature Reviews'),
 ("Skrivankova2021_STROBEMR",    'Skrivankova VW[Author] AND Mendelian Randomization'),
 ("Burgess2017_Egger",           'Burgess S[Author] AND MR-Egger method'),
 ("Bycroft2018_UKB",             'Bycroft C[Author] AND UK Biobank resource'),
 ("Lee2011_liability",           'Lee SH[Author] AND proportion of variation AND schizophrenia AND SNPs'),
 ("Steensma2018_CHIPclinical",   'Steensma DP[Author] AND clonal hematopoiesis'),
 ("Ulirsch2016_MPNgwas",         'Ulirsch JC[Author] AND myeloproliferative neoplasms AND genome-wide'),
 ("Young2018_AA",                'Young NS[Author] AND Aplastic Anemia AND New England Journal'),
 ("Zhou2022_saige",              'Zhou W[Author] AND whole-genome regression AND binary traits'),
]

def esearch(q):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=5&term="
           + urllib.parse.quote(q))
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["esearchresult"].get("idlist", [])

def esummary(ids):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
           + ",".join(ids))
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["result"]

rows = []
for tag, q in Q:
    try:
        ids = esearch(q)
    except Exception as e:
        print("[ERR]", tag, str(e)[:60]); continue
    time.sleep(0.4)
    if not ids:
        print("[MISS]", tag); continue
    res = esummary(ids[:3])
    # 取第一个有完整标题的
    picked = None
    for pid in ids[:3]:
        d = res.get(pid, {})
        if d.get("title"):
            au = d.get("authors") or []
            picked = {
                "tag": tag, "pmid": pid, "first_author": au[0]["name"] if au else "",
                "title": d.get("title", ""), "journal": d.get("fulljournalname", d.get("source", "")),
                "date": d.get("pubdate", ""), "vol": d.get("volume", ""), "pages": d.get("pages", ""),
                "doi": next((x["value"] for x in d.get("articleids", []) if x["idtype"] == "doi"), ""),
            }
            break
    if picked:
        rows.append(picked)
        print(f"{tag:28s} {picked['pmid']:>9} {picked['first_author'][:14]:14s} "
              f"{picked['journal'][:28]:28s} {picked['date'][:10]:10s} {picked['title'][:55]}")
    time.sleep(0.4)

with open(OUT, "w", encoding="utf-8") as f:
    f.write("tag\tpmid\tfirst_author\ttitle\tjournal\tdate\tvol\tpages\tdoi\n")
    for r in rows:
        f.write("\t".join([r["tag"], r["pmid"], r["first_author"], r["title"], r["journal"],
                           r["date"], r["vol"], r["pages"], r["doi"]]) + "\n")
print(f"\n[OK] {len(rows)} -> {OUT}")
