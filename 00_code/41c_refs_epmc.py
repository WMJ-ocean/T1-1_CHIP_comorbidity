#!/usr/bin/env python
# 41c_refs_epmc.py —— Europe PMC 精确检索题录（TITLE 字段）
import json, os, time, urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "03_meta", "refs_epmc.tsv")
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

Q = [
 ("Bick2020",      'TITLE:"Inherited causes of clonal haematopoiesis"'),
 ("Bao2020",       'TITLE:"Inherited myeloproliferative neoplasm risk"'),
 ("Olcaydu2009",   'TITLE:"JAK2 haplotype" AND TITLE:"myeloproliferative"'),
 ("Hinds2016",     'TITLE:"Germ line variants predispose to both JAK2"'),
 ("Hemani2018",    'TITLE:"MR-Base"'),
 ("Hemani2017",    'TITLE:"Orienting the causal relationship between imprecisely measured traits"'),
 ("Grinfeld2018",  'TITLE:"Classification and personalized prognosis in myeloproliferative neoplasms"'),
 ("Khoury2022",    'TITLE:"The 5th edition of the World Health Organization Classification of Haematolymphoid Tumours"'),
 ("Jaiswal2017",   'TITLE:"Clonal Hematopoiesis and Risk of Atherosclerotic Cardiovascular Disease"'),
 ("Abelson2018",   'TITLE:"Prediction of acute myeloid leukaemia risk in healthy elderly individuals"'),
 ("McKerrell2015", 'TITLE:"Leukemia-Associated Somatic Mutations Drive Distinct Patterns of Age-Related Clonal Hemopoiesis"'),
 ("Niroula2021",   'TITLE:"Distinction of lymphoid and myeloid clonal hematopoiesis"'),
 ("Zink2017",      'TITLE:"Clonal hematopoiesis, with and without candidate driver mutations"'),
 ("Bowman2018",    'TITLE:"Clonal Hematopoiesis and Evolution to Hematopoietic Malignancies"'),
 ("Mbatchou2021",  'TITLE:"Computationally efficient whole-genome regression for quantitative and binary traits"'),
 ("Vosa2021",      'TITLE:"Large-scale cis- and trans-eQTL analyses"'),
 ("Skrivankova2021", 'TITLE:"Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization"'),
 ("Burgess2017",   'TITLE:"Interpreting findings from Mendelian randomization using the MR-Egger method"'),
 ("Lee2011",       'TITLE:"Estimating the proportion of variation in susceptibility to schizophrenia"'),
 ("Steensma2018",  'TITLE:"Clinical implications of clonal hematopoiesis"'),
 ("Young2018",     'TITLE:"Aplastic Anemia" AND SRC:MED AND PUB_YEAR:2018'),
 ("Steensma2015",  'TITLE:"Clonal hematopoiesis of indeterminate potential"'),
 ("Weinstock2023", 'TITLE:"Clonal hematopoiesis" AND SRC:MED AND PUB_YEAR:2023 AND JOURNAL:"Nature reviews"'),
 ("Bick2020b",     'AUTH:"Bick AG" AND TITLE:"clonal haematopoiesis"'),
 ("Kar2022",       'TITLE:"clonal mosaicism" AND PUB_YEAR:2022'),
]

def search(q):
    url = API + "?query=" + urllib.parse.quote(q) + "&format=json&pageSize=3&resultType=core"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
    with urllib.request.urlopen(req, timeout=60) as r:
        js = json.load(r)
    return js.get("resultList", {}).get("result", [])

rows = []
for tag, q in Q:
    try:
        res = search(q)
    except Exception as e:
        print("[ERR]", tag, str(e)[:70]); continue
    time.sleep(0.3)
    if not res:
        print("[MISS]", tag); continue
    r = res[0]
    rows.append({
        "tag": tag, "pmid": r.get("pmid", ""), "doi": r.get("doi", ""),
        "title": (r.get("title", "") or "").rstrip("."),
        "journal": (r.get("journalInfo", {}).get("journal", {}) or {}).get("title", r.get("journalTitle", "")),
        "journal_abbr": r.get("journalTitle", ""),
        "year": r.get("pubYear", ""), "vol": r.get("journalInfo", {}).get("volume", ""),
        "issue": r.get("journalInfo", {}).get("issue", ""), "pages": r.get("pageInfo", ""),
        "authors": r.get("authorString", ""),
    })
    print(f"{tag:16s} PMID={r.get('pmid',''):>9} {r.get('pubYear','')} {r.get('journalTitle','')[:30]:30s} {(r.get('title','') or '')[:60]}")
    time.sleep(0.3)

with open(OUT, "w", encoding="utf-8") as f:
    f.write("tag\tpmid\tdoi\tauthors\ttitle\tjournal_abbr\tjournal\tyear\tvol\tissue\tpages\n")
    for r in rows:
        f.write("\t".join([r["tag"], r["pmid"], r["doi"], r["authors"], r["title"],
                           r["journal_abbr"], r["journal"], r["year"], r["vol"],
                           r["issue"], r["pages"]]) + "\n")
print(f"\n[OK] {len(rows)} -> {OUT}")
