#!/usr/bin/env python
import json, urllib.parse, urllib.request, time, os
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "03_meta", "refs_extra.tsv")
Q = [
 ("Bick2020Cell",    'AUTH:"Bick AG" AND JOURNAL:"Cell" AND PUB_YEAR:2020'),
 ("Hemani2018eLife", 'AUTH:"Hemani G" AND JOURNAL:"eLife" AND PUB_YEAR:2018'),
 ("STROBEMR2021",    'TITLE:"Strengthening the reporting of observational studies in epidemiology using Mendelian randomization" AND JOURNAL:"JAMA"'),
 ("Steensma2018",    'AUTH:"Steensma DP" AND PUB_YEAR:2018 AND TITLE:"clonal hematopoiesis"'),
 ("Abelson2018Nat",  'TITLE:"Prediction of acute myeloid leukaemia risk in healthy elderly individuals"'),
 ("Terao2020Nat",    'AUTH:"Terao C" AND TITLE:"Chromosomal alterations"'),
 ("Kar2022NatGen",   'TITLE:"Genome-wide analyses" AND TITLE:"clonal mosaicism"'),
 ("Bick2020alt",     'TITLE:"Inherited causes of clonal haematopoiesis in 97,691 whole genomes"'),
 ("Weinstock2023",   'AUTH:"Weinstock JS" AND TITLE:"clonal hematopoiesis"'),
]
rows = []
for tag, q in Q:
    try:
        url = API + "?query=" + urllib.parse.quote(q) + "&format=json&pageSize=5&resultType=core"
        r = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        js = json.load(urllib.request.urlopen(r, timeout=60))
        res = js.get("resultList", {}).get("result", [])
    except Exception as e:
        print("[ERR]", tag, str(e)[:60]); continue
    if not res:
        print("[MISS]", tag); continue
    for x in res[:3]:
        print(f"{tag:16s} PMID={x.get('pmid',''):>9} {x.get('pubYear','')} "
              f"{x.get('journalTitle','')[:28]:28s} {(x.get('title','') or '')[:66]}")
        rows.append({"tag": tag, "pmid": x.get("pmid", ""), "doi": x.get("doi", ""),
                     "authors": x.get("authorString", ""),
                     "title": (x.get("title", "") or "").rstrip("."),
                     "journal": x.get("journalTitle", ""), "year": x.get("pubYear", ""),
                     "vol": x.get("journalInfo", {}).get("volume", ""),
                     "pages": x.get("pageInfo", "")})
    time.sleep(0.3)
with open(OUT, "w", encoding="utf-8") as f:
    f.write("tag\tpmid\tdoi\tauthors\ttitle\tjournal\tyear\tvol\tpages\n")
    for r in rows:
        f.write("\t".join([r["tag"], r["pmid"], r["doi"], r["authors"], r["title"],
                           r["journal"], r["year"], r["vol"], r["pages"]]) + "\n")
print("\n[OK]", OUT)
