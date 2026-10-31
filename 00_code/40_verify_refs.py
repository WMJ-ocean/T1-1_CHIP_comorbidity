#!/usr/bin/env python
# 40_verify_refs.py —— 批量核对候选参考文献 PMID 与题录（NCBI esummary）
import json, os, time
import urllib.request, urllib.parse

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "03_meta", "refs_verified.txt")

PMIDS = """36450978 25426837 25426838 32888430 32541965 19278957 32807987 36653562
30024853 24830394 25642630 29686387 26050253 28581444 30405240 30405241
36646947 35835855 25764628 35732930 27325186 29212778 33986236 31697822
28827402 29686387 31792460 32060390 33139914 30643219 33299181 34594039
31367058 28916732 30699371 30108116 27989341 29058712 31015403 33753929
33149262 31217584 33139914 32214244 28327225 29058712"""
pmids = sorted(set(PMIDS.split()))
print("query", len(pmids), "pmids")

url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="
       + ",".join(pmids))
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
with urllib.request.urlopen(req, timeout=60) as r:
    js = json.load(r)

lines = []
res = js.get("result", {})
for p in res.get("uids", []):
    d = res.get(p, {})
    au = d.get("authors") or []
    first = au[0]["name"] if au else ""
    pubdate = d.get("pubdate", "")
    epub = d.get("epubdate", "")
    lines.append({
        "pmid": p,
        "title": d.get("title", ""),
        "journal": d.get("fulljournalname", d.get("source", "")),
        "date": pubdate,
        "vol": d.get("volume", ""),
        "pages": d.get("pages", ""),
        "first_author": first,
        "doi": next((x["value"] for x in d.get("articleids", []) if x["idtype"] == "doi"), ""),
    })
lines.sort(key=lambda x: x["date"], reverse=True)
with open(OUT, "w", encoding="utf-8") as f:
    for l in lines:
        f.write(f"{l['pmid']}\t{l['first_author']}\t{l['journal']}\t{l['date']}\t"
                f"{l['vol']}:{l['pages']}\t{l['doi']}\t{l['title']}\n")
print(f"[OK] {len(lines)} refs -> {OUT}")
for l in lines:
    print(f"{l['pmid']:>9}  {l['first_author'][:18]:18s} {l['journal'][:35]:35s} {l['date'][:12]:12s} {l['title'][:70]}")
