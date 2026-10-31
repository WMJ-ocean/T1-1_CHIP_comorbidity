# 53_probe_bbj_detail.py — 判定 BioBank Japan 能否用于东亚验证（阳性对照 + 表型清单）
import requests, json, re
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}
S = requests.Session(); S.headers.update(UA)

print("=== 1. 阳性对照：ALP rs9411471 (chr9:136120091 hg19) ===")
try:
    r = S.get("https://pheweb.jp/api/variant/9-136120091-G-A", timeout=40)
    print("status", r.status_code)
    j = r.json()
    if j:
        print(f"  命中 {len(j)} 个表型关联")
        top = sorted(j, key=lambda x: x.get("pval", 1))[:3]
        for t in top:
            print(f"    {t.get('phenocode')}  p={t.get('pval'):.3e}  beta={t.get('beta'):.3f}")
    else:
        print("  返回 null -> 坐标体系不对")
except Exception as e:
    print("ERR", e)

print("\n=== 2. BBJ 完整表型清单 ===")
try:
    r = S.get("https://pheweb.jp/api/phenotypes.json", timeout=60)
    ph = r.json()
    print(f"  共 {len(ph)} 个表型")
    kws = ["polycy", "myelo", "leuk", "lymph", "ane", "aemia", "emia", "thromb", "hema",
           "blood", "platelet", "hemoglobin", "white", "red", "cancer", "malig"]
    def show(kw_list, title):
        hits = [p for p in ph if any(k in (p.get("phenocode","")+" "+p.get("category","")+" "+str(p.get("phenostring",""))).lower() for k in kw_list)]
        print(f"\n  --{title}--  {len(hits)} 个")
        for p in hits[:25]:
            print(f"     {p.get('phenocode'):<14} N={p.get('num_samples'):>7}  Cat={p.get('category','')[:18]:<18} {str(p.get('phenostring',''))[:40]}")
    show(["polycy"], "真性红细胞增多症 PV")
    show(["myelo"], "骨髓增殖性肿瘤 MPN")
    show(["leuk"], "白血病")
    show(["thromb"], "血小板增多症")
    show(["platelet","hemoglobin","white","red"], "血液 DAG。-血细胞性状")
    show(["cancer","malig"], "恶性肿瘤")
except Exception as e:
    print("ERR", type(e).__name__, e)
print("\nDONE")
