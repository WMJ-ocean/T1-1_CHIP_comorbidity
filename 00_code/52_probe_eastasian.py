# 52_probe_eastasian.py — 东亚 GWAS / 公开数据源可达性实测
import requests, json, urllib3
urllib3.disable_warnings()
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research"}

def get(tag, url, **kw):
    try:
        r = requests.get(url, headers=UA, timeout=40, **kw)
        print(f"[{r.status_code}] {tag}\n         {url}")
        if r.status_code == 200:
            t = r.text[:220].replace("\n", " ")
            print(f"         -> {t}")
        return r
    except Exception as e:
        print(f"[ERR ] {tag}\n         {url}\n         -> {type(e).__name__}: {str(e)[:160]}")
        return None

print("========== 1. IEU OpenGWAS（自建 eQTL + 多祖先 GWAS 仓库） ==========")
get("status(verify=False)", "https://gwas-api.mrcieu.ac.uk/status", verify=False)
get("list datasets named eqtl", "https://gwas-api.mrcieu.ac.uk/list?start=0&size=3", verify=False)

print("\n========== 2. Pan-UKBB（含东亚祖源分层结果） ==========")
for u in [
    "https://pan.ukbb.broadinstitute.org/phenotypes",
    "https://pan.ukbb.broadinstitute.org/api/root/info/summary",
    "https://pan-ukbb-us-east-1.s3.amazonaws.com/?list-type=2&max-keys=5",
]:
    get("pan-ukbb", u)

print("\n========== 3. GWAS Atlas（中科院/国家基因组科学数据中心，含中国人群队列） ==========")
get("gwas-atlas traits",
    "https://bigd.big.ac.cn/gwas/api/v1/trait?page=1&pageSize=3")
get("gwas-atlas search MPN",
    "https://bigd.big.ac.cn/gwas/search?searchTerm=polycythemia")

print("\n========== 4. BioBank Japan PheWeb（换用坐标式 variant 端点） ==========")
# rs7705526 -> chr5:1285859 (hg38 约 1285744?)
get("pheweb.jp variant by pos", "https://pheweb.jp/api/variant/5-1285859-A-C")
get("pheweb.jp root", "https://pheweb.jp/")
get("pheweb.jp phenotypes.json", "https://pheweb.jp/api/phenotypes.json")

print("\n========== 5. GWAS Catalog v2 REST（现役接口） ==========")
get("catalog v2 search" ,
    "https://www.ebi.ac.uk/gwas/rest/api/v2/studies/search/findByEfoTrait?efoTrait=myeloproliferative%20neoplasm",
    headers={})
print("DONE")
