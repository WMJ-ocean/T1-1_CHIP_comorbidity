# Germline CHIP liability and myeloid neoplasms

Analysis code for the manuscript:

> **Germline CHIP liability causes myeloproliferative neoplasms independently of *JAK2*:
> Mendelian randomisation and colocalisation**

Mingjing Wang, Yi Chen, Richeng Quan, Lv Yan, Xiaomei Hu

Department of Hematology, Xiyuan Hospital, China Academy of Chinese Medical Sciences, Beijing, China

---

## What this repository contains

| Path | Content |
|---|---|
| `00_code/` | All analysis scripts (93 Python, 17 R, 1 shell), numbered in execution order |
| `05_iv/` | Instrument variable lists for the CHIP exposure |
| `07_mr/` | Mendelian randomisation outputs (per-outcome CSVs, harmonised data, leave-one-out, heterogeneity, pleiotropy) |
| `09_coloc/` | Colocalisation results (`eqtl_coloc_summary.csv`: 21 instruments x 511 gene-outcome pairs) |
| `09_ldsc/` | LD score regression genetic correlation matrix (`rg_matrix.csv`) |

**Not included**: raw summary statistics (17 GB), the source GWAS files, and the manuscript itself.
All input data are public and listed below.

---

## Data sources

| Data | Source | Accession |
|---|---|---|
| CHIP (aggregate) GWAS | Kessler et al. 2022, *Nature* | GCST90165267 |
| CHIP, *DNMT3A*-mutation subtype | Kessler et al. 2022, *Nature* | GCST90165271 |
| MPN, polycythaemia vera, MDS, ITP, aplastic anaemia | FinnGen release 12 | `MYELOPROF_NONCML_MPN`, `POLYCYTVERA_PV`, `MYELODYSPL`, `ITP`, `APLASTANAEMIA` |
| Blood cis-eQTL | eQTLGen consortium (34,672 samples) | https://www.eqtlgen.org |
| Plasma proteins (Olink inflammation panel) | Zhao et al. 2023, *Nat Immunol* | https://doi.org/10.1038/s41590-023-01588-w |
| Monocyte transcriptome | GEO | GSE283711 |
| LD scores | 1000 Genomes European, HapMap3 | https://zenodo.org/records/8182036 |

---

## Analysis pipeline

### 1. Instrument selection (`00_code/16_parse_chip.py` → `05_iv/`)

The CHIP GWAS distributes an **empty `standard_error` column**; the real standard errors are
retained only inside `additional_info` as `REGENIE_SE=`. Effect estimates and standard errors are
recovered by regular-expression parsing of that field. Variants in the *JAK2* 46/1 region
(chr9:4.8–5.4 Mb) were excluded so that the forward estimate cannot be driven by somatic V617F.

### 2. Two-sample MR (`00_code/`)

Inverse-variance weighted estimates with MR-Egger, weighted median and weighted mode as sensitivity
analyses; MR-PRESSO for outlier detection; Steiger filtering for direction. Reverse-direction MR
uses 22 independent MPN risk loci as instruments.

### 3. Colocalisation (`00_code/71_coloc_eqtl.R`)

`coloc.abf` with sdY = 1, N = 418,401 and p12 = 1e-5, run for every instrument × blood cis-eQTL pair.
eQTL z-scores are converted to betas as
`beta = z / sqrt(2f(1-f)(N + z^2))`, using the assessed-allele frequency.

### 4. Protein mediation (`00_code/69_protein_mr.py`, `61f_dl_pqtl_mt.py`)

Two-step MR across 20 Olink inflammation proteins: CHIP → protein (step 1), then cis-instrumented
protein → outcome (step 2), with the mediated proportion computed as a × b / c.

### 5. Subtype-restricted MR (`00_code/79_dnmt3a_mr2.py`)

Repeats the forward analysis using *DNMT3A*-mutation CHIP instruments. The subtype summary
statistics carry neither rsIDs nor beta/SE columns, so variants are matched to FinnGen by GRCh38
position and effect sizes recovered from `REGENIE_BETA`/`REGENIE_SE`. Variants with |beta| > 1.3
are excluded as somatic read-outs of the *DNMT3A* locus itself (the R882 hotspot reaches
beta ≈ +8), then pruned at 10 Mb.

### 6. LD score regression and transcriptome (`00_code/`, `09_ldsc/`)

SNP heritability and cross-trait genetic correlation using 1,290,028 HapMap3 variants and European
1000 Genomes reference scores, 200 block-jackknife standard errors.

---

## Reproducing

Requires R ≥ 4.5 (TwoSampleMR, MendelianRandomization, coloc, ieugwasr, gwasglue, limma, GSVA)
and Python ≥ 3.11 (numpy, pandas, requests).

```bash
# instruments
python 00_code/16_parse_chip.py

# forward and reverse MR
Rscript 00_code/20_mr_forward.R
Rscript 00_code/30_mr_reverse.R

# colocalisation
Rscript 00_code/71_coloc_eqtl.R

# subtype-restricted MR (downloads 3.5 GB; scan and fit run separately)
python 00_code/79_dnmt3a_mr2.py dl
python 00_code/79_dnmt3a_mr2.py scan
python 00_code/79_dnmt3a_mr2.py mr
```

Scripts use absolute paths in the `BASE` variable at the top of each file; adjust to your local
layout before running.

---

## Two technical notes for reusers

1. **The CHIP summary statistics have a silent defect.** `standard_error` is empty for every
   variant; the values exist only in `additional_info`. Analyses that read the standard column
   directly will fail without an error message.
2. **The GWAS Catalog files are multi-member gzip archives.** Decompressing with a single
   `gzip.open()` read returns only the first member (about 900 rows) and silently truncates the
   file. Iterate the decompressor over `unused_data` instead.

---

## Contact

Xiaomei Hu — huxiaomei_2@163.com

## Licence

Code released under the MIT licence. Derived tables are provided for research use.
