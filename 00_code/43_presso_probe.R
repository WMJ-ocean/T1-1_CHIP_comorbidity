#!/usr/bin/env Rscript
# 43_presso_probe.R —— 探查 mr_presso 返回对象的真实结构
suppressMessages({ library(TwoSampleMR); library(MRPRESSO) })
BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
norm <- function(d) {
  if (!"SNP" %in% names(d)) d$SNP <- if ("rsid" %in% names(d)) d$rsid else d[[1]]
  if ("sebeta" %in% names(d) && !"se" %in% names(d)) d$se <- d$sebeta
  if ("af" %in% names(d) && !"eaf" %in% names(d)) d$eaf <- d$af
  if ("alt" %in% names(d) && !"effect_allele" %in% names(d)) d$effect_allele <- d$alt
  if ("ref" %in% names(d) && !"other_allele" %in% names(d)) d$other_allele <- d$ref
  d
}
fmt <- function(d, ty, lab) {
  raw <- data.frame(SNP = d$SNP, beta = as.numeric(d$beta), se = as.numeric(d$se),
                    pval = as.numeric(d$pval), ea = d$effect_allele, oa = d$other_allele,
                    eaf = as.numeric(d$eaf), stringsAsFactors = FALSE)
  raw <- raw[is.finite(raw$beta) & is.finite(raw$se) & raw$se > 0, ]
  fd <- format_data(raw, type = ty, snp_col = "SNP", beta_col = "beta", se_col = "se",
                    pval_col = "pval", effect_allele_col = "ea", other_allele_col = "oa",
                    eaf_col = "eaf")
  if (ty == "exposure") { fd$id.exposure <- lab; fd$exposure <- lab }
  else                  { fd$id.outcome <- lab; fd$outcome <- lab }
  fd[!duplicated(fd$SNP), ]
}
iv <- norm(read.delim(file.path(BASE, "05_iv", "CHIP_germline_iv.tsv"), stringsAsFactors = FALSE))
o  <- norm(read.delim(file.path(BASE, "06_targets", "MPN_targets.tsv"), stringsAsFactors = FALSE))
h  <- harmonise_data(fmt(iv, "exposure", "e"), fmt(o, "outcome", "o"), action = 2)
h  <- h[h$mr_keep, ]
d  <- data.frame(SNP = h$SNP, BetaOutcome = h$beta.outcome, SdOutcome = h$se.outcome,
                 BetaExposure = h$beta.exposure, SdExposure = h$se.exposure)
cat("n =", nrow(d), "\n")
po <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure",
                SdOutcome = "SdOutcome", SdExposure = "SdExposure",
                data = d, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                NbDistribution = 1000, seed = 20261029)
cat("\n--- top names ---\n"); print(names(po))
cat("\n--- Main MR results ---\n"); print(po$`Main MR results`)
cat("\n--- MR-PRESSO results names ---\n"); print(names(po$`MR-PRESSO results`))
pr <- po$`MR-PRESSO results`
for (n in names(pr)) {
  cat("\n### ", n, "\n")
  x <- pr[[n]]
  if (is.data.frame(x)) { print(head(x, 3)); cat("  cols:", paste(names(x), collapse=", "), "\n") }
  else { print(str(x, max.level = 2)) }
}
