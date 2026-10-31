#!/usr/bin/env Rscript
# 38_presso.R —— MR-PRESSO 离群值校正（全部正向 + 反向组合）
# 输出: 07_mr/presso/presso_summary.csv
suppressMessages({
  library(TwoSampleMR); library(MRPRESSO); library(dplyr)
})

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
IVD  <- file.path(BASE, "05_iv")
TGT  <- file.path(BASE, "06_targets")
OUTD <- file.path(BASE, "07_mr", "presso")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

rd  <- function(p) read.delim(p, stringsAsFactors = FALSE)
norm <- function(d) {
  nm <- names(d)
  if (!"SNP" %in% nm) d$SNP <- if ("rsid" %in% nm) d$rsid else d[[1]]
  if ("sebeta" %in% nm && !"se" %in% nm) d$se <- d$sebeta
  if ("af" %in% nm && !"eaf" %in% nm) d$eaf <- d$af
  if ("alt" %in% nm && !"effect_allele" %in% nm) d$effect_allele <- d$alt
  if ("ref" %in% nm && !"other_allele" %in% nm) d$other_allele <- d$ref
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
  else                  { fd$id.outcome  <- lab; fd$outcome  <- lab }
  fd[!duplicated(fd$SNP), ]
}

PHE <- c("MPN", "PV", "MDS", "AA", "ITP")
jobs <- list()
for (p in PHE) {
  jobs[[length(jobs)+1]] <- list(
    tag = sprintf("CHIPg_to_%s", p),
    iv  = file.path(IVD, "CHIP_germline_iv.tsv"),
    out = file.path(TGT, sprintf("%s_targets.tsv", p)))
}
for (p in c("MPN", "PV")) {
  jobs[[length(jobs)+1]] <- list(
    tag = sprintf("%s_to_CHIP", p),
    iv  = file.path(IVD, sprintf("finngen_R12_%s_iv_gws.tsv",
            ifelse(p=="MPN","MYELOPROF_NONCML_MPN","POLYCYTVERA_PV"))),
    out = file.path(TGT, "CHIP_targets.tsv"))
}

res <- list()
for (j in jobs) {
  cat("\n==== ", j$tag, " ====\n")
  out <- try({
    iv <- norm(rd(j$iv)); o <- norm(rd(j$out))
    h  <- harmonise_data(fmt(iv, "exposure", j$tag), fmt(o, "outcome", j$tag), action = 2)
    h  <- h[h$mr_keep, ]
    if (nrow(h) < 3) stop("too few SNPs")
    d  <- data.frame(
      SNP = h$SNP,
      BetaOutcome  = h$beta.outcome,  SdOutcome  = h$se.outcome,
      BetaExposure = h$beta.exposure, SdExposure = h$se.exposure,
      stringsAsFactors = FALSE)
    d <- d[is.finite(d$BetaOutcome) & is.finite(d$BetaExposure), ]
    po <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure",
                    SdOutcome = "SdOutcome", SdExposure = "SdExposure",
                    data = d, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    NbDistribution = 1000, seed = 20261029)
    mrres <- po$`Main MR results`
    raw   <- mrres[mrres$`MR Analysis` == "Raw", ]
    corr  <- mrres[mrres$`MR Analysis` == "Outlier-corrected", ]
    pr    <- po$`MR-PRESSO results`
    gt    <- pr$`Global Test`
    glob  <- if (is.null(gt$Pvalue)) NA else as.character(gt$Pvalue)[1]
    dt    <- pr$`Distortion Test`
    ol    <- if (is.null(dt$`Outliers Indices`)) integer(0) else dt$`Outliers Indices`
    n_ol  <- length(ol)
    dist_p <- if (is.null(dt$Pvalue)) NA_real_ else suppressWarnings(as.numeric(dt$Pvalue))
    dist_coef <- if (is.null(dt$`Distortion Coefficient`)) NA_real_ else as.numeric(dt$`Distortion Coefficient`)[1]
    data.frame(
      analysis   = j$tag, nsnp = nrow(d),
      raw_b      = as.numeric(raw$`Causal Estimate`),
      raw_se     = as.numeric(raw$Sd),
      raw_p      = as.numeric(raw$`P-value`),
      global_p   = glob,
      n_outliers = n_ol,
      outliers   = paste(d$SNP[ol], collapse = ";"),
      dist_coef  = dist_coef,
      dist_p     = dist_p,
      corr_b     = if (nrow(corr) == 0) NA_real_ else as.numeric(corr$`Causal Estimate`),
      corr_se    = if (nrow(corr) == 0) NA_real_ else as.numeric(corr$Sd),
      corr_p     = if (nrow(corr) == 0) NA_real_ else as.numeric(corr$`P-value`),
      stringsAsFactors = FALSE)
  }, silent = TRUE)
  if (inherits(out, "try-error")) {
    cat("FAIL:", as.character(out), "\n")
    res[[j$tag]] <- data.frame(analysis = j$tag, note = substr(as.character(out), 1, 120))
  } else {
    print(out)
    res[[j$tag]] <- out
  }
}

all <- dplyr::bind_rows(res)
all$raw_or  <- exp(all$raw_b)
all$corr_or <- exp(all$corr_b)
write.csv(all, file.path(OUTD, "presso_summary.csv"), row.names = FALSE)
cat("\nwritten ->", file.path(OUTD, "presso_summary.csv"), "\n")
print(all[, c("analysis","nsnp","raw_or","raw_p","n_outliers","corr_or","corr_p","dist_p")])
