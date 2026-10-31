#!/usr/bin/env Rscript
# ============================================================
# T1-1 双向 MR: CHIP <-> {MPN, PV, MDS, AA, ITP}
# 用法: Rscript 20_mr.R <exposure_iv.tsv> <outcome_targets.tsv> <tag>
#   exposure_iv.tsv: SNP chr pos effect_allele other_allele beta se pval eaf
#   outcome_targets.tsv: chrom pos ref alt rsid pval beta sebeta af
#   (FinnGen: alt = effect allele, ref = other allele)
# ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("need: exposure_iv.tsv outcome_targets.tsv tag")
exp_f  <- args[1]
out_f  <- args[2]
tag    <- args[3]

OUTD <- file.path(dirname(dirname(exp_f)), "07_mr")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

# ---------- 读暴露 ----------
e <- read.delim(exp_f, stringsAsFactors = FALSE, check.names = FALSE)
e <- e[!is.na(e$beta) & !is.na(e$se) & e$se > 0, ]
exp_dat <- data.frame(
  SNP             = e$SNP,
  beta.exposure   = as.numeric(e$beta),
  se.exposure     = as.numeric(e$se),
  pval.exposure   = as.numeric(e$pval),
  effect_allele.exposure   = e$effect_allele,
  other_allele.exposure    = e$other_allele,
  eaf.exposure    = as.numeric(e$eaf),
  exposure        = tag,
  id.exposure     = tag,
  stringsAsFactors = FALSE
)
exp_dat <- exp_dat[!duplicated(exp_dat$SNP), ]

# ---------- 读结局 ----------
o <- read.delim(out_f, stringsAsFactors = FALSE, check.names = FALSE)
out_dat <- data.frame(
  SNP           = o$rsid,
  beta.outcome  = as.numeric(o$beta),
  se.outcome    = as.numeric(o$sebeta),
  pval.outcome  = as.numeric(o$pval),
  effect_allele.outcome  = o$alt,     # FinnGen: alt 为效应等位
  other_allele.outcome   = o$ref,
  eaf.outcome   = as.numeric(o$af),
  outcome       = basename(out_f),
  id.outcome    = basename(out_f),
  stringsAsFactors = FALSE
)
out_dat <- out_dat[!duplicated(out_dat$SNP), ]

cat(sprintf("exposure IV = %d ; outcome matched rows = %d\n", nrow(exp_dat), nrow(out_dat)))

# ---------- harmonise ----------
dat <- harmonise_data(exp_dat, out_dat, action = 2)
dat <- dat[dat$mr_keep, ]
cat(sprintf("after harmonise (mr_keep) = %d SNPs\n", nrow(dat)))
if (nrow(dat) < 3) {
  cat("WARNING: <3 SNPs after harmonise, MR unreliable\n")
}

# ---------- MR ----------
meth <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
if (nrow(dat) < 4) meth <- c("mr_ivw", "mr_wald_ratio")

res <- tryCatch(mr(dat, method_list = meth), error = function(e) NULL)

if (!is.null(res)) {
  res$OR      <- exp(res$b)
  res$OR_low  <- exp(res$b - 1.96 * res$se)
  res$OR_high <- exp(res$b + 1.96 * res$se)
  res$F_stat  <- (dat$beta.exposure / dat$se.exposure)^2 %>% mean(na.rm = TRUE)
  write.csv(res, file.path(OUTD, sprintf("%s_MR.csv", tag)), row.names = FALSE)

  cat("\n===== MR results =====\n")
  print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_low", "OR_high")])

  # 敏感性
  if (nrow(dat) >= 3) {
    het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    if (!is.null(het)) { write.csv(het, file.path(OUTD, sprintf("%s_heterogeneity.csv", tag)), row.names = FALSE)
      cat("\n-- heterogeneity --\n"); print(het) }
    plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    if (!is.null(plei)) { write.csv(plei, file.path(OUTD, sprintf("%s_pleiotropy.csv", tag)), row.names = FALSE)
      cat("\n-- pleiotropy (Egger intercept) --\n"); print(plei) }
    loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
    if (!is.null(loo)) write.csv(loo, file.path(OUTD, sprintf("%s_leaveoneout.csv", tag)), row.names = FALSE)
    ss <- tryCatch(mr_steiger(dat), error = function(e) NULL)
    if (!is.null(ss)) { write.csv(ss, file.path(OUTD, sprintf("%s_steiger.csv", tag)), row.names = FALSE) }
  }
  # MR-PRESSO
  mp <- tryCatch({
    if (requireNamespace("MRPRESSO", quietly = TRUE) && nrow(dat) >= 4) {
      r <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                               SdOutcome = "se.outcome", SdExposure = "se.exposure",
                               data = dat, NbDistribution = 3000, SignifThreshold = 0.05)
      saveRDS(r, file.path(OUTD, sprintf("%s_presso.rds", tag)))
      cat("\n-- MR-PRESSO global test --\n")
      print(r$`MR-PRESSO results`$`Global Test`)
      r
    } else NULL
  }, error = function(e) { cat("PRESSO skip:", conditionMessage(e), "\n"); NULL })

  # 散点图
  p1 <- tryCatch({
    g <- mr_scatter_plot(res, dat)[[1]] + ggtitle(sprintf("MR: %s", tag)) +
      theme_bw(base_size = 11)
    ggsave(file.path(OUTD, sprintf("%s_scatter.pdf", tag)), g, width = 5.2, height = 4.4)
    g
  }, error = function(e) NULL)
  # 森林图
  tryCatch({
    g <- mr_forest_plot(mr_singlesnp(dat))[[1]] + ggtitle(sprintf("Single-SNP: %s", tag)) +
      theme_bw(base_size = 11)
    ggsave(file.path(OUTD, sprintf("%s_forest.pdf", tag)), g, width = 6, height = max(3, 0.22 * nrow(dat) + 2))
  }, error = function(e) NULL)
  # 留一法
  tryCatch({
    g <- mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]] + ggtitle(sprintf("Leave-one-out: %s", tag)) +
      theme_bw(base_size = 11)
    ggsave(file.path(OUTD, sprintf("%s_loo.pdf", tag)), g, width = 6, height = max(3, 0.22 * nrow(dat) + 2))
  }, error = function(e) NULL)

  write.csv(dat, file.path(OUTD, sprintf("%s_harmonised.csv", tag)), row.names = FALSE)
  cat("\nwritten ->", OUTD, "\n")
} else {
  cat("MR failed\n")
}
