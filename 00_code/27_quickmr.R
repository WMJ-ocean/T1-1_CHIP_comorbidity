# 27_quickmr.R —— 敏感性分析批量 MR（跳过 PRESSO，秒级完成）
# 一次算完所有 IV × 表型组合，输出单一汇总表
suppressMessages({ library(TwoSampleMR); library(dplyr) })

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
SENS <- file.path(BASE, "05_iv", "sens")
TGT  <- file.path(BASE, "06_targets")
OUTD <- file.path(BASE, "07_mr", "sens")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) read.delim(p, stringsAsFactors = FALSE)

# 列名归一化：兼容 IV 格式(SNP/effect_allele/other_allele/se/eaf)
# 与 targets 格式(chrom/pos/ref/alt/rsid/pval/beta/sebeta/af, 其中 alt=效应等位)
norm <- function(d) {
  nm <- names(d)
  if (!"SNP" %in% nm) {
    if ("rsid" %in% nm) d$SNP <- d$rsid else d$SNP <- d[[1]]
  }
  if ("sebeta" %in% nm && !"se" %in% nm) d$se <- d$sebeta
  if ("af" %in% nm && !"eaf" %in% nm) d$eaf <- d$af
  if ("alt" %in% nm && !"effect_allele" %in% nm) d$effect_allele <- d$alt
  if ("ref" %in% nm && !"other_allele" %in% nm) d$other_allele <- d$ref
  d
}

PHE <- c("MPN", "PV", "MDS", "AA", "ITP")

jobs <- list()
for (s in c("noTERT", "noMHC")) {
  for (p in PHE) {
    jobs[[length(jobs) + 1]] <- list(
      tag = sprintf("CHIPg%s_to_%s", s, p),
      iv  = file.path(SENS, sprintf("CHIP_germline_%s.tsv", s)),
      out = file.path(TGT, sprintf("%s_targets.tsv", p)))
  }
}
for (s in c("noJAK2", "noJAK2TERT")) {
  for (p in c("MPN", "PV")) {
    jobs[[length(jobs) + 1]] <- list(
      tag = sprintf("%s%s_to_CHIP", p, s),
      iv  = file.path(SENS, sprintf("%s_%s.tsv", p, s)),
      out = file.path(TGT, "CHIP_targets.tsv"))
  }
}
# 基线重算（统一口径，便于对比）
for (p in PHE) {
  jobs[[length(jobs) + 1]] <- list(
    tag = sprintf("CHIPg_to_%s_BASE", p),
    iv  = file.path(BASE, "05_iv", "CHIP_germline_iv.tsv"),
    out = file.path(TGT, sprintf("%s_targets.tsv", p)))
}
for (p in c("MPN", "PV")) {
  jobs[[length(jobs) + 1]] <- list(
    tag = sprintf("%s_to_CHIP_BASE", p),
    iv  = file.path(BASE, "05_iv",
                    sprintf("finngen_R12_%s_iv_gws.tsv",
                            ifelse(p == "MPN", "MYELOPROF_NONCML_MPN", "POLYCYTVERA_PV"))),
    out = file.path(TGT, "CHIP_targets.tsv"))
}

res <- list()
for (j in jobs) {
  if (!file.exists(j$iv) || !file.exists(j$out)) {
    cat(sprintf("%-28s SKIP (missing input)\n", j$tag)); next
  }
  e <- try({
    iv  <- norm(rd(j$iv))
    iv$Phenotype <- "exposure"
    o   <- norm(rd(j$out))
    o$Phenotype <- "outcome"
    fmt <- function(d, ty, lab) {
      raw <- data.frame(SNP = d$SNP, beta = as.numeric(d$beta), se = as.numeric(d$se),
                        pval = as.numeric(d$pval),
                        ea = d$effect_allele, oa = d$other_allele,
                        eaf = as.numeric(d$eaf), stringsAsFactors = FALSE)
      raw <- raw[is.finite(raw$beta) & is.finite(raw$se) & raw$se > 0, ]
      fd <- format_data(raw, type = ty, snp_col = "SNP",
                        beta_col = "beta", se_col = "se", pval_col = "pval",
                        effect_allele_col = "ea", other_allele_col = "oa",
                        eaf_col = "eaf")
      if (ty == "exposure") { fd$id.exposure <- lab; fd$exposure <- lab }
      else                  { fd$id.outcome  <- lab; fd$outcome  <- lab }
      fd[!duplicated(fd$SNP), ]
    }
    h <- harmonise_data(fmt(iv, "exposure", j$tag), fmt(o, "outcome", j$tag), action = 2)
    h <- h[h$mr_keep, ]
    if (nrow(h) < 3) stop("too few SNPs")
    m <- mr(h, method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode"))
    q <- mr_heterogeneity(h)
    pl <- mr_pleiotropy_test(h)
    m$tag <- j$tag; m$nsnp_h <- nrow(h)
    m$Q_pval <- ifelse(m$method == "Inverse variance weighted",
                       q$Q_pval[q$method == "Inverse variance weighted"][1], NA)
    m$egger_p <- pl$pval[1]
    m
  }, silent = TRUE)
  if (inherits(e, "try-error")) {
    cat(sprintf("%-28s FAIL: %s\n", j$tag, as.character(e)))
  } else {
    res[[j$tag]] <- e
    ivw <- e[e$method == "Inverse variance weighted", ]
    cat(sprintf("%-28s n=%2d  IVW OR=%.2f (%.2f-%.2f) P=%.2g  EggerP=%.2f\n",
                j$tag, ivw$nsnp_h, exp(ivw$b), exp(ivw$b - 1.96 * ivw$se),
                exp(ivw$b + 1.96 * ivw$se), ivw$pval, ivw$egger_p))
  }
}

if (length(res) > 0) {
  all <- do.call(rbind, res)
  rownames(all) <- NULL
  write.csv(all, file.path(OUTD, "sens_MR_all.csv"), row.names = FALSE)
  cat("\nwritten ->", file.path(OUTD, "sens_MR_all.csv"), "\n")
}
