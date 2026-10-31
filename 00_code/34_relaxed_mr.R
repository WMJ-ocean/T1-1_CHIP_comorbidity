# 34_relaxed_mr.R —— P1: 放宽阈值(p<5e-6)的双向 MR
#   A) 正向: CHIP 胚系 85 IV (p<5e-6) -> 五病   [检验 OR 1.34 等在更多工具下是否稳健]
#   B) 反向: MDS / ITP / AA 各自 p<5e-6 的 IV -> CHIP  [补齐无 GWS 位点的三病双向证据]
suppressMessages({ library(TwoSampleMR); library(dplyr) })

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
SENS <- file.path(BASE, "05_iv", "sens")
IVD  <- file.path(BASE, "05_iv")
TGT  <- file.path(BASE, "06_targets")
RLX  <- file.path(TGT, "relaxed")
OUTD <- file.path(BASE, "07_mr", "relaxed")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) read.delim(p, stringsAsFactors = FALSE)
norm <- function(d) {
  nm <- names(d)
  if (!"SNP" %in% nm) d$SNP <- if ("rsid" %in% nm) d$rsid else d[[1]]
  if ("sebeta" %in% nm && !"se" %in% nm) d$se <- d$sebeta
  if ("af" %in% nm && !"eaf" %in% nm) d$eaf <- d$af
  if ("alt" %in% nm && !"effect_allele" %in% nm) d$effect_allele <- d$alt
  if ("ref" %in% nm && !"other_allele" %in% nm) d$other_allele <- d$ref
  d
}

jobs <- list()
## --- A) 正向：放宽 CHIP IV ---
for (p in c("PV", "MPN", "MDS", "ITP", "AA")) {
  jobs[[length(jobs) + 1]] <- list(
    tag = sprintf("CHIPgP5E6_to_%s", p),
    iv  = file.path(SENS, "CHIP_germline_p5e6.tsv"),
    out = file.path(RLX, sprintf("%s_targets_relaxed.tsv", p)))
}
## --- B) 反向：MDS / ITP / AA (p<5e-6) ---
rev_files <- c(
  MDS = "finngen_R12_C3_MYELODYSP_SYNDR_EXALLC_MDS_iv.tsv",
  ITP = "finngen_R12_D3_ITP_ITP_iv.tsv",
  AA  = "finngen_R12_D3_OTHERAPLASTICANAEMIA_AA_iv.tsv")
for (p in names(rev_files)) {
  jobs[[length(jobs) + 1]] <- list(
    tag = sprintf("%s_P5E6_to_CHIP", p),
    iv  = file.path(IVD, rev_files[[p]]),
    out = file.path(TGT, "CHIP_targets.tsv"))
}

res <- list()
for (j in jobs) {
  if (!file.exists(j$iv) || !file.exists(j$out)) {
    cat(sprintf("%-26s SKIP (missing)\n", j$tag)); next
  }
  e <- try({
    iv <- norm(rd(j$iv)); o <- norm(rd(j$out))
    # 注意: F>10 的弱工具过滤只能作用于 exposure。
    # 若也用于 outcome, 就变成"按结局强度挑选 IV", 会造成严重选择偏倚
    # (初次误用时 PV 只剩 10 个 SNP, OR 虚高到 7.10)。
    fmt <- function(d, ty, lab) {
      raw <- data.frame(SNP = d$SNP, beta = as.numeric(d$beta), se = as.numeric(d$se),
                        pval = as.numeric(d$pval),
                        ea = d$effect_allele, oa = d$other_allele,
                        eaf = as.numeric(d$eaf), stringsAsFactors = FALSE)
      raw <- raw[is.finite(raw$beta) & is.finite(raw$se) & raw$se > 0, ]
      raw$F <- raw$beta^2 / raw$se^2
      fd <- format_data(raw, type = ty, snp_col = "SNP",
                        beta_col = "beta", se_col = "se", pval_col = "pval",
                        effect_allele_col = "ea", other_allele_col = "oa",
                        eaf_col = "eaf")
      if (ty == "exposure") { fd$id.exposure <- lab; fd$exposure <- lab }
      else                  { fd$id.outcome  <- lab; fd$outcome  <- lab }
      fd <- fd[!duplicated(fd$SNP), ]
      if (ty == "exposure") {
        keep <- raw$F[match(fd$SNP, raw$SNP)] > 10
        fd <- fd[!is.na(keep) & keep, ]
      }
      fd
    }
    h <- harmonise_data(fmt(iv, "exposure", j$tag), fmt(o, "outcome", j$tag), action = 2)
    h <- h[h$mr_keep, ]
    if (nrow(h) < 3) stop(sprintf("too few SNPs (%d)", nrow(h)))
    m <- mr(h, method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode"))
    q <- mr_heterogeneity(h); pl <- mr_pleiotropy_test(h)
    m$tag <- j$tag; m$nsnp_h <- nrow(h)
    m$Q_pval <- ifelse(m$method == "Inverse variance weighted",
                       q$Q_pval[q$method == "Inverse variance weighted"][1], NA)
    m$egger_p <- pl$pval[1]
    m
  }, silent = TRUE)
  if (inherits(e, "try-error")) {
    cat(sprintf("%-26s FAIL: %s\n", j$tag, trimws(as.character(e))))
  } else {
    res[[j$tag]] <- e
    ivw <- e[e$method == "Inverse variance weighted", ]
    wm  <- e[e$method == "Weighted median", ]
    cat(sprintf("%-26s n=%2d  IVW OR=%.2f (%.2f-%.2f) P=%.2g | WM OR=%.2f P=%.2g | EggerP=%.2f\n",
                j$tag, ivw$nsnp_h, exp(ivw$b), exp(ivw$b - 1.96 * ivw$se),
                exp(ivw$b + 1.96 * ivw$se), ivw$pval,
                exp(wm$b), wm$pval, ivw$egger_p))
  }
}

if (length(res) > 0) {
  all <- do.call(rbind, res); rownames(all) <- NULL
  write.csv(all, file.path(OUTD, "relaxed_MR_all.csv"), row.names = FALSE)
  cat("\nwritten ->", file.path(OUTD, "relaxed_MR_all.csv"), "\n")
}
