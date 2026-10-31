# 26_coloc.R —— 区域级共定位 (coloc.abf)
# 用法: Rscript 26_coloc.R
suppressMessages(library(coloc))

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
IND  <- file.path(BASE, "09_coloc", "input")
OUTD <- file.path(BASE, "09_coloc")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

files <- list.files(IND, pattern = "^.*_CHIP_vs_.*\\.tsv$", full.names = FALSE)
if (length(files) == 0) { cat("no input files\n"); quit(status = 0) }

res <- list()
for (fn in files) {
  d <- read.delim(file.path(IND, fn), stringsAsFactors = FALSE)
  d <- d[!duplicated(d$SNP), ]
  d <- d[order(d$pos), ]
  # 过滤极端 MAF 与缺失
  d <- d[is.finite(d$b1) & is.finite(d$b2) & is.finite(d$se1) & is.finite(d$se2), ]
  d <- d[d$se1 > 0 & d$se2 > 0, ]
  if (nrow(d) < 50) { cat(sprintf("%-34s only %d SNPs, skip\n", fn, nrow(d))); next }

  ds1 <- list(snp = d$SNP, position = d$pos, beta = d$b1, varbeta = d$se1^2,
              type = "quant", sdY = 1, MAF = d$maf1, N = 418401)
  ds2 <- list(snp = d$SNP, position = d$pos, beta = d$b2, varbeta = d$se2^2,
              type = "quant", sdY = 1, MAF = d$maf2, N = 500000)

  r <- try(suppressWarnings(coloc.abf(dataset1 = ds1, dataset2 = ds2, p12 = 1e-5)),
           silent = TRUE)
  if (inherits(r, "try-error")) {
    cat(sprintf("%-34s FAILED: %s\n", fn, as.character(r)))
    next
  }
  s <- r$summary
  # lead SNP（各自最小 p）
  lead1 <- d$SNP[which.min(d$p1)]
  lead2 <- d$SNP[which.min(d$p2)]
  best  <- r$results$snp[which.max(r$results$SNP.PP.H4)]

  res[[fn]] <- data.frame(
    pair = sub("\\.tsv$", "", fn),
    nSNP = nrow(d),
    PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
    PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
    lead_CHIP = lead1, lead_pheno = lead2, same_lead = lead1 == lead2,
    top_H4_snp = best, top_H4_pp = max(r$results$SNP.PP.H4),
    stringsAsFactors = FALSE)
  cat(sprintf("%-34s n=%5d  H3=%.3f  H4=%.3f  same_lead=%s\n",
              sub("\\.tsv$", "", fn), nrow(d), s[["PP.H3.abf"]], s[["PP.H4.abf"]],
              ifelse(lead1 == lead2, "YES", "no")))
}

if (length(res) > 0) {
  out <- do.call(rbind, res)
  write.csv(out, file.path(OUTD, "coloc_summary.csv"), row.names = FALSE)
  cat("\n=== summary ===\n")
  print(out[, c("pair","nSNP","PP.H3","PP.H4","same_lead","top_H4_snp")], row.names = FALSE)
  cat("\nwritten ->", file.path(OUTD, "coloc_summary.csv"), "\n")
}
