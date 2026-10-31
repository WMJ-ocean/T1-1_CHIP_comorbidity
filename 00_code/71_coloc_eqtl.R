## 71_coloc_eqtl.R —— eQTL coloc：CHIP 胚系信号 vs 血液 cis-eQTL（eQTLGen 全量）
## 输入：
##   04_clean/region/CHIPIV_<rs>.tsv        CHIP GWAS 各 IV ±1Mb 全量 sumstats
##   06_targets/eqtlgen_windows_full.tsv.gz eQTLGen 全量 cis-eQTL 同窗口全部 SNP-基因对
##   02_raw/eqtlgen/eQTLGen_SNP_AF_AlleleB.txt.gz  等位基因频率（Z -> beta/se 必需）
## 输出：09_coloc/eqtl_coloc_summary.csv

suppressPackageStartupMessages({ library(data.table); library(coloc) })

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
REG  <- file.path(BASE, "04_clean", "region")
EQ   <- file.path(BASE, "06_targets", "eqtlgen_windows_full.tsv.gz")
AF   <- file.path(BASE, "02_raw", "eqtlgen", "eQTLGen_SNP_AF_AlleleB.txt.gz")
OUT  <- file.path(BASE, "09_coloc")
IVF  <- file.path(BASE, "05_iv", "CHIP_germline_iv.tsv")
WIN  <- 1e6

iv <- fread(IVF)
setnames(iv, c("SNP", "chr", "pos", "ea", "oa", "beta", "se", "pval", "eaf"))
cat("[setup] IV", nrow(iv), "个\n")

## ---- AF ----
af <- fread(AF, select = c("SNP", "AlleleA", "AlleleB", "AlleleB_all"))
setkey(af, SNP)
cat("[af] 位点频率", nrow(af), "\n")

## ---- eQTL 全量窗口 ----
eq <- fread(EQ)
cat("[eqtl] 行", nrow(eq), "\n")
eq[, f := as.numeric(NA)]

comp <- c(A = "T", T = "A", C = "G", G = "C")

res <- list()
for (i in seq_len(nrow(iv))) {
  rs <- iv$SNP[i]; ch <- as.character(iv$chr[i]); p0 <- iv$pos[i]
  fc <- file.path(REG, paste0("CHIPIV_", rs, ".tsv"))
  if (!file.exists(fc)) { cat("  [skip] 无 CHIP 区域文件:", rs, "\n"); next }
  chip <- fread(fc)
  if (nrow(chip) < 100) { cat("  [skip] CHIP 区域 SNP 太少:", rs, "\n"); next }
  chip[, pos := as.integer(pos)]

  sub <- eq[SNPChr == ch & SNPPos >= p0 - WIN & SNPPos <= p0 + WIN]
  if (nrow(sub) == 0) { cat("  [skip] 无 eQTL:", rs, "\n"); next }
  ## 基因本体需落在窗口内（cis 才有解释力）
  genes <- sub[GenePos >= p0 - WIN & GenePos <= p0 + WIN,
               .(minFDR = min(FDR, na.rm = TRUE)), by = GeneSymbol]
  genes <- genes[minFDR < 0.05][order(minFDR)]
  if (nrow(genes) == 0) { cat("  [", rs, "] 窗口内无 FDR<0.05 的基因\n"); next }
  cat(sprintf("[%s] chr%s:%d  窗口 SNP-基因对 %d  显著基因 %d\n",
              rs, ch, p0, nrow(sub), nrow(genes)))

  ## 频率
  sub <- merge(sub, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
  sub[, f := fifelse(AssessedAllele == AlleleB, AlleleB_all,
             fifelse(AssessedAllele == AlleleA, 1 - AlleleB_all, NA_real_))]
  sub[, k := 2 * f * (1 - f) * (NrSamples + Zscore^2)]
  sub[, b_e := Zscore / sqrt(k)]
  sub[, se_e := 1 / sqrt(k)]

  for (g in genes$GeneSymbol) {
    d <- sub[GeneSymbol == g & is.finite(b_e) & is.finite(se_e) & f > 0.01 & f < 0.99]
    if (nrow(d) == 0) next
    d <- d[!duplicated(SNP)]
    m <- merge(d[, .(SNP, ea = AssessedAllele, oa = OtherAllele, b_e, se_e, f, NrSamples)],
               chip[, .(SNP, c_ea = ea, c_oa = oa, c_b = beta, c_se = se, pos)],
               by = "SNP")
    if (nrow(m) < 50) next
    ## 把 CHIP 的 beta 对齐到 eQTL 的 AssessedAllele
    m[, c_b2 := fifelse(ea == c_ea & oa == c_oa, c_b,
                fifelse(ea == c_oa & oa == c_ea, -c_b,
                fifelse(ea == comp[c_oa] & oa == comp[c_ea], c_b,
                fifelse(ea == comp[c_ea] & oa == comp[c_oa], -c_b, NA_real_))))]
    m <- m[is.finite(c_b2) & is.finite(c_se) & c_se > 0]
    if (nrow(m) < 50) next
    m <- m[order(pos)]

    d1 <- list(snp = m$SNP, beta = m$b_e, varbeta = m$se_e^2,
               type = "quant", sdY = 1, MAF = pmin(m$f, 1 - m$f),
               N = as.integer(median(m$NrSamples, na.rm = TRUE)))
    d2 <- list(snp = m$SNP, beta = m$c_b2, varbeta = m$c_se^2,
               type = "quant", sdY = 1, MAF = pmin(m$f, 1 - m$f), N = 418401)
    r <- try(suppressWarnings(coloc.abf(d1, d2, p12 = 1e-5)), silent = TRUE)
    if (inherits(r, "try-error")) {
      cat(sprintf("   %-12s %-14s FAILED\n", rs, g)); next
    }
    s <- r$summary
    lead1 <- m$SNP[which.max(m$b_e^2 / m$se_e^2)][1]
    lead2 <- m$SNP[which.max(m$c_b2^2 / m$c_se^2)][1]
    best  <- r$results$snp[which.max(r$results$SNP.PP.H4)]
    res[[paste(rs, g)]] <- data.table(
      IV = rs, gene = g, nSNP = nrow(m),
      PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
      lead_eqtl = lead1, lead_chip = lead2, same_lead = lead1 == lead2,
      z_eqtl_lead = m$b_e[which.max(m$b_e^2 / m$se_e^2)][1] / m$se_e[which.max(m$b_e^2 / m$se_e^2)][1],
      z_chip_lead = m$c_b2[which.max(m$c_b2^2 / m$c_se^2)][1] / m$c_se[which.max(m$c_b2^2 / m$c_se^2)][1],
      top_H4_snp = best, top_H4_pp = max(r$results$SNP.PP.H4))
    cat(sprintf("   %-12s %-14s n=%5d  H3=%.3f  H4=%.3f  same=%s\n",
                rs, g, nrow(m), s[["PP.H3.abf"]], s[["PP.H4.abf"]],
                ifelse(lead1 == lead2, "YES", "no")))
  }
}

if (length(res)) {
  out <- rbindlist(res)
  setorder(out, -PP.H4)
  fwrite(out, file.path(OUT, "eqtl_coloc_summary.csv"))
  cat("\n[done]", nrow(out), "对 ->", file.path(OUT, "eqtl_coloc_summary.csv"), "\n")
  print(head(out[, .(IV, gene, nSNP, PP.H4, same_lead)], 15))
  cat("\nPP.H4 >= 0.8 :", sum(out$PP.H4 >= 0.8),
      "  0.5-0.8 :", sum(out$PP.H4 >= 0.5 & out$PP.H4 < 0.8), "\n")
} else {
  cat("\n[warn] 没有可运行的对\n")
}
