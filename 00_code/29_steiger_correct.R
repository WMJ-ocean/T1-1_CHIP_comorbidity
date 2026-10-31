# 29_steiger_correct.R —— 正确的方向性检验（case-control 专用 get_r_from_lor）
# 修正 28 号脚本用线性近似导致极罕见结局端 R² 被放大的问题
suppressMessages(library(TwoSampleMR))

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
OUTD <- file.path(BASE, "07_mr")

# 精确样本量（FinnGen R12 manifest + Kessler 2022）
N <- list(
  CHIP = c(ncase = 25657, ncontrol = 392744, prev = 0.061),
  MPN  = c(ncase =  2791, ncontrol = 497256, prev = 0.0003),
  PV   = c(ncase =  1205, ncontrol = 378022, prev = 0.00006),
  MDS  = c(ncase =   979, ncontrol = 378296, prev = 0.00005),
  AA   = c(ncase =   374, ncontrol = 479581, prev = 0.00004),
  ITP  = c(ncase =  1064, ncontrol = 492275, prev = 0.0001)
)

rd <- function(p) { d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE); names(d)[1] <- "SNP"; d }

pairs <- list(
  list(tag = "CHIPg_to_MPN", f = "CHIPg_to_MPN_harmonised.csv", e = "CHIP", o = "MPN"),
  list(tag = "CHIPg_to_PV",  f = "CHIPg_to_PV_harmonised.csv",  e = "CHIP", o = "PV"),
  list(tag = "CHIPg_to_MDS", f = "CHIPg_to_MDS_harmonised.csv", e = "CHIP", o = "MDS"),
  list(tag = "MPN_to_CHIP",  f = "MPN_to_CHIP_harmonised.csv",  e = "MPN",  o = "CHIP"),
  list(tag = "PV_to_CHIP",   f = "PV_to_CHIP_harmonised.csv",   e = "PV",   o = "CHIP")
)

res <- list()
for (p in pairs) {
  fp <- file.path(OUTD, p$f)
  if (!file.exists(fp)) next
  d <- rd(fp); d <- d[d$mr_keep %in% c(TRUE, "TRUE"), ]
  if (nrow(d) < 3) next
  ne <- N[[p$e]]; no <- N[[p$o]]

  re <- get_r_from_lor(lor = d$beta.exposure, af = d$eaf.exposure,
                       ncase = ne[["ncase"]], ncontrol = ne[["ncontrol"]],
                       prevalence = ne[["prev"]], model = "logit")
  ro <- get_r_from_lor(lor = d$beta.outcome, af = d$eaf.outcome,
                       ncase = no[["ncase"]], ncontrol = no[["ncontrol"]],
                       prevalence = no[["prev"]], model = "logit")

  r2e <- sum(re^2, na.rm = TRUE); r2o <- sum(ro^2, na.rm = TRUE)
  res[[p$tag]] <- data.frame(
    pair = p$tag, exposure = p$e, outcome = p$o, nSNP = nrow(d),
    R2_exp = r2e, R2_out = r2o, ratio = r2e / r2o,
    steiger_dir = ifelse(r2e > r2o, "exposure -> outcome", "outcome -> exposure"),
    # Steiger 显著性检验需要独立样本的 r，此处仅报告 R2 比值方向
    stringsAsFactors = FALSE)
}

if (length(res) > 0) {
  out <- do.call(rbind, res); rownames(out) <- NULL
  print(out, row.names = FALSE, digits = 4)
  write.csv(out, file.path(OUTD, "steiger_corrected.csv"), row.names = FALSE)
  cat("\nwritten -> steiger_corrected.csv\n")
}
