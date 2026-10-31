# 28_steiger_r2.R —— 方向性辅助检验：比较双向的遗传解释方差 R²
# Steiger 的核心逻辑：真方向暴露端的 R² 应显著大于结局端
suppressMessages(library(dplyr))

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
OUTD <- file.path(BASE, "07_mr")

rd <- function(p) { d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE); names(d)[1] <- "SNP"; d }

pairs <- list(
  list(tag = "CHIPg_to_MPN", f = "CHIPg_to_MPN_harmonised.csv", dir = "CHIP -> MPN"),
  list(tag = "CHIPg_to_PV",  f = "CHIPg_to_PV_harmonised.csv",  dir = "CHIP -> PV"),
  list(tag = "MPN_to_CHIP",  f = "MPN_to_CHIP_harmonised.csv",  dir = "MPN -> CHIP"),
  list(tag = "PV_to_CHIP",   f = "PV_to_CHIP_harmonised.csv",   dir = "PV -> CHIP")
)

res <- list()
for (p in pairs) {
  fp <- file.path(OUTD, p$f)
  if (!file.exists(fp)) next
  d <- rd(fp)
  if (!"mr_keep" %in% names(d)) next
  d <- d[d$mr_keep, ]
  # R² = sum(2*eaf*(1-eaf)*beta^2)   （单位方差标准化近似）
  r2e <- sum(2 * d$eaf.exposure * (1 - d$eaf.exposure) * d$beta.exposure^2, na.rm = TRUE)
  r2o <- sum(2 * d$eaf.outcome  * (1 - d$eaf.outcome)  * d$beta.outcome^2,  na.rm = TRUE)
  res[[p$tag]] <- data.frame(pair = p$tag, direction = p$dir, nSNP = nrow(d),
                             R2_exposure = r2e, R2_outcome = r2o,
                             ratio_exp_over_out = r2e / r2o,
                             supports = ifelse(r2e > r2o, "stated direction", "reverse"),
                             stringsAsFactors = FALSE)
}

if (length(res) > 0) {
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  print(out, row.names = FALSE, digits = 4)
  write.csv(out, file.path(OUTD, "steiger_R2_summary.csv"), row.names = FALSE)
  cat("\nwritten ->", file.path(OUTD, "steiger_R2_summary.csv"), "\n")
}
