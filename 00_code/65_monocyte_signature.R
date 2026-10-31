## 65_monocyte_signature.R
## P2：GSE283711 CD14+ 单核细胞 (NBM 14 / PV 8 / ET 15 / MF 31) 验证 CHIP 因果基因集
## 分析：
##   A) limma-voom 差异表达：MF vs NBM、PV vs NBM、ET vs NBM、MF vs PV
##   B) 69 个 eQTLGen 靶基因的方向加权复合评分（权重 = exp_dir_chip）四组趋势
##   C) GSVA + MSigDB hallmark 通路活性
## 原则：效应量优先；不做事后挑基因；方向来自 CHIP MR（exp_dir_chip），不是表达数据本身。

suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(GSVA)
})

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
DAT  <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/analysis/00_data"
OUT  <- file.path(BASE, "07_p2_transcriptome")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------------- 1. 读入 ----------------
cnt <- fread(file.path(DAT, "GSE283711_CD14_rawcounts.csv.gz"), sep = ",", header = TRUE)
setnames(cnt, 1, "gene")
cnt <- cnt[!is.na(gene) & gene != ""]
## 同名基因取均值（避免重复行干扰）
m <- as.matrix(cnt[, -1])
rownames(m) <- cnt$gene
if (any(duplicated(rownames(m)))) {
  m <- rowsum(m, group = rownames(m))
  m <- m / as.vector(table(rownames(m)))
}
storage.mode(m) <- "double"

grp <- sub("^([A-Za-z]+)[0-9]*_.*$", "\\1", colnames(m))   # MF12_CD14 -> MF
grp <- factor(grp, levels = c("NBM", "ET", "PV", "MF"))
cat("[data] 基因 x 样本 =", nrow(m), "x", ncol(m), "\n")
print(table(grp))

## ---------------- 2. 差异表达 limma-voom ----------------
keep <- filterByExpr(m, group = grp)
cat("[filter] 保留基因", sum(keep), "/", nrow(m), "\n")
dge <- DGEList(counts = m[keep, ], group = grp)
dge <- calcNormFactors(dge, method = "TMM")
v   <- voom(dge, design = model.matrix(~ 0 + grp), plot = FALSE)
fit <- lmFit(v, design = model.matrix(~ 0 + grp))
cm  <- makeContrasts(MF_NBM = grpMF - grpNBM,
                     PV_NBM = grpPV - grpNBM,
                     ET_NBM = grpET - grpNBM,
                     MF_PV  = grpMF - grpPV,
                     levels = fit$design)
fit2 <- eBayes(contrasts.fit(fit, cm))

de <- list()
for (k in colnames(cm)) {
  tt <- topTable(fit2, coef = k, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  de[[k]] <- tt
  write.table(tt, file.path(OUT, paste0("DE_", k, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[DE] %-6s  FDR<0.05 = %d  (up %d / down %d)\n", k,
              sum(tt$adj.P.Val < 0.05),
              sum(tt$adj.P.Val < 0.05 & tt$logFC > 0),
              sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)))
}

logcpm <- cpm(dge, log = TRUE, prior.count = 2)

## ---------------- 3. CHIP 因果基因集复合评分 ----------------
mr <- fread(file.path(BASE, "06_targets/gene_mr_wald_corrected.tsv"), sep = "\t")
mr <- mr[!is.na(exp_dir_chip) & exp_dir_chip %in% c("up", "down")]
## 每个基因只保留一条：优先 MPN 结局，其次 PV/其他
mr[, prio := fifelse(outcome == "MPN", 1L, 2L)]
setorder(mr, GeneSymbol, prio)
gs <- mr[, .SD[1], by = GeneSymbol]
gs <- gs[GeneSymbol %in% rownames(logcpm)]
cat("[signature] 可用基因", nrow(gs), "个 (up", sum(gs$exp_dir_chip == "up"),
    "/ down", sum(gs$exp_dir_chip == "down"), ")\n")

z <- t(scale(t(logcpm[gs$GeneSymbol, , drop = FALSE])))
z[!is.finite(z)] <- 0
w <- ifelse(gs$exp_dir_chip == "up", 1, -1)
score <- colMeans(z * w)

## 分组汇总 + 趋势检验（组编码 NBM=0, ET=1, PV=2, MF=3）
gcode <- as.integer(grp) - 1L
kt <- kruskal.test(score ~ grp)
sp <- cor.test(gcode, score, method = "spearman")
sm <- data.table(sample = names(score), group = grp, score = as.numeric(score))
summ <- sm[, .(n = .N, median = median(score), mean = mean(score),
               sd = sd(score),
               q1 = quantile(score, .25), q3 = quantile(score, .75)), by = group]
setorder(summ, group)
write.table(summ, file.path(OUT, "signature_score_by_group.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sm, file.path(OUT, "signature_score_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[sig] Kruskal p =", signif(kt$p.value, 3),
    " Spearman trend rho =", signif(sp$estimate, 3),
    " p =", signif(sp$p.value, 3), "\n")
print(summ)

## 亚组：只看 CHIP 上调 / 下调基因，分别评分（方向一致性检验）
for (d in c("up", "down")) {
  idx <- gs$GeneSymbol[gs$exp_dir_chip == d]
  idx <- idx[idx %in% rownames(logcpm)]
  if (length(idx) >= 5) {
    sc <- colMeans(z[idx, , drop = FALSE])
    s2 <- data.table(sample = names(sc), group = grp, score = as.numeric(sc))
    sm2 <- s2[, .(median = median(score)), by = group]
    cat(sprintf("[sig-%s] n=%d  ", d, length(idx)))
    cat(paste(sprintf("%s=%.3f", sm2$group, sm2$median), collapse = "  "),
        "  trend p =", signif(cor.test(gcode, as.numeric(sc), method = "spearman")$p.value, 3), "\n")
    fwrite(s2, file.path(OUT, paste0("signature_", d, "_per_sample.tsv")), sep = "\t")
  }
}

## 单个基因的方向一致性：MF vs NBM 的 logFC 与 exp_dir_chip 是否同号
dn <- de[["MF_NBM"]]
dn <- dn[dn$gene %in% gs$GeneSymbol, c("gene", "logFC", "P.Value", "adj.P.Val")]
dn <- merge(dn, gs[, .(gene = GeneSymbol, dir = exp_dir_chip)], by = "gene")
dn <- as.data.table(dn)
dn[, concordant := sign(logFC) == ifelse(dir == "up", 1, -1)]
write.table(dn, file.path(OUT, "signature_gene_level_MF_vs_NBM.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
tb <- table(dn$dir, dn$concordant)
cat("[concordance] MF vs NBM 方向一致率 =",
    round(100 * sum(dn$concordant) / nrow(dn), 1), "%\n")
print(tb)
if (nrow(tb) == 2 && ncol(tb) == 2) {
  cat("[concordance] Fisher p =", signif(fisher.test(tb)$p.value, 3), "\n")
}

## ---------------- 4. GSVA hallmark ----------------
gmt <- readLines(file.path(DAT, "msigdb_v2024.1_Hs.gmt"))
sets <- list()
for (ln in gmt) {
  f <- strsplit(ln, "\t")[[1]]
  if (length(f) > 2 && grepl("^HALLMARK_", f[1])) sets[[f[1]]] <- f[-c(1, 2)]
}
cat("[gsva] hallmark 集", length(sets), "个\n")
if (length(sets) > 0) {
  par <- gsvaParam(logcpm, sets, kcdf = "Gaussian", minSize = 10, maxSize = 500)
  gsv <- gsva(par, verbose = FALSE)
  write.table(data.frame(pathway = rownames(gsv), gsv, check.names = FALSE),
              file.path(OUT, "gsva_hallmark_scores.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  gfit <- lmFit(gsv, design = model.matrix(~ 0 + grp))
  gfit2 <- eBayes(contrasts.fit(gfit, cm))
  gt <- topTable(gfit2, coef = "MF_NBM", number = Inf, sort.by = "P")
  gt$pathway <- rownames(gt)
  write.table(gt, file.path(OUT, "gsva_hallmark_MF_vs_NBM.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[gsva] MF vs NBM  FDR<0.05 =", sum(gt$adj.P.Val < 0.05), "个通路\n")
  print(head(gt[, c("pathway", "logFC", "P.Value", "adj.P.Val")], 10))
}

## ---------------- 5. 图 ----------------
png(file.path(OUT, "FigP2_signature_boxplot.png"), width = 1600, height = 1200, res = 200)
boxplot(score ~ grp, col = c("#9ecae1", "#a1d99b", "#fdd0a2", "#fc9272"),
        xlab = "CD14+ monocyte group", ylab = "CHIP-direction weighted signature score",
        main = sprintf("GSE283711 (NBM14/ET15/PV8/MF31)  Kruskal p=%.2g, trend p=%.2g",
                       kt$p.value, sp$p.value),
        outline = FALSE)
stripchart(score ~ grp, vertical = TRUE, method = "jitter", pch = 19,
           col = rgb(0, 0, 0, .35), add = TRUE, cex = .7)
dev.off()

cat("[done] 输出目录:", OUT, "\n")
