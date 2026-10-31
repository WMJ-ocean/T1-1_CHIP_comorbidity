#!/usr/bin/env Rscript
# 22 个 MPN 显著位点在各病中的效应方向（Figure 1 候选）
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(cowplot); library(ggrepel)
})

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
TD   <- file.path(BASE, "06_targets")
OUTD <- file.path(BASE, "08_figures"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

rd <- function(tag) {
  d <- read.delim(file.path(TD, sprintf("%s_targets.tsv", tag)), stringsAsFactors = FALSE)
  data.frame(SNP = d$rsid, beta = as.numeric(d$beta), p = as.numeric(d$pval),
             chr = d$chrom, pos = as.numeric(d$pos))
}
tags <- c("MPN", "PV", "MDS", "AA", "ITP")
L <- lapply(tags, rd); names(L) <- tags
snps <- Reduce(intersect, lapply(L, function(x) x$SNP))
cat("common SNPs:", length(snps), "\n")

mpn <- L$MPN[match(snps, L$MPN$SNP), ]

JAK2 <- "rs1327494"
make_panel <- function(tag, lbl) {
  y <- L[[tag]][match(snps, L[[tag]]$SNP), ]
  df <- data.frame(SNP = snps, x = mpn$beta, y = y$beta, p = y$p)
  df$lab <- ifelse(df$SNP == JAK2, "JAK2 (46/1)", "")
  # 显著位点加标签（y 方向）
  df$lab[is.na(df$lab)] <- ""
  r  <- cor(df$x, df$y)
  pr <- cor.test(df$x, df$y)$p.value
  same <- sum(sign(df$x) == sign(df$y))
  n <- nrow(df)
  sub <- sprintf("r = %.2f (P = %.3g)\nsame direction: %d/%d", r, pr, same, n)
  ggplot(df, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#b2182b", fill = "#f4a582", alpha = .25, linewidth = .7) +
    geom_point(aes(colour = p < 5e-8), size = 2.2, alpha = .85) +
    scale_colour_manual(values = c("FALSE" = "#2166ac", "TRUE" = "#b2182b"),
                        labels = c("p \u2265 5e-8", "p < 5e-8"), name = NULL) +
    geom_text_repel(aes(label = lab), size = 3, max.overlaps = 20, min.segment.length = 0,
                    box.padding = .35, seed = 1) +
    labs(x = "Effect on MPN (log OR)", y = sprintf("Effect on %s (log OR)", lbl),
         title = lbl, subtitle = sub) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8.5, colour = "grey25"),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
}

labs <- c(PV = "PV", MDS = "MDS", AA = "Aplastic anaemia", ITP = "ITP")
ps <- lapply(names(labs), function(t) make_panel(t, labs[[t]]))
names(ps) <- names(labs)

p <- plot_grid(plotlist = ps, nrow = 2, labels = "AUTO", align = "hv")
ggsave(file.path(OUTD, "Fig1_direction_consistency.pdf"), p, width = 9, height = 7.6)
ggsave(file.path(OUTD, "Fig1_direction_consistency.png"), p, width = 9, height = 7.6, dpi = 300)

# 单独输出一张只含 MPN vs AA 的放大图（核心发现）
ggsave(file.path(OUTD, "Fig1B_MPN_vs_AA.png"), ps[["AA"]], width = 4.6, height = 4.4, dpi = 300)

cat("written ->", OUTD, "\n")
for (t in names(ps)) {
  d <- data.frame(x = mpn$beta, y = L[[t]][match(snps, L[[t]]$SNP), ]$beta)
  cat(sprintf("%-4s  r=%+.3f  P=%.4g  same=%d/%d\n", t,
              cor(d$x, d$y), cor.test(d$x, d$y)$p.value,
              sum(sign(d$x) == sign(d$y)), nrow(d)))
}
