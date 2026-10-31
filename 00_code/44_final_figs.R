#!/usr/bin/env Rscript
# 44_final_figs.R —— 生成投稿规格 Figure 1 与 Figure 3（Leukemia: >=300 dpi, PDF + PNG）
# Figure 2 由 32_fig3_model.R 生成的机制图直接复用
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
FIG  <- file.path(BASE, "10_manuscript", "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

TH <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.position = "none")

# ============================================================
# Figure 1: A 正向 MR 森林图 | B coloc 后验概率
# ============================================================
mr <- data.frame(
  pheno = factor(c("PV", "MPN", "MDS", "ITP", "AA"), levels = c("PV", "MPN", "MDS", "ITP", "AA")),
  OR  = c(3.48, 2.93, 1.34, 0.95, 0.84),
  low = c(2.26, 2.02, 1.02, 0.74, 0.58),
  high= c(5.37, 4.25, 1.76, 1.23, 1.21),
  pn  = c(1.6e-08, 1.5e-08, 0.037, 0.70, 0.35))
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
mr$star <- sig(mr$pn)

p1A <- ggplot(mr, aes(x = OR, y = pheno, colour = pheno)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.20, linewidth = 0.6) +
  geom_point(size = 3.0) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f) %s", OR, low, high, star)),
            hjust = 0, nudge_x = 0.45, size = 2.5, colour = "grey20") +
  scale_x_log10(limits = c(0.42, 60), breaks = c(0.5, 1, 2, 4, 8)) +
  scale_colour_manual(values = c(PV = "#B03A2E", MPN = "#CA6F1E", MDS = "#7D3C98",
                                 ITP = "#1E8449", AA = "#2471A3")) +
  labs(x = "Odds ratio (95% CI) per unit CHIP liability", y = NULL,
       title = "A   CHIP liability -> haematological disease") + TH

coloc <- data.frame(
  region = factor(c("TERT", "RUNX1\n(21q22)", "ATM\n(11q22)", "TP53\n(17p13)", "JAK2\n(9p24)", "3q25"),
                  levels = c("TERT", "RUNX1\n(21q22)", "ATM\n(11q22)", "TP53\n(17p13)", "JAK2\n(9p24)", "3q25")),
  trait  = rep(c("MPN", "Polycythaemia vera"), each = 6),
  H4     = c(0.998, 0.984, 0.969, 0.875, 0.000, 0.000,
             0.046, 0.984, 0.767, 0.461, 0.000, 0.000),
  H3     = c(0.002, 0.016, 0.031, 0.125, 1.000, 1.000,
             0.954, 0.016, 0.233, 0.539, 1.000, 1.000))
coloc$lab  <- ifelse(coloc$H3 > 0.8, "", sprintf("%.3f", coloc$H4))
coloc$inside <- coloc$H4 >= 0.30
coloc$ypos <- ifelse(coloc$inside, coloc$H4 - 0.055, coloc$H4 + 0.045)
coloc$tcol <- ifelse(coloc$inside, "white", "grey25")
coloc_note <- data.frame(region = factor(c("JAK2\n(9p24)", "3q25"),
                                         levels = levels(coloc$region)),
                         y = 0.05, lab = "H3 = 1.00")

p1B <- ggplot(coloc, aes(x = region, y = H4, fill = trait)) +
  geom_col(position = position_dodge(0.75), width = 0.68, colour = "grey25", linewidth = 0.25) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#B03A2E", linewidth = 0.35) +
  geom_text(aes(y = ypos, label = lab, group = trait, colour = tcol),
            size = 2.2, position = position_dodge(0.75), show.legend = FALSE) +
  scale_colour_identity() +
  geom_text(data = coloc_note, aes(x = region, y = y, label = lab),
            inherit.aes = FALSE, size = 2.1, colour = "grey35", angle = 90) +
  scale_y_continuous(limits = c(0, 1.06), expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_fill_manual(values = c("MPN" = "#7D3C98", "Polycythaemia vera" = "#AEB6BF")) +
  labs(x = NULL, y = "PP.H4 (shared causal variant)",
       title = "B   Colocalisation at CHIP loci", fill = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        axis.title = element_text(size = 8.5), axis.text = element_text(size = 6.8),
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.position = c(0.72, 0.88), legend.background = element_blank(),
        legend.key.size = unit(3, "mm"))

fig1 <- (p1A | p1B) + plot_layout(widths = c(1, 1.05))
ggsave(file.path(FIG, "Figure1_MR_coloc.pdf"), fig1, width = 180, height = 78, units = "mm")
ggsave(file.path(FIG, "Figure1_MR_coloc.png"), fig1, width = 180, height = 78, units = "mm", dpi = 600)

# ============================================================
# Figure 3: A 敏感性 | B 反向 MR | C 方向一致性 | D 双向不对称
# ============================================================
sens <- data.frame(
  model = factor(c("Baseline\n(20 IV)", "Excl. TERT\n(19 IV)", "Excl. MHC\n(19 IV)",
                   "Relaxed p<5e-6\n(79 IV)", "MR-PRESSO\ncorrected"),
                 levels = c("Baseline\n(20 IV)", "Excl. TERT\n(19 IV)", "Excl. MHC\n(19 IV)",
                            "Relaxed p<5e-6\n(79 IV)", "MR-PRESSO\ncorrected")),
  trait  = rep(c("MPN", "PV"), each = 5),
  OR   = c(2.93, 2.70, 2.96, 2.38, 3.45, 3.48, 3.19, 3.52, 2.74, 3.62),
  low  = c(2.02, 1.72, 2.02, 1.93, 2.66, 2.26, 1.88, 2.25, 2.10, 2.71),
  high = c(4.25, 4.24, 4.34, 2.94, 4.46, 5.37, 5.41, 5.50, 3.58, 4.84))

p3A <- ggplot(sens, aes(x = OR, y = model, colour = trait, shape = trait)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.18, linewidth = 0.5,
                 position = position_dodge(0.55)) +
  geom_point(size = 2.4, position = position_dodge(0.55)) +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 4, 8)) +
  scale_colour_manual(values = c("MPN" = "#CA6F1E", "PV" = "#B03A2E")) +
  labs(x = "Odds ratio (95% CI)", y = NULL, title = "A   Sensitivity analyses") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        axis.title = element_text(size = 8.5), axis.text = element_text(size = 6.8),
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.position = "bottom", legend.background = element_blank(),
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0))

rev <- data.frame(
  exposure = factor(c("MPN", "PV", "MDS", "ITP", "AA"),
                    levels = c("MPN", "PV", "MDS", "ITP", "AA")),
  OR  = c(1.27, 1.19, 0.99, 1.01, 1.02),
  low = c(1.11, 1.00, 0.96, 1.00, 1.00),
  high= c(1.46, 1.43, 1.01, 1.03, 1.03),
  pn  = c(7.4e-04, 0.056, 0.20, 0.098, 0.068))
rev$star <- sig(rev$pn)
p3B <- ggplot(rev, aes(x = OR, y = exposure)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.20, linewidth = 0.6, colour = "#2874A6") +
  geom_point(size = 3.0, colour = "#2874A6") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f) %s", OR, low, high, star)),
            hjust = 0, nudge_x = 0.035, size = 2.3, colour = "grey20") +
  scale_x_continuous(limits = c(0.93, 2.15), breaks = c(1.0, 1.2, 1.4)) +
  labs(x = "Odds ratio (95% CI) for CHIP", y = NULL,
       title = "B   Reverse direction") + TH

dir <- data.frame(
  pheno = factor(c("PV", "MDS", "ITP", "AA"), levels = c("PV", "MDS", "ITP", "AA")),
  same  = c(22, 16, 12, 5), opp = c(0, 6, 10, 17),
  r     = c(0.956, -0.085, -0.260, -0.496),
  rp    = c("4.0e-12", "0.71", "0.24", "0.011"))
d2 <- tidyr::pivot_longer(dir, cols = c(same, opp), names_to = "dir", values_to = "n")
d2$dir <- factor(d2$dir, levels = c("same", "opp"), labels = c("Concordant", "Opposite"))
p3C <- ggplot(d2, aes(x = pheno, y = n, fill = dir)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.6, colour = "white") +
  geom_text(data = dir, aes(x = pheno, y = 23.5,
                            label = sprintf("r = %+.2f\nP = %s", r, rp)),
            inherit.aes = FALSE, size = 2.3, colour = "grey25", lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 27), expand = c(0, 0)) +
  scale_fill_manual(values = c("Concordant" = "#7D3C98", "Opposite" = "#2471A3")) +
  labs(x = NULL, y = "MPN risk loci (n = 22)", title = "C   Direction concordance") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        axis.title = element_text(size = 8.5), axis.text = element_text(size = 7),
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.position = "bottom", legend.background = element_blank(),
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0))

asym <- data.frame(
  lbl  = factor(c("PV  forward", "PV  reverse", "MPN  forward", "MPN  reverse"),
                levels = c("PV  forward", "PV  reverse", "MPN  forward", "MPN  reverse")),
  OR   = c(3.48, 1.19, 2.93, 1.27),
  low  = c(2.26, 1.00, 2.02, 1.11),
  high = c(5.37, 1.43, 4.25, 1.46),
  grp  = c("PV", "PV", "MPN", "MPN"))
p3D <- ggplot(asym, aes(x = lbl, y = OR, fill = grp)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.25) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.13, linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", OR), y = high), vjust = -0.6, size = 2.7, colour = "grey20") +
  scale_y_continuous(limits = c(0, 6.2), expand = c(0, 0)) +
  scale_fill_manual(values = c("PV" = "#B03A2E", "MPN" = "#CA6F1E")) +
  labs(x = NULL, y = "Odds ratio (95% CI)",
       title = "D   Effect asymmetry") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(size = 6.2, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 7), legend.position = "none")

fig3 <- (p3A | p3B) / (p3C | p3D) + plot_layout(heights = c(1, 1))
ggsave(file.path(FIG, "Figure3_sensitivity.pdf"), fig3, width = 180, height = 150, units = "mm")
ggsave(file.path(FIG, "Figure3_sensitivity.png"), fig3, width = 180, height = 150, units = "mm", dpi = 600)

cat("Figure 1 & 3 written to", FIG, "\n")
