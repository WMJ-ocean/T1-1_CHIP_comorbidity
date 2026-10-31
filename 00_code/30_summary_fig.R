# 30_summary_fig.R —— Figure 2: MR 效应梯度 + coloc 共定位证据
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
FIG  <- file.path(BASE, "08_figures")
dir.create(FIG, showWarnings = FALSE)

# ---------- Panel A: 正向 MR 五种血液病 ----------
mr <- data.frame(
  pheno = factor(c("PV", "MPN", "MDS", "ITP", "AA"),
                 levels = c("PV", "MPN", "MDS", "ITP", "AA")),
  OR  = c(3.48, 2.93, 1.34, 0.95, 0.84),
  low = c(2.26, 2.02, 1.02, 0.74, 0.58),
  high= c(5.37, 4.25, 1.76, 1.23, 1.21),
  p   = c("1.6e-08", "1.5e-08", "0.037", "0.70", "0.35"),
  nsnp= c(20, 20, 20, 20, 20)
)
mr$lab <- sprintf("%s (n=%d)", mr$pheno, mr$nsnp)

pA <- ggplot(mr, aes(x = OR, y = pheno, colour = pheno)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.22, linewidth = 0.7) +
  geom_point(size = 3.4) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", OR, low, high)),
            hjust = 0, nudge_x = 0.30, size = 2.8, colour = "grey25") +
  scale_x_log10(limits = c(0.4, 30), breaks = c(0.5, 1, 2, 4, 8)) +
  scale_colour_manual(values = c(PV = "#C0392B", MPN = "#E67E22", MDS = "#8E44AD",
                                 ITP = "#27AE60", AA = "#2980B9"), guide = "none") +
  labs(x = "OR (95% CI) per genetically predicted CHIP", y = NULL,
       title = "A  CHIP (germline liability) → haematological disease") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

# ---------- Panel B: 双向对比 ----------
bidir <- data.frame(
  grp = factor(c("MPN | CHIP → disease", "MPN | disease → CHIP",
                 "PV | CHIP → disease",  "PV | disease → CHIP"),
               levels = rev(c("MPN | CHIP → disease", "MPN | disease → CHIP",
                              "PV | CHIP → disease",  "PV | disease → CHIP"))),
  dir   = factor(c("CHIP → disease", "disease → CHIP", "CHIP → disease", "disease → CHIP"),
                 levels = c("CHIP → disease", "disease → CHIP")),
  OR    = c(2.93, 1.27, 3.48, 1.19),
  low   = c(2.02, 1.11, 2.26, 1.00),
  high  = c(4.25, 1.46, 5.37, 1.43)
)

pB <- ggplot(bidir, aes(x = OR, y = grp, fill = dir)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.2, linewidth = 0.65) +
  geom_point(shape = 21, size = 3.2, stroke = 0.5) +
  scale_x_log10(limits = c(0.9, 7.5), breaks = c(1, 2, 4)) +
  scale_fill_manual(values = c("CHIP → disease" = "#34495E",
                               "disease → CHIP" = "#95A5A6")) +
  labs(x = "OR (95% CI)", y = NULL, fill = NULL,
       title = "B  Bidirectional MR: forward effect dominates") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey93"),
        plot.title = element_text(face = "bold", size = 10))

# ---------- Panel C: coloc 共定位热图 ----------
cl <- read.csv(file.path(BASE, "09_coloc", "coloc_summary.csv"))
cl <- cl[grep("CHIP_vs_", cl$pair), ]
cl$region <- sub("_CHIP_vs_.*", "", cl$pair)
cl$pheno  <- sub(".*_CHIP_vs_", "", cl$pair)
cl <- cl[cl$pheno %in% c("MPN", "PV"), ]
cl$region <- factor(cl$region, levels = rev(c("TERT", "C21q22", "C11q22", "C17p13", "JAK2", "C3q25")))
cl$pheno  <- factor(cl$pheno, levels = c("MPN", "PV"))
# 关键修正：H3 优势时不能把 H3 数值画在 H4 色阶上
cl$state  <- ifelse(cl$PP.H4 >= 0.8, "shared",
              ifelse(cl$PP.H3 >= 0.8, "distinct", "unclear"))
cl$fillv  <- ifelse(cl$state == "shared", cl$PP.H4, NA)
cl$lab    <- ifelse(cl$state == "shared", sprintf("%.3f", cl$PP.H4),
              ifelse(cl$state == "distinct", "distinct", sprintf("%.2f?", cl$PP.H4)))

pC <- ggplot(cl, aes(x = pheno, y = region)) +
  geom_tile(data = cl[cl$state == "distinct", ], fill = "#D5D8DC", colour = "white", linewidth = 1.2) +
  geom_tile(data = cl[cl$state == "unclear",  ], fill = "#F4F6F7", colour = "white", linewidth = 1.2) +
  geom_tile(data = cl[cl$state == "shared",   ], aes(fill = fillv), colour = "white", linewidth = 1.2) +
  geom_text(data = cl[cl$state == "distinct", ], aes(label = lab), size = 2.7, colour = "grey35") +
  geom_text(data = cl[cl$state == "unclear",  ], aes(label = lab), size = 2.9, colour = "grey45") +
  geom_text(data = cl[cl$state == "shared",   ], aes(label = lab), size = 3.2, colour = "white",
            fontface = "bold") +
  scale_fill_gradient2(low = "#F6DDCC", mid = "#E59866", high = "#922B21",
                       midpoint = 0.85, limits = c(0.8, 1),
                       name = "PP.H4\n(shared\nvariant)", na.value = "#D5D8DC") +
  labs(x = NULL, y = NULL,
       title = "C  Colocalisation with CHIP",
       subtitle = "grey: distinct variants (PP.H3>0.8);  ? = inconclusive") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, colour = "grey40"),
        axis.text = element_text(colour = "grey15"))

# ---------- Panel D: 敏感性稳定性 ----------
sens <- read.csv(file.path(BASE, "07_mr", "sens", "sens_MR_all.csv"))
sens <- sens[sens$method %in% c("Inverse variance weighted") &
             grepl("^CHIPg", sens$tag) & grepl("(MPN|PV)$|(MPN|PV)_BASE", sens$tag), ]
sens$set <- ifelse(grepl("BASE", sens$tag), "Base",
            ifelse(grepl("noTERT", sens$tag), "No TERT",
            ifelse(grepl("noMHC", sens$tag), "No MHC", "other")))
sens$pheno <- ifelse(grepl("MPN", sens$tag), "MPN", "PV")
sens$set <- factor(sens$set, levels = c("Base", "No TERT", "No MHC"))
sens <- sens[sens$set %in% c("Base", "No TERT", "No MHC"), ]
sens$lo <- exp(sens$b - 1.96 * sens$se); sens$hi <- exp(sens$b + 1.96 * sens$se)

pD <- ggplot(sens, aes(x = set, y = exp(b), colour = pheno, group = pheno)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.65) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c(MPN = "#E67E22", PV = "#C0392B")) +
  labs(x = NULL, y = "OR (IVW, 95% CI)", colour = NULL,
       title = "D  Sensitivity to pleiotropic loci") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))

fig <- (pA / pB) | (pC / pD)
fig <- fig + plot_layout(widths = c(1.25, 1))

ggsave(file.path(FIG, "Fig2_MR_coloc.png"), fig, width = 9.2, height = 7.4, dpi = 400)
ggsave(file.path(FIG, "Fig2_MR_coloc.pdf"), fig, width = 9.2, height = 7.4)
cat("written ->", file.path(FIG, "Fig2_MR_coloc.png"), "\n")
