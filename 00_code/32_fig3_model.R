## ============================================================
## Fig 3  Mechanism model
##   A  three convergent germline pathways -> CHIP   (coloc PP.H4)
##   B  bifurcation: one germline cause -> four phenotypes  (IVW OR)
##   C  JAK2 as an independent somatic route
## Sources: 07_mr/, 09_coloc/, 05_iv/CHIP_germline_iv_annotated.tsv
## ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid)
})

BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
OUT  <- file.path(BASE, "08_figures")

COL  <- c(PV = "#C0392B", MPN = "#E67E22", MDS = "#8E44AD",
          ITP = "#27AE60", AA = "#2980B9")
cTEL <- "#1F6F8B"; cDDR <- "#6C3483"; cHSC <- "#B9770E"; cCHIP <- "#2C3E50"

bx <- function(x1, x2, y1, y2, fill = "white", col = "grey35", lw = 0.5,
               alpha = 1, lty = "solid")
  annotate("rect", xmin = x1, xmax = x2, ymin = y1, ymax = y2,
           fill = fill, colour = col, linewidth = lw, alpha = alpha, linetype = lty)

ar <- function(x, xe, y, ye, col = "grey30", lw = 1.4, asz = 2.6)
  annotate("segment", x = x, xend = xe, y = y, yend = ye, colour = col,
           linewidth = lw, arrow = arrow(length = unit(asz, "mm"), type = "closed"))

crv <- function(x, xe, y, ye, col = "grey40", lw = 1.1, asz = 2.4, curv = 0.22)
  annotate("curve", x = x, xend = xe, y = y, yend = ye, colour = col,
           linewidth = lw, curvature = curv, angle = 90, ncp = 60,
           arrow = arrow(length = unit(asz, "mm"), type = "closed"))

tx <- function(x, y, lab, size = 2.6, col = "grey15", face = 1,
               hjust = 0.5, vjust = 0.5, lh = 0.95)
  annotate("text", x = x, y = y, label = lab, size = size, colour = col,
           fontface = face, hjust = hjust, vjust = vjust, lineheight = lh)

canv <- function(p) p +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2),
        plot.background = element_rect(fill = "white", colour = NA))

## ============================================================
## A  convergent germline architecture
## ============================================================
pnl <- data.frame(
  y   = c(72, 46, 20),
  col = c(cTEL, cDDR, cHSC),
  nm  = c("Telomere maintenance", "DNA damage response", "HSC transcriptional control"),
  sub = c("replicative lifespan  \u00b7  longevity",
          "PARP1\u2013ATM\u2013p53  \u00b7  genome stability",
          "self-renewal set-point  \u00b7  lineage priming"),
  ttl = c("TERT", "PARP1 \u00b7 ATM \u00b7 TP53", "RUNX1 \u00b7 MSI2 \u00b7 CD164"),
  gns = c("5p15.33", "1q42 \u00b7 11q22 \u00b7 17p13", "21q22 \u00b7 17q22 \u00b7 6q21"),
  ph4 = c("PP.H4 = 0.998", "ATM 0.969  \u00b7  TP53 0.875", "PP.H4 = 0.984 (RUNX1)"),
  stringsAsFactors = FALSE
)

pA <- ggplot()
for (i in seq_len(nrow(pnl))) {
  r <- pnl[i, ]
  pA <- pA +
    bx(0.5, 34.5, r$y - 8, r$y + 8, col = r$col, lw = 0.6) +
    tx(17.5, r$y + 4.9, r$ttl, size = 3.1, face = 2, col = r$col) +
    tx(17.5, r$y - 1.3, r$gns, size = 2.25, col = "grey45") +
    tx(17.5, r$y - 5.2, r$ph4, size = 2.35, col = "grey20", face = 3) +
    ar(36, 62.5, r$y, r$y, col = r$col, lw = 1.9, asz = 3.0) +
    tx(47.6, r$y + 6.2, r$nm,  size = 3.15, face = 2, col = r$col) +
    tx(47.6, r$y - 5.6, r$sub, size = 2.15, col = "grey40")
}
pA <- pA +
  annotate("segment", x = 63.1, xend = 63.1, y = 20, yend = 72,
           colour = "grey65", linewidth = 0.9) +
  ar(63.3, 64.7, 46, 46, col = cCHIP, lw = 2.4, asz = 3.2) +
  bx(65, 95, 8, 84, fill = cCHIP, col = cCHIP, lw = 0.6) +
  tx(80, 74, "CHIP", size = 6.2, face = 2, col = "white") +
  tx(80, 65, "clonal expansion of a", size = 2.65, col = "grey90") +
  tx(80, 61, "mutated HSC", size = 2.65, col = "grey90") +
  annotate("segment", x = 69, xend = 91, y = 53, yend = 53, colour = "grey70", linewidth = 0.35) +
  tx(80, 46, "UK Biobank", size = 2.35, col = "grey80") +
  tx(80, 41, "418,401 adults", size = 2.35, col = "grey80") +
  tx(80, 36, "25,657 cases (6.1%)", size = 2.35, col = "grey80") +
  annotate("segment", x = 69, xend = 91, y = 28, yend = 28, colour = "grey70", linewidth = 0.35) +
  tx(80, 21, "Kessler et al.", size = 2.25, col = "grey75") +
  tx(80, 17, "Nature 2022", size = 2.25, col = "grey75") +
  tx(0.5, 97, "A   Convergent germline architecture linking CHIP to myeloproliferation",
     size = 4.3, face = 2, col = "grey10", hjust = 0) +
  tx(0.5, 3.0,
     paste0("PP.H4 = posterior probability of a shared causal variant with MPN (coloc, \u00b1 500 kb); ",
            "each pathway carries at least one locus with PP.H4 > 0.8."),
     size = 2.15, col = "grey45", hjust = 0)
pA <- canv(pA)

## ============================================================
## B  bifurcation into four phenotypes
## ============================================================
br <- data.frame(
  y   = c(80, 57.5, 35, 12.5),
  key = c("PV", "MDS", "ITP", "AA"),
  nm  = c("MPN  \u00b7  PV / ET", "MDS", "ITP", "AA"),
  mech= c("somatic driver + proliferative signalling",
          "ineffective erythropoiesis  \u00b7  spliceosome lesions",
          "peripheral platelet destruction, no HSC lesion",
          "telomere attrition + immune destruction"),
  or  = c("OR 3.48 (PV) \u00b7 2.93 (MPN)", "OR 1.34  (P = 0.037)",
          "OR 0.95  (P = 0.70)", "OR 0.84  (P = 0.35)  reversed"),
  stringsAsFactors = FALSE
)

pB <- ggplot() +
  bx(0.5, 18, 34, 66, fill = cCHIP, col = cCHIP, lw = 0.6) +
  tx(9.2, 53, "CHIP", size = 4.6, face = 2, col = "white") +
  tx(9.2, 45, "clone in", size = 2.3, col = "grey85") +
  tx(9.2, 41, "aged HSC", size = 2.3, col = "grey85")
for (i in seq_len(nrow(br))) {
  r <- br[i, ]
  cv <- ifelse(r$y > 50, -0.30, 0.30)
  pB <- pB +
    crv(18.6, 29.4, 50, r$y, col = COL[[r$key]], lw = 1.2, asz = 2.3, curv = cv) +
    bx(30, 99, r$y - 7.5, r$y + 7.5, col = COL[[r$key]], lw = 0.55) +
    tx(32.2, r$y + 3.6, r$nm, size = 3.2, face = 2, col = COL[[r$key]], hjust = 0) +
    tx(32.2, r$y - 3.4, r$mech, size = 2.0, col = "grey45", hjust = 0, face = 3) +
    tx(97.2, r$y + 3.6, r$or, size = 2.25, col = COL[[r$key]], hjust = 1, face = 2)
}
pB <- pB +
  tx(0.5, 99, "B   One germline cause, divergent phenotypes",
     size = 3.5, face = 2, col = "grey10", hjust = 0) +
  tx(30.2, 1.2,
     "MPN risk loci reverse direction in AA (17 / 22; r = \u22120.50, P = 0.011)",
     size = 2.0, col = "grey50", hjust = 0, face = 3)
pB <- canv(pB)

## ============================================================
## C  JAK2 parallel route
## ============================================================
pC <- ggplot() +
  tx(0.5, 99, "C   JAK2: a parallel route, not the shared cause",
     size = 3.15, face = 2, col = "grey10", hjust = 0) +
  bx(1, 99, 78, 92, fill = "#EAF2F8", col = cTEL, lw = 0.55) +
  tx(50, 87, "Germline JAK2 46/1 haplotype", size = 2.5, face = 2, col = cTEL) +
  tx(50, 81.5, "rs1327494  \u00b7  MPN P = 1.1e-107", size = 2.0, col = "grey40") +
  ar(50, 50, 77, 70.5, col = "grey45", lw = 1.2, asz = 2.2) +
  tx(52.5, 74.0, "second hit", size = 1.9, col = "grey50", hjust = 0) +
  bx(1, 99, 56, 70, fill = "#FDEDEC", col = COL[["PV"]], lw = 0.55) +
  tx(50, 65, "Somatic JAK2 V617F", size = 2.5, face = 2, col = COL[["PV"]]) +
  tx(50, 59.5, "rs77375493  \u00b7  acquired in HSC", size = 2.0, col = "grey40") +
  ar(50, 50, 55, 50.5, col = COL[["PV"]], lw = 1.4, asz = 2.4) +
  bx(1, 99, 37, 50, fill = COL[["PV"]], col = COL[["PV"]], lw = 0.55) +
  tx(50, 45, "PV / ET", size = 3.3, face = 2, col = "white") +
  tx(50, 39.5, "clonal myeloproliferation", size = 2.1, col = "white") +
  bx(1, 99, 21, 33, fill = "#F4F6F7", col = "grey55", lw = 0.5) +
  tx(50, 29, "coloc at JAK2 (\u00b1 500 kb): PP.H3 = 1.00", size = 2.0, col = "grey20") +
  tx(50, 24.5, "CHIP (V617F readout) and MPN (46/1)\nare distinct causal variants",
     size = 2.0, col = "grey40", lh = 1.15) +
  annotate("rect", xmin = 1, xmax = 99, ymin = 3, ymax = 17,
           fill = "#FDF2E9", colour = cHSC, linewidth = 0.5, linetype = "21") +
  tx(50, 12.5, "0 of 21 CHIP instruments", size = 2.9, face = 2, col = "#7E5109") +
  tx(50, 7.2, "fall in chr9:4.8\u20135.4 Mb \u2192 no 46/1 pleiotropy",
     size = 1.95, col = "#7E5109", face = 3)
pC <- canv(pC)

## ============================================================
fig <- pA / (pB | pC) + plot_layout(heights = c(1.0, 1.12), widths = c(1.42, 1))

ggsave(file.path(OUT, "Fig3_mechanism_model.png"), fig,
       width = 180, height = 182, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(OUT, "Fig3_mechanism_model.pdf"), fig,
       width = 180, height = 182, units = "mm", device = cairo_pdf, bg = "white")
cat("written\n")
