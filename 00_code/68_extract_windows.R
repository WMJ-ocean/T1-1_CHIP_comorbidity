# 68_extract_windows.R
# 从 eQTLGen 全量 cis-eQTL (3.88 GB) 中，抽取 21 个胚系 IV 周围 ±1 Mb 窗口内
# 的全部 SNP-基因对（含零效应者）—— 这是做标准 coloc.abf 所必需的。
# 注意：该文件是多成员 gzip，必须顺序读完；用 data.table::fread + cmd 流式(stderr 命令) 处理会失败，
#       这里用纯 R 的 gzfile() 连接逐块读取（R 内部用 zlib，能正确处理多成员）。
suppressMessages(library(data.table))

BASE <- dirname(dirname(parent.frame(1)$ofile %||% normalizePath(".")))
if (is.null(BASE) || is.na(BASE)) {
  BASE <- "C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity"
}
f_full  <- file.path(BASE, "02_raw", "eqtlgen", "eQTLGen_cis_eQTL_full.txt.gz")
f_iv    <- file.path(BASE, "05_iv", "CHIP_germline_iv_annotated.tsv")
f_out   <- file.path(BASE, "06_targets", "eqtlgen_windows_full.tsv")
f_log   <- file.path(BASE, "01_log", "68_extract_windows.txt")

cat("输入:", f_full, "\n")
cat("大小:", round(file.info(f_full)$size / 1e6, 1), "MB\n\n")

IV <- fread(f_iv, sep = "\t", header = TRUE)
cat("IV 表列:", paste(names(IV), collapse = ", "), "\n")
cat("IV 行數:", nrow(IV), "\n\n")
stopifnot("SNP" %in% names(IV))
IVS <- unique(IV$SNP)
cat("唯一 IV:", length(IVS), "个\n")

con <- gzfile(f_full, "rt")
hdr <- readLines(con, n = 1)
cat("表头:", hdr, "\n")
nm  <- strsplit(hdr, "\t")[[1]]
cat("列数:", length(nm), "\n")
i_snp  <- which(nm == "SNP")[1]
i_gene <- which(nm == "GeneSymbol")[1]
i_chr  <- which(nm == "SNPChr")[1]
i_pos  <- which(nm == "SNPPos")[1]
i_p    <- which(nm == "Pvalue")[1]
i_z    <- which(nm == "Zscore")[1]
i_gt   <- which(nm == "GeneChr")[1]
cat(sprintf("索引 SNP=%d GeneSymbol=%d Chr=%d Pos=%d P=%d Z=%d\n",
            i_snp, i_gene, i_chr, i_pos, i_p, i_z))

# 第一遍：定位每个 IV 的位置
cat("\n--- 第一遍：扫描全表定位 IV 坐标 ---\n")
pos_tab <- list()
nrow_total <- 0L
buf_lines <- 200000L
repeat {
  ln <- readLines(con, n = buf_lines)
  if (length(ln) == 0) break
  nrow_total <- nrow_total + length(ln)
  parts <- strsplit(ln, "\t", fixed = TRUE)
  snp   <- vapply(parts, function(p) if (length(p) >= max(c(i_snp, i_chr, i_pos))) p[i_snp] else NA_character_,
                  FUN.VALUE = character(1))
  idx <- which(snp %in% IVS)
  if (length(idx)) {
    for (k in idx) {
      p <- parts[[k]]
      pos_tab[[p[i_snp]]] <- list(chr = suppressWarnings(as.integer(p[i_chr])),
                                  pos = suppressWarnings(as.numeric(p[i_pos])))
    }
  }
  if (nrow_total %% 20000000 == 0)
    cat("  已扫", format(nrow_total, big.mark = ","), "行 ...\n")
}
close(con)
cat("总行数:", format(nrow_total, big.mark = ","), "\n")
cat("定位到 IV:", length(pos_tab), "个\n")

if (length(pos_tab) == 0) stop("未定位到任何 IV，请检查 SNP 列")

PT <- rbindlist(lapply(names(pos_tab), function(s)
  data.table(SNP = s, chr = pos_tab[[s]]$chr, pos = pos_tab[[s]]$pos)))
cat("\n已定位 IV 坐标:\n"); print(PT[order(chr, pos)])

# 第二遍：按窗口抽取
cat("\n--- 第二遍：抽取 ±1 Mb 窗口 ---\n")
WIN <- 1e6
res <- list(); nmade <- 0L
con <- gzfile(f_full, "rt")
invisible(readLines(con, n = 1))
n_row <- 0L
repeat {
  ln <- readLines(con, n = buf_lines)
  if (length(ln) == 0) break
  n_row <- n_row + length(ln)
  parts <- strsplit(ln, "\t", fixed = TRUE)
  keep_idx <- integer(0)
  for (j in seq_len(nrow(PT))) {
    want_chr <- PT$chr[j]; lo <- PT$pos[j] - WIN; hi <- PT$pos[j] + WIN
    hit <- vapply(parts, function(p) {
      if (length(p) < max(c(i_chr, i_pos))) return(FALSE)
      suppressWarnings(as.integer(p[i_chr]) == want_chr &&
                         as.numeric(p[i_pos]) >= lo && as.numeric(p[i_pos]) <= hi)
    }, FUN.VALUE = logical(1))
    hit[is.na(hit)] <- FALSE
    keep_idx <- union(keep_idx, which(hit))
  }
  if (length(keep_idx)) { nmade <- nmade + length(keep_idx); res[[length(res) + 1]] <- ln[keep_idx] }
  if (n_row %% 20000000 == 0) cat("  已扫", format(n_row, big.mark = ","), "行，累计抽出", nmade, "\n")
}
close(con)
cat("第二遍完成，抽取行数:", nmade, "\n")

if (nmade > 0) {
  out <- rbindlist(lapply(unlist(res, use.names = FALSE), function(L) {
    as.data.table(t(strsplit(L, "\t", fixed = TRUE)[[1]]))
  }), fill = TRUE)
  names(out) <- nm[seq_len(min(ncol(out), length(nm)))]
  fwrite(out, f_out, sep = "\t")
  cat("\n写出:", f_out, "\n行数:", nrow(out), "基因数:", uniqueN(out$GeneSymbol), "\n")
} else cat("\n未抽到任何行，请检查窗口或坐标。\n")

sink(f_log)
cat("rows_scanned:", nrow_total, "\nrows_extracted:", nmade, "\n")
sink()
