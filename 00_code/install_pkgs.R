options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
pk <- c("coloc", "data.table", "ggrepel", "cowplot")
for (p in pk) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("installing", p, "...\n")
    tryCatch(install.packages(p, dependencies = TRUE), error = function(e) cat("FAIL", p, conditionMessage(e), "\n"))
  } else {
    cat("present:", p, "\n")
  }
}
for (p in pk) {
  cat(sprintf("%-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
}
