cat("=== R version ===\n"); print(R.version.string)
pk <- c("TwoSampleMR","ieugwasr","MendelianRandomization","gwasglue","data.table","dplyr","MRPRESSO","RadialMR","coloc","locuscomparer")
for (p in pk) {
  v <- tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED")
  cat(sprintf("%-24s %s\n", p, v))
}
cat("\n=== OPENGWAS_JWT set? ===\n")
cat(Sys.getenv("OPENGWAS_JWT") == "" , "\n")
cat("\n=== TwoSampleMR available_outcomes test (offline?) ===\n")
r <- tryCatch({
  library(ieugwasr)
  x <- api_query("gwasinfo/ieu-b-4965", access_token=ieugwasr::get_opengwas_jwt())
  head(x, 2)
}, error=function(e) paste("ERR:", conditionMessage(e)))
print(r)
cat("\nDONE\n")
