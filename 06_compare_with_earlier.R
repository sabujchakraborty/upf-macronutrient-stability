# =============================================================================
# 06_compare_with_earlier.R  --  reproducibility check
# Compares the numbers this pipeline produced (OUT_DIR) with the results CSVs from the original
# step-by-step analysis (WORK_DIR). Differences above the tolerance are listed.
# =============================================================================
source("00_setup.R")
cat("\n== 06 Compare with earlier results ==\n")

files <- c("results_fndds_paired_tests_3sets.csv", "results_ilr_paired.csv", "results_fdc_two_part.csv",
           "results_fdc_two_part_off_confirmed.csv", "results_fdc_primary_cluster_robust.csv",
           "results_classifier_validation_off.csv", "results_fdc_gdsn_zero_patterns.csv",
           "results_fdc_primary_cluster_robust_off_confirmed.csv")

cmp_one <- function(f, tol = 1e-6) {
  old_f <- file.path(WORK_DIR, f); new_f <- out_path(f)
  if (!file.exists(old_f)) return(data.frame(file = f, status = "no earlier file to compare"))
  if (!file.exists(new_f)) return(data.frame(file = f, status = "NOT produced by pipeline"))
  o <- read.csv(old_f, check.names = FALSE); n <- read.csv(new_f, check.names = FALSE)
  common <- intersect(names(o), names(n)); common <- common[vapply(common, function(cn) is.numeric(o[[cn]]) && is.numeric(n[[cn]]), logical(1))]
  if (nrow(o) != nrow(n)) return(data.frame(file = f, status = sprintf("row count differs (%d vs %d)", nrow(o), nrow(n))))
  worst <- 0; where <- ""
  for (cn in common) {
    d <- abs(sort(o[[cn]]) - sort(n[[cn]]))
    d <- d / pmax(1, abs(sort(o[[cn]])))                       # relative for large values, absolute for small
    if (any(!is.na(d)) && max(d, na.rm = TRUE) > worst) { worst <- max(d, na.rm = TRUE); where <- cn }
  }
  data.frame(file = f, status = if (worst <= tol) "matches" else paste0("DIFFERS in ", where),
             max_rel_diff = signif(worst, 3))
}
res <- do.call(rbind, lapply(files, cmp_one))
print(res, row.names = FALSE)
save_csv(res, "reproducibility_check.csv")
cat("06_compare_with_earlier.R done.\n")
