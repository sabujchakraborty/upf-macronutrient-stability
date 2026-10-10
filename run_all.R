# =============================================================================
# run_all.R  --  run the whole pipeline in order
# 1. Edit the three paths at the top of 00_setup.R.
# 2. Put the scripts in WORK_DIR (or setwd() to the folder that holds them) and run:  source("run_all.R")
# Each step writes its results to WORK_DIR/pipeline_output. Steps 02 and 03 are slow (large files, first OFF download);
# you can re-run later steps on their own because every step saves its output.
# =============================================================================
steps <- c("01_fndds.R", "02_fdc_extract.R", "03_off_validation.R", "04_fdc_analysis.R",
           "05_tables_figures.R", "06_compare_with_earlier.R")
script_dir <- getwd()                       # folder holding the scripts
for (stp in steps) {
  t0 <- Sys.time()
  cat("\n################ ", stp, " ################\n")
  setwd(script_dir)
  source(file.path(script_dir, stp), local = FALSE)
  setwd(script_dir)
  cat(sprintf("%s finished in %.1f min\n", stp, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
sink(file.path(OUT_DIR, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("\nAll steps finished. Results are in", OUT_DIR, "\n")
