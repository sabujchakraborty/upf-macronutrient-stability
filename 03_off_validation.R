# =============================================================================
# 03_off_validation.R  --  validate the ingredient-keyword UPF classifier against Open Food Facts NOVA
# Inputs : fdc_pairs_raw.rds (from 02), Open Food Facts Parquet export (read remotely; only 2 columns are fetched)
# Outputs: off_nova.rds (cache), results_classifier_validation_off.csv, results_classifier_crosstab_off.csv,
#          results_off_counts.csv, fdc_dd_final.rds (analysis table used by 04 and 05)
# Needs internet on the first run. OFF NOVA is itself derived from ingredient/additive lists, so agreement is
# convergent validity, not ground truth.
# =============================================================================
source("00_setup.R")
cat("\n== 03 OFF validation ==\n")
dat <- readRDS(out_path("fdc_pairs_raw.rds"))

# ---- 1. OFF barcodes with a NOVA group (cached after the first download) --------------------------------------------
off_cache <- out_path("off_nova.rds")
if (file.exists(off_cache)) {
  off <- readRDS(off_cache)
} else {
  con2 <- dbConnect(duckdb())
  dbExecute(con2, "INSTALL httpfs"); dbExecute(con2, "LOAD httpfs")
  cnt <- dbGetQuery(con2, sprintf("SELECT COUNT(*) AS n_products, COUNT(nova_group) AS n_with_nova
                                   FROM read_parquet('%s')", OFF_URL))
  print(cnt); save_csv(cnt, "results_off_counts.csv")
  off <- dbGetQuery(con2, sprintf("SELECT code, nova_group FROM read_parquet('%s')
                                   WHERE nova_group IS NOT NULL AND code IS NOT NULL", OFF_URL))
  dbDisconnect(con2, shutdown = TRUE)
  saveRDS(off, off_cache)
}
cat("OFF rows with NOVA:", nrow(off), "\n")
off$key <- norm_barcode(off$code)
off <- off[!duplicated(off$key), ]

# ---- 2. classifier agreement (all pairs, before the scale screen, as in the original analysis) ----------------------
v <- dat[dat$cohort %in% c("UPF_both", "nonUPF_both"), c("gtin_upc", "cohort")]
v$key <- norm_barcode(v$gtin_upc)
v$off_nova <- off$nova_group[match(v$key, off$key)]
m <- v[!is.na(v$off_nova), ]
cat("\nMatch rate (%):\n"); print(round(100 * table(m$cohort) / table(v$cohort), 1))
xt <- table(m$cohort, m$off_nova)
cat("\nOur cohort x OFF NOVA group:\n"); print(xt)
save_csv(cbind(our_cohort = rownames(xt), as.data.frame.matrix(xt)), "results_classifier_crosstab_off.csv")

ref  <- m$off_nova == 4
pred <- m$cohort == "UPF_both"
tp <- as.numeric(sum(pred & ref)); fp <- as.numeric(sum(pred & !ref))
fn <- as.numeric(sum(!pred & ref)); tn <- as.numeric(sum(!pred & !ref))
n <- tp + fp + fn + tn; po <- (tp + tn) / n
pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
metrics <- data.frame(n = n, sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
                      PPV = tp / (tp + fp), NPV = tn / (tn + fn), accuracy = po, kappa = (po - pe) / (1 - pe))
print(metrics, digits = 3)
check("matched products", n, 153393); check("kappa", round(metrics$kappa, 3), 0.845)
k <- m$off_nova %in% c(1, 4)
cat("NOVA 1 vs 4 only: n =", sum(k), "agreement =", round(mean((m$cohort[k] == "UPF_both") == (m$off_nova[k] == 4)), 3), "\n")
save_csv(metrics, "results_classifier_validation_off.csv")

# ---- 3. analysis table: drop scale-flagged pairs, keep the two stable cohorts, attach OFF NOVA -----------------------
dd <- dat[!dat$scale_flag & dat$cohort %in% c("UPF_both", "nonUPF_both"), ]
dd$key <- norm_barcode(dd$gtin_upc)
dd$off_nova <- off$nova_group[match(dd$key, off$key)]
cat("analysis pairs:", nrow(dd), "| with OFF NOVA:", sum(!is.na(dd$off_nova)), "\n")
print(table(dd$cohort, dd$src_path))
saveRDS(dd, out_path("fdc_dd_final.rds"))
cat("03_off_validation.R done.\n")
