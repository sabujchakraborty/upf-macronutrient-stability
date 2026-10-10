# =============================================================================
# 02_fdc_extract.R  --  FoodData Central Branded Foods: label versions, first-vs-last pairs, UPF flag
# Uses DuckDB so the 1.4 GB nutrient file is streamed from disk (works on an 8 GB laptop).
# Inputs : branded_food.csv, food.csv, food_nutrient.csv, nutrient.csv  (in FDC_DIR)
# Outputs: fdc_pairs_raw.rds (one row per barcode pair), results_fdc_counts.csv,
#          results_fdc_batch_months_owner.csv, results_fdc_batch_months_category.csv
# Run time: a few minutes (three scans of large CSV files).
# =============================================================================
source("00_setup.R")
cat("\n== 02 FDC extract ==\n")

p   <- function(f) file.path(FDC_DIR, paste0(f, ".csv"))
con <- dbConnect(duckdb()); invisible(dbExecute(con, "SET memory_limit='3GB'"))
q   <- function(sql) dbGetQuery(con, sql)

# ---- 1. small working table: one row per FDC record ---------------------------------------------------------------
dbExecute(con, sprintf("CREATE OR REPLACE TABLE bf AS
  SELECT b.fdc_id, b.gtin_upc, b.branded_food_category, b.market_country, b.data_source,
         b.modified_date, b.available_date, b.discontinued_date,
         (b.ingredients IS NOT NULL AND length(b.ingredients) > 0) AS has_ingredients,
         f.publication_date
  FROM read_csv_auto('%s') b JOIN read_csv_auto('%s') f USING (fdc_id)", p("branded_food"), p("food")))
dbExecute(con, sprintf("CREATE OR REPLACE TABLE bo AS
  SELECT fdc_id, brand_owner FROM read_csv_auto('%s')", p("branded_food")))
n_all <- q("SELECT COUNT(*) AS n_records, COUNT(DISTINCT gtin_upc) AS n_barcodes FROM bf")
print(n_all)
check("FDC records", n_all$n_records, 1993975)
check("distinct barcodes", n_all$n_barcodes, 462636)
print(q("SELECT data_source, COUNT(*) AS n FROM bf GROUP BY 1 ORDER BY n DESC LIMIT 5"))

# ---- 2. one row per label version (barcode + modified_date), US products ----------------------------------------
dbExecute(con, "CREATE OR REPLACE TABLE ver AS
  SELECT * EXCLUDE (rn) FROM (
    SELECT fdc_id, gtin_upc, modified_date, available_date, branded_food_category,
           ROW_NUMBER() OVER (PARTITION BY gtin_upc, modified_date
                              ORDER BY available_date DESC, fdc_id DESC) AS rn
    FROM bf WHERE modified_date IS NOT NULL AND market_country IN ('United States','US'))
  WHERE rn = 1")
n_ver <- q("SELECT COUNT(*) AS n_label_versions, COUNT(DISTINCT gtin_upc) AS n_barcodes FROM ver")
print(n_ver)

# ---- 3. nutrient values for those versions (scans food_nutrient.csv once) -----------------------------------------
dbExecute(con, sprintf("CREATE OR REPLACE TABLE vn AS
  SELECT v.fdc_id,
    MAX(CASE WHEN n.nutrient_id = 1008 THEN n.amount END) AS kcal,
    MAX(CASE WHEN n.nutrient_id = 1003 THEN n.amount END) AS protein,
    MAX(CASE WHEN n.nutrient_id = 1004 THEN n.amount END) AS fat,
    MAX(CASE WHEN n.nutrient_id = 1005 THEN n.amount END) AS carb,
    MAX(CASE WHEN n.nutrient_id = 2000 THEN n.amount END) AS sugars,
    MAX(CASE WHEN n.nutrient_id = 1079 THEN n.amount END) AS fibre,
    MAX(CASE WHEN n.nutrient_id = 1258 THEN n.amount END) AS sat_fat,
    MAX(CASE WHEN n.nutrient_id = 1093 THEN n.amount END) AS sodium
  FROM read_csv_auto('%s') n JOIN ver v ON n.fdc_id = v.fdc_id
  WHERE n.nutrient_id IN (1003,1004,1005,1008,1079,1093,1258,2000)
  GROUP BY v.fdc_id", p("food_nutrient")))

# ---- 4. first vs last label version per barcode (complete energy + macronutrients at both ends) -------------------
dbExecute(con, "CREATE OR REPLACE TABLE pairs AS
  WITH ranked AS (
    SELECT v.gtin_upc, v.fdc_id, v.modified_date, v.branded_food_category,
           ROW_NUMBER() OVER (PARTITION BY v.gtin_upc ORDER BY v.modified_date ASC)  AS rk_first,
           ROW_NUMBER() OVER (PARTITION BY v.gtin_upc ORDER BY v.modified_date DESC) AS rk_last
    FROM ver v JOIN vn ON vn.fdc_id = v.fdc_id
    WHERE vn.kcal IS NOT NULL AND vn.protein IS NOT NULL AND vn.fat IS NOT NULL AND vn.carb IS NOT NULL)
  SELECT f.gtin_upc, f.branded_food_category, f.fdc_id AS fdc_first, l.fdc_id AS fdc_last,
         f.modified_date AS date_first, l.modified_date AS date_last
  FROM ranked f JOIN ranked l ON f.gtin_upc = l.gtin_upc AND f.rk_first = 1 AND l.rk_last = 1
  WHERE f.modified_date < l.modified_date")
n_pairs <- q("SELECT COUNT(*) AS n_pairs, quantile_cont(date_diff('day', date_first, date_last), 0.5) AS med_span_days FROM pairs")
print(n_pairs)

# ---- 5. ultra-processed flag from the ingredient text (six marker groups) -------------------------------------------
dbExecute(con, sprintf("CREATE OR REPLACE TABLE ing AS
  SELECT fdc_id, lower(ingredients) AS ing
  FROM read_csv_auto('%s')
  WHERE fdc_id IN (SELECT fdc_first FROM pairs UNION SELECT fdc_last FROM pairs)", p("branded_food")))
dbExecute(con, "CREATE OR REPLACE TABLE upf AS
  SELECT *, (flavour OR emulsifier_gum OR sweetener_syrup OR modified_ingr OR colour OR preservative_other) AS upf_any
  FROM (SELECT fdc_id,
    regexp_matches(ing, 'flavor|flavour') AS flavour,
    regexp_matches(ing, 'lecithin|diglyceride|polysorbate|cellulose gum|carboxymethylcellulose|xanthan|guar gum|carrageenan|gellan|datem|polyglycerol|sodium stearoyl|locust bean|gum arabic') AS emulsifier_gum,
    regexp_matches(ing, 'high.fructose|corn syrup|glucose syrup|fructose syrup|maltodextrin|dextrose|invert sugar|sucralose|aspartame|acesulfame|saccharin|sorbitol|maltitol|xylitol|erythritol|polydextrose|isomaltulose') AS sweetener_syrup,
    regexp_matches(ing, 'modified (food |corn |tapioca |potato |wheat |rice )?starch|hydroly[sz]ed|protein isolate|protein concentrate|hydrogenated|interesterified|mechanically separated') AS modified_ingr,
    regexp_matches(ing, 'caramel colou?r|artificial colou?r|fd&c|red (#)?40|yellow (#)?[56]|blue (#)?[12]|titanium dioxide|carmine|colou?r added') AS colour,
    regexp_matches(ing, 'sodium benzoate|potassium sorbate|sorbic acid|\\b(bha|bht|tbhq|msg|edta)\\b|sodium nitrite|sodium nitrate|calcium propionate|sodium propionate|sodium phosphate|monosodium glutamate|yeast extract|disodium (inosinate|guanylate)|sodium erythorbate|sulfite|sulphite|propylene glycol|sodium acid pyrophosphate|artificial') AS preservative_other
  FROM ing WHERE ing IS NOT NULL)")
marker_share <- q("SELECT COUNT(*) AS n, ROUND(AVG(flavour::INT),3) AS flavour, ROUND(AVG(emulsifier_gum::INT),3) AS emulsifier_gum,
         ROUND(AVG(sweetener_syrup::INT),3) AS sweetener_syrup, ROUND(AVG(modified_ingr::INT),3) AS modified_ingr,
         ROUND(AVG(colour::INT),3) AS colour, ROUND(AVG(preservative_other::INT),3) AS preservative_other,
         ROUND(AVG(upf_any::INT),3) AS upf_any FROM upf")
print(marker_share)
save_csv(marker_share, "results_fdc_marker_shares.csv")

# ---- 6. analysis table: one row per barcode pair ------------------------------------------------------------------
sel1 <- paste0("x.", NUTS, " AS ", NUTS, "_1", collapse = ", ")
sel2 <- paste0("y.", NUTS, " AS ", NUTS, "_2", collapse = ", ")
dat <- q(sprintf("SELECT p.gtin_upc, p.branded_food_category AS category, p.date_first, p.date_last,
    date_diff('day', p.date_first, p.date_last) / 365.25 AS span_years,
    a.upf_any AS upf_1, b.upf_any AS upf_2,
    (a.flavour OR a.emulsifier_gum OR a.sweetener_syrup OR a.modified_ingr OR a.colour) AS strict_1,
    (b.flavour OR b.emulsifier_gum OR b.sweetener_syrup OR b.modified_ingr OR b.colour) AS strict_2,
    sa.data_source AS src_1, sb.data_source AS src_2, bo.brand_owner AS brand_owner,
    %s, %s
  FROM pairs p
  JOIN upf a ON a.fdc_id = p.fdc_first JOIN upf b ON b.fdc_id = p.fdc_last
  JOIN vn  x ON x.fdc_id = p.fdc_first JOIN vn  y ON y.fdc_id = p.fdc_last
  JOIN bf sa ON sa.fdc_id = p.fdc_first JOIN bf sb ON sb.fdc_id = p.fdc_last
  LEFT JOIN bo ON bo.fdc_id = p.fdc_last", sel1, sel2))

dat$cohort <- ifelse(dat$upf_1 & dat$upf_2, "UPF_both", ifelse(!dat$upf_1 & !dat$upf_2, "nonUPF_both", "switched"))
dat$cohort_strict <- ifelse(dat$strict_1 & dat$strict_2, "UPF_both",
                     ifelse(!dat$strict_1 & !dat$strict_2, "nonUPF_both", "switched"))
dat$src_path <- paste(dat$src_1, "->", dat$src_2)
dat$cl <- ifelse(is.na(dat$brand_owner) | dat$brand_owner == "",
                 paste0("NA_", dat$gtin_upc), toupper(trimws(dat$brand_owner)))      # cluster id = brand owner
cat("pairs:", nrow(dat), "\n"); print(table(dat$cohort))
check("pairs analysed", nrow(dat), 208590)
check("UPF_both / nonUPF_both / switched", c(sum(dat$cohort == "UPF_both"), sum(dat$cohort == "nonUPF_both"),
                                           sum(dat$cohort == "switched")), c(145337, 61635, 1618))

# ---- 7. scale-error screen -----------------------------------------------------------------------------------------
# A pair is flagged when >= 3 of the 4 macronutrient ratios (version 2 / version 1) agree within 3%
# (max/min < 1.03) while energy changed by more than 10% on the log scale.
ratio <- function(a, b) ifelse(a > 0 & b > 0, b / a, NA_real_)
r <- cbind(kcal = ratio(dat$kcal_1, dat$kcal_2), protein = ratio(dat$protein_1, dat$protein_2),
           fat  = ratio(dat$fat_1, dat$fat_2),   carb    = ratio(dat$carb_1, dat$carb_2))
n_ok <- rowSums(!is.na(r))
rmax <- suppressWarnings(apply(r, 1, max, na.rm = TRUE))
rmin <- suppressWarnings(apply(r, 1, min, na.rm = TRUE))
dat$scale_flag <- ((n_ok >= 3) & (rmax / rmin < 1.03) & (abs(log(r[, "kcal"])) > 0.10)) %in% TRUE
cat("flagged pairs:", sum(dat$scale_flag), "of", nrow(dat), "\n")
check("scale-flagged pairs", sum(dat$scale_flag), 374)
print(table(dat$cohort, dat$scale_flag))
flag_cat <- as.data.frame(sort(table(dat$category[dat$scale_flag]), decreasing = TRUE), stringsAsFactors = FALSE)
names(flag_cat) <- c("category", "n_flagged")
print(head(flag_cat, 12))
save_csv(flag_cat, "results_fdc_scale_flag_categories.csv")
saveRDS(dat, out_path("fdc_pairs_raw.rds"))

# ---- 8. counts table for the Methods / supplement -------------------------------------------------------------------
counts <- data.frame(
  item = c("FDC branded records", "distinct barcodes", "US label versions", "barcodes with >= 2 label versions and complete macros (pairs)",
           "pairs with ingredient text (analysed)", "UPF_both", "nonUPF_both", "switched (excluded)", "scale-flagged pairs",
           "median span between versions (days)"),
  value = c(n_all$n_records, n_all$n_barcodes, n_ver$n_label_versions, n_pairs$n_pairs, nrow(dat),
            sum(dat$cohort == "UPF_both"), sum(dat$cohort == "nonUPF_both"), sum(dat$cohort == "switched"),
            sum(dat$scale_flag), n_pairs$med_span_days))
print(counts); save_csv(counts, "results_fdc_counts.csv")

# ---- 9. artefact catalogue: which brand owners / categories drive the batch-update months --------------------------
spike_months <- "DATE '2021-05-01', DATE '2022-02-01', DATE '2022-10-01', DATE '2025-02-01', DATE '2025-07-01'"
dbWriteTable(con, "upf_gtin", data.frame(gtin_upc = unique(dat$gtin_upc[dat$cohort == "UPF_both"])), overwrite = TRUE)
dbExecute(con, sprintf("CREATE OR REPLACE TABLE spike AS
  SELECT * FROM (
    SELECT v.fdc_id, v.gtin_upc, date_trunc('month', v.modified_date) AS month, v.branded_food_category AS category,
           LN(n.kcal)    - AVG(LN(n.kcal))    OVER (PARTITION BY v.gtin_upc) AS kcal_l,
           LN(n.protein) - AVG(LN(n.protein)) OVER (PARTITION BY v.gtin_upc) AS protein_l
    FROM ver v JOIN vn n ON n.fdc_id = v.fdc_id
    WHERE v.gtin_upc IN (SELECT gtin_upc FROM upf_gtin)
      AND n.kcal > 0 AND n.protein > 0 AND n.fat > 0 AND n.carb > 0)
  WHERE month IN (%s)", spike_months))
dbExecute(con, "CREATE OR REPLACE TABLE spike_b AS SELECT s.*, bo.brand_owner FROM spike s JOIN bo USING (fdc_id)")
by_owner <- q("SELECT month, brand_owner, COUNT(*) AS n_versions, ROUND(AVG(kcal_l) * 100, 1) AS mean_kcal_logpts,
               ROUND(AVG(protein_l) * 100, 1) AS mean_protein_logpts
        FROM spike_b GROUP BY month, brand_owner
        QUALIFY ROW_NUMBER() OVER (PARTITION BY month ORDER BY SUM(kcal_l) ASC) <= 4
        ORDER BY month, SUM(kcal_l)")
by_cat <- q("SELECT month, category, COUNT(*) AS n_versions, ROUND(AVG(kcal_l) * 100, 1) AS mean_kcal_logpts
        FROM spike_b GROUP BY month, category
        QUALIFY ROW_NUMBER() OVER (PARTITION BY month ORDER BY SUM(kcal_l) ASC) <= 3
        ORDER BY month, SUM(kcal_l)")
print(by_owner); print(by_cat)
save_csv(by_owner, "results_fdc_batch_months_owner.csv")
save_csv(by_cat,   "results_fdc_batch_months_category.csv")

dbDisconnect(con, shutdown = TRUE)
cat("02_fdc_extract.R done.\n")
