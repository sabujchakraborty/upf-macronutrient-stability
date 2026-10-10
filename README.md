# Reproducible pipeline: macronutrient stability in ultra-processed (NOVA 4) foods

R scripts that rebuild every table, figure and statistic in the manuscript from the raw files.
They consolidate the step-by-step blocks run during the original analysis.

## Before you run

1. Put all files from this folder into your work folder
   (`E:/MSc Docs/Mini Project/KNIME output+work`), next to the FNDDS exports and the NOVA lookup.
2. Open `00_setup.R` and check the three paths (`WORK_DIR`, `FDC_DIR`, `OUT_DIR`).
3. Packages (install once): `install.packages(c("readr","dplyr","readxl","stringr","DBI","duckdb"), repos = "https://cloud.r-project.org")`
4. Inputs needed
   - `fndds_foods.csv`, `fndds_nutrients.csv` (KNIME exports), `nova_lookup_martinez_steele_2017_2018.xlsx`
   - FoodData Central Branded Foods CSVs (`branded_food.csv`, `food.csv`, `food_nutrient.csv`, `nutrient.csv`) in `FDC_DIR`
   - Internet on the first run of step 03 (Open Food Facts; only two columns are fetched, then cached)

## Run

`source("run_all.R")` runs everything. Or run the scripts one at a time, in order; each saves its output.

| Step | Script | What it does | Main outputs (in `pipeline_output/`) |
|---|---|---|---|
| 0 | `00_setup.R` | paths, packages, helper functions (SPD, two-part estimand, cluster-robust inference, ilr) | none |
| 1 | `01_fndds.R` | FNDDS NOVA 4 pairing, sets A/B/C, paired / Welch / Wilcoxon tests, ilr + Hotelling, zero sensitivity | `results_fndds_paired_tests_3sets.csv`, `results_ilr_*.csv` |
| 2 | `02_fdc_extract.R` | Branded Foods label versions, first-vs-last pairs, UPF flag, scale-error screen, batch-month catalogue (DuckDB) | `fdc_pairs_raw.rds`, `results_fdc_counts.csv`, `results_fdc_batch_months_*.csv` |
| 3 | `03_off_validation.R` | classifier validation against Open Food Facts NOVA; builds the analysis table | `results_classifier_validation_off.csv`, `fdc_dd_final.rds` |
| 4 | `04_fdc_analysis.R` | source strata, two-part estimand, brand-clustered primary analysis, OFF-confirmed and GDSN sensitivity | `results_fdc_primary_cluster_robust*.csv`, `results_fdc_two_part*.csv` |
| 5 | `05_tables_figures.R` | formatted manuscript tables, Figure 1 (PDF + 600 dpi PNG) | `table_*.csv`, `fig_fdc_forest.*` |
| 6 | `06_compare_with_earlier.R` | compares pipeline output with the earlier results files | `reproducibility_check.csv` |

Slow steps: 02 (about three scans of the 1.4 GB nutrient file and the ingredient file) and the first run of 03 (download).

## Built-in checks

Each script prints `[OK]` / `[CHECK]` lines comparing key counts with the values from the original analysis
(for example 12,707 FNDDS rows, 1,118 / 1,101 / 1,060 pairs per set, 208,590 FDC pairs, 374 scale-flagged pairs,
kappa 0.845). Any `[CHECK]` means something differs and should be looked at before using the results.

## Analysis plan (frozen after exploration; disclosed in the paper)

- Primary: FNDDS paired tests on set A; Branded Foods LI -> LI, full cohorts, brand-owner cluster-robust, equivalence +/-1% SPD.
- Sensitivity: FNDDS sets B and C and the ilr coordinates; the GDSN -> GDSN stratum; OFF-confirmed cohorts; margins +/-0.5% and +/-2%.
- Decided after seeing data: the scale-error screen, stratification by data source, the two-part (zero status + magnitude) estimand.
- Not part of inference: Mann-Kendall / Sen's slope / Bai-Perron on monthly indices (driven by batch updates), and the capped
  percent-deviation and winsorised relative-change indices (asymmetric). They are intentionally absent from this pipeline.

## Notes

- The Open Food Facts NOVA label is itself derived from ingredient and additive lists: agreement is convergent validity, not ground truth.
- `NUT` abbreviations: kcal = energy, carb = carbohydrate by difference, sat_fat = saturated fat, SPD = symmetric percent difference.
- FDC nutrient IDs: 1008 energy, 1003 protein, 1004 fat, 1005 carbohydrate, 1079 fibre, 1093 sodium, 1258 saturated fat, 2000 total sugars.

## Data sources (not redistributed here)

Download these yourself and place them as described above.

- FNDDS 2017-2018 and 2019-2020 and FoodData Central Branded Foods (December 2025 release): https://fdc.nal.usda.gov/download-datasets
- Open Food Facts Parquet export (accessed 9 October 2026): https://world.openfoodfacts.org/data
- NOVA lookup of Martinez Steele et al. (2017-2018): obtain from the authors' publication or supplement and respect its terms.

## Citing

Use the Zenodo DOI of the tagged release (see `CITATION.cff`). Software: R 4.6.1, duckdb 1.5.6, DBI 1.3.0, dplyr 1.2.1, readr 2.2.0, readxl 1.5.0.1, stringr 1.6.0.
