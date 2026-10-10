# =============================================================================
# 01_fndds.R  --  FNDDS 2017-18 vs 2019-20: NOVA 4 pairing, paired tests, ilr
# Inputs : fndds_foods.csv, fndds_nutrients.csv (exported from KNIME), NOVA lookup xlsx
# Outputs: results_fndds_paired_tests_3sets.csv, results_ilr_paired.csv,
#          results_ilr_hotelling.csv, results_ilr_zero_sensitivity.csv
# =============================================================================
source("00_setup.R")
cat("\n== 01 FNDDS ==\n")

# ---- 1. load (read as text so food codes are never altered) -------------------------------------------------------
read_flex <- function(f) {                      # true CSV, or an Excel workbook saved with a .csv extension
  x <- tryCatch(read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) x <- read_excel(f, col_types = "text")
  x
}
foods  <- read_flex("fndds_foods.csv")
nutr   <- read_flex("fndds_nutrients.csv")
lookup <- read_excel(LOOKUP_FILE, sheet = "nova_lookup")
nutr   <- nutr %>% mutate(across(-(1:6), as.numeric))
cat("Checks against the values obtained in the original analysis:\n")
check("foods rows x cols", dim(foods),  c(12707, 13))
check("nutrients rows x cols", dim(nutr), c(12707, 71))
check("lookup rows x cols", dim(lookup), c(1861, 5))

# ---- 2. join, NOVA assignment ------------------------------------------------------------------------------------
master <- nutr %>%
  left_join(foods %>% select(`Food code`, cycle, final_food_code), by = c("Food code", "cycle"))
names(master) <- gsub("\n", " ", names(master))
check("rows in master", nrow(master), 12707)
check("missing final_food_code", sum(is.na(master$final_food_code)), 0)
check("rows 2017-2018", sum(master$cycle == "2017-2018"), 7083)
check("rows 2019-2020", sum(master$cycle == "2019-2020"), 5624)

lk <- lookup %>% select(match_key, nova_group, nova_subgroup) %>% distinct(match_key, .keep_all = TRUE)
master <- master %>%
  mutate(match_key = str_trim(str_to_lower(`Main food description`))) %>%
  left_join(lk, by = "match_key") %>%
  mutate(nova_group = as.integer(nova_group),
         nova_group_final = case_when(
           str_detect(match_key, "^infant formula") & is.na(nova_group) ~ 4L,   # unlabelled infant formula -> NOVA 4
           TRUE ~ nova_group))
check("rows after NOVA join", nrow(master), 12707)

# ---- 3. group status per cycle and final_food_code ---------------------------------------------------------------
groups <- master %>%
  group_by(cycle, final_food_code) %>%
  summarise(n_foods   = n(),
            n_missing = sum(is.na(nova_group_final)),
            nmin = ifelse(all(is.na(nova_group_final)), NA_integer_, min(nova_group_final, na.rm = TRUE)),
            nmax = ifelse(all(is.na(nova_group_final)), NA_integer_, max(nova_group_final, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(group_status = case_when(is.na(nmin)   ~ "unclassified",
                                  nmin != nmax  ~ "mixed",
                                  n_missing > 0 ~ "partial",
                                  TRUE          ~ "full"))
print(table(groups$cycle, groups$group_status))

nova4  <- groups %>% filter(nmin == 4, group_status != "mixed")
paired <- intersect(nova4$final_food_code[nova4$cycle == "2017-2018"],
                    nova4$final_food_code[nova4$cycle == "2019-2020"])
check("NOVA 4 codes paired across cycles", length(paired), 1118)

# ---- 4. collapse merged codes and line the cycles up ----------------------------------------------------------------
macros <- c("Energy (kcal)", "Protein (g)", "Carbohydrate (g)", "Total Fat (g)", "Sugars, total (g)",
            "Fiber, total dietary (g)", "Fatty acids, total saturated (g)", "Sodium (mg)", "Water (g)")
collapsed <- master %>%
  filter(final_food_code %in% paired) %>%
  group_by(cycle, final_food_code) %>%                                  # one row per cycle and code: mean of members
  summarise(across(all_of(macros), ~ mean(.x, na.rm = TRUE)), n_foods = n(), .groups = "drop") %>%
  left_join(groups %>% select(cycle, final_food_code, group_status), by = c("cycle", "final_food_code"))
a <- collapsed %>% filter(cycle == "2017-2018") %>% arrange(final_food_code)
b <- collapsed %>% filter(cycle == "2019-2020") %>% arrange(final_food_code)
stopifnot(identical(a$final_food_code, b$final_food_code))

strict     <- a$group_status == "full" & b$group_status == "full"
one_to_one <- a$n_foods == 1 & b$n_foods == 1
sets <- list(A_primary_full = strict,                 # primary: full in both cycles
             B_with_partial = rep(TRUE, nrow(a)),     # sensitivity: adds partially classified codes
             C_one_to_one   = strict & one_to_one)    # sensitivity: no averaging across different foods
cat("pairs per set:\n"); print(sapply(sets, sum))
check("set A / B / C", unname(sapply(sets, sum)), c(1101, 1118, 1060))

# ---- 5. paired tests -------------------------------------------------------------------------------------------------
res3 <- do.call(rbind, lapply(names(sets), function(s) run_tests(a[sets[[s]], ], b[sets[[s]], ], macros, s)))
base_mean <- sapply(seq_len(nrow(res3)), function(i) mean(a[sets[[res3$set[i]]], ][[res3$nutrient[i]]]))
res3$pct_change <- 100 * res3$mean_change / base_mean
res3$q_welch    <- ave(res3$p_welch,    res3$set, FUN = bh)
res3$q_paired   <- ave(res3$p_paired,   res3$set, FUN = bh)
res3$q_wilcoxon <- ave(res3$p_wilcoxon, res3$set, FUN = bh)
print(res3[, c("set", "nutrient", "n", "mean_change", "pct_change", "prop_changed", "q_paired", "q_wilcoxon")], digits = 3)
save_csv(res3, "results_fndds_paired_tests_3sets.csv")

# ---- 6. compositional (ilr) analysis -------------------------------------------------------------------------------
make_parts <- function(df) {
  P <- df$`Protein (g)`; Ch <- df$`Carbohydrate (g)`; Ft <- df$`Total Fat (g)`
  cbind(P = P, Ch = Ch, Ft = Ft, Rm = 100 - (P + Ch + Ft))
}
parts_a <- make_parts(a); parts_b <- make_parts(b)
cat("remainder <= 0 (2017-18, 2019-20):", sum(parts_a[, 4] <= 0), sum(parts_b[, 4] <= 0), "\n")
cat("rows with a zero part:", sum(rowSums(parts_a <= 0) > 0), sum(rowSums(parts_b <= 0) > 0), "\n")
Za <- ilr3(mult_replace(parts_a)); Zb <- ilr3(mult_replace(parts_b))

ilr_res <- list(); hot_res <- list()
for (s in names(sets)) {
  D <- Zb[sets[[s]], ] - Za[sets[[s]], ]
  ilr_res[[s]] <- do.call(rbind, lapply(colnames(D), function(cn) {
    t <- t.test(D[, cn])
    data.frame(set = s, coord = cn, n = nrow(D), mean_diff = unname(t$estimate),
               ci_low = t$conf.int[1], ci_high = t$conf.int[2], p = t$p.value)
  }))
  ilr_res[[s]]$q <- bh(ilr_res[[s]]$p)
  hot_res[[s]] <- data.frame(set = s, n = nrow(D), t(hotelling_paired(D)))
}
ilr_tab <- do.call(rbind, ilr_res); hot_tab <- do.call(rbind, hot_res)
rownames(ilr_tab) <- NULL; rownames(hot_tab) <- NULL
print(ilr_tab, digits = 3); print(hot_tab, digits = 3)
save_csv(ilr_tab, "results_ilr_paired.csv")
save_csv(hot_tab, "results_ilr_hotelling.csv")

# zero-replacement sensitivity: other delta values, and dropping pairs with any zero part
sens <- do.call(rbind, lapply(c(0.001, 0.005, 0.05), function(dl) {
  Za2 <- ilr3(mult_replace(parts_a, dl)); Zb2 <- ilr3(mult_replace(parts_b, dl))
  do.call(rbind, lapply(names(sets), function(s) {
    D <- Zb2[sets[[s]], ] - Za2[sets[[s]], ]
    data.frame(variant = paste0("delta=", dl), set = s, n = nrow(D), t(hotelling_paired(D)))
  }))
}))
nz <- rowSums(parts_a <= 0) == 0 & rowSums(parts_b <= 0) == 0
sens_nz <- do.call(rbind, lapply(names(sets), function(s) {
  D <- Zb[sets[[s]] & nz, ] - Za[sets[[s]] & nz, ]
  data.frame(variant = "no zero parts", set = s, n = nrow(D), t(hotelling_paired(D)))
}))
sens_all <- rbind(sens, sens_nz); rownames(sens_all) <- NULL
print(sens_all, digits = 3)
save_csv(sens_all, "results_ilr_zero_sensitivity.csv")

cat("01_fndds.R done.\n")
