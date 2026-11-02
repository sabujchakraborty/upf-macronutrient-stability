# =============================================================================
# 07_supplementary.R  --  numbers quoted in the paper that steps 01-05 do not write out
#   (1) Table 2b and Section 3.3: effect of the scale-error screen
#       (percent change of the mean = mean change / mean first value; unclustered paired t; BH within sample)
#   (2) Limitations: size of the classifier-error bias in the UPF - nonUPF contrast
#       (OFF NOVA taken as correct; errors assumed unrelated to nutrient change)
# Inputs : fdc_pairs_raw.rds (02), results_classifier_crosstab_off.csv (03),
#          results_fdc_primary_cluster_robust.csv and results_fdc_primary_cluster_robust_off_confirmed.csv (04)
# Outputs: results_fdc_scale_screen_effect.csv, results_classifier_attenuation.csv
# Added in v1.0.1. Steps 00-06 are unchanged, so every other result is identical to v1.0.0.
# =============================================================================
source("00_setup.R")
cat("\n== 07 Supplementary ==\n")

# ---- 1. scale-error screen: before vs after, UPF and non-UPF cohorts (all source paths) ----------------------------
dat <- readRDS(out_path("fdc_pairs_raw.rds"))

change_of_mean <- function(d, nn, label) {
  x1 <- d[[paste0(nn, "_1")]]; x2 <- d[[paste0(nn, "_2")]]
  ok <- !is.na(x1) & !is.na(x2); x1 <- x1[ok]; x2 <- x2[ok]
  tt <- t.test(x2, x1, paired = TRUE)
  data.frame(sample = label, nutrient = nn, n = length(x1), mean_first = mean(x1),
             mean_change = mean(x2 - x1), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
             pct_change = 100 * mean(x2 - x1) / mean(x1), p = tt$p.value)
}
samples <- list(
  UPF_before    = dat[dat$cohort == "UPF_both", ],
  UPF_after     = dat[dat$cohort == "UPF_both" & !dat$scale_flag, ],
  nonUPF_before = dat[dat$cohort == "nonUPF_both", ],
  nonUPF_after  = dat[dat$cohort == "nonUPF_both" & !dat$scale_flag, ])
screen <- do.call(rbind, lapply(names(samples), function(s)
  do.call(rbind, lapply(NUTS, function(nn) change_of_mean(samples[[s]], nn, s)))))
screen$q <- ave(screen$p, screen$sample, FUN = bh)
rownames(screen) <- NULL
print(screen[, c("sample", "nutrient", "n", "mean_change", "ci_low", "ci_high", "pct_change", "q")], digits = 3)
save_csv(screen, "results_fdc_scale_screen_effect.csv")

g <- function(s, nn, col) screen[screen$sample == s & screen$nutrient == nn, col]
cat("Checks against the values quoted in the paper (Table 2b, Section 3.3):\n")
check("UPF pairs before / after screen", c(nrow(samples$UPF_before), nrow(samples$UPF_after)), c(145337, 145057))
check("Table 2b before, % change (kcal, protein, fat, carb, sugars, fibre, sat_fat, sodium)",
      round(sapply(NUTS, function(nn) g("UPF_before", nn, "pct_change")), 2),
      c(-0.91, -1.79, -0.92, -0.77, -1.21, -0.71, -1.00, 3.52))
check("Table 2b after, % change",
      round(sapply(NUTS, function(nn) g("UPF_after", nn, "pct_change")), 2),
      c(-0.17, -0.33, -0.16, -0.17, -0.48, -0.39, -0.36, 3.84))
check("UPF energy before: mean change and 95% CI (kcal)",
      round(c(g("UPF_before", "kcal", "mean_change"), g("UPF_before", "kcal", "ci_low"), g("UPF_before", "kcal", "ci_high")), 2),
      c(-2.47, -3.37, -1.57))
check("UPF energy after: mean change (kcal)", round(g("UPF_after", "kcal", "mean_change"), 2), -0.47)
check("UPF protein before: mean change and 95% CI (g)",
      round(c(g("UPF_before", "protein", "mean_change"), g("UPF_before", "protein", "ci_low"), g("UPF_before", "protein", "ci_high")), 3),
      c(-0.125, -0.163, -0.088))
check("q, UPF energy and protein before", signif(c(g("UPF_before", "kcal", "q"), g("UPF_before", "protein", "q")), 2), c(1.9e-7, 4.2e-10))
cat(sprintf("  [info] sodium %% change of the mean: UPF %.2f / %.2f, non-UPF %.2f / %.2f (before / after); min q = %.2g\n",
            g("UPF_before", "sodium", "pct_change"), g("UPF_after", "sodium", "pct_change"),
            g("nonUPF_before", "sodium", "pct_change"), g("nonUPF_after", "sodium", "pct_change"),
            min(screen$q[screen$nutrient == "sodium"])))

# ---- 2. classifier error: attenuation of the UPF - nonUPF contrast ---------------------------------------------------
xt <- read.csv(out_path("results_classifier_crosstab_off.csv"), check.names = FALSE)
upf_row <- xt[xt$our_cohort == "UPF_both", c("1", "2", "3", "4")]
non_row <- xt[xt$our_cohort == "nonUPF_both", c("1", "2", "3", "4")]
a <- sum(upf_row[c("1", "2", "3")]) / sum(upf_row)      # share of our UPF cohort that OFF places in NOVA 1-3
b <- non_row[["4"]] / sum(non_row)                       # share of our non-UPF cohort that OFF places in NOVA 4
f <- 1 - a - b                                           # a true difference D is observed as f * D

res  <- read.csv(out_path("results_fdc_primary_cluster_robust.csv"))
conf <- read.csv(out_path("results_fdc_primary_cluster_robust_off_confirmed.csv"))
d  <- res[res$what == "UPF - nonUPF", ];  d  <- d[match(NUTS, d$nutrient), ]
u  <- res[res$what == "UPF_both", ];      u  <- u[match(NUTS, u$nutrient), ]
dc <- conf[conf$what == "UPF - nonUPF", ]; dc <- dc[match(NUTS, dc$nutrient), ]
att <- data.frame(nutrient = NUTS, share_upf_misassigned = a, share_nonupf_misassigned = b, attenuation = f,
                  diff_observed = d$est, diff_corrected = d$est / f, diff_off_confirmed = dc$est,
                  upf_mean_observed = u$est, upf_mean_corrected = u$est + a * d$est / f)
att$upf_mean_shift <- att$upf_mean_corrected - att$upf_mean_observed
print(att, digits = 3)
save_csv(att, "results_classifier_attenuation.csv")
cat("Checks against the values quoted in the paper (Limitations):\n")
check("misassigned shares (UPF, non-UPF) and attenuation factor", round(c(a, b, f), 3), c(0.044, 0.112, 0.844))
check("corrected differences, fibre and sugars (pp)",
      round(att$diff_corrected[match(c("fibre", "sugars"), NUTS)], 2), c(-0.61, 0.25))
check("largest shift in a UPF cohort mean <= 0.03 pp", max(abs(att$upf_mean_shift)) <= 0.03, TRUE)

cat("07_supplementary.R done.\n")
