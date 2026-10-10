# =============================================================================
# 04_fdc_analysis.R  --  Branded Foods: within-product comparison (all primary and sensitivity analyses)
# Input  : fdc_dd_final.rds (from 03)
# Outputs: see the save_csv() calls below (all prefixed results_fdc_)
# Strata : LI -> LI (primary, 92% of pairs; last versions through 2023) and GDSN -> GDSN (separate stratum).
# =============================================================================
source("00_setup.R")
cat("\n== 04 FDC analysis ==\n")
dd <- readRDS(out_path("fdc_dd_final.rds"))

li  <- dd[dd$src_path == "LI -> LI"     & dd$cohort %in% c("UPF_both", "nonUPF_both"), ]
gd  <- dd[dd$src_path == "GDSN -> GDSN" & dd$cohort %in% c("UPF_both", "nonUPF_both"), ]
both_strata <- rbind(li, gd)
confirmed <- function(d) d[!is.na(d$off_nova) &
                           ((d$cohort == "UPF_both" & d$off_nova == 4) |
                            (d$cohort == "nonUPF_both" & d$off_nova %in% 1:3)), ]
li_conf <- confirmed(li)
cat("LI->LI pairs:", nrow(li), "| GDSN->GDSN pairs:", nrow(gd), "| OFF-confirmed LI->LI:", nrow(li_conf), "\n")
check("LI->LI pairs", nrow(li), 191818); check("OFF-confirmed LI->LI", nrow(li_conf), 135867)

# ---- 1. mean SPD by cohort and source path (all pairs; zeros included) -- shows the GDSN zero artefact ---------------
src_tab <- do.call(rbind, lapply(NUTS, function(nn) {
  s <- spd(dd, nn)
  do.call(rbind, lapply(split(seq_len(nrow(dd)), paste(dd$cohort, dd$src_path)), function(ix) {
    z <- s[ix]; z <- z[!is.na(z)]
    if (length(z) < 300) return(NULL)
    tt <- t.test(z)
    data.frame(nutrient = nn, cohort = dd$cohort[ix[1]], src_path = dd$src_path[ix[1]], n = length(z),
               mean_spd = mean(z), lo = tt$conf.int[1], hi = tt$conf.int[2])
  }))
}))
rownames(src_tab) <- NULL
print(src_tab[src_tab$nutrient %in% c("sat_fat", "fibre", "fat"), ], digits = 3)
save_csv(src_tab, "results_fdc_by_source_path.csv")

# ---- 2. why GDSN differs: timing of last versions and zero prevalence ---------------------------------------------
gdsn_li <- both_strata
gdsn_li$yr_last <- format(as.Date(gdsn_li$date_last), "%Y")
yr_tab <- as.data.frame.matrix(round(100 * prop.table(table(gdsn_li$yr_last, gdsn_li$src_path), 2), 1))
cat("\n% of pairs by year of last version:\n"); print(yr_tab)
save_csv(cbind(year_last = rownames(yr_tab), yr_tab), "results_fdc_last_year_by_source.csv")

zero_tab <- do.call(rbind, lapply(c("sat_fat", "fibre", "sugars", "fat"), function(nn) {
  x1 <- gdsn_li[[paste0(nn, "_1")]]; x2 <- gdsn_li[[paste0(nn, "_2")]]
  do.call(rbind, lapply(split(seq_len(nrow(gdsn_li)), gdsn_li$src_path), function(ix)
    data.frame(nutrient = nn, path = gdsn_li$src_path[ix[1]], n = length(ix),
               pct_na_1 = 100 * mean(is.na(x1[ix])), pct_na_2 = 100 * mean(is.na(x2[ix])),
               pct_zero_1 = 100 * mean(x1[ix] == 0, na.rm = TRUE), pct_zero_2 = 100 * mean(x2[ix] == 0, na.rm = TRUE))))
}))
rownames(zero_tab) <- NULL
print(zero_tab, digits = 3)
save_csv(zero_tab, "results_fdc_gdsn_zero_patterns.csv")

# ---- 3. two-part estimand (unclustered), both strata, full cohorts and OFF-confirmed -------------------------------
tp_full <- run_two_part(both_strata)
tp_conf <- run_two_part(confirmed(both_strata))
print(tp_full[order(tp_full$group, tp_full$nutrient),
              c("group", "nutrient", "n_pos", "pct_to_zero", "pct_from_zero", "q_zero", "mean_spd_pos", "lo", "hi",
                "pct_identical", "q_tost", "equivalent_1pct")], digits = 3)
save_csv(tp_full, "results_fdc_two_part.csv")
save_csv(tp_conf, "results_fdc_two_part_off_confirmed.csv")

# ---- 4. PRIMARY: LI -> LI, brand-owner cluster-robust; equivalence at +/-0.5, 1, 2 ---------------------------------
res_li <- run_cluster(li)
print(res_li[order(res_li$what, res_li$nutrient),
             c("nutrient", "what", "n", "G", "est", "lo", "hi", "q", "q_tost", "equiv_0.5", "equiv_1", "equiv_2")], digits = 3)
save_csv(res_li, "results_fdc_primary_cluster_robust.csv")

# ---- 5. sensitivity: OFF-confirmed cohorts (classifier-error check) and the GDSN stratum ------------------------------
res_conf <- run_cluster(li_conf)
save_csv(res_conf, "results_fdc_primary_cluster_robust_off_confirmed.csv")
res_gd <- run_cluster(gd)
save_csv(res_gd, "results_fdc_gdsn_cluster_robust.csv")

a <- res_li[res_li$what == "UPF - nonUPF", ]; b <- res_conf[res_conf$what == "UPF - nonUPF", ]
cmp <- data.frame(nutrient = a$nutrient,
                  full = sprintf("%.2f (%.2f, %.2f)", a$est, a$lo, a$hi), q_full = signif(a$q, 2),
                  confirmed = sprintf("%.2f (%.2f, %.2f)", b$est, b$lo, b$hi), q_confirmed = signif(b$q, 2))
cat("\nUPF - nonUPF difference (percentage points), brand-clustered:\n"); print(cmp, row.names = FALSE)
save_csv(cmp, "results_fdc_diff_clustered_full_vs_confirmed.csv")

cat("04_fdc_analysis.R done.\n")
