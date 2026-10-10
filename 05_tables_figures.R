# =============================================================================
# 05_tables_figures.R  --  formatted tables and Figure 1 for the manuscript
# Inputs : result CSVs written by 01 and 04 (in OUT_DIR)
# Outputs: table_fndds_setA.csv, table_fdc_cohort_means.csv, table_fdc_differences.csv,
#          fig_fdc_forest.pdf (vector) and fig_fdc_forest.png (600 dpi)
# =============================================================================
source("00_setup.R")
cat("\n== 05 Tables and figures ==\n")

fmt_ci  <- function(e, l, h, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), e, l, h)
fmt_q   <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))
fmt_n   <- function(x) format(x, big.mark = ",", trim = TRUE)

# ---- Table 1: FNDDS, primary set A ---------------------------------------------------------------------------------
f <- read.csv(out_path("results_fndds_paired_tests_3sets.csv"), check.names = FALSE)
fa <- f[f$set == "A_primary_full", ]
t1 <- data.frame(Nutrient = fa$nutrient,
                 `Mean change (95% CI)` = fmt_ci(fa$mean_change, fa$ci_low, fa$ci_high, 3),
                 `% change` = sprintf("%+.2f", fa$pct_change),
                 `Codes changed (%)` = sprintf("%.1f", 100 * fa$prop_changed),
                 `q (paired t)` = fmt_q(fa$q_paired), `q (Wilcoxon)` = fmt_q(fa$q_wilcoxon),
                 check.names = FALSE)
print(t1, row.names = FALSE); save_csv(t1, "table_fndds_setA.csv")

# ---- Table 3 and 4: Branded Foods LI -> LI -----------------------------------------------------------------------
res  <- read.csv(out_path("results_fdc_primary_cluster_robust.csv"))
conf <- read.csv(out_path("results_fdc_primary_cluster_robust_off_confirmed.csv"))
u <- res[res$what == "UPF_both", ];    u <- u[match(NUTS, u$nutrient), ]
v <- res[res$what == "nonUPF_both", ]; v <- v[match(NUTS, v$nutrient), ]
t3 <- data.frame(Nutrient = NUT_LABEL[NUTS],
                 `n (UPF / non-UPF)` = paste(fmt_n(u$n), "/", fmt_n(v$n)),
                 `UPF mean SPD (95% CI)` = fmt_ci(u$est, u$lo, u$hi), `q (UPF)` = fmt_q(u$q),
                 `Non-UPF mean SPD (95% CI)` = fmt_ci(v$est, v$lo, v$hi), `q (non-UPF)` = fmt_q(v$q),
                 `Equivalent at +/-0.5%` = ifelse(u$equiv_0.5 & v$equiv_0.5, "both", ifelse(u$equiv_0.5, "UPF only", ifelse(v$equiv_0.5, "non-UPF only", "neither"))),
                 `Equivalent at +/-1%`   = ifelse(u$equiv_1 & v$equiv_1, "both", ifelse(u$equiv_1, "UPF only", ifelse(v$equiv_1, "non-UPF only", "neither"))),
                 `Equivalent at +/-2%`   = ifelse(u$equiv_2 & v$equiv_2, "both", ifelse(u$equiv_2, "UPF only", ifelse(v$equiv_2, "non-UPF only", "neither"))),
                 check.names = FALSE)
print(t3, row.names = FALSE); save_csv(t3, "table_fdc_cohort_means.csv")

dfull <- res[res$what == "UPF - nonUPF", ];  dfull <- dfull[match(NUTS, dfull$nutrient), ]
dconf <- conf[conf$what == "UPF - nonUPF", ]; dconf <- dconf[match(NUTS, dconf$nutrient), ]
t4 <- data.frame(Nutrient = NUT_LABEL[NUTS],
                 `Full cohorts (95% CI)` = fmt_ci(dfull$est, dfull$lo, dfull$hi), `q` = fmt_q(dfull$q),
                 `OFF-confirmed (95% CI)` = fmt_ci(dconf$est, dconf$lo, dconf$hi), `q (confirmed)` = fmt_q(dconf$q),
                 check.names = FALSE)
print(t4, row.names = FALSE); save_csv(t4, "table_fdc_differences.csv")

# ---- Figure 1: forest plot --------------------------------------------------------------------------------------------
ord  <- NUTS
cols <- c(UPF_both = "#D55E00", nonUPF_both = "#0072B2")     # colour-blind-safe pair

draw_fig <- function() {
  par(mfrow = c(1, 2), mar = c(6.2, 7, 2.5, 0.8), cex.axis = 0.8, cex.lab = 0.85)
  yl <- c(-0.7, 8.6)
  plot(NA, xlim = c(-1.6, 1.6), ylim = yl, yaxt = "n", ylab = "",
       xlab = "Mean symmetric % difference", main = "A. Within-cohort change", cex.main = 0.95)
  rect(-1, yl[1], 1, yl[2], col = "grey92", border = NA); abline(v = 0, lty = 2)
  for (i in seq_along(ord)) for (co in names(cols)) {
    r <- res[res$nutrient == ord[i] & res$what == co, ]
    yy <- (9 - i) + ifelse(co == "UPF_both", 0.17, -0.17)
    segments(r$lo, yy, r$hi, yy, col = cols[co], lwd = 2)
    points(r$est, yy, pch = 19, col = cols[co], cex = 0.9)
  }
  axis(2, at = 8:1, labels = NUT_LABEL[ord], las = 1, tick = FALSE, cex.axis = 0.8)
  legend("bottom", c("UPF", "Non-UPF"), col = cols, pch = 19, lwd = 2, bty = "n", horiz = TRUE, cex = 0.8)
  mtext("Shaded band: +/-1% equivalence margin", side = 1, line = 4.9, cex = 0.7); box()

  plot(NA, xlim = c(-1.1, 1.1), ylim = yl, yaxt = "n", ylab = "",
       xlab = "UPF minus non-UPF (percentage points)", main = "B. Difference between cohorts", cex.main = 0.95)
  abline(v = 0, lty = 2)
  for (i in seq_along(ord)) {
    r <- res[res$nutrient == ord[i] & res$what == "UPF - nonUPF", ]
    yy <- 9 - i
    segments(r$lo, yy, r$hi, yy, lwd = 2, col = "grey30")
    points(r$est, yy, pch = ifelse(r$q < 0.05, 19, 1), col = "grey30", cex = 1)
  }
  axis(2, at = 8:1, labels = NUT_LABEL[ord], las = 1, tick = FALSE, cex.axis = 0.8)
  legend("bottom", c("BH q < 0.05", "BH q >= 0.05"), pch = c(19, 1), bty = "n", horiz = TRUE, cex = 0.8)
  mtext("Bars: brand-owner cluster-robust 95% CI", side = 1, line = 4.9, cex = 0.7); box()
}
png(out_path("fig_fdc_forest.png"), width = 5220, height = 3000, res = 600); draw_fig(); invisible(dev.off())
pdf(out_path("fig_fdc_forest.pdf"), width = 8.7, height = 5);                draw_fig(); invisible(dev.off())

cat("Files in", OUT_DIR, ":\n"); print(list.files(OUT_DIR))
cat("05_tables_figures.R done.\n")
