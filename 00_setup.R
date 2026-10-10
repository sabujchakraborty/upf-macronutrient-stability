# =============================================================================
# 00_setup.R  --  configuration, packages and helper functions
# Study: macronutrient stability in ultra-processed (NOVA 4) foods
#        (USDA FNDDS 2017-18 vs 2019-20; FoodData Central Branded Foods;
#         Open Food Facts used only to validate the UPF classifier)
#
# EDIT THE THREE PATHS BELOW, then run run_all.R (or each script in order).
# =============================================================================

WORK_DIR <- "E:/MSc Docs/Mini Project/KNIME output+work"  # fndds_foods.csv, fndds_nutrients.csv, NOVA lookup
FDC_DIR  <- "E:/MSc Docs/Mini Project/FDC_branded"         # branded_food.csv, food.csv, food_nutrient.csv, nutrient.csv
OUT_DIR  <- file.path(WORK_DIR, "pipeline_output")         # every result of this pipeline is written here

OFF_URL  <- "https://huggingface.co/datasets/openfoodfacts/product-database/resolve/main/food.parquet"
LOOKUP_FILE <- "nova_lookup_martinez_steele_2017_2018.xlsx"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORK_DIR)
options(width = 200, stringsAsFactors = FALSE)

# ---- packages (install once by hand: install.packages(c("readr","dplyr","readxl","stringr","DBI","duckdb"),
#                                                         repos = "https://cloud.r-project.org")) ------------------
need <- c("readr", "dplyr", "readxl", "stringr", "DBI", "duckdb")
missing_pkgs <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs))
  stop("Install these packages first: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(readxl); library(stringr); library(DBI); library(duckdb)
})

# ---- constants --------------------------------------------------------------------------------------------------
NUTS <- c("kcal", "protein", "fat", "carb", "sugars", "fibre", "sat_fat", "sodium")
NUT_LABEL <- c(kcal = "Energy", protein = "Protein", fat = "Total fat", carb = "Carbohydrate",
               sugars = "Total sugars", fibre = "Dietary fibre", sat_fat = "Saturated fat", sodium = "Sodium")
EQUIV_MARGINS <- c(0.5, 1, 2)   # equivalence margins in % SPD; 1 is primary

# ---- small utilities -------------------------------------------------------------------------------------------
bh        <- function(p) p.adjust(p, method = "BH")
out_path  <- function(f) file.path(OUT_DIR, f)
save_csv  <- function(x, f) { write.csv(x, out_path(f), row.names = FALSE); invisible(x) }
norm_barcode <- function(x) sub("^0+", "", trimws(as.character(x)))   # strip leading zeros (FDC vs OFF formats differ)
check <- function(label, got, expected) {
  ok <- isTRUE(all.equal(as.numeric(got), as.numeric(expected)))
  cat(sprintf("  [%s] %s: got %s (expected %s)\n", if (ok) "OK" else "CHECK", label,
              paste(got, collapse = "/"), paste(expected, collapse = "/")))
  invisible(ok)
}

# Symmetric percent difference (bounded, zero-safe)
spd <- function(d, n) {
  x1 <- d[[paste0(n, "_1")]]; x2 <- d[[paste0(n, "_2")]]
  ifelse(!is.na(x1) & !is.na(x2) & (x1 + x2) > 0, 100 * (x2 - x1) / ((x1 + x2) / 2), NA_real_)
}

# ---- FNDDS: paired tests ---------------------------------------------------------------------------------------
run_tests <- function(a, b, macros, label) {
  out <- lapply(macros, function(m) {
    x <- a[[m]]; y <- b[[m]]; d <- y - x
    w  <- t.test(y, x)                                   # Welch, unpaired (as specified in the protocol)
    p  <- t.test(y, x, paired = TRUE)                    # paired t-test (inference of record)
    wl <- suppressWarnings(wilcox.test(y, x, paired = TRUE, exact = FALSE))
    data.frame(set = label, nutrient = m, n = length(d),
               mean_change = mean(d), ci_low = p$conf.int[1], ci_high = p$conf.int[2],
               dz = mean(d) / sd(d), prop_changed = mean(d != 0),
               p_welch = w$p.value, p_paired = p$p.value, p_wilcoxon = wl$p.value)
  })
  do.call(rbind, out)
}

# ---- compositional data (ilr) ----------------------------------------------------------------------------------
mult_replace <- function(x, delta = 0.005) {            # multiplicative zero replacement
  k <- rowSums(x <= 0); out <- x * (1 - k * delta / 100); out[x <= 0] <- delta; out
}
ilr3 <- function(x) {                                    # sequential binary partition; parts = P, Ch, Ft, Rm
  l <- log(x)
  cbind(z1 = sqrt(3/4) * ((l[, 1] + l[, 2] + l[, 3]) / 3 - l[, 4]),
        z2 = sqrt(2/3) * (l[, 3] - 0.5 * (l[, 1] + l[, 2])),
        z3 = sqrt(1/2) * (l[, 1] - l[, 2]))
}
hotelling_paired <- function(D) {
  n <- nrow(D); p <- ncol(D); dbar <- colMeans(D)
  T2 <- n * drop(t(dbar) %*% solve(cov(D)) %*% dbar)
  Fs <- (n - p) / (p * (n - 1)) * T2
  c(T2 = T2, F = Fs, p = pf(Fs, p, n - p, lower.tail = FALSE))
}

# ---- Branded Foods: two-part estimand --------------------------------------------------------------------------
# Part A: share of pairs moving to / from zero (exact McNemar). Part B: SPD among pairs with both values > 0.
two_part <- function(d, n, label, margin = 1) {
  x1 <- d[[paste0(n, "_1")]]; x2 <- d[[paste0(n, "_2")]]
  ok <- !is.na(x1) & !is.na(x2); x1 <- x1[ok]; x2 <- x2[ok]
  to0   <- sum(x1 > 0 & x2 == 0)
  from0 <- sum(x1 == 0 & x2 > 0)
  p_zero <- if (to0 + from0 > 0) binom.test(to0, to0 + from0, 0.5)$p.value else NA_real_
  pos <- x1 > 0 & x2 > 0
  s <- 100 * (x2[pos] - x1[pos]) / ((x1[pos] + x2[pos]) / 2)
  m <- mean(s); se <- sd(s) / sqrt(length(s)); df <- length(s) - 1
  ci <- m + c(-1, 1) * qt(0.975, df) * se
  p_tost <- max(pt((m + margin) / se, df, lower.tail = FALSE), pt((m - margin) / se, df))
  data.frame(group = label, nutrient = n, n_pairs = length(x1), n_pos = length(s),
             pct_to_zero = 100 * to0 / length(x1), pct_from_zero = 100 * from0 / length(x1),
             p_zero = p_zero, mean_spd_pos = m, lo = ci[1], hi = ci[2],
             median_spd_pos = median(s), pct_identical = 100 * mean(s == 0), p_tost_1pct = p_tost)
}

run_two_part <- function(g) {
  res <- do.call(rbind, lapply(NUTS, function(nn)
    do.call(rbind, lapply(split(seq_len(nrow(g)), paste(g$cohort, g$src_path)), function(ix)
      two_part(g[ix, ], nn, paste(g$cohort[ix[1]], g$src_path[ix[1]]))))))
  rownames(res) <- NULL
  res$q_zero <- bh(res$p_zero)
  res$q_tost <- bh(res$p_tost_1pct)
  res$equivalent_1pct <- res$q_tost < 0.05
  res
}

# ---- Branded Foods: brand-owner cluster-robust inference -------------------------------------------------------
# OLS with CR1 cluster-robust variance; returns coefficient j with a t(G-1) interval
clu_ols <- function(y, X, cl, j) {
  fit <- lm.fit(X, y); e <- fit$residuals
  U <- rowsum(X * e, cl)
  G <- nrow(U); n <- length(y); k <- ncol(X)
  XtXi <- solve(crossprod(X))
  V <- XtXi %*% crossprod(U) %*% XtXi * (G / (G - 1)) * ((n - 1) / (n - k))
  b <- unname(fit$coefficients[j]); se <- sqrt(V[j, j])
  c(est = b, se = se, lo = b - qt(.975, G - 1) * se, hi = b + qt(.975, G - 1) * se,
    p = 2 * pt(-abs(b / se), G - 1), G = G, n = n)
}
tost_p <- function(est, se, G, m)
  max(pt((est + m) / se, G - 1, lower.tail = FALSE), pt((est - m) / se, G - 1))

# Cohort means (UPF, non-UPF) with equivalence tests at several margins, plus the UPF - nonUPF difference.
# `L` must contain cohort, cl (cluster id) and the <nutrient>_1 / _2 columns; only pairs with both values > 0 are used.
run_cluster <- function(L) {
  rows <- list()
  for (nn in NUTS) {
    x1 <- L[[paste0(nn, "_1")]]; x2 <- L[[paste0(nn, "_2")]]
    ok <- !is.na(x1) & !is.na(x2) & x1 > 0 & x2 > 0
    s <- 100 * (x2[ok] - x1[ok]) / ((x1[ok] + x2[ok]) / 2)
    cl <- L$cl[ok]; upf <- as.numeric(L$cohort[ok] == "UPF_both")
    for (co in c("UPF_both", "nonUPF_both")) {
      i <- (upf == 1) == (co == "UPF_both")
      r <- clu_ols(s[i], matrix(1, sum(i), 1), cl[i], 1)
      row <- data.frame(nutrient = nn, what = co, n = unname(r["n"]), G = unname(r["G"]), est = unname(r["est"]),
                        se = unname(r["se"]), lo = unname(r["lo"]), hi = unname(r["hi"]), p = unname(r["p"]))
      for (m in EQUIV_MARGINS) row[[paste0("p_tost_", m)]] <- tost_p(r["est"], r["se"], r["G"], m)
      rows[[length(rows) + 1]] <- row
    }
    r <- clu_ols(s, cbind(1, upf), cl, 2)
    row <- data.frame(nutrient = nn, what = "UPF - nonUPF", n = unname(r["n"]), G = unname(r["G"]), est = unname(r["est"]),
                      se = unname(r["se"]), lo = unname(r["lo"]), hi = unname(r["hi"]), p = unname(r["p"]))
    for (m in EQUIV_MARGINS) row[[paste0("p_tost_", m)]] <- NA_real_
    rows[[length(rows) + 1]] <- row
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out$q <- ave(out$p, out$what, FUN = bh)
  for (m in EQUIV_MARGINS) {
    pc <- paste0("p_tost_", m)
    out[[paste0("q_tost_", m)]]  <- ave(out[[pc]], out$what, FUN = bh)
    out[[paste0("equiv_", m)]]   <- out[[paste0("q_tost_", m)]] < 0.05
  }
  out$p_tost <- out$p_tost_1; out$q_tost <- out$q_tost_1; out$equivalent_1pct <- out$equiv_1
  out
}

cat("00_setup.R loaded. Work dir:", getwd(), "| outputs:", OUT_DIR, "\n")
