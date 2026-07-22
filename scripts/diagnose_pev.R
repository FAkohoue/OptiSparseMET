# diagnose_pev.R
# ------------------------------------------------------------------------------
# Does met_evaluate_alpha_efficiency() (random, IID) return the correct
# design-based mean PEV? We compare it against an INDEPENDENT dense mixed-model-
# equations (MME) solve using the SAME fixed variance components. That dense MME
# is exactly what a mixed-model engine (e.g. sommer) returns when its variance
# components are fixed at these values -- with no estimation noise.
#
# Model (random-treatment, IID entries), matching the evaluator:
#   y = 1*mu + Zrep u_rep + Zib u_ib + Zrow u_row + Zcol u_col + Zg g + e
#   g ~ N(0, sigma_g2 I),  e ~ N(0, sigma_e2 I),  nuisance ~ N(0, their variances)
#   Fixed part = intercept only (no checks); IID residual.
#
# Coefficient matrix:  C = [[X'QX, X'QZ], [Z'QX, Z'QZ + Ginv]],  Q = R^{-1}
# PEV of genotype g = diagonal of C^{-1} on the genotype block.
# ------------------------------------------------------------------------------

library(OptiSparseMET)

vc <- list(sigma_e2 = 2, sigma_g2 = 1,
           sigma_rep2 = 1, sigma_ib2 = 1, sigma_r2 = 1, sigma_c2 = 1)

## Independent dense-MME mean genotype PEV -------------------------------------
dense_mean_pev <- function(fb, vc) {
  # One 0/1 column per unique level (handles single-level factors, matching the
  # evaluator's .make_sparse_incidence).
  inc <- function(x) {
    x  <- as.character(x)
    lv <- unique(x)
    m  <- matrix(0, length(x), length(lv))
    for (j in seq_along(lv)) m[x == lv[j], j] <- 1
    m
  }
  nn  <- nrow(fb)
  X   <- matrix(1, nn, 1)                                  # intercept (fixed)
  Zrep <- inc(fb$Rep); Zb <- inc(fb$IBlock)
  Zr   <- inc(fb$Row); Zc <- inc(fb$Column); Zg <- inc(fb$Treatment)
  Z   <- cbind(Zrep, Zb, Zr, Zc, Zg)
  ng  <- c(ncol(Zrep), ncol(Zb), ncol(Zr), ncol(Zc), ncol(Zg))
  Q   <- diag(nn) / vc$sigma_e2
  Ginv <- diag(c(rep(1 / vc$sigma_rep2, ng[1]), rep(1 / vc$sigma_ib2, ng[2]),
                 rep(1 / vc$sigma_r2,  ng[3]), rep(1 / vc$sigma_c2, ng[4]),
                 rep(1 / vc$sigma_g2,  ng[5])))
  C <- rbind(cbind(t(X) %*% Q %*% X, t(X) %*% Q %*% Z),
             cbind(t(Z) %*% Q %*% X, t(Z) %*% Q %*% Z + Ginv))
  Cinv <- solve(C)
  gpos <- (1 + sum(ng[1:4]) + 1):(1 + sum(ng))
  mean(diag(Cinv)[gpos])
}

## ---- Case A: the single-cell design that produced 0.66 ----------------------
fbA <- data.frame(
  Plot = 1:12, Row = 1L, Column = 1L, Rep = 1L, IBlock = 1L, BlockInRep = 1L,
  Treatment = rep(paste0("G", 1:4), each = 3), Check = FALSE,
  stringsAsFactors = FALSE)

effA <- met_evaluate_alpha_efficiency(
  field_book = fbA, n_rows = 1, n_cols = 1, check_treatments = character(0),
  treatment_effect = "random", prediction_type = "IID",
  varcomp = vc, residual_structure = "IID")

cat("== CASE A: single-cell (all nuisance factors single-level) ==\n")
cat("  evaluator mean_PEV :", round(effA$mean_PEV, 5), "\n")
cat("  dense-MME mean_PEV :", round(dense_mean_pev(fbA, vc), 5), "\n\n")

## ---- Case B: non-degenerate design (nuisance factors have >= 2 levels) ------
fbB <- data.frame(
  Plot = 1:12,
  Row = rep(1:3, times = 4), Column = rep(1:4, each = 3),
  Rep = rep(1:2, each = 6), IBlock = rep(1:4, each = 3),
  BlockInRep = rep(1:2, times = 6),
  Treatment = c(paste0("G", 1:6), paste0("G", 1:6)), Check = FALSE,
  stringsAsFactors = FALSE)

effB <- met_evaluate_alpha_efficiency(
  field_book = fbB, n_rows = 3, n_cols = 4, check_treatments = character(0),
  treatment_effect = "random", prediction_type = "IID",
  varcomp = vc, residual_structure = "IID")

cat("== CASE B: non-degenerate ==\n")
cat("  evaluator mean_PEV :", round(effB$mean_PEV, 5), "\n")
cat("  dense-MME mean_PEV :", round(dense_mean_pev(fbB, vc), 5), "\n\n")

cat("INTERPRETATION\n")
cat("  If evaluator == dense-MME in both cases -> evaluator is correct.\n")
cat("  If they differ, the gap and the case pinpoint the issue.\n\n")

## ---- Optional: sommer::mmes fit on Case B (estimates VC; approximate) --------
## Current sommer API: mmes() with vsm(ism(...)) (mmer()/vsr() are deprecated).
## A fit ESTIMATES the variance components from only 12 observations, so its PEV
## will NOT equal the fixed-VC design PEV above -- it is a rough sanity check
## only. The base-R dense-MME comparison above (with fixed VC) is the definitive
## one. (mmes has no init/theta argument to fix VC at arbitrary known values.)
if (requireNamespace("sommer", quietly = TRUE)) {
  datB <- fbB
  datB$Geno <- factor(datB$Treatment)
  set.seed(1)
  G0 <- diag(nlevels(datB$Geno))
  dimnames(G0) <- list(levels(datB$Geno), levels(datB$Geno))
  g  <- rnorm(nlevels(datB$Geno), 0, sqrt(vc$sigma_g2))
  datB$y <- g[datB$Geno] + rnorm(nrow(datB), 0, sqrt(vc$sigma_e2))

  fit <- try(sommer::mmes(
    y ~ 1,
    random = ~ sommer::vsm(sommer::ism(Geno), Gu = G0),
    rcov   = ~ units, data = datB,
    getPEV = TRUE, verbose = FALSE, dateWarning = FALSE), silent = TRUE)

  if (!inherits(fit, "try-error")) {
    # Genotype PEV: prefer uPevList; fall back to the diagonal of Ci (inverse of
    # the coefficient matrix) at the genotype partition.
    pev <- tryCatch(as.numeric(unlist(fit$uPevList)), error = function(e) NULL)
    if (is.null(pev) || !length(pev)) {
      pev <- tryCatch({
        part <- fit$partitions[[1]]                 # start/end of Geno block
        diag(as.matrix(fit$Ci))[part[1, 1]:part[1, 2]]
      }, error = function(e) NA_real_)
    }
    cat("== sommer::mmes (VC ESTIMATED, approximate) ==\n")
    cat("  sommer mean PEV(Geno):", round(mean(pev), 5), "\n")
    cat("  (VC estimated from 12 rows => not equal to the fixed-VC design PEV;\n")
    cat("   the base-R dense MME above is the exact fixed-VC comparison.)\n")
  } else {
    cat("sommer::mmes fit failed:\n  ", as.character(fit), "\n")
  }
} else {
  cat("sommer not installed; skipping optional fit.\n")
}
