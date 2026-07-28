# Tests for the Mothukuri et al. (2025) additions: breeder outcome metrics,
# the 8-criterion training-set menu, factor-analytic Sigma_E, and the
# design-criterion validation harness.

mk <- function(J = 40L, E = 4L, seed = 11) {
  set.seed(seed)
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  SigE <- matrix(c(1, .6, .2, -.1,  .6, 1, .3, 0,
                   .2, .3, 1, .5,  -.1, 0, .5, 1), 4, 4)
  list(G = G, SigE = SigE, J = J, E = E)
}

# ---- M1: breeder outcome metrics in simulate_met ----------------------------

test_that("simulate_met returns valid common-selected and average-rank metrics", {
  d <- mk()
  full <- matrix(1L, d$J, d$E,
                 dimnames = list(rownames(d$G), paste0("E", seq_len(d$E))))
  res <- simulate_met(full, G = d$G, Sigma_E = d$SigE, n_sim = 30,
                      select_fraction = 0.1, seed = 1)
  expect_true(res$common_selected_mean >= 0 && res$common_selected_mean <= 1)
  n_sel <- res$n_selected
  # average rank lies between the ideal (n_sel+1)/2 and the population mean rank
  expect_gte(res$avg_rank_mean, (n_sel + 1) / 2 - 1e-8)
  expect_lte(res$avg_rank_mean, (d$J + 1) / 2 + 5)
})

test_that("a balanced design selects more truly-best lines than a sparse one", {
  d <- mk()
  ln <- rownames(d$G); en <- paste0("E", seq_len(d$E))
  full <- matrix(1L, d$J, d$E, dimnames = list(ln, en))
  set.seed(3)
  sparse <- matrix(0L, d$J, d$E, dimnames = list(ln, en))
  for (j in seq_len(d$J)) sparse[j, sample(d$E, 2)] <- 1L
  cf <- simulate_met(full,   G = d$G, Sigma_E = d$SigE, n_sim = 40, seed = 5)$common_selected_mean
  cs <- simulate_met(sparse, G = d$G, Sigma_E = d$SigE, n_sim = 40, seed = 5)$common_selected_mean
  expect_gte(cf, cs)
})

# ---- M3: 8-criterion training-set selection ---------------------------------

test_that("select_individuals reports all eight relationship measurements", {
  d <- mk()
  res <- select_individuals(d$G, n_train = 12, criterion = "cdmean",
                            method = "exchange", iters = 200, seed = 5)
  expect_named(res$measures,
               c("cdmean", "cdmax", "pevmean", "pevmax",
                 "aopt", "dopt", "goptpev", "neg_dist"),
               ignore.order = TRUE)
  expect_equal(res$measures$aopt, res$measures$pevmean * d$J, tolerance = 1e-6)
  expect_true(res$measures$cdmean >= 0 && res$measures$cdmean <= 1)
})

test_that("optimising a criterion improves it over a random set", {
  d <- mk()
  for (crit in c("cdmean", "pevmean", "aopt", "dopt")) {
    ex  <- select_individuals(d$G, n_train = 12, criterion = crit,
                              method = "exchange", iters = 300, seed = 7)
    rnd <- select_individuals(d$G, n_train = 12, criterion = crit,
                              method = "random", seed = 7)
    # `score` is reported in the native direction; the exchange maximises an
    # internally-signed version, so compare via that sign.
    maximise <- crit %in% c("cdmean", "cdmax", "neg_dist")
    if (maximise) expect_gte(ex$score, rnd$score, label = crit)
    else          expect_lte(ex$score, rnd$score, label = crit)
  }
})

test_that("the cdmean_exchange alias and mean_CD field still work", {
  d <- mk()
  res <- select_individuals(d$G, n_train = 10, method = "cdmean_exchange",
                            iters = 100, seed = 5)
  expect_true(is.numeric(res$mean_CD))
  expect_equal(res$method, "exchange")
})

# ---- M4: factor-analytic Sigma_E --------------------------------------------

test_that("fa_sigma_e approximates a covariance preserving its diagonal and PD", {
  d <- mk()
  for (k in 1:3) {
    Sfa <- fa_sigma_e(d$SigE, n_factors = k)
    expect_equal(diag(Sfa), diag(d$SigE), tolerance = 1e-8, info = k)
    expect_true(all(eigen(Sfa, symmetric = TRUE, only.values = TRUE)$values > 0),
                info = k)
  }
})

test_that("fa_sigma_e builds a covariance from loadings + psi", {
  L <- matrix(c(1, 0.8, 0.5, -0.2), 4, 1)
  S <- fa_sigma_e(loadings = L, psi = 0.3, n_factors = 1)
  expect_equal(dim(S), c(4L, 4L))
  expect_equal(S, L %*% t(L) + diag(0.3, 4), tolerance = 1e-10)
})

test_that("an FA-approximated Sigma_E is accepted by met_information", {
  d <- mk()
  full <- matrix(1L, d$J, d$E,
                 dimnames = list(rownames(d$G), paste0("E", seq_len(d$E))))
  Sfa <- fa_sigma_e(d$SigE, n_factors = 2)
  info <- met_information(full, G = d$G, Sigma_E = Sfa)
  expect_true(is.finite(info$mean_PEV))
})

# ---- M2: validate_design_criteria harness -----------------------------------

test_that("validate_design_criteria correlates criteria with outcomes", {
  d <- mk()
  ln <- rownames(d$G); en <- paste0("E", seq_len(d$E))
  full <- matrix(1L, d$J, d$E, dimnames = list(ln, en))
  set.seed(2)
  sp2 <- matrix(0L, d$J, d$E, dimnames = list(ln, en))
  for (j in seq_len(d$J)) sp2[j, sample(d$E, 2)] <- 1L
  sp3 <- matrix(0L, d$J, d$E, dimnames = list(ln, en))
  for (j in seq_len(d$J)) sp3[j, sample(d$E, 3)] <- 1L

  vd <- validate_design_criteria(
    list(full = full, sparse2 = sp2, sparse3 = sp3),
    G = d$G, Sigma_E = d$SigE, n_sim = 20, seed = 4)

  expect_equal(nrow(vd$table), 3L)
  expect_true(all(c("mean_PEV", "CDmean", "accuracy", "gain",
                    "common_selected", "avg_rank") %in% names(vd$table)))
  # lower design PEV should track higher realized accuracy across these designs
  r <- vd$correlations$correlation[
    vd$correlations$criterion == "mean_PEV" & vd$correlations$outcome == "accuracy"]
  expect_lt(r, 0)
})
