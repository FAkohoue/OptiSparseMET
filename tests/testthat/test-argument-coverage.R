# Argument- and branch-coverage tests. These deliberately exercise the
# arguments, alternative methods, and error paths of the exported functions
# that the topic-focused test files touch on only lightly, so that every
# documented argument is executed at least once.

mk <- function(J = 40L, E = 4L, seed = 11) {
  set.seed(seed)
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  SigE <- matrix(c(1, .8, .1, 0,  .8, 1, 0, .1,
                   .1, 0, 1, .8,   0, .1, .8, 1), 4, 4)
  dimnames(SigE) <- list(paste0("E", 1:4), paste0("E", 1:4))
  list(G = G, SigE = SigE, J = J, E = E)
}

# full 4-environment allocation (every entry in every environment)
full_alloc <- function(d) {
  matrix(1L, d$J, d$E,
         dimnames = list(rownames(d$G), colnames(d$SigE)))
}

# =============================================================================
# build_relationship_matrix: type dispatch, error paths, omega/tau
# =============================================================================

test_that("build_relationship_matrix rejects missing / malformed inputs", {
  expect_error(build_relationship_matrix(type = "genomic"),
               "markers.*required")
  expect_error(build_relationship_matrix(type = "pedigree"),
               "pedigree.*required")
  M <- matrix(sample(0:2, 20, TRUE), 4, 5)
  expect_error(build_relationship_matrix(markers = M, type = "hmatrix"),
               "both required")
  bad_ped <- data.frame(a = 1:3, b = 0:2, c = 0:2)
  expect_error(build_relationship_matrix(pedigree = bad_ped, type = "pedigree"),
               "columns id, sire, dam")
})

test_that("omega blends more pedigree into H as it grows", {
  set.seed(7)
  ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 3, 0),
                    dam = c(0, 0, 2, 3, 4, 0))
  M <- matrix(sample(0:2, 3 * 120, replace = TRUE), 3, 120,
              dimnames = list(c("3", "4", "5"), NULL))
  H_lo <- build_relationship_matrix(markers = M, pedigree = ped,
                                    type = "hmatrix", omega = 0.05)
  H_hi <- build_relationship_matrix(markers = M, pedigree = ped,
                                    type = "hmatrix", omega = 0.50)
  # both valid H matrices, but a different blend => different genotyped block
  expect_false(isTRUE(all.equal(H_lo, H_hi)))
  expect_equal(H_hi, t(H_hi), tolerance = 1e-8)
})

# =============================================================================
# build_enviromic_covariates: standardize, direct soil, error paths
# =============================================================================

test_that("standardize = TRUE centres and scales numeric covariates", {
  sites <- data.frame(environment = c("E1", "E2", "E3", "E4"))
  wx <- data.frame(environment = sites$environment,
                   mean_temp = c(18, 22, 26, 30),
                   precip    = c(500, 400, 300, 250))
  X <- build_enviromic_covariates(sites, weather = wx, standardize = TRUE)
  expect_equal(unname(colMeans(X)), c(0, 0), tolerance = 1e-8)
  expect_equal(rownames(X), sites$environment)
})

test_that("directly supplied soil is merged and no-source input errors", {
  sites <- data.frame(environment = c("E1", "E2", "E3"))
  soil  <- data.frame(environment = sites$environment,
                      clay = c(30, 25, 40), ph = c(6.1, 6.4, 5.8))
  X <- build_enviromic_covariates(sites, soil = soil, standardize = FALSE)
  expect_true(all(c("clay", "ph") %in% colnames(X)))
  expect_error(build_enviromic_covariates(sites), "No covariate sources")
  expect_error(build_enviromic_covariates(data.frame(x = 1)),
               "environment.*column")
})

# =============================================================================
# cluster_environments: kmeans method, supplied Sigma_E, bounds
# =============================================================================

test_that("kmeans clustering recovers the two blocks", {
  D <- matrix(c(1, .9, .1, .1,  .9, 1, .1, .1,
                .1, .1, 1, .9,  .1, .1, .9, 1), 4, 4,
              dimnames = list(paste0("E", 1:4), paste0("E", 1:4)))
  cl <- cluster_environments(D, n_clusters = 2, method = "kmeans")
  expect_equal(cl$membership[["E1"]], cl$membership[["E2"]])
  expect_false(cl$membership[["E1"]] == cl$membership[["E3"]])
})

test_that("a supplied Sigma_E is subset per cluster and bounds are enforced", {
  D <- matrix(c(1, .9, .1, .1,  .9, 1, .1, .1,
                .1, .1, 1, .9,  .1, .1, .9, 1), 4, 4,
              dimnames = list(paste0("E", 1:4), paste0("E", 1:4)))
  S <- diag(c(2, 3, 4, 5)); dimnames(S) <- dimnames(D)
  cl <- cluster_environments(D, n_clusters = 2, method = "hclust", Sigma_E = S)
  # each cluster's covariance is the matching 2x2 diagonal sub-block of S
  vals <- unname(sort(unlist(lapply(cl$Sigma_E, diag))))
  expect_equal(vals, c(2, 3, 4, 5))
  expect_error(cluster_environments(D, n_clusters = 5), "between 1 and")
})

# =============================================================================
# select_individuals: every criterion, target set, validation, reproducibility
# =============================================================================

test_that("all eight criteria return a valid selection and full measures", {
  d <- mk(J = 30L)
  for (crit in c("cdmean", "cdmax", "pevmean", "pevmax",
                 "aopt", "dopt", "goptpev", "neg_dist")) {
    res <- select_individuals(d$G, n_train = 10, criterion = crit,
                              method = "exchange", iters = 120, seed = 3)
    expect_length(res$selected, 10L)
    expect_named(res$measures,
                 c("cdmean", "cdmax", "pevmean", "pevmax",
                   "aopt", "dopt", "goptpev", "neg_dist"),
                 ignore.order = TRUE)
    expect_true(is.finite(res$score))
  }
})

test_that("a target set restricts the CD evaluation and n_train is validated", {
  d <- mk(J = 30L)
  tgt <- rownames(d$G)[1:6]
  res <- select_individuals(d$G, n_train = 10, target = tgt,
                            criterion = "cdmean", method = "exchange",
                            iters = 150, seed = 2)
  expect_length(res$selected, 10L)
  expect_true(res$measures$cdmean >= 0 && res$measures$cdmean <= 1)
  expect_error(select_individuals(d$G, n_train = 30), "1 <= n_train")
})

test_that("select_individuals is reproducible under a fixed seed", {
  d <- mk(J = 30L)
  a <- select_individuals(d$G, n_train = 8, method = "exchange", seed = 99)$selected
  b <- select_individuals(d$G, n_train = 8, method = "exchange", seed = 99)$selected
  expect_identical(a, b)
})

# =============================================================================
# select_environments: out-of-range n is rejected
# =============================================================================

test_that("select_environments rejects n greater than the pool", {
  D <- diag(4); dimnames(D) <- list(paste0("E", 1:4), paste0("E", 1:4))
  expect_error(select_environments(D, n = 5))
})

# =============================================================================
# met_information: environment_specific target branch
# =============================================================================

test_that("met_information computes a per-cell (environment_specific) summary", {
  d <- mk()
  M <- full_alloc(d)
  info <- met_information(M, G = d$G, Sigma_E = d$SigE,
                          target = "environment_specific")
  expect_equal(info$target, "environment_specific")
  expect_true(is.finite(info$mean_PEV))
})

# =============================================================================
# sensitivity_varcomp: custom lambda grid, structure, monotone PEV
# =============================================================================

test_that("sensitivity_varcomp returns sorted unique lambdas with finite PEV", {
  d <- mk()
  M <- full_alloc(d)
  tab <- sensitivity_varcomp(M, G = d$G, Sigma_E = d$SigE,
                             lambda_grid = c(2, 0.5, 0.5, 1))
  expect_named(tab, c("lambda", "mean_PEV", "CDmean"))
  expect_equal(tab$lambda, c(0.5, 1, 2))          # sorted + de-duplicated
  expect_true(all(is.finite(tab$mean_PEV)))
  expect_true(tab$mean_PEV[3] >= tab$mean_PEV[1]) # more resid. variance -> more PEV
})

# =============================================================================
# sparsity_grid: plot_budget filter and best-first ordering
# =============================================================================

test_that("sparsity_grid honours plot_budget and returns ranked rows", {
  d <- mk()
  tab <- sparsity_grid(G = d$G, Sigma_E = d$SigE,
                       ia_values = c(20, 40), ipf_values = c(10, 20),
                       noe_values = c(2, 4), plot_budget = 120,
                       n_sim = 10, seed = 1)
  expect_true(all(tab$total_plots <= 120))
  expect_equal(tab$accuracy_mean, sort(tab$accuracy_mean, decreasing = TRUE))
  expect_error(
    sparsity_grid(G = d$G, Sigma_E = d$SigE, ia_values = 20,
                  ipf_values = 10, noe_values = 2, plot_budget = 1),
    "No feasible grid cells")
})

# =============================================================================
# recommend_replication: custom levels + seed-budget feasibility flag
# =============================================================================

test_that("recommend_replication reports the requested levels and feasibility", {
  d <- mk()
  M <- full_alloc(d)
  rr <- recommend_replication(M, G = d$G, Sigma_E = d$SigE,
                              replication_levels = c(1, 2, 3),
                              n_sim = 12, seed = 1)
  expect_equal(rr$table$replication, c(1, 2, 3))
  expect_true("feasible" %in% names(rr$table))
  expect_true(is.na(rr$recommended) || rr$recommended %in% rr$table$replication)
})

# =============================================================================
# optimize_allocation_gxe: reproducibility and margin preservation
# =============================================================================

test_that("optimize_allocation_gxe is reproducible and preserves column sums", {
  d <- mk()
  set.seed(5)
  M <- matrix(0L, d$J, d$E, dimnames = list(rownames(d$G), colnames(d$SigE)))
  for (j in seq_len(d$E)) M[sample(d$J, 20), j] <- 1L
  o1 <- optimize_allocation_gxe(M, G = d$G, Sigma_E = d$SigE, iter = 50, seed = 8)
  o2 <- optimize_allocation_gxe(M, G = d$G, Sigma_E = d$SigE, iter = 50, seed = 8)
  expect_identical(o1$allocation_matrix, o2$allocation_matrix)
  expect_equal(colSums(o1$allocation_matrix), colSums(M))
})

# =============================================================================
# run_design_strategy: evaluate toggle, recommend_reps, common treatments
# =============================================================================

test_that("run_design_strategy skips evaluation when evaluate = FALSE", {
  d <- mk()
  out <- run_design_strategy(
    treatments = rownames(d$G),
    environments = colnames(d$SigE),
    n_test_entries_per_environment = 20,
    G = d$G, Sigma_E = d$SigE,
    allocation_method = "random_balanced",
    evaluate = FALSE, seed = 1)
  expect_null(out$evaluation)
  expect_equal(dim(out$allocation_matrix), c(d$J, 4L))
})

test_that("run_design_strategy runs replication advice and keeps common entries", {
  d <- mk()
  common <- rownames(d$G)[1:4]
  out <- run_design_strategy(
    treatments = rownames(d$G),
    environments = colnames(d$SigE),
    n_test_entries_per_environment = 20,
    G = d$G, Sigma_E = d$SigE,
    allocation_method = "random_balanced",
    common_treatments = common,
    recommend_reps = TRUE, replication_levels = c(1, 2),
    evaluate = TRUE, n_sim = 10, seed = 1)
  expect_false(is.null(out$replication))
  # common treatments appear in every environment
  expect_true(all(rowSums(out$allocation_matrix[common, , drop = FALSE]) == 4))
})

test_that("run_design_strategy also resolves environments from a relationship", {
  d <- mk()
  set.seed(3)
  GxE <- matrix(rnorm(d$J * d$E), d$J, d$E,
                dimnames = list(NULL, colnames(d$SigE)))
  Denv <- build_environment_relationship(GxE, source = "genetic_correlation")
  out <- run_design_strategy(
    treatments = rownames(d$G),
    n_test_entries_per_environment = 20,
    env_relationship = Denv, n_environments = 3,
    env_select_method = "representative",
    G = d$G, Sigma_E = d$SigE,
    allocation_method = "random_balanced",
    evaluate = FALSE, seed = 1)
  expect_length(out$decisions$environments, 3L)
})

# =============================================================================
# fa_sigma_e: build-from-loadings validation and no-input error
# =============================================================================

test_that("fa_sigma_e validates its inputs", {
  expect_error(fa_sigma_e(n_factors = 0), "positive integer")
  expect_error(fa_sigma_e(), "Sigma_E.*or.*loadings")
  L <- matrix(c(1, .8, .6, .2), ncol = 1)
  S <- fa_sigma_e(loadings = L, psi = 0.1)
  expect_equal(dim(S), c(4L, 4L))
  expect_true(all(eigen(S, symmetric = TRUE, only.values = TRUE)$values > 0))
})

# =============================================================================
# compare_spatial_models: famoptg evaluator + P-spline model
# =============================================================================

test_that("compare_spatial_models runs the famoptg evaluator including pspline", {
  n_rows <- 6L; n_cols <- 6L; N <- n_rows * n_cols
  fb <- data.frame(
    Plot      = seq_len(N),
    Row       = rep(seq_len(n_rows), times = n_cols),
    Column    = rep(seq_len(n_cols), each = n_rows),
    Block     = rep(1:6, each = N / 6),
    Treatment = rep(paste0("G", 1:9), length.out = N),
    Family    = "F1",
    stringsAsFactors = FALSE)
  vc <- list(sigma_e2 = 2, sigma_g2 = 1, sigma_b2 = 0.4,
             sigma_r2 = 0.2, sigma_c2 = 0.2)
  res <- compare_spatial_models(
    field_book = fb, n_rows = n_rows, n_cols = n_cols,
    check_treatments = character(0),
    models = c("IID", "AR1xAR1_nugget", "pspline"),
    evaluator = "famoptg",
    spline_lambda_row = 1e6, spline_lambda_col = 1e6,
    treatment_effect = "random", prediction_type = "IID", varcomp = vc,
    spline_knots_row = 4, spline_knots_col = 4)
  expect_equal(nrow(res), 3L)
  expect_true(all(c("model", "A_criterion", "CDmean", "mean_PEV") %in% names(res)))
  expect_true(all(is.na(res$error)))
})
