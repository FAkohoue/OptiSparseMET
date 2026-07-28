# Tests for the three field-capacity / limiting-variable features:
#  1. variable selection + weighting in build_environment_relationship()
#  2. named / variable per-site capacity in allocate_sparse_met()
#  3. suggest_site_capacity() for an unconstrained site

# ---- Feature 1: limiting-variable selection and weighting -------------------

# Temperature/precip group {E1,E2}; clay groups {E1,E3}. So the choice of
# limiting variable(s) genuinely flips which environment E1 is closest to.
mk_cov <- function() {
  matrix(c(20, 400, 30,   # E1
           21, 390, 55,   # E2
           28, 260, 32,   # E3
           29, 250, 54),  # E4
         nrow = 4, byrow = TRUE,
         dimnames = list(paste0("E", 1:4),
                         c("mean_temp", "precip", "clay")))
}

test_that("variables restricts the kinship to the chosen covariates", {
  X <- mk_cov()
  D_all  <- build_environment_relationship(X, source = "enviromic",
                                           kernel = "gaussian")
  D_temp <- build_environment_relationship(X, source = "enviromic",
                                           kernel = "gaussian",
                                           variables = "mean_temp")
  expect_equal(dim(D_temp), c(4L, 4L))
  # using only temperature changes the similarity structure
  expect_false(isTRUE(all.equal(D_all, D_temp)))
  # by temperature, E1 (20) is closest to E2 (21), not to E3 (28)
  expect_gt(D_temp["E1", "E2"], D_temp["E1", "E3"])
})

test_that("unknown variable names and bad weights are rejected", {
  X <- mk_cov()
  expect_error(
    build_environment_relationship(X, source = "enviromic",
                                   variables = "rainfall"),
    "not found")
  expect_error(
    build_environment_relationship(X, source = "enviromic",
                                   weights = c(1, 2)),   # wrong length (3 cols)
    "one value per")
  expect_error(
    build_environment_relationship(X, source = "enviromic",
                                   weights = c(mean_temp = 1, precip = -1,
                                               clay = 1)),
    "non-negative")
})

test_that("weights emphasise the chosen limiting factor", {
  X <- mk_cov()
  # Heavily weight clay: E1 (30) should be most similar to E3 (32).
  D_w <- build_environment_relationship(
    X, source = "enviromic", kernel = "gaussian",
    weights = c(mean_temp = 0.01, precip = 0.01, clay = 100))
  expect_gt(D_w["E1", "E3"], D_w["E1", "E2"])
  expect_equal(diag(D_w), rep(1, 4), tolerance = 1e-8, ignore_attr = TRUE)
})

# ---- Feature 2: named / variable per-site capacity --------------------------

test_that("a named capacity vector is aligned to environment order", {
  set.seed(1)
  trt <- paste0("L", 1:40)
  envs <- c("Loc_A", "Loc_B", "Loc_C")
  cap <- c(Loc_C = 12, Loc_A = 20, Loc_B = 8)   # deliberately out of order
  al <- allocate_sparse_met(
    treatments = trt, environments = envs,
    allocation_method = "random_balanced",
    n_test_entries_per_environment = cap, seed = 1)
  sizes <- colSums(al$allocation_matrix)
  expect_equal(sizes[["Loc_A"]], 20)
  expect_equal(sizes[["Loc_B"]], 8)
  expect_equal(sizes[["Loc_C"]], 12)
})

test_that("a missing environment in a named capacity vector errors", {
  trt <- paste0("L", 1:30); envs <- c("A", "B", "C")
  expect_error(
    allocate_sparse_met(treatments = trt, environments = envs,
                        allocation_method = "random_balanced",
                        n_test_entries_per_environment = c(A = 10, B = 8)),
    "missing capacities")
})

test_that("equireplicate keeps EQUAL REPLICATION with UNEQUAL site sizes (M4)", {
  set.seed(3)
  trt <- paste0("L", 1:40); envs <- c("A", "B", "C", "D")
  cap <- c(20, 20, 16, 24)                 # unequal, sum = 80 = 40 * 2
  al <- allocate_sparse_met(
    treatments = trt, environments = envs,
    allocation_method = "equireplicate",
    n_test_entries_per_environment = cap,
    target_replications = 2, seed = 3)
  M <- al$allocation_matrix
  expect_equal(unname(colSums(M)), cap)    # site-specific (unequal) sizes honoured
  expect_true(all(rowSums(M) == 2))        # equal replication preserved
})

test_that("equireplicate still errors when slot totals are infeasible", {
  trt <- paste0("L", 1:40); envs <- c("A", "B", "C", "D")
  expect_error(
    allocate_sparse_met(treatments = trt, environments = envs,
                        allocation_method = "equireplicate",
                        n_test_entries_per_environment = c(20, 20, 16, 25),
                        target_replications = 2),
    "infeasible")                          # sum = 81 != 40 * 2
})

test_that(".equireplicate_degree_feasible implements the Gale-Ryser condition", {
  expect_true(.equireplicate_degree_feasible(c(20, 20, 16, 24), r = 2, n_sparse = 40))
  expect_false(.equireplicate_degree_feasible(c(20, 20, 16, 24), r = 3, n_sparse = 40)) # sum
  expect_false(.equireplicate_degree_feasible(c(6, 1, 1), r = 2, n_sparse = 4))         # k>n
})

test_that("positional variable capacities produce unequal site sizes", {
  set.seed(2)
  trt <- paste0("L", 1:50); envs <- c("A", "B", "C")
  al <- allocate_sparse_met(
    treatments = trt, environments = envs,
    allocation_method = "random_balanced",
    n_test_entries_per_environment = c(25, 10, 15), seed = 2)
  expect_equal(unname(colSums(al$allocation_matrix)), c(25, 10, 15))
})

# ---- Feature 3: suggest_site_capacity ---------------------------------------

mk_gs <- function(J = 60L, seed = 3) {
  set.seed(seed)
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  SigE <- diag(4) * 0.4 + 0.6
  dimnames(SigE) <- list(paste0("E", 1:4), paste0("E", 1:4))
  list(G = G, SigE = SigE)
}

test_that("suggest_site_capacity returns a ranked table and a recommendation", {
  d <- mk_gs()
  res <- suggest_site_capacity(
    d$G, d$SigE, candidate_plots = c(10, 20, 40, 60),
    other_n_env = 2, other_capacity = 20, n_sim = 10, seed = 1)
  expect_true(all(c("plots", "total_plots", "accuracy_mean", "gain_mean",
                    "marginal_accuracy_per_plot") %in% names(res$table)))
  expect_equal(res$table$plots, c(10, 20, 40, 60))     # sorted
  expect_true(res$recommended_plots %in% res$table$plots)
  # accuracy should not fall as the focal site tests more of the TPG
  expect_gte(res$table$accuracy_mean[nrow(res$table)],
             res$table$accuracy_mean[1] - 0.05)
})

test_that("a strict min_gain recommends fewer plots than a lax one", {
  d <- mk_gs()
  strict <- suggest_site_capacity(d$G, d$SigE, candidate_plots = c(10, 20, 40, 60),
                                  other_capacity = 20, n_sim = 12,
                                  select = "diminishing",
                                  min_gain = 1, seed = 1)$recommended_plots
  lax    <- suggest_site_capacity(d$G, d$SigE, candidate_plots = c(10, 20, 40, 60),
                                  other_capacity = 20, n_sim = 12,
                                  select = "diminishing",
                                  min_gain = 0, seed = 1)$recommended_plots
  expect_lte(strict, lax)
})

test_that("cost_per_plot adds a gain-per-cost column and inputs are validated", {
  d <- mk_gs()
  res <- suggest_site_capacity(d$G, d$SigE, candidate_plots = c(10, 20),
                               other_capacity = 15, cost_per_plot = 5,
                               n_sim = 8, seed = 1)
  expect_true("gain_per_cost" %in% names(res$table))
  expect_error(
    suggest_site_capacity(d$G, d$SigE, candidate_plots = c(10, 200),
                          other_capacity = 15, n_sim = 5),
    "fewer genotypes")
})

test_that("scope = 'all' optimises every site and scales total plots", {
  d <- mk_gs()
  res <- suggest_site_capacity(d$G, d$SigE, candidate_plots = c(10, 20, 40),
                               scope = "all", n_env = 3, n_sim = 10, seed = 1)
  expect_equal(res$table$plots, c(10, 20, 40))
  expect_equal(res$table$total_plots, c(10, 20, 40) * 3)   # applied to all sites
  expect_true(res$recommended_plots %in% c(10, 20, 40))
})

test_that("select rule picks within the candidate interval", {
  d <- mk_gs()
  cp <- c(20, 30, 40, 50)                       # <= nrow(G) = 60
  rmax <- suggest_site_capacity(d$G, d$SigE, candidate_plots = cp, scope = "all",
                                n_env = 3, select = "max", n_sim = 10, seed = 1)
  # "max" exploits capacity -> the candidate with the highest accuracy
  expect_equal(rmax$recommended_plots,
               rmax$table$plots[which.max(rmax$table$accuracy_mean)])
  expect_true(rmax$recommended_plots %in% cp)
  # a very strict knee threshold falls back to the smallest candidate
  rknee <- suggest_site_capacity(d$G, d$SigE, candidate_plots = cp, scope = "all",
                                 n_env = 3, select = "diminishing", min_gain = 1,
                                 n_sim = 10, seed = 1)
  expect_equal(rknee$recommended_plots, min(cp))
})

test_that("scope = 'all' validates n_env against Sigma_E", {
  d <- mk_gs()   # Sigma_E is 4 x 4
  expect_error(
    suggest_site_capacity(d$G, d$SigE, candidate_plots = c(10, 20),
                          scope = "all", n_env = 9, n_sim = 5),
    "at least")
})

test_that("scope = 'subset' holds partner capacities and counts real plots", {
  d <- mk_gs()
  caps <- c(E1 = 10, E2 = 10, E3 = 20, E4 = 25)
  res <- suggest_site_capacity(
    d$G, d$SigE, candidate_plots = c(10, 20),
    scope = "subset", focal_envs = c("E1", "E2"),
    site_capacities = caps,
    plots_per_entry = c(E1 = 2, E2 = 2, E3 = 1, E4 = 1),
    check_plots_per_site = 4, n_sim = 6, seed = 3)
  expect_equal(res$table$total_plots, c(101, 141))
  expect_equal(diff(res$table$total_plots), 40)
  expect_equal(res$table$marginal_accuracy_per_plot[2],
               diff(res$table$accuracy_mean) / 40)
})

test_that("capacity selection uses maximin covariance uncertainty", {
  d <- mk_gs()
  independent <- diag(4)
  dimnames(independent) <- dimnames(d$SigE)
  scenarios <- robust_scenarios(
    sigma_g2 = 1, sigma_e2 = 1, sigmaE_shrink = 0,
    Sigma_E_candidates = list(central = d$SigE,
                              interaction = independent)
  )
  robust <- suggest_site_capacity(
    d$G, d$SigE, candidate_plots = c(10, 20),
    other_n_env = 3, other_capacity = 15,
    robust = scenarios, robust_aggregate = "min",
    n_sim = 4, seed = 7
  )
  central <- suggest_site_capacity(
    d$G, d$SigE, candidate_plots = c(10, 20),
    other_n_env = 3, other_capacity = 15, n_sim = 4, seed = 7
  )
  interaction <- suggest_site_capacity(
    d$G, independent, candidate_plots = c(10, 20),
    other_n_env = 3, other_capacity = 15, n_sim = 4, seed = 7
  )
  expect_true(robust$robust)
  expect_equal(robust$robust_aggregate, "min")
  expect_equal(robust$n_scenarios, 2L)
  expect_equal(
    robust$table$accuracy_mean,
    pmin(central$table$accuracy_mean, interaction$table$accuracy_mean)
  )
})

# ---- Feature 1 (cont.): variable catalogue ----------------------------------

test_that("enviromic_variable_catalog documents fetchable names", {
  cat_all <- enviromic_variable_catalog("all")
  expect_true(all(c("variable", "source", "description", "units") %in%
                    names(cat_all)))
  expect_true(all(c("mean_temp", "total_precip", "clay", "phh2o", "bdod") %in%
                    cat_all$variable))
  expect_setequal(unique(enviromic_variable_catalog("soil")$source), "soil")
  expect_gte(nrow(enviromic_variable_catalog("weather")), 6L)
})

test_that("an unknown variable name lists the available columns", {
  X <- mk_cov()   # columns: mean_temp, precip, clay
  expect_error(
    build_environment_relationship(X, source = "enviromic",
                                   variables = "phh2o"),
    "Available columns")
})
