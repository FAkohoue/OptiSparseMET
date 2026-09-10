small_information_fixture <- function() {
  ids <- paste0("L", 1:12); envs <- paste0("E", 1:3)
  G <- diag(12); dimnames(G) <- list(ids, ids)
  S <- diag(3) * 0.4 + 0.6; dimnames(S) <- list(envs, envs)
  M <- matrix(0L, 12, 3, dimnames = list(ids, envs))
  M[1:8, 1] <- 1L; M[3:10, 2] <- 1L; M[5:12, 3] <- 1L
  list(M = M, G = G, S = S)
}

test_that("TPE weights, heterogeneous residuals, and local information are coherent", {
  d <- small_information_fixture()
  a <- met_information(d$M, d$G, d$S)
  b <- met_information(d$M, d$G, d$S,
                       tpe_weights = c(E3 = 1, E1 = 1, E2 = 1))
  expect_equal(a$mean_PEV, b$mean_PEV, tolerance = 1e-10)

  h1 <- met_information(d$M, d$G, d$S,
                        sigma_e2 = c(E1 = 1, E2 = 2, E3 = 4))
  h2 <- met_information(d$M, d$G, d$S,
                        env_efficiency = c(E1 = 1, E2 = .5, E3 = .25))
  expect_equal(h1$C_uu, h2$C_uu, tolerance = 1e-10)

  J <- nrow(d$M)
  blocks <- lapply(seq_len(ncol(d$M)), function(e) {
    idx <- (e - 1L) * J + seq_len(J)
    a$data_information[idx, idx, drop = FALSE]
  })
  names(blocks) <- colnames(d$M)
  blocks <- Map(function(L, e) {
    dimnames(L) <- list(rownames(d$M), rownames(d$M)); L
  }, blocks, seq_along(blocks))
  c <- met_information(d$M, d$G, d$S, local_information = blocks)
  expect_equal(a$mean_PEV, c$mean_PEV, tolerance = 1e-10)
  expect_true(c$local_information_used)
})

test_that("matrix-free information agrees with the dense engine", {
  d <- small_information_fixture()
  dense <- met_information(d$M, d$G, d$S, solver = "dense")
  pcg <- met_information(d$M, d$G, d$S, solver = "pcg",
                         trace_samples = 250L)
  expect_equal(pcg$solver_used, "pcg")
  expect_true(pcg$approximation)
  expect_equal(pcg$mean_PEV, dense$mean_PEV, tolerance = 0.04)
  expect_equal(pcg$CDmean, dense$CDmean, tolerance = 0.06)
  expect_equal(pcg$solver_diagnostics$converged_fraction, 1)
})

test_that("real field layouts produce contrast information", {
  fb <- data.frame(
    Treatment = rep(c("A", "B", "C"), 2),
    Block = rep(1:2, each = 3), Row = rep(1:2, each = 3),
    Column = rep(1:3, 2))
  L <- local_treatment_information(fb, fixed_effects = "Block")
  expect_equal(dim(L), c(3L, 3L))
  expect_equal(unname(rowSums(L)), rep(0, 3), tolerance = 1e-8)
  expect_true(min(eigen(L, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
})

test_that("optimizer consumes candidate replication and reports diagnostics", {
  d <- small_information_fixture()
  evaluator <- function(M) {
    reps <- M
    first <- which(M == 1L)[1]
    reps[first] <- 2L
    list(reps = reps, feasible = TRUE)
  }
  out <- optimize_design(d$M, d$G, d$S, design_evaluator = evaluator,
                         preserve = "margins", n_starts = 1L, iters = 8L,
                         seed = 1)
  expect_s3_class(out$design, "sparse_met_design")
  expect_equal(out$components$plots, sum(out$design$reps))
  expect_equal(out$diagnostics$proposals, 8L)
  expect_gt(out$diagnostics$unique_design_evaluations, 0L)
  expect_true(validate_sparse_met_design(out$design))
})

test_that("fieldbook evaluator builds the exact supplied allocation", {
  engine <- function(treatments, ...) {
    n <- length(treatments)
    list(field_book = data.frame(
      Treatment = treatments, Family = rep("F", n), Block = rep(1L, n),
      Plot = seq_len(n), Row = seq_len(n), Column = rep(1L, n)))
  }
  specs <- list(E1 = list(design = engine), E2 = list(design = engine))
  ev <- fieldbook_design_evaluator(specs, seed = 1L)
  M <- matrix(c(1, 1, 0, 1, 0, 1), 3, 2,
              dimnames = list(c("A", "B", "C"), c("E1", "E2")))
  ans <- ev(M)
  expect_equal(ans$reps, M)
  expect_named(ans$local_information, c("E1", "E2"))
  expect_equal(names(ans$fieldbooks), c("E1", "E2"))
  expect_equal(ans$metadata$plan_summary$allocation_method, "provided")
})

test_that("design schema catches fieldbook-replication disagreement", {
  d <- small_information_fixture()
  x <- sparse_met_design(d$M, G = d$G, Sigma_E = d$S)
  expect_true(validate_sparse_met_design(x))
  x$fieldbooks <- setNames(lapply(colnames(d$M), function(e)
    data.frame(Treatment = rownames(d$M)[d$M[, e] == 1L])), colnames(d$M))
  x$fieldbooks$E1 <- x$fieldbooks$E1[-1, , drop = FALSE]
  audit <- validate_sparse_met_design(x, error = FALSE)
  expect_false(audit$valid)
  expect_match(paste(audit$issues, collapse = " "), "does not reconcile")
})

test_that("historical REML fit returns a valid covariance", {
  set.seed(22)
  E <- 3L; J <- 45L; env <- paste0("E", seq_len(E))
  S <- matrix(c(1, .6, .25, .6, 1, .35, .25, .35, 1), E)
  U <- matrix(rnorm(J * E), J, E) %*% chol(S)
  dat <- expand.grid(Genotype = paste0("G", seq_len(J)),
                     Environment = env, Rep = 1:2,
                     KEEP.OUT.ATTRS = FALSE)
  idx <- cbind(match(dat$Genotype, paste0("G", seq_len(J))),
               match(dat$Environment, env))
  dat$Value <- U[idx] + rnorm(nrow(dat), sd = .45)
  fit <- fit_historical_met(dat, model = "unstructured")
  expect_s3_class(fit, "historical_met_fit")
  expect_equal(fit$convergence, 0)
  expect_true(all(is.finite(fit$Sigma_E)))
  expect_true(min(eigen(fit$Sigma_E, symmetric = TRUE,
                        only.values = TRUE)$values) > 0)
  expect_gt(fit$correlation[1, 2], 0.2)
})

test_that("replication recommendations are plantable integers", {
  d <- small_information_fixture()
  rr <- recommend_replication(d$M, d$G, d$S,
                              replication_levels = c(1, 1.5),
                              n_sim = 4L, seed = 3)
  expect_true(all(rr$replication_matrices[[2]] ==
                    round(rr$replication_matrices[[2]])))
  expect_equal(rr$table$achieved_replication[2], 1.5)
})
