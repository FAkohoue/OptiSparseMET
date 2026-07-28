# Tests for desired-gain multi-trait index weights (Pesek-Baker built-in and the
# DGQGSI-result bridge). Reference values verified independently in numpy.

Gt <- matrix(c(1, .3, -.2, .3, .8, .1, -.2, .1, .5), 3,
             dimnames = list(c("YLD", "DIS", "MOI"), c("YLD", "DIS", "MOI")))

test_that("Pesek-Baker weights make realized gains proportional to desired gains", {
  dg <- c(YLD = 1.5, DIS = 1.0, MOI = 1.0)
  w <- desired_gain_weights(dg, gen_cov = Gt, lower_is_better = "DIS")
  # numpy reference b for d_signed = c(1.5, -1, 1) (DIS sign-flipped)
  expect_equal(unname(w$weights), c(3.139535, -2.906977, 3.837209),
               tolerance = 1e-5)
  # G b must equal the signed desired gains exactly
  expect_equal(as.numeric(Gt %*% w$weights), c(1.5, -1.0, 1.0), tolerance = 1e-8)
  expect_equal(w$sigma_index, 3.3843, tolerance = 1e-4)
})

test_that("sign convention: lower_is_better negates that trait's desired gain", {
  dg <- c(YLD = 1, DIS = 1, MOI = 1)
  w_none <- desired_gain_weights(dg, gen_cov = Gt)
  w_dis  <- desired_gain_weights(dg, gen_cov = Gt, lower_is_better = "DIS")
  expect_equal(unname(w_none$desired_gain), c(1, 1, 1))
  expect_equal(unname(w_dis$desired_gain),  c(1, -1, 1))
  expect_false(isTRUE(all.equal(w_none$weights, w_dis$weights)))
})

test_that("gen_cov is aligned to the desired-gain trait order", {
  dg <- c(MOI = 1, YLD = 2, DIS = 1)                 # deliberately reordered
  w <- desired_gain_weights(dg, gen_cov = Gt)
  expect_equal(names(w$weights), c("MOI", "YLD", "DIS"))
  # weights solve the reordered system
  Gre <- Gt[names(dg), names(dg)]
  expect_equal(as.numeric(Gre %*% w$weights), unname(w$desired_gain),
               tolerance = 1e-8)
})

test_that("desired-gain weights plug into the multi-trait design objective", {
  dg <- c(YLD = 1.5, DIS = 1.0, MOI = 1.0)
  w <- desired_gain_weights(dg, gen_cov = Gt, lower_is_better = "DIS")

  set.seed(1)
  G <- crossprod(matrix(rnorm(20 * 50), 50, 20)) / 50 + diag(20) * 0.3
  dimnames(G) <- list(paste0("L", 1:20), paste0("L", 1:20))
  SigE <- diag(3) * 0.5 + 0.5
  dimnames(SigE) <- list(paste0("E", 1:3), paste0("E", 1:3))
  M <- matrix(0L, 20, 3, dimnames = list(rownames(G), colnames(SigE)))
  for (e in 1:3) M[sample(20, 10), e] <- 1L

  gg <- expected_genetic_gain(reliability = 0.5,
                              trait_weights = w$weights, trait_gencov = w$gen_cov)
  expect_equal(gg$sigma, w$sigma_index, tolerance = 1e-6)   # index SD carried through

  o <- design_objective(M, G, SigE,
                        trait_weights = w$weights, trait_gencov = w$gen_cov)
  expect_true(is.finite(o$gain) && o$gain > 0)
})

test_that("the DGQGSI bridge extracts weights and errors helpfully", {
  # mock a DGQGSI result shape (nested dg_result with a weight vector)
  res <- list(dg_result = list(weights = c(0.5, -0.2, 0.7),
                               best_q = 0.03, mean_replicate_cor = 0.95))
  w <- desired_gain_weights(c(YLD = 1, DIS = 1, MOI = 1), gen_cov = Gt,
                            method = "dgqgsi", dgqgsi_result = res)
  expect_equal(unname(w$weights), c(0.5, -0.2, 0.7))

  bad <- list(dg_result = list(best_q = 0.03))    # no weight vector
  expect_error(
    desired_gain_weights(c(YLD = 1, DIS = 1, MOI = 1), gen_cov = Gt,
                         method = "dgqgsi", dgqgsi_result = bad),
    "Could not find")
  expect_error(
    desired_gain_weights(c(YLD = 1, DIS = 1, MOI = 1), gen_cov = Gt,
                         method = "dgqgsi"),
    "needs `dgqgsi_result`")
})
