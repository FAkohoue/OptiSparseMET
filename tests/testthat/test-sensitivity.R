# Test for sensitivity_varcomp() (P3.5).

test_that("mean PEV increases with the residual/genetic variance ratio", {
  set.seed(4)
  J <- 20L; E <- 3L
  B <- matrix(rnorm(J * J), J, J); G <- tcrossprod(B) / J + diag(J) * 0.4
  dimnames(G) <- list(paste0("L", seq_len(J)), paste0("L", seq_len(J)))
  M <- matrix(1L, J, E, dimnames = list(rownames(G), paste0("E", seq_len(E))))

  tab <- sensitivity_varcomp(M, G = G, lambda_grid = c(0.5, 1, 2, 4))
  expect_equal(nrow(tab), 4L)
  # PEV is monotone increasing in lambda; CDmean monotone decreasing.
  expect_true(all(diff(tab$mean_PEV) > 0))
  expect_true(all(diff(tab$CDmean) < 0))
})
