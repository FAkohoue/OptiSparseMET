# Tests for construct_exact_bibd(): a true BIBD is only buildable in the small-J
# regime where the necessary conditions hold. The classic case is the symmetric
# Fano-plane design v = b = 7, r = k = 3, lambda = 1.

test_that("the Fano-plane BIBD (v=b=7, r=k=3, lambda=1) is constructed and verified", {
  skip_if_not_installed("crossdes")
  out <- construct_exact_bibd(
    treatments   = paste0("L", 1:7),
    environments = paste0("E", 1:7),
    k            = 3
  )
  M <- out$allocation_matrix
  expect_equal(dim(M), c(7L, 7L))
  expect_equal(out$r, 3)
  expect_equal(out$lambda, 1)
  expect_true(out$is_bibd)

  # Structural guarantees of a BIBD, checked directly from the incidence matrix.
  expect_true(all(rowSums(M) == 3))          # equal replication
  expect_true(all(colSums(M) == 3))          # equal block size
  conc <- tcrossprod(M)                      # pairwise concurrence
  expect_true(all(conc[upper.tri(conc)] == 1))   # constant lambda = 1
})

test_that("Fisher's inequality is enforced", {
  skip_if_not_installed("crossdes")
  # J = 16, I = 8, k = 6 satisfies BOTH counting conditions -- r = I*k/J = 3 and
  # lambda = r(k-1)/(J-1) = 1 are integers -- so Fisher's inequality (I >= J) is
  # the only condition left to fail. This isolates the Fisher check; using
  # parameters that also break the counting conditions would trip an earlier
  # error instead.
  expect_error(
    construct_exact_bibd(paste0("L", 1:16), paste0("E", 1:8), k = 6),
    regexp = "Fisher"
  )
})

test_that("non-integer lambda is rejected", {
  skip_if_not_installed("crossdes")
  # v = 10, b = 10, k = 5 -> r = 5, lambda = 5*4/9 = 20/9, not an integer.
  expect_error(
    construct_exact_bibd(paste0("L", 1:10), paste0("E", 1:10), k = 5),
    regexp = "lambda"
  )
})

test_that("indivisible replication is rejected", {
  skip_if_not_installed("crossdes")
  # I*k = 5*3 = 15 is not divisible by J = 7, so r is not an integer.
  expect_error(
    construct_exact_bibd(paste0("L", 1:7), paste0("E", 1:5), k = 3),
    regexp = "divisible"
  )
})

test_that("degenerate dimensions are rejected", {
  skip_if_not_installed("crossdes")
  expect_error(construct_exact_bibd(paste0("L", 1:7), paste0("E", 1:7), k = 7))
  expect_error(construct_exact_bibd("L1", c("E1", "E2"), k = 1))
})

test_that("check_equireplicate_feasibility agrees on when a BIBD is possible", {
  # The Fano parameters: feasible AND a strict BIBD.
  f_ok <- check_equireplicate_feasibility(
    n_treatments_total = 7, n_environments = 7,
    n_test_entries_per_environment = 3, target_replications = 3)
  expect_true(f_ok$strict_bibd_possible)

  # A realistic sparse-testing case: slot identity holds, BIBD impossible.
  f_no <- check_equireplicate_feasibility(
    n_treatments_total = 110, n_environments = 4,
    n_test_entries_per_environment = 55, target_replications = 2)
  expect_true(f_no$feasible)
  expect_false(f_no$strict_bibd_possible)
})
