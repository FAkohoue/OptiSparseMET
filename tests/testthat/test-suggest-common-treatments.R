# Tests for the auto common-treatment suggester (seed feasibility + family /
# genetic-group representativeness + connectivity).

test_that("seed feasibility filters the pool and all families are covered", {
  set.seed(1)
  J <- 60L
  G <- crossprod(matrix(rnorm(J * 80), 80, J)) / 80 + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", 1:J), paste0("L", 1:J))
  ti <- data.frame(Treatment = rownames(G),
                   Family = paste0("F", (seq_len(J) %% 6) + 1),
                   stringsAsFactors = FALSE)
  si <- data.frame(Treatment = rownames(G),
                   SeedAvailable = c(rep(300, 40), rep(20, 20)))   # 20 low-seed
  res <- suggest_common_treatments(G, ti, si, seed_required_per_plot = c(10, 16, 20),
                                   reps_per_site = 2, min_per_family = 1, seed = 1)
  expect_equal(res$seed_need, (10 + 16 + 20) * 2)          # 92 g
  expect_equal(res$n_feasible, 40L)                        # only the high-seed lines
  expect_true(all(res$selected %in% rownames(G)[1:40]))    # never picks infeasible
  fams_sel <- ti$Family[match(res$selected, ti$Treatment)]
  expect_true(all(unique(ti$Family) %in% fams_sel))        # every family represented
})

test_that("a user-supplied n_common is honoured", {
  set.seed(2)
  J <- 40L
  G <- crossprod(matrix(rnorm(J * 60), 60, J)) / 60 + diag(J) * 0.3
  dimnames(G) <- list(paste0("L", 1:J), paste0("L", 1:J))
  ti <- data.frame(Treatment = rownames(G),
                   Family = paste0("F", (seq_len(J) %% 4) + 1), stringsAsFactors = FALSE)
  si <- data.frame(Treatment = rownames(G), SeedAvailable = rep(300, J))
  res <- suggest_common_treatments(G, ti, si, seed_required_per_plot = c(10, 10),
                                   n_common = 12, reps_per_site = 1, seed = 1)
  expect_equal(res$n_common, 12L)
  expect_true(isTRUE(res$rationale$user_supplied))
})

test_that("errors when no line has enough seed for all sites", {
  G <- diag(5); dimnames(G) <- list(paste0("L", 1:5), paste0("L", 1:5))
  ti <- data.frame(Treatment = rownames(G), Family = "F1", stringsAsFactors = FALSE)
  si <- data.frame(Treatment = rownames(G), SeedAvailable = rep(5, 5))
  expect_error(
    suggest_common_treatments(G, ti, si, seed_required_per_plot = c(50, 50)),
    "enough seed")
})
