make_block_environment_kernel <- function(within = 0.9, between = 0.05) {
  env <- paste0("E", 1:8)
  group <- rep(1:2, each = 4)
  K <- outer(group, group, function(i, j) ifelse(i == j, within, between))
  diag(K) <- 1
  dimnames(K) <- list(env, env)
  K
}

test_that("mega-environment count is inferred from stable structure", {
  D <- make_block_environment_kernel()
  D2 <- D
  D2[upper.tri(D2)] <- D2[upper.tri(D2)] +
    rep(c(-0.015, 0.015), length.out = sum(upper.tri(D2)))
  D2[lower.tri(D2)] <- t(D2)[lower.tri(D2)]
  diag(D2) <- 1

  out <- infer_mega_environments(
    D,
    relationships = list(weather_year_1 = D, weather_year_2 = D2),
    n_boot = 20L, seed = 9L
  )

  expect_equal(out$n_clusters, 2L)
  expect_equal(out$status, "stable")
  expect_true(out$hard_groups)
  expect_equal(length(unique(out$membership[1:4])), 1L)
  expect_equal(length(unique(out$membership[5:8])), 1L)
  expect_false(out$membership[1] == out$membership[5])
  expect_true(all(c("silhouette", "stability_lower",
                    "separation_pass", "stability_pass") %in%
                  names(out$diagnostics)))
})

test_that("one relationship supports only provisional grouping", {
  D <- make_block_environment_kernel()
  out <- infer_mega_environments(
    D, relationships = list(central = D), n_boot = 0L, seed = 2L
  )
  expect_equal(out$n_clusters, 2L)
  expect_equal(out$status, "provisional")
  expect_false(out$hard_groups)
})

test_that("repeated years remain one evidence block", {
  D <- make_block_environment_kernel()
  rel <- list(weather_2022 = D, weather_2023 = D, weather_2024 = D)
  groups <- stats::setNames(rep("weather", length(rel)), names(rel))
  out <- infer_mega_environments(
    D, relationships = rel, relationship_groups = groups,
    n_boot = 5L, seed = 3L
  )
  expect_equal(out$n_clusters, 2L)
  expect_equal(out$status, "provisional")
  expect_equal(out$resampling$n_relationship_blocks, 1L)
  expect_equal(out$resampling$n_boot_used, 5L)
})

test_that("unsupported environmental grouping falls back safely", {
  D <- diag(8)
  dimnames(D) <- list(paste0("E", 1:8), paste0("E", 1:8))
  out <- NULL
  expect_no_warning(
    out <- infer_mega_environments(
      D, relationships = list(year_1 = D, year_2 = D),
      n_boot = 10L, seed = 4L
    )
  )
  expect_equal(out$n_clusters, 1L)
  expect_equal(out$status, "unstable")
  expect_false(out$hard_groups)
  expect_equal(unique(out$membership), 1L)
  expect_match(out$reason, "No partition passed")
})

test_that("fixed-k clustering remains available for external evidence", {
  D <- make_block_environment_kernel()
  out <- cluster_environments(D, n_clusters = 2L, method = "kmeans",
                              seed = 7L)
  expect_equal(length(out$clusters), 2L)
  expect_equal(unname(sort(vapply(out$clusters, length, integer(1)))),
               c(4L, 4L))
  expect_equal(names(out$membership), rownames(D))
})

test_that("finite-descent k-means converges with non-empty groups", {
  D <- diag(8)
  dimnames(D) <- list(paste0("E", 1:8), paste0("E", 1:8))
  cl <- NULL
  expect_no_warning(
    cl <- .finite_descent_kmeans(
      .mega_coordinates(D), k = 3L, nstart = 10L, seed = 19L
    )
  )
  expect_length(cl, 8L)
  expect_equal(sort(unique(cl)), 1:3)
  expect_true(all(tabulate(cl, nbins = 3L) > 0L))
})

test_that("mega-environment inference preserves the caller RNG state", {
  D <- make_block_environment_kernel()
  set.seed(77)
  before <- .Random.seed
  infer_mega_environments(
    D, relationships = list(year_1 = D, year_2 = D),
    n_boot = 5L, seed = 12L
  )
  expect_equal(.Random.seed, before)
})
