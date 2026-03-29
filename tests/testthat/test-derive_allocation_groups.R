# ==============================================================================
# test-derive_allocation_groups.R
# ==============================================================================

make_trt <- function(n = 12) paste0("L", sprintf("%03d", seq_len(n)))

make_tinfo <- function(trt) {
  data.frame(
    Treatment = trt,
    Family    = rep(c("F1", "F2", "F3"), length.out = length(trt)),
    stringsAsFactors = FALSE
  )
}

# -- None mode ----------------------------------------------------------------

test_that("'none' assigns all treatments to 'ALL'", {
  trt <- make_trt()
  res <- derive_allocation_groups(trt, allocation_group_source = "none")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(trt))
  expect_true(all(res$AllocationGroup == "ALL"))
  expect_equal(names(res), c("Treatment", "AllocationGroup"))
})

# -- Family mode --------------------------------------------------------------

test_that("'Family' reads labels from treatment_info", {
  trt   <- make_trt()
  tinfo <- make_tinfo(trt)
  res   <- derive_allocation_groups(trt, "Family", treatment_info = tinfo)
  expect_equal(nrow(res), length(trt))
  expect_equal(sort(unique(res$AllocationGroup)), c("F1", "F2", "F3"))
  # Order matches treatment order
  expect_equal(res$Treatment, trt)
})

test_that("'Family' stops when treatment missing from treatment_info", {
  trt   <- make_trt(12)
  tinfo <- make_tinfo(trt[1:10])   # missing last two
  expect_error(
    derive_allocation_groups(trt, "Family", treatment_info = tinfo),
    regexp = "no family"
  )
})

test_that("'Family' stops when treatment_info not supplied", {
  expect_error(
    derive_allocation_groups(make_trt(), "Family"),
    regexp = "treatment_info"
  )
})

# -- GRM mode -----------------------------------------------------------------

make_grm <- function(trt, seed = 1L) {
  n   <- length(trt)
  K   <- make_sym_pd(n, seed = seed)
  rownames(K) <- colnames(K) <- trt
  K
}

test_that("'GRM' returns one row per treatment with GRP_G prefix", {
  trt   <- make_trt()
  tinfo <- make_tinfo(trt)
  GRM   <- make_grm(trt)
  res   <- derive_allocation_groups(
    trt, "GRM", GRM = GRM, treatment_info = tinfo,
    group_seed = 1L, group_attempts = 5L
  )
  expect_equal(nrow(res), length(trt))
  expect_true(all(grepl("^GRP_G", res$AllocationGroup)))
  expect_equal(res$Treatment, trt)
})

test_that("'GRM' anchors cluster count to unique families when supplied", {
  trt   <- make_trt()
  tinfo <- make_tinfo(trt)  # 3 families
  GRM   <- make_grm(trt)
  res   <- derive_allocation_groups(
    trt, "GRM", GRM = GRM, treatment_info = tinfo,
    group_seed = 1L, group_attempts = 5L
  )
  n_groups <- length(unique(res$AllocationGroup))
  expect_equal(n_groups, 3L)
})

test_that("'A' prefix used for pedigree matrix", {
  trt <- make_trt()
  A   <- make_grm(trt, seed = 7L)
  res <- derive_allocation_groups(
    trt, "A", A = A, group_seed = 1L, group_attempts = 5L
  )
  expect_true(all(grepl("^GRP_A", res$AllocationGroup)))
})

test_that("deduplication of treatments before grouping", {
  trt   <- c(make_trt(6), make_trt(3))  # 9 with 3 duplicates
  tinfo <- make_tinfo(unique(trt))
  res   <- derive_allocation_groups(unique(trt), "Family", treatment_info = tinfo)
  expect_equal(nrow(res), 6L)
})

test_that("stops when GRM required but NULL", {
  expect_error(
    derive_allocation_groups(make_trt(), "GRM"),
    regexp = "Relationship matrix"
  )
})

test_that("stops when GRM lacks rownames", {
  trt <- make_trt()
  K   <- make_sym_pd(length(trt))
  expect_error(
    derive_allocation_groups(trt, "GRM", GRM = K),
    regexp = "row names"
  )
})
