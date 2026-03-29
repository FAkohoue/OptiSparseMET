# ==============================================================================
# test-combine_met_fieldbooks.R
# ==============================================================================

make_fb <- function(env_id, treatments, extra_col = NULL) {
  fb <- data.frame(
    Treatment = treatments,
    Family    = rep("F1", length(treatments)),
    Block     = 1L,
    Plot      = seq_along(treatments),
    Row       = seq_along(treatments),
    Column    = 1L,
    stringsAsFactors = FALSE
  )
  if (!is.null(extra_col)) fb[[extra_col]] <- rnorm(nrow(fb))
  fb
}

# -- Return structure ---------------------------------------------------------

test_that("returns a data frame with correct row count", {
  fb_list <- list(
    E1 = make_fb("E1", c("T1", "T2", "CHK")),
    E2 = make_fb("E2", c("T3", "T4", "CHK"))
  )
  res <- combine_met_fieldbooks(fb_list)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 6L)
})

test_that("Environment column populated from list names", {
  fb_list <- list(
    SiteA = make_fb("A", c("T1", "T2")),
    SiteB = make_fb("B", c("T3", "T4"))
  )
  res <- combine_met_fieldbooks(fb_list)
  expect_equal(sort(unique(res$Environment)), c("SiteA", "SiteB"))
})

test_that("standard columns always present in output", {
  fb_list <- list(E1 = make_fb("E1", "T1"))
  res     <- combine_met_fieldbooks(fb_list)
  required <- c("Environment", "LocalDesign", "ReplicationMode", "SparseMethod",
                "IsCommonTreatment", "Treatment", "Family", "Gcluster",
                "Block", "Plot", "Row", "Column")
  expect_true(all(required %in% names(res)))
})

test_that("standard columns appear at the front", {
  fb_list <- list(E1 = make_fb("E1", "T1"))
  res     <- combine_met_fieldbooks(fb_list)
  front   <- names(res)[1:12]
  expect_equal(front[1], "Environment")
  expect_equal(front[6], "Treatment")
})

# -- Metadata columns ---------------------------------------------------------

test_that("LocalDesign populated from local_designs argument", {
  fb_list <- list(
    E1 = make_fb("E1", "T1"),
    E2 = make_fb("E2", "T2")
  )
  res <- combine_met_fieldbooks(
    fb_list,
    local_designs = c(E1 = "met_prep_famoptg", E2 = "met_alpha_rc_stream")
  )
  expect_equal(res$LocalDesign[res$Environment == "E1"], "met_prep_famoptg")
  expect_equal(res$LocalDesign[res$Environment == "E2"], "met_alpha_rc_stream")
})

test_that("SparseMethod applied uniformly to all rows", {
  fb_list <- list(E1 = make_fb("E1", c("T1","T2")))
  res <- combine_met_fieldbooks(fb_list, sparse_method = "random_balanced")
  expect_true(all(res$SparseMethod == "random_balanced"))
})

test_that("IsCommonTreatment flags correctly", {
  fb_list <- list(
    E1 = make_fb("E1", c("T1", "CHK1")),
    E2 = make_fb("E2", c("T2", "CHK1"))
  )
  res <- combine_met_fieldbooks(fb_list, common_treatments = "CHK1")
  expect_equal(sum(res$IsCommonTreatment), 2L)
  expect_true(all(res$IsCommonTreatment[res$Treatment == "CHK1"]))
})

test_that("NULL common_treatments -> IsCommonTreatment all FALSE", {
  fb_list <- list(E1 = make_fb("E1", c("T1", "T2")))
  res <- combine_met_fieldbooks(fb_list)
  expect_true(all(!res$IsCommonTreatment))
})

# -- Heterogeneous columns ----------------------------------------------------

test_that("missing columns filled with NA for environments lacking them", {
  fb_list <- list(
    E1 = make_fb("E1", c("T1", "T2")),
    E2 = make_fb("E2", c("T3", "T4"), extra_col = "SpatialResidual")
  )
  res <- combine_met_fieldbooks(fb_list)
  expect_true("SpatialResidual" %in% names(res))
  expect_true(all(is.na(res$SpatialResidual[res$Environment == "E1"])))
  expect_false(any(is.na(res$SpatialResidual[res$Environment == "E2"])))
})

test_that("rownames are clean sequential integers in output", {
  # combine_met_fieldbooks() calls rownames(out) <- NULL which resets row names
  # to the default sequential integer sequence ("1", "2", ..., n).
  # This ensures rbind() does not leave environment-prefixed names like "E1.1".
  fb_list <- list(E1 = make_fb("E1", c("T1","T2","T3")))
  res     <- combine_met_fieldbooks(fb_list)
  expect_equal(rownames(res), as.character(seq_len(nrow(res))))
})

# -- Validation ---------------------------------------------------------------

test_that("stops when field_books is not a named list", {
  expect_error(combine_met_fieldbooks(list(make_fb("E1", "T1"))))
})

test_that("stops when an element is not a data frame", {
  expect_error(
    combine_met_fieldbooks(list(E1 = "not a data frame")),
    regexp = "data frame"
  )
})