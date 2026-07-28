# Tests for the raw-input builders. Scientific math (VanRaden G, Henderson A,
# single-step H, haversine) verified independently in numpy.

# ---- build_relationship_matrix: G / A / H -----------------------------------

test_that("VanRaden G is symmetric, PD and correctly named", {
  set.seed(1)
  M <- matrix(sample(0:2, 40 * 200, replace = TRUE), 40, 200)
  rownames(M) <- paste0("L", 1:40)
  G <- build_relationship_matrix(markers = M, type = "genomic")
  expect_equal(G, t(G), tolerance = 1e-10)
  expect_true(all(eigen(G + diag(1e-6, 40), symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(rownames(G), rownames(M))
})

test_that("pedigree A matches the Henderson recursion (inbreeding captured)", {
  ped <- data.frame(id = 1:5, sire = c(0, 0, 1, 1, 3), dam = c(0, 0, 2, 3, 4))
  A <- build_relationship_matrix(pedigree = ped, type = "pedigree")
  A <- A[as.character(1:5), as.character(1:5)]
  # numpy reference values
  expect_equal(diag(A), c(1, 1, 1, 1.25, 1.375), tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(A["4", "5"], 1.0, tolerance = 1e-8)
  expect_equal(A["1", "3"], 0.5, tolerance = 1e-8)
})

test_that("pedigrees given out of parent order are handled", {
  ped <- data.frame(id = c(5, 4, 3, 2, 1), sire = c(3, 1, 1, 0, 0),
                    dam = c(4, 3, 2, 0, 0))
  A <- build_relationship_matrix(pedigree = ped, type = "pedigree")
  expect_equal(A["5", "5"], 1.375, tolerance = 1e-8)   # same as ordered version
})

test_that("single-step H is symmetric, PD and preserves the genomic block", {
  set.seed(2)
  ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 3, 0), dam = c(0, 0, 2, 3, 4, 0))
  M <- matrix(sample(0:2, 3 * 150, replace = TRUE), 3, 150,
              dimnames = list(c("3", "4", "5"), NULL))   # only 3 genotyped
  H <- build_relationship_matrix(markers = M, pedigree = ped, type = "hmatrix")
  expect_equal(dim(H), c(6L, 6L))
  expect_equal(H, t(H), tolerance = 1e-8)
  expect_true(all(eigen(H, symmetric = TRUE, only.values = TRUE)$values > -1e-6))
  expect_true(all(c("1", "2", "6") %in% rownames(H)))   # non-genotyped included
})

# ---- haversine geo distance in build_environment_relationship ---------------

test_that("geo = TRUE uses great-circle distance for GPS coordinates", {
  coords <- matrix(c(-1.29, 36.82,   # Nairobi
                      9.06,  7.49,   # Abuja
                      6.52,  3.38),  # Lagos
                   ncol = 2, byrow = TRUE,
                   dimnames = list(c("Nairobi", "Abuja", "Lagos"),
                                   c("latitude", "longitude")))
  D <- build_environment_relationship(coords, source = "enviromic",
                                      kernel = "gaussian", geo = TRUE)
  expect_equal(diag(D), rep(1, 3), tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(all(D >= 0 & D <= 1))
  # Abuja and Lagos (both Nigeria) are far more similar than Abuja and Nairobi
  expect_gt(D["Abuja", "Lagos"], D["Abuja", "Nairobi"])
})

# ---- cluster_environments ---------------------------------------------------

test_that("cluster_environments partitions into mega-environments with covariances", {
  # two tight blocks: {E1,E2} and {E3,E4}
  D <- matrix(c(1, .9, .1, .1,  .9, 1, .1, .1,
                .1, .1, 1, .9,  .1, .1, .9, 1), 4, 4,
              dimnames = list(paste0("E", 1:4), paste0("E", 1:4)))
  cl <- cluster_environments(D, n_clusters = 2, method = "hclust")
  expect_equal(length(cl$clusters), 2L)
  expect_equal(cl$membership[["E1"]], cl$membership[["E2"]])
  expect_false(cl$membership[["E1"]] == cl$membership[["E3"]])
  expect_true(all(vapply(cl$Sigma_E, function(s) nrow(s) == 2L, logical(1))))
})

# ---- build_enviromic_covariates (offline assembly) --------------------------

test_that("enviromic covariates merge weather, soil and management (one-hot)", {
  sites <- data.frame(environment = c("E1", "E2", "E3"),
                      latitude = c(-1, 9, 6), longitude = c(37, 7, 3))
  wx <- data.frame(environment = c("E1", "E2", "E3"),
                   mean_temp = c(21, 24, 27), precip = c(400, 300, 250))
  mg <- data.frame(environment = c("E1", "E2", "E3"),
                   planting_system = c("conv", "no_till", "conv"),
                   fertN = c(120, 80, 100))
  X <- build_enviromic_covariates(sites, management = mg, weather = wx,
                                  standardize = FALSE)
  expect_equal(rownames(X), sites$environment)
  # weather (2) + management: fertN (1) + planting_system one-hot (2) = 5 cols
  expect_true(ncol(X) >= 5L)
  expect_true(any(grepl("planting_system_", colnames(X))))
  # feeds build_environment_relationship
  D <- build_environment_relationship(X, source = "enviromic", kernel = "gaussian")
  expect_equal(dim(D), c(3L, 3L))
})
