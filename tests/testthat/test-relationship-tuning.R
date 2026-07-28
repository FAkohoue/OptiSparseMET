# Tests for AGHmatrix / ASRgenomics integration: G build method, G.tuneup
# conditioning, and kinship PCA. Assertions hold whether or not the optional
# packages are installed (they exercise the built-in fallbacks otherwise).

mk_markers <- function(n = 20L, m = 300L, seed = 1) {
  set.seed(seed)
  M <- matrix(sample(0:2, n * m, replace = TRUE), n, m)
  rownames(M) <- paste0("L", seq_len(n)); M
}

test_that("method = 'AGHmatrix' builds a symmetric G (or falls back cleanly)", {
  M <- mk_markers()
  G <- build_relationship_matrix(markers = M, type = "genomic",
                                 method = "AGHmatrix", ploidy = 2)
  expect_equal(dim(G), c(nrow(M), nrow(M)))
  expect_equal(G, t(G), tolerance = 1e-8)
  expect_equal(rownames(G), rownames(M))
})

test_that("ploidy feeds the built-in VanRaden estimator", {
  M <- mk_markers()
  g2 <- build_relationship_matrix(markers = M, type = "genomic", ploidy = 2)
  expect_equal(dim(g2), c(nrow(M), nrow(M)))
  expect_true(all(is.finite(g2)))
})

test_that("dominance relationship is symmetric and distinct from additive", {
  set.seed(4)
  M <- matrix(sample(0:2, 25 * 300, TRUE), 25, 300); rownames(M) <- paste0("H", 1:25)
  G_A <- build_relationship_matrix(markers = M, type = "genomic",
                                   relationship = "additive")
  G_D <- build_relationship_matrix(markers = M, type = "genomic",
                                   relationship = "dominance")
  expect_equal(dim(G_D), dim(G_A))
  expect_equal(G_D, t(G_D), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(G_A, G_D)))          # captures different info
  expect_error(build_relationship_matrix(pedigree = data.frame(id = 1, sire = 0,
                 dam = 0), type = "pedigree", relationship = "dominance"),
               "only available for type")
})

test_that("additive is the default and dominance remains explicit", {
  M <- mk_markers()
  G_default <- build_relationship_matrix(markers = M, type = "genomic")
  G_additive <- build_relationship_matrix(
    markers = M, type = "genomic", relationship = "additive")
  expect_equal(G_default, G_additive, tolerance = 1e-12)
})

test_that("combine_relationship_matrices makes a variance-weighted blend", {
  A <- matrix(c(2, 0, 0, 2), 2, dimnames = list(c("H1","H2"), c("H1","H2")))
  D <- matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(c("H1","H2"), c("H1","H2")))
  G <- combine_relationship_matrices(list(additive = A, dominance = D),
                                     weights = c(3, 1))
  expect_equal(G["H1","H1"], 0.75 * 2 + 0.25 * 1)    # weights normalised to 3:1
  expect_equal(G["H1","H2"], 0.75 * 0 + 0.25 * 0.5)
  G_named <- combine_relationship_matrices(
    list(additive = A, dominance = D),
    weights = c(dominance = 1, additive = 3))
  expect_equal(G_named, G)
  expect_error(combine_relationship_matrices(
    list(additive = A, dominance = D),
    weights = c(additive = 3, other = 1)), "matching")
  expect_error(combine_relationship_matrices(list()), "empty")
})

test_that("native 'bend' makes a matrix positive definite without ASRgenomics", {
  # a rank-deficient (singular) matrix
  set.seed(5); B <- matrix(rnorm(4 * 2), 4, 2); S <- tcrossprod(B)
  dimnames(S) <- list(paste0("L", 1:4), paste0("L", 1:4))
  Sb <- tune_relationship_matrix(S, method = "bend")
  ev <- eigen(Sb, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))                            # PD
  expect_false(inherits(try(solve(Sb), silent = TRUE), "try-error"))
})

test_that("tune_relationship_matrix returns an invertible, same-shape matrix", {
  M <- mk_markers()
  G <- build_relationship_matrix(markers = M, type = "genomic")
  Gb <- tune_relationship_matrix(G, blend = TRUE, rcn = TRUE)
  expect_equal(dim(Gb), dim(G))
  expect_equal(rownames(Gb), rownames(G))
  expect_false(inherits(try(solve(Gb), silent = TRUE), "try-error"))  # invertible
})

test_that("tuneup = TRUE inside build_relationship_matrix yields an invertible G", {
  M <- mk_markers()
  G <- build_relationship_matrix(markers = M, type = "genomic", tuneup = TRUE)
  expect_false(inherits(try(solve(G), silent = TRUE), "try-error"))
})

test_that("build_hybrid_relationship equals VanRaden G on parent-mean dosages", {
  set.seed(3)
  Mpar <- matrix(sample(0:2, 8 * 400, TRUE), 8, 400,
                 dimnames = list(c("A","B","C","D","E","F","T1","T2"), NULL))
  Gp <- build_relationship_matrix(markers = Mpar, type = "genomic")
  crosses <- data.frame(p1 = c("A","B","C","D","E","F"),
                        p2 = c("T1","T1","T2","T2","T1","T2"))
  ids <- paste0("H", 1:6)
  Ghyb <- build_hybrid_relationship(Gp, crosses, hybrid_ids = ids)
  # direct: VanRaden on the mean-dosage hybrid markers, centred on the SAME
  # (parental) allele frequencies -- the correct centring for hybrids.
  p <- colMeans(Mpar) / 2
  Zh <- sweep(0.5 * (Mpar[crosses$p1, ] + Mpar[crosses$p2, ]), 2, 2 * p)
  Gdirect <- tcrossprod(Zh) / (2 * sum(p * (1 - p)))
  expect_equal(unname(Ghyb), unname(Gdirect), tolerance = 1e-8)
  expect_equal(rownames(Ghyb), ids)
  expect_error(build_hybrid_relationship(Gp, data.frame(p1 = "A", p2 = "ZZ")),
               "not found")
})

test_that("hybrid relationship creates stable automatic IDs", {
  Mpar <- mk_markers(n = 6L, m = 100L)
  rownames(Mpar) <- c("A", "B", "C", "T1", "T2", "T3")
  Gp <- build_relationship_matrix(Mpar, type = "genomic")
  crosses <- data.frame(parent1 = c("A", "B", "C"),
                        parent2 = c("T1", "T2", "T3"))
  Gh <- build_hybrid_relationship(Gp, crosses)
  expect_equal(rownames(Gh), c("A_x_T1", "B_x_T2", "C_x_T3"))
  expect_error(build_hybrid_relationship(
    Gp, crosses, hybrid_ids = c("H1", "H1", "H3")), "unique")
})

test_that("kinship_pca returns PCs (ASRgenomics object or eigen fallback)", {
  M <- mk_markers()
  G <- tune_relationship_matrix(build_relationship_matrix(markers = M, type = "genomic"))
  pcs <- kinship_pca(G, ncp = 5)
  expect_true(is.list(pcs))
  ## the built-in fallback exposes scores with the requested number of PCs
  if (!requireNamespace("ASRgenomics", quietly = TRUE)) {
    expect_equal(ncol(pcs$scores), 5L)
    expect_equal(nrow(pcs$scores), nrow(G))
  }
})
