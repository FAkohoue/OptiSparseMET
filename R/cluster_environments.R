#' Partition environments into a specified number of mega-environments
#'
#' `cluster_environments()` is the low-level function for situations where the
#' number of groups is fixed by external evidence. In routine use,
#' [infer_mega_environments()] should be preferred because it estimates the
#' number of supported groups and refuses to impose unstable hard groupings.
#'
#' @param D Square, finite, named environment relationship/similarity matrix.
#' @param n_clusters Number of mega-environments fixed by external evidence.
#' @param method `"hclust"` (default) or `"kmeans"`. Hierarchical clustering
#'   operates on kernel-induced distances. K-means operates on positive-eigen
#'   kernel principal coordinates, not directly on raw similarities, and uses
#'   a finite-descent relocation algorithm that always terminates.
#' @param Sigma_E Optional environment covariance to subset per cluster
#'   (defaults to `D`).
#' @param seed Optional integer used by k-means.
#'
#' @return A list with `membership` (a named integer cluster vector), `clusters`
#'   (a list of environment-name vectors), and `Sigma_E` (a list of
#'   within-cluster covariance matrices).
#'
#' @seealso [infer_mega_environments()],
#'   [build_environment_relationship()], [simulate_met()].
#' @export
cluster_environments <- function(
    D, n_clusters, method = c("hclust", "kmeans"),
    Sigma_E = NULL, seed = NULL) {
  method <- match.arg(method)
  D <- .validate_mega_similarity(D, "D")
  E <- nrow(D)
  if (!is.numeric(n_clusters) || length(n_clusters) != 1L ||
      !is.finite(n_clusters) ||
      abs(n_clusters - round(n_clusters)) > 1e-8 ||
      n_clusters < 1L || n_clusters > E)
    stop("`n_clusters` must be one integer between 1 and nrow(D).")
  n_clusters <- as.integer(round(n_clusters))

  if (is.null(Sigma_E)) Sigma_E <- D
  Sigma_E <- .validate_mega_similarity(Sigma_E, "Sigma_E")
  envs <- rownames(D)
  if (!setequal(rownames(Sigma_E), envs))
    stop("`Sigma_E` must contain exactly the environments in `D`.")
  Sigma_E <- Sigma_E[envs, envs, drop = FALSE]

  cl <- .cluster_mega_similarity(D, n_clusters, method, seed)
  if (is.null(cl))
    stop("The requested clustering is not estimable from `D`; reduce ",
         "`n_clusters` or revise the environment relationship matrix.")
  .assemble_mega_result(cl, Sigma_E)
}


#' Infer stable mega-environments without specifying their number
#'
#' `infer_mega_environments()` chooses both the number of mega-environments and
#' the clustering algorithm from the environmental evidence. It compares
#' hierarchical clustering with k-means across a defensible candidate range,
#' requires adequate average silhouette separation, rejects clusters smaller
#' than `min_cluster_size`, and evaluates reproducibility across supplied
#' modality-, year-, or bootstrap-specific relationship matrices. Bootstrap
#' draws resample whole relationship matrices and recompute their rank-median
#' consensus; they do not estimate or apply arbitrary modality weights.
#'
#' Candidate solutions first have to pass the separation and stability
#' thresholds. Selection is then lexicographic: largest silhouette, largest
#' lower stability quantile, fewer groups, and finally method order. Thus no
#' user-defined weighted score determines the answer. With only one
#' relationship matrix a well-separated, algorithm-reproducible solution is
#' labelled `"provisional"`. With multiple independent evidence blocks it can
#' be labelled `"stable"`. If no candidate is supported, the function returns one
#' unpartitioned group with status `"unstable"` rather than forcing a
#' biologically fragile classification.
#'
#' K-means uses a finite-descent Hartigan relocation algorithm. It starts from
#' non-empty partitions and accepts only moves that strictly reduce
#' within-cluster sum of squares. Because the number of partitions is finite,
#' every start terminates at a valid local optimum without an iteration-limit
#' failure. Multiple starts and deterministic seed handling are retained.
#'
#' These groups are descriptive environmental strata. They become genetic
#' mega-environments only when historical MET response data support the
#' corresponding genotype-by-environment structure.
#'
#' @param D Central square, finite, named environment similarity matrix.
#' @param relationships Optional named list of aligned environment relationship
#'   matrices representing distinct modalities, years, or defensible
#'   uncertainty draws. These are used to assess grouping stability, not to
#'   estimate genetic covariance or user weights.
#' @param relationship_groups Optional named character vector assigning each
#'   relationship matrix to an evidence block, such as all yearly weather
#'   matrices to `"weather"`. Agreement is first aggregated within blocks, and
#'   bootstrap sampling draws within every block, so a modality with many years
#'   cannot dominate modalities represented by one matrix.
#' @param Sigma_E Optional covariance matrix returned within inferred groups.
#'   It does not affect group selection and defaults to `D`.
#' @param k_range Optional integer candidate counts. The default is
#'   `2:min(8, floor(n_environment / min_cluster_size))`.
#' @param methods Any non-empty subset of `"hclust"` and `"kmeans"`.
#' @param min_cluster_size Smallest admissible group size. The default of two
#'   prevents an unsupported singleton from being called a mega-environment.
#' @param min_silhouette Minimum average silhouette required for a candidate.
#' @param min_stability Minimum lower stability quantile required. Stability is
#'   measured by the adjusted Rand index against alternative methods,
#'   relationship matrices, and relationship-bootstrap consensuses.
#' @param stability_quantile Lower quantile used for the stability requirement.
#' @param n_boot Number of relationship-level bootstrap draws. Bootstrap is
#'   performed when at least one evidence block contains repeated matrices.
#' @param seed Integer seed for reproducible k-means and bootstrap sampling.
#'
#' @return A list containing `membership`, `clusters`, `Sigma_E`, inferred
#'   `n_clusters`, selected `method`, `status` (`"stable"`, `"provisional"`, or
#'   `"unstable"`), `hard_groups`, a human-readable `reason`, candidate
#'   `diagnostics`, and resampling metadata. When status is `"unstable"`,
#'   `membership` contains one group and `hard_groups` is `FALSE`.
#'
#' @seealso [cluster_environments()], [consensus_environment_kernels()],
#'   [assess_envirotype_stability()].
#' @export
infer_mega_environments <- function(
    D, relationships = NULL, relationship_groups = NULL,
    Sigma_E = NULL, k_range = NULL,
    methods = c("hclust", "kmeans"), min_cluster_size = 2L,
    min_silhouette = 0.25, min_stability = 0.70,
    stability_quantile = 0.10, n_boot = 200L, seed = 1L) {
  D <- .validate_mega_similarity(D, "D")
  envs <- rownames(D)
  E <- length(envs)

  if (is.null(Sigma_E)) Sigma_E <- D
  Sigma_E <- .validate_mega_similarity(Sigma_E, "Sigma_E")
  if (!setequal(rownames(Sigma_E), envs))
    stop("`Sigma_E` must contain exactly the environments in `D`.")
  Sigma_E <- Sigma_E[envs, envs, drop = FALSE]

  methods <- unique(match.arg(methods, c("hclust", "kmeans"),
                              several.ok = TRUE))
  if (!length(methods))
    stop("`methods` must contain at least one supported clustering method.")
  .check_mega_count(min_cluster_size, "min_cluster_size", minimum = 1L)
  min_cluster_size <- as.integer(round(min_cluster_size))
  .check_mega_probability(min_silhouette, "min_silhouette", -1, 1)
  .check_mega_probability(min_stability, "min_stability", 0, 1)
  .check_mega_probability(stability_quantile, "stability_quantile", 0, 1)
  .check_mega_count(n_boot, "n_boot", minimum = 0L)
  n_boot <- as.integer(round(n_boot))
  if (!is.null(seed)) {
    .check_mega_count(seed, "seed", minimum = 0L)
    seed <- as.integer(round(seed))
  }

  max_k <- min(8L, floor(E / min_cluster_size))
  if (is.null(k_range)) {
    k_range <- if (max_k >= 2L) seq.int(2L, max_k) else integer(0)
  } else {
    if (!is.numeric(k_range) || any(!is.finite(k_range)) ||
        any(abs(k_range - round(k_range)) > 1e-8))
      stop("`k_range` must contain finite integers.")
    k_range <- sort(unique(as.integer(round(k_range))))
    k_range <- k_range[k_range >= 2L & k_range <= E]
    if (any(k_range * min_cluster_size > E))
      stop("Every `k_range` value must permit `min_cluster_size` ",
           "environments per group.")
  }

  rel <- .validate_mega_relationships(relationships, envs)
  n_relationships <- length(rel)
  if (is.null(relationship_groups)) {
    relationship_groups <- stats::setNames(names(rel), names(rel))
  } else {
    if (!is.character(relationship_groups) ||
        is.null(names(relationship_groups)) ||
        anyDuplicated(names(relationship_groups)) ||
        !setequal(names(relationship_groups), names(rel)) ||
        anyNA(relationship_groups) || any(!nzchar(relationship_groups)))
      stop("`relationship_groups` must be a named, non-missing character ",
           "vector assigning every `relationships` matrix to one block.")
    relationship_groups <- relationship_groups[names(rel)]
  }
  n_relationship_blocks <- length(unique(relationship_groups))

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (old_seed_exists)
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)

  fallback <- function(reason, diagnostics = .empty_mega_diagnostics()) {
    cl <- stats::setNames(rep.int(1L, E), envs)
    out <- .assemble_mega_result(cl, Sigma_E)
    c(out, list(
      n_clusters = 1L, method = NA_character_, status = "unstable",
      hard_groups = FALSE, reason = reason, diagnostics = diagnostics,
      resampling = list(
        n_relationships = n_relationships,
        n_relationship_blocks = n_relationship_blocks,
        n_boot_requested = n_boot,
        n_boot_used = 0L
      )
    ))
  }

  if (!length(k_range))
    return(fallback(
      "Too few environments support two groups at the requested minimum group size."
    ))

  # Precompute alternative partitions for each k. Relationship matrices are
  # kept separate; bootstrap consensuses resample whole evidence blocks.
  alternatives <- stats::setNames(vector("list", length(k_range)),
                                  as.character(k_range))
  for (k in k_range) {
    alt_k <- list()
    for (m in methods) {
      z <- .cluster_mega_similarity(D, k, m, seed)
      if (!is.null(z))
        alt_k[[paste0("central:", m)]] <- z
    }
    if (n_relationships) {
      for (r in seq_along(rel)) {
        for (m in methods) {
          z <- .cluster_mega_similarity(
            rel[[r]], k, m,
            if (is.null(seed)) NULL else seed + 1000L + 31L * r + k
          )
          if (!is.null(z))
            alt_k[[paste0("relationship:", names(rel)[r], ":", m)]] <- z
        }
      }
    }
    alternatives[[as.character(k)]] <- alt_k
  }

  boot_used <- 0L
  replicated_blocks <- split(seq_along(rel), relationship_groups)
  can_bootstrap <- any(lengths(replicated_blocks) > 1L)
  if (can_bootstrap && n_boot > 0L) {
    for (b in seq_len(n_boot)) {
      # One draw per evidence block prevents temporal replication within one
      # modality from becoming an implicit modality weight.
      take <- vapply(replicated_blocks, function(ii)
        sample(ii, 1L), integer(1))
      boot_rel <- rel[take]
      names(boot_rel) <- paste0("draw", seq_along(boot_rel))
      Db <- consensus_environment_kernels(boot_rel)
      for (k in k_range) {
        for (m in methods) {
          z <- .cluster_mega_similarity(
            Db, k, m,
            if (is.null(seed)) NULL else seed + 100000L + 101L * b + k
          )
          if (!is.null(z)) {
            alternatives[[as.character(k)]][[
              paste0("bootstrap:", b, ":", m)
            ]] <- z
          }
        }
      }
      boot_used <- boot_used + 1L
    }
  }

  dst <- .mega_distance(D)
  rows <- list()
  row_i <- 0L
  for (k in k_range) {
    alt_k <- alternatives[[as.character(k)]]
    central_names <- paste0("central:", methods)
    for (nm in central_names) {
      if (!nm %in% names(alt_k)) {
        row_i <- row_i + 1L
        rows[[row_i]] <- data.frame(
          k = k, method = sub("^central:", "", nm), estimable = FALSE,
          min_cluster_size = NA_integer_, max_cluster_size = NA_integer_,
          silhouette = NA_real_, stability_mean = NA_real_,
          stability_lower = NA_real_, n_stability_comparisons = 0L,
          n_stability_blocks = 0L, separation_pass = FALSE,
          size_pass = FALSE, stability_pass = FALSE,
          stringsAsFactors = FALSE
        )
        next
      }
      candidate <- alt_k[[nm]]
      method <- sub("^central:", "", nm)
      sizes <- tabulate(candidate, nbins = k)
      silhouette <- .mean_mega_silhouette(candidate, dst)
      comparisons <- alt_k[setdiff(names(alt_k), nm)]
      stability_raw <- if (length(comparisons))
        vapply(comparisons, .adjusted_rand_index, numeric(1),
               y = candidate) else numeric(0)
      stability_raw <- pmax(0, pmin(1, stability_raw))
      comparison_names <- names(comparisons)
      comparison_blocks <- vapply(comparison_names, function(z) {
        if (startsWith(z, "central:")) return("algorithm")
        if (startsWith(z, "bootstrap:")) return("bootstrap")
        relationship_name <- sub(
          ":(hclust|kmeans)$", "",
          sub("^relationship:", "", z)
        )
        paste0("relationship:", relationship_groups[[relationship_name]])
      }, character(1))
      stability <- if (length(stability_raw)) {
        stability_split <- split(stability_raw, comparison_blocks)
        vapply(names(stability_split), function(block) {
          z <- stability_split[[block]]
          if (identical(block, "bootstrap"))
            unname(stats::quantile(z, stability_quantile,
                                   names = FALSE, type = 8)) else mean(z)
        }, numeric(1))
      } else numeric(0)
      stability_mean <- if (length(stability)) mean(stability) else NA_real_
      stability_lower <- if (length(stability))
        unname(stats::quantile(stability, stability_quantile,
                               names = FALSE, type = 8)) else NA_real_
      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        k = k, method = method, estimable = TRUE,
        min_cluster_size = min(sizes), max_cluster_size = max(sizes),
        silhouette = silhouette,
        stability_mean = stability_mean,
        stability_lower = stability_lower,
        n_stability_comparisons = length(stability_raw),
        n_stability_blocks = length(stability),
        separation_pass = is.finite(silhouette) &&
          silhouette >= min_silhouette,
        size_pass = min(sizes) >= min_cluster_size,
        stability_pass = is.finite(stability_lower) &&
          stability_lower >= min_stability,
        stringsAsFactors = FALSE
      )
    }
  }
  diagnostics <- if (length(rows)) do.call(rbind, rows) else
    .empty_mega_diagnostics()
  if (!nrow(diagnostics) || !any(diagnostics$estimable))
    return(fallback("No candidate partition was estimable.", diagnostics))

  eligible <- diagnostics$separation_pass & diagnostics$size_pass &
    diagnostics$stability_pass
  if (!any(eligible))
    return(fallback(
      paste0(
        "No partition passed minimum separation, group-size, and stability ",
        "requirements; environments are retained as one unpartitioned set."
      ),
      diagnostics
    ))

  pool <- diagnostics[eligible, , drop = FALSE]
  method_rank <- match(pool$method, methods)
  ord <- order(-pool$silhouette, -pool$stability_lower,
               pool$k, method_rank)
  selected <- pool[ord[1L], , drop = FALSE]
  cl <- alternatives[[as.character(selected$k)]][[
    paste0("central:", selected$method)
  ]]
  out <- .assemble_mega_result(cl, Sigma_E)
  status <- if (n_relationship_blocks >= 2L) "stable" else "provisional"
  reason <- if (status == "stable") {
    paste0(
      "Partition passed separation, minimum-size, and cross-relationship ",
      "stability requirements."
    )
  } else {
    paste0(
      "Partition is separated and reproducible across algorithms, but only ",
      "one or no independent evidence block was supplied."
    )
  }
  c(out, list(
    n_clusters = as.integer(selected$k),
    method = as.character(selected$method),
    status = status,
    hard_groups = identical(status, "stable"),
    reason = reason,
    diagnostics = diagnostics,
    resampling = list(
      n_relationships = n_relationships,
      n_relationship_blocks = n_relationship_blocks,
      n_boot_requested = n_boot,
      n_boot_used = boot_used
    )
  ))
}


.validate_mega_similarity <- function(x, label) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || nrow(x) != ncol(x) ||
      nrow(x) < 1L)
    stop("`", label, "` must be a non-empty square numeric matrix.")
  if (any(!is.finite(x)))
    stop("`", label, "` must contain only finite values.")
  rn <- rownames(x)
  cn <- colnames(x)
  if (is.null(rn) && is.null(cn)) {
    rn <- cn <- paste0("E", seq_len(nrow(x)))
    dimnames(x) <- list(rn, cn)
  } else if (is.null(rn) || is.null(cn) || any(!nzchar(rn)) ||
             any(!nzchar(cn)) || anyDuplicated(rn) || anyDuplicated(cn) ||
             !setequal(rn, cn)) {
    stop("`", label, "` must have unique matching row and column names.")
  } else {
    x <- x[rn, rn, drop = FALSE]
  }
  asym <- max(abs(x - t(x)))
  scale_ref <- max(1, max(abs(x)))
  if (asym > 1e-7 * scale_ref)
    stop("`", label, "` must be symmetric.")
  x <- (x + t(x)) / 2
  dimnames(x) <- list(rn, rn)
  x
}

.validate_mega_relationships <- function(relationships, envs) {
  if (is.null(relationships)) return(list())
  if (!is.list(relationships) || !length(relationships))
    stop("`relationships` must be NULL or a non-empty list of matrices.")
  if (is.null(names(relationships)) || any(!nzchar(names(relationships))) ||
      anyDuplicated(names(relationships)))
    stop("`relationships` must have unique, non-empty names.")
  out <- lapply(seq_along(relationships), function(i) {
    z <- .validate_mega_similarity(
      relationships[[i]], paste0("relationships$", names(relationships)[i])
    )
    if (!setequal(rownames(z), envs))
      stop("Every `relationships` matrix must contain exactly the ",
           "environments in `D`.")
    z[envs, envs, drop = FALSE]
  })
  names(out) <- names(relationships)
  out
}

.check_mega_count <- function(x, label, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      abs(x - round(x)) > 1e-8 || x < minimum)
    stop("`", label, "` must be one integer >= ", minimum, ".")
  invisible(TRUE)
}

.check_mega_probability <- function(x, label, lower, upper) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper)
    stop("`", label, "` must be one finite number in [", lower,
         ", ", upper, "].")
  invisible(TRUE)
}

.mega_distance <- function(D) {
  d2 <- outer(diag(D), diag(D), `+`) - 2 * D
  d2[d2 < 0 & d2 > -1e-8 * max(1, max(abs(d2)))] <- 0
  d2 <- pmax(d2, 0)
  diag(d2) <- 0
  sqrt(d2)
}

.mega_coordinates <- function(D) {
  n <- nrow(D)
  if (n == 1L) return(matrix(0, 1L, 1L))
  H <- diag(n) - matrix(1 / n, n, n)
  B <- (H %*% D %*% H + t(H %*% D %*% H)) / 2
  ee <- eigen(B, symmetric = TRUE)
  tol <- max(1, max(abs(ee$values))) * 1e-9
  keep <- which(ee$values > tol)
  if (!length(keep)) return(matrix(0, n, 1L))
  sweep(ee$vectors[, keep, drop = FALSE], 2L,
        sqrt(ee$values[keep]), `*`)
}

.finite_descent_kmeans <- function(X, k, nstart = 50L, seed = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < k || k < 1L || any(!is.finite(X))) return(NULL)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  # Hierarchical initialisation supplies one deterministic, non-empty start.
  # Remaining starts are balanced random partitions, so every group contains
  # at least one environment before finite descent begins.
  starts <- vector("list", nstart)
  starts[[1L]] <- stats::cutree(
    stats::hclust(stats::dist(X), method = "average"), k = k
  )
  if (nstart > 1L) {
    for (s in 2:nstart) {
      cl <- integer(n)
      cl[sample.int(n)] <- rep(seq_len(k), length.out = n)
      starts[[s]] <- cl
    }
  }

  objective_scale <- max(1, sum(X * X))
  improvement_tolerance <- 1e-12 * objective_scale
  summarise_partition <- function(cl) {
    sizes <- tabulate(cl, nbins = k)
    centers <- do.call(rbind, lapply(seq_len(k), function(g)
      colMeans(X[cl == g, , drop = FALSE])))
    wss <- sum(vapply(seq_len(n), function(i)
      sum((X[i, ] - centers[cl[i], ])^2), numeric(1)))
    list(sizes = sizes, centers = centers, withinss = wss)
  }
  optimise_start <- function(cl) {
    cl <- as.integer(cl)
    if (length(unique(cl)) != k) return(NULL)
    repeat {
      current <- summarise_partition(cl)
      sizes <- current$sizes
      centers <- current$centers
      best_delta <- 0
      best_i <- NA_integer_
      best_group <- NA_integer_
      for (i in seq_len(n)) {
        from <- cl[i]
        if (sizes[from] <= 1L) next
        removal <- sizes[from] / (sizes[from] - 1L) *
          sum((X[i, ] - centers[from, ])^2)
        for (to in seq_len(k)) {
          if (to == from) next
          addition <- sizes[to] / (sizes[to] + 1L) *
            sum((X[i, ] - centers[to, ])^2)
          delta <- addition - removal
          if (is.finite(delta) &&
              delta < best_delta - improvement_tolerance) {
            best_delta <- delta
            best_i <- i
            best_group <- to
          }
        }
      }
      if (is.na(best_i)) break
      proposed <- cl
      proposed[best_i] <- best_group
      proposed_summary <- summarise_partition(proposed)
      if (!is.finite(proposed_summary$withinss) ||
          proposed_summary$withinss >=
            current$withinss - improvement_tolerance)
        break
      cl <- proposed
    }
    final <- summarise_partition(cl)
    list(cluster = cl, withinss = final$withinss)
  }

  fits <- lapply(starts, optimise_start)
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) return(NULL)
  objective <- vapply(fits, `[[`, numeric(1), "withinss")
  fits[[which.min(objective)]]$cluster
}

.cluster_mega_similarity <- function(D, k, method, seed = NULL) {
  n <- nrow(D)
  if (k == 1L)
    return(stats::setNames(rep.int(1L, n), rownames(D)))
  if (n < k) return(NULL)
  z <- tryCatch({
    if (method == "hclust") {
      stats::cutree(
        stats::hclust(stats::as.dist(.mega_distance(D)), method = "average"),
        k = k
      )
    } else {
      X <- .mega_coordinates(D)
      if (nrow(unique(round(X, 12))) < k) return(NULL)
      .finite_descent_kmeans(X, k, nstart = 50L, seed = seed)
    }
  }, error = function(e) NULL)
  if (is.null(z) || length(unique(z)) != k) return(NULL)
  # Deterministic relabelling by first environment makes printed output stable.
  first <- vapply(split(seq_along(z), z), min, integer(1))
  map <- stats::setNames(seq_along(order(first)), names(first)[order(first)])
  z <- unname(map[as.character(z)])
  stats::setNames(as.integer(z), rownames(D))
}

.mean_mega_silhouette <- function(cl, distance) {
  groups <- unique(cl)
  if (length(groups) < 2L) return(NA_real_)
  s <- numeric(length(cl))
  for (i in seq_along(cl)) {
    same <- which(cl == cl[i] & seq_along(cl) != i)
    if (!length(same)) {
      s[i] <- 0
      next
    }
    a <- mean(distance[i, same])
    other <- groups[groups != cl[i]]
    b <- min(vapply(other, function(g)
      mean(distance[i, cl == g]), numeric(1)))
    den <- max(a, b)
    s[i] <- if (den <= .Machine$double.eps) 0 else (b - a) / den
  }
  mean(s)
}

.adjusted_rand_index <- function(x, y) {
  if (length(x) != length(y))
    stop("Partitions must have equal lengths.")
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  if (n < 2L) return(1)
  index <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  expected <- row_pairs * col_pairs / total_pairs
  maximum <- (row_pairs + col_pairs) / 2
  den <- maximum - expected
  if (abs(den) <= .Machine$double.eps)
    return(if (identical(as.integer(x), as.integer(y))) 1 else 0)
  (index - expected) / den
}

.assemble_mega_result <- function(cl, Sigma_E) {
  envs <- names(cl)
  groups <- sort(unique(cl))
  clusters <- lapply(groups, function(g) envs[cl == g])
  sig <- lapply(clusters, function(members)
    Sigma_E[members, members, drop = FALSE])
  names(clusters) <- names(sig) <- paste0("mega", groups)
  list(membership = cl, clusters = clusters, Sigma_E = sig)
}

.empty_mega_diagnostics <- function() {
  data.frame(
    k = integer(0), method = character(0), estimable = logical(0),
    min_cluster_size = integer(0), max_cluster_size = integer(0),
    silhouette = numeric(0), stability_mean = numeric(0),
    stability_lower = numeric(0), n_stability_comparisons = integer(0),
    n_stability_blocks = integer(0),
    separation_pass = logical(0), size_pass = logical(0),
    stability_pass = logical(0), stringsAsFactors = FALSE
  )
}
