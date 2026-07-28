#' Enumerate candidate interaction terms among covariates
#'
#' @description
#' Lists every possible interaction among management/environmental covariate
#' columns so a breeder can inspect the options before defining a short,
#' pre-specified set with [build_variable_interaction_kernels()] or adding
#' explicit covariate columns with [add_management_interactions()]. Each candidate is
#' labelled by type (categorical, numeric, or mixed) and by the number of columns
#' it would generate -- the dimensionality cost -- so cheap, low-risk terms are
#' obvious and the explosion of high-order terms is visible. Enumeration is not a
#' recommendation: with few environments, most interactions add noise, so pair
#' this with agronomic judgement ([agronomic_interaction_hypotheses()]) or the
#' data-driven [screen_interactions()]. Exhaustive enumeration is not a
#' recommendation to construct every listed kernel.
#'
#' @param x A data frame of covariates (column classes are used), or a named list
#'   giving, per column, either `"numeric"` or the character vector of factor
#'   levels.
#' @param order Maximum interaction order (2 = pairwise, 3 = three-way, ...).
#'   Default 2.
#' @param exclude Column names to exclude (e.g. an `environment` id, or doses
#'   already handled by [nest_dose_within_type()]). `"environment"` is always
#'   excluded.
#' @return A data frame with one row per candidate: `term` (colon-joined
#'   variables), `order`, `type` (`"numeric"`, `"categorical"`, or `"mixed"`),
#'   and `est_columns` (columns it would generate), sorted cheapest first.
#' @seealso [build_variable_interaction_kernels()],
#'   [add_management_interactions()], [screen_interactions()],
#'   [agronomic_interaction_hypotheses()].
#' @examples
#' df <- data.frame(temp = 1:3, rain = 4:6,
#'                  planting = c("a", "b", "a"), irrigation = c("x", "y", "x"))
#' list_candidate_interactions(df)
#' @export
list_candidate_interactions <- function(x, order = 2L, exclude = character(0)) {
  # column -> number of levels (numeric = 1) and type
  levels_of <- function(nm, col) {
    if (is.list(x) && !is.data.frame(x)) {
      v <- x[[nm]]
      if (is.character(v) && length(v) == 1L && v == "numeric")
        return(list(k = 1L, type = "numeric"))
      return(list(k = length(v), type = "categorical"))
    }
    if (is.numeric(col)) list(k = 1L, type = "numeric")
    else list(k = length(unique(col[!is.na(col)])), type = "categorical")
  }

  nms <- if (is.data.frame(x)) names(x) else names(x)
  nms <- setdiff(nms, unique(c("environment", exclude)))
  if (length(nms) < 2L) stop("Need at least two covariates to form interactions.")
  info <- lapply(nms, function(nm)
    levels_of(nm, if (is.data.frame(x)) x[[nm]] else NULL))
  names(info) <- nms

  rows <- list()
  for (o in 2:max(2L, order)) {
    if (o > length(nms)) break
    combos <- utils::combn(nms, o, simplify = FALSE)
    for (cb in combos) {
      ks    <- vapply(cb, function(nm) info[[nm]]$k, integer(1))
      types <- vapply(cb, function(nm) info[[nm]]$type, character(1))
      type  <- if (all(types == "numeric")) "numeric"
               else if (all(types == "categorical")) "categorical" else "mixed"
      rows[[length(rows) + 1L]] <- data.frame(
        term = paste(cb, collapse = ":"), order = o, type = type,
        est_columns = prod(ks), stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  tab[order(tab$est_columns, tab$order), , drop = FALSE]
}


#' Screen environmental covariates/interactions by variance explained in G x E
#'
#' @description
#' The principled way to decide which covariates or interactions to keep is to
#' ask which of them explain the observed genotype-by-environment interaction
#' (G x E) in historical multi-environment data. `screen_interactions()` extracts
#' the environment-level G x E signal from a genotype-by-environment matrix
#' (double-centring to remove genotype and environment main effects, then an SVD
#' to get the leading environment G x E scores) and ranks each environmental
#' covariate by how well it explains those scores. Add only the top few -- the
#' rest are noise, which is what overfits a kinship built from many interactions.
#'
#' @param gxe Genotype-by-environment numeric matrix of adjusted means
#'   (BLUEs/BLUPs); row names genotypes, column names environments.
#' @param env_covariates Environment-by-covariate matrix (row names
#'   environments), typically including candidate interaction columns already
#'   built with [add_management_interactions()].
#' @param method `"anova"` (default; R^2 of the covariate against the leading
#'   G x E scores, no extra dependency) or `"rf"` (random-forest importance, needs
#'   \pkg{ranger}).
#' @param n_scores Number of leading environment G x E scores to explain.
#'   Default 2.
#' @param top Optional number of top covariates to return; default all.
#' @return A data frame sorted by `importance` (descending): `covariate`,
#'   `importance`, and `method`.
#' @seealso [list_candidate_interactions()],
#'   [build_variable_interaction_kernels()],
#'   [assess_variable_interactions()], [add_management_interactions()],
#'   [build_environment_relationship()].
#' @examples
#' set.seed(1)
#' envs <- paste0("E", 1:6); x <- c(-2, -1, 0, 1, 2, 3)
#' slopes <- rnorm(20)
#' gxe <- outer(slopes, x) + matrix(rnorm(120, 0, 0.1), 20)
#' dimnames(gxe) <- list(paste0("L", 1:20), envs)
#' cov <- cbind(x = x, noise = rnorm(6)); rownames(cov) <- envs
#' screen_interactions(gxe, cov)
#' @export
screen_interactions <- function(gxe, env_covariates,
                                method = c("anova", "rf"),
                                n_scores = 2L, top = NULL) {
  method <- match.arg(method)
  G <- as.matrix(gxe)
  C <- as.matrix(env_covariates)

  # align environments
  common <- intersect(colnames(G), rownames(C))
  if (length(common) < 3L)
    stop("Need at least 3 shared environments between `gxe` columns and ",
         "`env_covariates` rows.")
  G <- G[, common, drop = FALSE]; C <- C[common, , drop = FALSE]

  # residual G x E: remove genotype and environment main effects (double-centre)
  rm <- rowMeans(G, na.rm = TRUE); cm <- colMeans(G, na.rm = TRUE); gm <- mean(G, na.rm = TRUE)
  R <- G - outer(rm, rep(1, ncol(G))) - outer(rep(1, nrow(G)), cm) + gm
  R[is.na(R)] <- 0
  sv <- svd(R)
  k  <- min(n_scores, length(sv$d))
  scores <- sv$v[, seq_len(k), drop = FALSE]          # env x k
  wts <- sv$d[seq_len(k)]^2; wts <- wts / sum(wts)

  covs <- colnames(C)
  if (method == "rf") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      warning("Package 'ranger' not installed; falling back to method = 'anova'.",
              call. = FALSE)
      method <- "anova"
    }
  }

  if (method == "rf") {
    imp <- rowSums(vapply(seq_len(k), function(j) {
      dat <- data.frame(y = scores[, j], C)
      fit <- ranger::ranger(y ~ ., data = dat, importance = "impurity",
                            num.trees = 500)
      iv <- fit$variable.importance
      iv <- iv[covs]; iv[is.na(iv)] <- 0
      wts[j] * iv / max(sum(iv), 1e-12)
    }, numeric(length(covs))))
    names(imp) <- covs
  } else {
    imp <- vapply(covs, function(cv) {
      sum(vapply(seq_len(k), function(j) {
        v <- C[, cv]
        if (stats::sd(v) == 0) return(0)
        wts[j] * summary(stats::lm(scores[, j] ~ v))$r.squared
      }, numeric(1)))
    }, numeric(1))
  }

  out <- data.frame(covariate = covs, importance = as.numeric(imp),
                    method = method, stringsAsFactors = FALSE)
  out <- out[order(-out$importance), , drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(top)) out <- utils::head(out, top)
  out
}


#' Curated agronomic interaction hypotheses
#'
#' @description
#' The primary basis for adding an interaction is agronomic knowledge: a
#' combination that physiology says should drive genotype-by-environment
#' interaction. `agronomic_interaction_hypotheses()` returns a curated shortlist
#' of physiologically plausible environmental and management interactions with
#' their rationale, as a hypothesis-driven starting point. Map the generic
#' roles to programme-specific post-quality-control column names, then construct
#' dedicated kernels with [build_variable_interaction_kernels()]. Use
#' [add_management_interactions()] only when explicit covariate columns are
#' required, and [nest_dose_within_type()] when doses must be nested within
#' product type.
#'
#' @return A data frame with `domain`, `factor_1`, `factor_2`, and `rationale`.
#' @seealso [build_variable_interaction_kernels()],
#'   [assess_variable_interactions()], [add_management_interactions()],
#'   [nest_dose_within_type()], [screen_interactions()].
#' @examples
#' agronomic_interaction_hypotheses()
#' @export
agronomic_interaction_hypotheses <- function() {
  data.frame(
    domain = c("abiotic stress", "evaporative demand", "soil compaction",
               "resource use",
               "nutrient x water", "canopy", "establishment", "phenology",
               "soil water", "nutrient availability", "disease"),
    factor_1 = c("temperature", "temperature", "temperature", "radiation",
                 "nitrogen",
                 "planting density", "planting mode", "sowing date",
                 "soil clay content", "soil pH", "humidity"),
    factor_2 = c("water availability", "vapour pressure deficit / humidity",
                 "root-zone bulk density", "water availability",
                 "water availability", "nitrogen",
                 "irrigation method", "temperature (photothermal time)",
                 "rainfall", "phosphorus / micronutrients", "temperature"),
    rationale = c(
      "Heat under drought is a compounding stress; ranking flips vs heat under moisture.",
      "High temperature with high VPD drives transpirational stress beyond either alone.",
      "Bulk density constrains rooting and water extraction; its consequence can intensify under high evaporative demand.",
      "Radiation-use efficiency depends on water status; light only helps if water allows.",
      "Nitrogen uptake and response require moisture; N x water is a classic crossover.",
      "Optimal density depends on N supply through competition for light and nutrients.",
      "Establishment regime (transplanting vs direct seeding) interacts with water delivery.",
      "Photothermal time: development responds to temperature conditional on sowing timing.",
      "Water-holding capacity (clay) modulates how rainfall translates to available water.",
      "pH controls availability of P and micronutrients, changing fertiliser response.",
      "Many diseases require warm, humid combinations; pressure is a temperature x humidity effect."),
    stringsAsFactors = FALSE)
}
