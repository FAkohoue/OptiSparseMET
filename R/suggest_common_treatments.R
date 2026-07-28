#' Fast heuristic screen for common (connectivity) treatments
#'
#' @description
#' Common treatments are tested at **every** site to connect the sparse trials so
#' that G x E, G x E x M and between-site genetic correlations can be estimated.
#' This function is a lightweight screening rule for early planning. It assumes
#' that every common treatment receives the replication supplied in
#' `reps_per_site`; it does not optimise replication.
#'
#' For the statistically robust production workflow, use
#' [optimize_common_treatments()]. That function compares common-set sizes under
#' variance uncertainty and returns binary common presence. Local replication
#' is assigned afterwards from the remaining seed inventory and never increases
#' connectivity.
#'
#' If the breeder supplies `n_common`, that number is used. Otherwise the package
#' suggests it from three requirements a good check set must meet:
#' \enumerate{
#'   \item **Seed feasibility** -- a common must have enough seed to be sown at
#'     all sites (\eqn{\sum_s q_s r_s} grams); only those enter the candidate pool.
#'   \item **Representativeness** -- every family *and* every genetic group is
#'     included, so the connectivity anchor spans the whole target population.
#'   \item **Estimability** -- enough commons to estimate the site correlations
#'     and G x E with acceptable precision (a connectivity target from the
#'     environment correlation, with a floor scaled to the number of sites).
#' }
#' The chosen set is seed-feasible, covers all families/groups, and is filled out
#' with the most genetically *central* (representative) candidates.
#'
#' @param G Genomic relationship matrix over the treatments (row/col names).
#' @param treatment_info Data frame with `Treatment` and `Family`.
#' @param seed_info Data frame with `Treatment` and `SeedAvailable`.
#' @param seed_required_per_plot Data frame (`Environment`,
#'   `SeedRequiredPerPlot`) or a numeric vector of per-site seed-per-plot.
#' @param n_common If supplied, the number to select (representative,
#'   seed-feasible). If `NULL`, the number is suggested (see Description).
#' @param reps_per_site Replicates of a common at each site (scalar or per-site
#'   vector) -- used for the seed requirement.
#' @param Sigma_E Optional environment covariance; sets the connectivity target
#'   for the suggested number.
#' @param n_test_entries Typical entries per site (for the connectivity calc);
#'   defaults to `n_treatments / n_sites`.
#' @param target_se Target standard error of an estimated site correlation.
#'   Default 0.15.
#' @param min_per_family Minimum commons per family. Default 1.
#' @param n_groups Number of genetic groups to guarantee coverage of (defaults to
#'   the number of families).
#' @param gxe_floor Minimum commons for G x E estimability (defaults to
#'   `max(20, 2 * n_sites)`).
#' @param seed Optional RNG seed (group clustering).
#' @return A list with `selected` (the common treatment ids), `n_common`,
#'   `seed_need`, `n_feasible`, `n_families`, `n_groups`, and `rationale`.
#' @seealso [optimize_common_treatments()], [recommend_common_treatments()],
#'   [tune_common_treatments()], [select_individuals()].
#' @export
suggest_common_treatments <- function(G, treatment_info, seed_info,
                                      seed_required_per_plot, n_common = NULL,
                                      reps_per_site = 1, Sigma_E = NULL,
                                      n_test_entries = NULL, target_se = 0.15,
                                      min_per_family = 1, n_groups = NULL,
                                      gxe_floor = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G  <- as.matrix(G)
  ti <- as.data.frame(treatment_info, stringsAsFactors = FALSE)
  si <- as.data.frame(seed_info, stringsAsFactors = FALSE)
  if (!all(c("Treatment", "Family") %in% names(ti)))
    stop("`treatment_info` needs columns Treatment, Family.")
  if (!all(c("Treatment", "SeedAvailable") %in% names(si)))
    stop("`seed_info` needs columns Treatment, SeedAvailable.")
  trts   <- as.character(ti$Treatment)
  fam    <- stats::setNames(as.character(ti$Family), trts)
  seedav <- stats::setNames(as.numeric(si$SeedAvailable), as.character(si$Treatment))

  ## 1. seed need to be a common at ALL sites, and the feasible pool
  req  <- if (is.data.frame(seed_required_per_plot))
    as.numeric(seed_required_per_plot$SeedRequiredPerPlot) else
    as.numeric(seed_required_per_plot)
  n_sites <- length(req)
  reps <- if (length(reps_per_site) == 1L) rep(reps_per_site, n_sites) else reps_per_site
  seed_need <- sum(req * reps)
  feasible <- intersect(trts, names(seedav)[seedav[trts] >= seed_need])
  if (!length(feasible))
    stop("No treatment has enough seed (", seed_need,
         ") to be a common at all sites.")

  ## genetic groups over the feasible pool (for representativeness coverage)
  n_fam <- length(unique(fam))
  k_grp <- if (is.null(n_groups)) n_fam else as.integer(n_groups)
  grp   <- .cluster_G_groups(G[feasible, feasible, drop = FALSE], k_grp)

  ## 2. suggested number (only when the user did not fix it)
  if (is.null(gxe_floor)) gxe_floor <- max(20L, 2L * n_sites)
  n_conn <- 0L
  if (!is.null(Sigma_E)) {
    nte <- if (is.null(n_test_entries)) max(10L, length(trts) %/% max(1L, n_sites))
           else n_test_entries
    rr <- suggest_common_range(n_test_entries_per_environment = nte,
            n_treatments = length(trts), n_environments = n_sites,
            rho = .mean_offdiag_cor(Sigma_E), target_se = target_se)
    n_conn <- max(rr$common_values)
  }
  n_cover <- max(min_per_family * n_fam, length(unique(grp)))
  target <- if (is.null(n_common))
    min(length(feasible), max(n_cover, n_conn, gxe_floor)) else
    min(as.integer(n_common), length(feasible))

  ## 3. select: families first, then genetic groups, then fill by centrality
  chosen <- character(0)
  for (f in unique(fam[feasible])) {
    if (length(chosen) >= target) break
    cand <- feasible[fam[feasible] == f]
    chosen <- unique(c(chosen, utils::head(.rank_representative(G, cand),
                                           min_per_family)))
  }
  chosen <- utils::head(chosen, target)
  for (g in unique(grp)) {
    if (length(chosen) >= target) break
    if (!any(grp[intersect(chosen, names(grp))] == g)) {
      cand <- setdiff(names(grp)[grp == g], chosen)
      if (length(cand)) chosen <- c(chosen, .rank_representative(G, cand)[1])
    }
  }
  chosen <- unique(chosen)
  if (length(chosen) < target) {
    remain <- setdiff(feasible, chosen)
    chosen <- c(chosen, utils::head(.rank_representative(G, remain),
                                    target - length(chosen)))
  }
  chosen <- unique(chosen)

  list(selected = chosen, n_common = length(chosen), seed_need = seed_need,
       n_feasible = length(feasible), n_families = n_fam,
       n_groups = length(unique(grp)),
       rationale = list(seed_feasible = length(feasible),
                        family_coverage = min_per_family * n_fam,
                        group_coverage = length(unique(grp)),
                        connectivity = n_conn, gxe_floor = gxe_floor,
                        user_supplied = !is.null(n_common)))
}


# ---- helpers ----------------------------------------------------------------

# Rank ids by genetic centrality (mean relationship to the whole population):
# central lines are the most representative connectivity anchors.
.rank_representative <- function(G, ids) {
  ids <- intersect(ids, rownames(G))
  if (!length(ids)) return(character(0))
  cent <- rowMeans(G[ids, , drop = FALSE])
  ids[order(-cent)]
}

# Cluster the relationship matrix into k genetic groups (k-means on its leading
# eigenvectors).
.cluster_G_groups <- function(Gsub, k) {
  n <- nrow(Gsub)
  k <- max(1L, min(as.integer(k), n - 1L))
  if (k <= 1L) return(stats::setNames(rep(1L, n), rownames(Gsub)))
  ev <- eigen((Gsub + t(Gsub)) / 2, symmetric = TRUE)
  pc <- ev$vectors[, seq_len(min(5L, ncol(Gsub))), drop = FALSE]
  cl <- tryCatch(stats::kmeans(pc, centers = k, nstart = 5)$cluster,
                 error = function(e) rep(1L, n))
  stats::setNames(cl, rownames(Gsub))
}

.mean_offdiag_cor <- function(Sigma_E) {
  R <- tryCatch(stats::cov2cor(as.matrix(Sigma_E)), error = function(e) as.matrix(Sigma_E))
  m <- mean(R[upper.tri(R)]); if (!is.finite(m)) 0.5 else max(0, min(0.999, m))
}
