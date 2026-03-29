# ==============================================================================
# met_prep_famoptg.R
#
# OptiSparseMET version of prep_famoptg() from OptiDesign.
# The met_ prefix prevents namespace conflicts when both packages are loaded.
# Function body is identical to prep_famoptg(). Only the name changes.
# ==============================================================================

#' Construct a repeated-check block design with flexible replication
#'
#' @description
#' `met_prep_famoptg()` is the OptiSparseMET version of `prep_famoptg()` from
#' the OptiDesign package. The `met_` prefix avoids namespace conflicts when
#' both packages are loaded simultaneously. All arguments, return values, and
#' internal logic are identical to `prep_famoptg()` in OptiDesign.
#'
#' `met_prep_famoptg()` builds a repeated-check block design for plant breeding
#' and agronomic experiments. Check treatments are included in every block,
#' while non-check treatments are allocated according to user-specified
#' replication levels. Depending on the treatment structure supplied, the
#' same function generates three design classes:
#'
#' - **Augmented repeated-check design**: checks repeated in every block,
#'   all test entries unreplicated (`p_rep_treatments = NULL`).
#' - **Partially replicated (p-rep) repeated-check design**: a mixture of
#'   replicated and unreplicated non-check entries.
#' - **RCBD-type repeated-check design**: all non-check entries given the
#'   same replication across blocks.
#'
#' **Design evaluation has been separated from construction.** This function
#' returns the field book, layout matrix, and seed only. To compute A, D, and
#' CDmean optimality criteria call [met_evaluate_famoptg_efficiency()] on the
#' returned `field_book`. To search for an optimised design across many
#' randomisations call [met_optimize_famoptg()].
#'
#' @details
#' ## Design structure and allocation rules
#'
#' Let:
#' - \eqn{c} = number of check treatments
#' - \eqn{b} = `n_blocks`
#' - \eqn{v_p} = number of p-rep treatments with replication counts
#'   \eqn{r_1, r_2, \ldots, r_{v_p}}
#' - \eqn{v_u} = number of unreplicated treatments
#'
#' Total plots required:
#' \deqn{\text{total} = b \times c + \sum_{i=1}^{v_p} r_i + v_u}
#'
#' Allocation rules enforced by construction:
#' \describe{
#'   \item{Checks}{Appear in every block. Each check contributes \eqn{b}
#'     plots total.}
#'   \item{P-rep treatments}{Each treatment \eqn{i} appears in exactly
#'     \eqn{r_i} blocks, always in distinct blocks - never twice in the
#'     same block. This is the core p-rep constraint.}
#'   \item{Unreplicated treatments}{Appear exactly once in the full design,
#'     distributed as evenly as possible across blocks.}
#' }
#'
#' ## Three design classes
#'
#' **Augmented repeated-check design** - obtained when `p_rep_treatments`
#' is empty and all non-check entries are supplied via
#' `unreplicated_treatments`. Checks are the only repeated treatments.
#'
#' **P-rep repeated-check design** - obtained when some non-check entries
#' in `p_rep_treatments` have `p_rep_reps > 1` and others are unreplicated.
#' The most general use case: a mixture of replicated candidate entries and
#' single-plot entries alongside repeated checks.
#'
#' **RCBD-type repeated-check design** - obtained when all non-check
#' treatments are supplied via `p_rep_treatments` with a common replication
#' count greater than 1. If that count equals `n_blocks`, every non-check
#' treatment appears once in every block - the closest repeated-check
#' analogue of a classical RCBD.
#'
#' ## Field dimension handling
#'
#' If the supplied `n_rows x n_cols` field does not exactly match the
#' required total plots, the function either stops (when
#' `warn_and_correct = FALSE`) or silently adjusts one dimension
#' (`warn_and_correct = TRUE`). When adjusting, `fix_rows = TRUE` holds
#' `n_rows` fixed and expands `n_cols`; `fix_rows = FALSE` does the reverse.
#' Surplus cells remain `NA` in the layout matrix.
#'
#' ## Within-block arrangement
#'
#' Treatments are shuffled within each block over `attempts` random
#' permutations to minimise the count of adjacent same-group treatment pairs
#' in the 1D block ordering. The group label used for adjacency scoring is
#' determined by `cluster_source` (see below). Checks are inserted at
#' positions controlled by `check_placement`:
#' \describe{
#'   \item{`"random"`}{Check positions drawn uniformly.}
#'   \item{`"systematic"`}{Check positions spaced evenly through the block.}
#'   \item{`"optimal"`}{Runs `check_opt_attempts` candidate layouts and
#'     keeps the one maximising the mean nearest-neighbour distance between
#'     check plots in the field - encouraging spatial spread of checks.}
#' }
#'
#' ## Grouping and adjacency scoring
#'
#' `cluster_source` determines the group label used to penalise adjacent
#' same-group treatments within blocks:
#' \describe{
#'   \item{`"Family"`}{Uses `p_rep_families` and `unreplicated_families`
#'     directly. `Gcluster` column in the field book is `NA` for all plots.}
#'   \item{`"GRM"`}{Derives clusters from PCA of the genomic relationship
#'     matrix `GRM`. Number of clusters = `length(unique(p_rep_families))`.
#'     Cluster labels are prefixed `"G"`.}
#'   \item{`"A"`}{Derives clusters from PCA of the pedigree relationship
#'     matrix `A`. Cluster labels are prefixed `"A"`.}
#' }
#' Check treatments always use `check_families` for grouping and always
#' have `Gcluster = NA`.
#'
#' ## Dispersion optimisation
#'
#' When `use_dispersion = TRUE`, a local swap search runs after layout
#' construction. At each of `dispersion_iters` iterations two non-check plots
#' are selected at random; the swap is accepted if it reduces the total
#' neighbourhood relatedness score - the sum of relationship matrix values
#' between non-check neighbour pairs within Chebyshev distance
#' `dispersion_radius`. The matrix used is selected by `dispersion_source`.
#'
#' @param check_treatments Character vector of check treatment identifiers.
#'   Checks appear in every block. Must not overlap with `p_rep_treatments`
#'   or `unreplicated_treatments` and must not contain duplicates.
#'
#' @param check_families Character vector of the same length as
#'   `check_treatments`. Family labels for checks, used for adjacency
#'   scoring. Always stored in the `Family` column; `Gcluster` is always
#'   `NA` for check plots.
#'
#' @param p_rep_treatments Character vector of replicated non-check treatment
#'   identifiers. Each entry appears in exactly `p_rep_reps[i]` distinct
#'   blocks. May be `NULL` or `character(0)` for augmented designs.
#'
#' @param p_rep_reps Integer vector of the same length as `p_rep_treatments`.
#'   Number of blocks each p-rep treatment appears in. Must satisfy
#'   `1 <= p_rep_reps[i] <= n_blocks` for all \eqn{i}.
#'
#' @param p_rep_families Character vector of the same length as
#'   `p_rep_treatments`. Family labels for p-rep entries.
#'
#' @param unreplicated_treatments Character vector of treatment identifiers
#'   that appear exactly once in the full design. May be `NULL` or
#'   `character(0)`.
#'
#' @param unreplicated_families Character vector of the same length as
#'   `unreplicated_treatments`. Family labels for unreplicated entries.
#'
#' @param n_blocks Positive integer. Number of experimental blocks.
#'
#' @param n_rows Positive integer. Number of field rows.
#'
#' @param n_cols Positive integer. Number of field columns.
#'
#' @param order Character. Field traversal direction: `"column"` fills
#'   column by column; `"row"` fills row by row.
#'
#' @param serpentine Logical. If `TRUE`, alternate columns (when
#'   `order = "column"`) or alternate rows (when `order = "row"`) traverse
#'   in reverse direction, producing a boustrophedon planting path.
#'
#' @param seed Optional integer. If `NULL`, a seed is drawn from
#'   `1:.Machine$integer.max` and returned as `seed_used`. Controls all
#'   randomised steps: p-rep block assignment, unreplicated entry
#'   distribution, within-block shuffling, check placement, and dispersion
#'   optimisation.
#'
#' @param attempts Positive integer. Number of within-block shuffle attempts
#'   to minimise adjacent same-group pairs. Default 1000.
#'
#' @param warn_and_correct Logical. If `TRUE`, silently adjust field
#'   dimensions when `n_rows * n_cols` does not match the required total
#'   plots. If `FALSE`, stop with an error.
#'
#' @param fix_rows Logical. Used only when `warn_and_correct = TRUE`. If
#'   `TRUE`, hold `n_rows` fixed and expand `n_cols`; otherwise hold
#'   `n_cols` fixed and expand `n_rows`.
#'
#' @param cluster_source Character. Grouping source for adjacency scoring:
#'   `"Family"`, `"GRM"`, or `"A"`.
#'
#' @param GRM Optional square numeric matrix with rownames and colnames equal
#'   to line IDs. Required when `cluster_source = "GRM"` or
#'   `dispersion_source = "GRM"`.
#'
#' @param A Optional square numeric matrix. Required when
#'   `cluster_source = "A"` or `dispersion_source = "A"`.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'   Required when treatment labels do not match relationship matrix rownames.
#'
#' @param cluster_method Character. Clustering algorithm: `"kmeans"` or
#'   `"hclust"` (Ward's D2 linkage).
#'
#' @param cluster_seed Integer. Seed for k-means initialisation, run in an
#'   isolated RNG scope.
#'
#' @param cluster_attempts Positive integer. Number of k-means random
#'   restarts (`nstart`).
#'
#' @param n_pcs_use Positive number or `Inf`. Leading PCs retained for
#'   matrix-based clustering. Actual number used is
#'   `min(n_pcs_use, n_positive_eigenvalues - 1)`.
#'
#' @param check_placement Character. Check position rule within blocks:
#' \describe{
#'   \item{`"random"`}{Check positions drawn uniformly.}
#'   \item{`"systematic"`}{Check positions spaced evenly through the block.}
#'   \item{`"optimal"`}{Keeps the layout that maximises mean nearest-neighbour
#'     distance between check plots in the field across `check_opt_attempts`
#'     candidates.}
#' }
#'
#' @param check_opt_attempts Positive integer. Number of candidate layouts
#'   evaluated when `check_placement = "optimal"`. Default 200.
#'
#' @param use_dispersion Logical. If `TRUE`, apply a post-layout swap-based
#'   dispersion optimisation to reduce pairwise relatedness among neighbouring
#'   non-check plots.
#'
#' @param dispersion_source Character. Matrix used for dispersion scoring:
#'   `"K"`, `"A"`, or `"GRM"`. The corresponding argument must be non-`NULL`.
#'
#' @param dispersion_radius Positive integer. Chebyshev distance radius
#'   defining the neighbourhood for dispersion scoring. Default 1.
#'
#' @param dispersion_iters Non-negative integer. Number of swap proposals for
#'   dispersion optimisation. Default 2000.
#'
#' @param dispersion_seed Integer. Seed for the dispersion step, run in an
#'   isolated RNG scope.
#'
#' @param K Optional square numeric matrix. Used for dispersion scoring when
#'   `dispersion_source = "K"`.
#'
#' @param line_id_map Optional data frame with columns `Treatment` and
#'   `LineID`. Required when treatment labels do not match `rownames(K)` or
#'   the dispersion matrix.
#'
#' @param verbose Logical. If `TRUE`, prints a one-line summary of field
#'   dimensions, block count, and treatment counts.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{`layout_matrix`}{Character matrix of dimension `n_rows x n_cols`.
#'     Used cells contain treatment IDs; unused cells (when field is larger
#'     than required) are `NA`.}
#'   \item{`field_book`}{Data frame with one row per assigned plot. Columns:
#'     `Treatment`, `Family` (original user-supplied label), `Gcluster` (`NA`
#'     when `cluster_source = "Family"` or for check plots; otherwise the
#'     matrix-derived cluster label prefixed `"G"` or `"A"`), `Block`
#'     (integer 1 to `n_blocks`), `Plot` (sequential integer), `Row`,
#'     `Column`.}
#'   \item{`seed_used`}{Integer. The random seed used.}
#' }
#'
#' @seealso [met_evaluate_famoptg_efficiency()] to compute A, D, and CDmean
#'   optimality criteria on the returned `field_book`.
#'   [met_optimize_famoptg()] to search for a criterion-optimal design using
#'   Random Restart.
#'
#' @examples
#' ## Augmented repeated-check design: checks repeated, entries unreplicated
#' design_aug <- met_prep_famoptg(
#'   check_treatments        = c("CHK1", "CHK2"),
#'   check_families          = c("CHECK", "CHECK"),
#'   p_rep_treatments        = character(0),
#'   p_rep_reps              = integer(0),
#'   p_rep_families          = character(0),
#'   unreplicated_treatments = paste0("E", 1:80),
#'   unreplicated_families   = rep(paste0("F", 1:4), 20),
#'   n_blocks = 5, n_rows = 10, n_cols = 20
#' )
#' dim(design_aug$layout_matrix)
#' head(design_aug$field_book)
#'
#' ## P-rep design: some entries replicated twice, others once
#' design_prep <- met_prep_famoptg(
#'   check_treatments        = c("CHK1", "CHK2"),
#'   check_families          = c("CHECK", "CHECK"),
#'   p_rep_treatments        = paste0("P", 1:20),
#'   p_rep_reps              = rep(2L, 20),
#'   p_rep_families          = rep(paste0("F", 1:4), 5),
#'   unreplicated_treatments = paste0("U", 1:60),
#'   unreplicated_families   = rep(paste0("F", 1:4), 15),
#'   n_blocks = 5, n_rows = 15, n_cols = 20,
#'   check_placement = "systematic"
#' )
#' table(design_prep$field_book$Block)
#'
#' ## Evaluate separately
#' eff <- met_evaluate_famoptg_efficiency(
#'   field_book         = design_prep$field_book,
#'   n_rows = 15, n_cols = 20,
#'   check_treatments   = c("CHK1", "CHK2"),
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row = 0.10, rho_col = 0.10
#' )
#' eff$A_criterion
#'
#' @importFrom stats dist cutree hclust kmeans setNames
#' @importFrom utils head
#' @importFrom pracma mod
#' @export
met_prep_famoptg <- function(
    check_treatments,
    check_families,
    p_rep_treatments,
    p_rep_reps,
    p_rep_families,
    unreplicated_treatments,
    unreplicated_families,
    n_blocks,
    n_rows,
    n_cols,
    order              = "column",
    serpentine         = FALSE,
    seed               = NULL,
    attempts           = 1000,
    warn_and_correct   = TRUE,
    fix_rows           = TRUE,
    cluster_source     = c("Family", "GRM", "A"),
    GRM                = NULL,
    A                  = NULL,
    id_map             = NULL,
    cluster_method     = c("kmeans", "hclust"),
    cluster_seed       = 1,
    cluster_attempts   = 25,
    n_pcs_use          = Inf,
    check_placement    = c("random", "systematic", "optimal"),
    check_opt_attempts = 200,
    use_dispersion     = FALSE,
    dispersion_source  = c("K", "A", "GRM"),
    dispersion_radius  = 1,
    dispersion_iters   = 2000,
    dispersion_seed    = 1,
    K                  = NULL,
    line_id_map        = NULL,
    verbose            = TRUE
) {

  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required.")

  # ============================================================
  # 0. RNG
  # ============================================================
  seed_used <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed
  set.seed(seed_used)

  # ============================================================
  # 1. VALIDATION AND NORMALISATION
  # ============================================================
  cluster_source    <- match.arg(cluster_source)
  cluster_method    <- match.arg(cluster_method)
  check_placement   <- match.arg(check_placement)
  dispersion_source <- match.arg(dispersion_source)

  if (is.null(p_rep_treatments))        p_rep_treatments        <- character(0)
  if (is.null(p_rep_reps))              p_rep_reps               <- integer(0)
  if (is.null(p_rep_families))          p_rep_families           <- character(0)
  if (is.null(unreplicated_treatments)) unreplicated_treatments  <- character(0)
  if (is.null(unreplicated_families))   unreplicated_families    <- character(0)

  if (length(check_treatments) != length(check_families))
    stop("Length of check_families must match length of check_treatments.")
  if (length(p_rep_treatments) != length(p_rep_reps) ||
      length(p_rep_treatments) != length(p_rep_families))
    stop("Lengths of p_rep_treatments, p_rep_reps, and p_rep_families must all match.")
  if (length(unreplicated_treatments) != length(unreplicated_families))
    stop("Length of unreplicated_families must match length of unreplicated_treatments.")

  if (anyDuplicated(check_treatments))        stop("Duplicate check_treatments found.")
  if (anyDuplicated(p_rep_treatments))        stop("Duplicate p_rep_treatments found.")
  if (anyDuplicated(unreplicated_treatments)) stop("Duplicate unreplicated_treatments found.")

  if (length(intersect(check_treatments, p_rep_treatments)) > 0 ||
      length(intersect(check_treatments, unreplicated_treatments)) > 0 ||
      length(intersect(p_rep_treatments, unreplicated_treatments)) > 0)
    stop("A treatment cannot appear in more than one treatment class.")

  if (n_blocks < 1)           stop("n_blocks must be at least 1.")
  if (any(p_rep_reps < 1))    stop("All p_rep_reps must be >= 1.")
  if (any(p_rep_reps > n_blocks))
    stop("Each p-rep treatment replication count must not exceed n_blocks.")
  if (!order %in% c("column", "row"))
    stop("Invalid 'order'. Use 'row' or 'column'.")
  if (!(is.numeric(n_pcs_use) && length(n_pcs_use) == 1 &&
        (is.finite(n_pcs_use) || is.infinite(n_pcs_use)) && n_pcs_use > 0))
    stop("n_pcs_use must be a single positive number or Inf.")
  if (dispersion_radius < 1)  stop("dispersion_radius must be >= 1.")
  if (dispersion_iters  < 0)  stop("dispersion_iters must be >= 0.")
  if (check_opt_attempts < 1) stop("check_opt_attempts must be >= 1.")

  # ============================================================
  # 2. FIELD DIMENSION ACCOUNTING
  # ============================================================
  total_checks   <- n_blocks * length(check_treatments)
  total_prep     <- sum(p_rep_reps)
  total_unrep    <- length(unreplicated_treatments)
  total_required <- total_checks + total_prep + total_unrep
  field_size     <- n_rows * n_cols

  if (field_size != total_required) {
    if (warn_and_correct) {
      warning(paste0(
        "Field size (", n_rows, " x ", n_cols, " = ", field_size,
        ") does not match required plots (", total_required,
        "). Adjusting dimensions."
      ))
      if (fix_rows) {
        n_cols <- ceiling(total_required / n_rows)
      } else {
        n_rows <- ceiling(total_required / n_cols)
      }
      field_size <- n_rows * n_cols
    } else {
      stop(paste0(
        "Field size (", field_size,
        ") does not match required plots (", total_required,
        "). Adjust n_rows/n_cols or enable warn_and_correct."
      ))
    }
  }

  if (verbose) message(paste0(
    "Field = ", n_rows, " x ", n_cols,
    " | total plots = ", field_size,
    " | checks/block = ", length(check_treatments),
    " | p-rep entries = ", length(p_rep_treatments),
    " | unreplicated = ", length(unreplicated_treatments),
    " | n_blocks = ", n_blocks
  ))

  # ============================================================
  # 3. LOOKUPS AND GROUPING
  # ============================================================
  family_lookup <- setNames(
    c(check_families, p_rep_families, unreplicated_families),
    c(check_treatments, p_rep_treatments, unreplicated_treatments)
  )

  noncheck_trt <- c(p_rep_treatments, unreplicated_treatments)

  gcluster_lookup <- setNames(
    rep(NA_character_, length(c(check_treatments, noncheck_trt))),
    c(check_treatments, noncheck_trt)
  )

  if (cluster_source %in% c("GRM", "A")) {
    Kc <- if (cluster_source == "GRM") GRM else A

    if (is.null(Kc))
      stop(paste0("cluster_source = '", cluster_source, "' but matrix is NULL."))
    if (is.null(rownames(Kc)) || is.null(colnames(Kc)))
      stop("GRM/A must have rownames and colnames.")

    noncheck_fams <- c(p_rep_families, unreplicated_families)
    k_clusters    <- length(unique(noncheck_fams))
    if (k_clusters < 2)
      stop("Non-check treatments have < 2 unique families; clustering is not meaningful.")

    if (is.null(id_map)) {
      line_ids <- setNames(noncheck_trt, noncheck_trt)
    } else {
      if (!is.data.frame(id_map) || !all(c("Treatment", "LineID") %in% names(id_map)))
        stop("id_map must be a data.frame with columns: Treatment, LineID.")
      line_ids <- setNames(id_map$LineID, id_map$Treatment)
    }

    noncheck_line_ids <- unname(line_ids[noncheck_trt])
    missing_ids <- setdiff(noncheck_line_ids, rownames(Kc))
    if (length(missing_ids) > 0)
      stop("Some non-check LineIDs are not found in GRM/A rownames.")

    Ksub <- Kc[noncheck_line_ids, noncheck_line_ids, drop = FALSE]
    eg   <- eigen(Ksub, symmetric = TRUE)
    pos  <- which(eg$values > 1e-10)
    if (length(pos) < 2)
      stop("GRM/A has too few positive eigenvalues for PCA clustering.")

    max_possible_pcs <- min(length(pos), nrow(Ksub) - 1)
    n_pcs <- if (is.infinite(n_pcs_use)) {
      max_possible_pcs
    } else {
      n_pcs_req <- as.integer(n_pcs_use)
      if (n_pcs_req > max_possible_pcs) {
        warning(paste0(
          "Requested n_pcs_use = ", n_pcs_use,
          " but only ", max_possible_pcs, " PCs available; using ", max_possible_pcs, "."
        ))
      }
      min(n_pcs_req, max_possible_pcs)
    }
    if (n_pcs < 2) stop("Fewer than 2 PCs available for clustering.")

    pcs <- eg$vectors[, pos[seq_len(n_pcs)], drop = FALSE]
    pcs <- sweep(pcs, 2, sqrt(eg$values[pos[seq_len(n_pcs)]]), `*`)

    if (cluster_method == "kmeans") {
      clust <- .with_local_seed(cluster_seed,
        stats::kmeans(pcs, centers = k_clusters, nstart = cluster_attempts)$cluster
      )
    } else {
      hc    <- stats::hclust(stats::dist(pcs), method = "ward.D2")
      clust <- stats::cutree(hc, k = k_clusters)
    }

    prefix <- if (cluster_source == "GRM") "G" else "A"
    gcluster_lookup[noncheck_trt] <- paste0(prefix, clust)
  }

  get_adj_group <- function(trt) {
    if (cluster_source == "Family") {
      family_lookup[trt]
    } else if (trt %in% check_treatments) {
      family_lookup[trt]
    } else {
      gcluster_lookup[trt]
    }
  }

  # ============================================================
  # 4. FIELD POSITION BUILDER
  # ============================================================
  build_positions <- function(nr, nc, ord, serp) {
    out <- vector("list", nr * nc); kk <- 1L
    if (ord == "row") {
      for (r in seq_len(nr)) {
        cols <- seq_len(nc)
        if (serp && mod(r, 2) == 0) cols <- rev(cols)
        for (c in cols) { out[[kk]] <- c(Row = r, Column = c); kk <- kk + 1L }
      }
    } else {
      for (c in seq_len(nc)) {
        rows <- seq_len(nr)
        if (serp && mod(c, 2) == 0) rows <- rev(rows)
        for (r in rows) { out[[kk]] <- c(Row = r, Column = c); kk <- kk + 1L }
      }
    }
    ans        <- as.data.frame(do.call(rbind, out))
    ans$Row    <- as.integer(ans$Row)
    ans$Column <- as.integer(ans$Column)
    ans
  }

  # ============================================================
  # 5. BLOCK ASSIGNMENT
  # ============================================================
  p_rep_assignments <- setNames(vector("list", length(p_rep_treatments)), p_rep_treatments)
  available_blocks  <- lapply(seq_along(p_rep_treatments), function(i) sample(seq_len(n_blocks)))

  for (i in seq_along(p_rep_treatments)) {
    trt_i  <- p_rep_treatments[i]
    reps_i <- p_rep_reps[i]
    if (reps_i > length(available_blocks[[i]]))
      stop(paste0("Not enough unique blocks for p-rep treatment '", trt_i, "'."))
    p_rep_assignments[[trt_i]] <- sample(available_blocks[[i]], reps_i)
  }

  if (length(unreplicated_treatments) > 0) {
    unrep_shuffled   <- sample(unreplicated_treatments)
    block_unrep_list <- split(
      unrep_shuffled,
      rep(seq_len(n_blocks), length.out = length(unrep_shuffled))
    )
  } else {
    block_unrep_list <- vector("list", n_blocks)
  }

  # ============================================================
  # 6. WITHIN-BLOCK ARRANGEMENT
  # ============================================================
  place_checks <- function(checks, others, mode) {
    total_len <- length(checks) + length(others)
    if (mode == "random") return(sample(c(checks, others)))
    pos <- round(seq(1, total_len, length.out = length(checks) + 2L))[2:(length(checks) + 1L)]
    pos <- pmin(pmax(pos, 1L), total_len)
    out <- rep(NA_character_, total_len)
    out[pos] <- checks
    if (length(others) > 0) out[is.na(out)] <- sample(others)
    out
  }

  score_adj_block <- function(trt_vec) {
    if (length(trt_vec) <= 1L) return(0L)
    grp <- vapply(trt_vec, get_adj_group, character(1))
    sum(grp[-1L] == grp[-length(grp)], na.rm = TRUE)
  }

  build_blocks_once <- function(placement_mode) {
    blocks <- vector("list", n_blocks)
    for (b in seq_len(n_blocks)) {
      others <- character(0)
      for (trt in names(p_rep_assignments)) {
        if (b %in% p_rep_assignments[[trt]]) others <- c(others, trt)
      }
      if (!is.null(block_unrep_list[[b]]) && length(block_unrep_list[[b]]) > 0)
        others <- c(others, block_unrep_list[[b]])

      block_trt <- place_checks(check_treatments, others, placement_mode)

      best    <- block_trt
      best_sc <- score_adj_block(block_trt)
      if (best_sc > 0 && length(block_trt) > 1L) {
        for (aa in seq_len(attempts)) {
          cand <- sample(block_trt)
          sc   <- score_adj_block(cand)
          if (sc < best_sc) { best_sc <- sc; best <- cand; if (best_sc == 0L) break }
        }
      }

      blocks[[b]] <- data.frame(
        Treatment = best,
        Family    = unname(family_lookup[best]),
        Gcluster  = unname(gcluster_lookup[best]),
        Block     = b,
        stringsAsFactors = FALSE
      )
    }
    blocks
  }

  # ============================================================
  # 7. CHECK PLACEMENT STRATEGY
  # ============================================================
  if (check_placement != "optimal") {
    blocks <- build_blocks_once(check_placement)
  } else {
    best_score  <- -Inf
    best_blocks <- NULL
    pos_mat     <- build_positions(n_rows, n_cols, order, serpentine)

    for (k in seq_len(check_opt_attempts)) {
      cand_blocks <- build_blocks_once("random")
      fb_cand     <- do.call(rbind, cand_blocks)
      n_cand      <- nrow(fb_cand)
      fb_cand$Row    <- pos_mat$Row[seq_len(n_cand)]
      fb_cand$Column <- pos_mat$Column[seq_len(n_cand)]

      is_chk <- fb_cand$Treatment %in% check_treatments
      chk_rc <- unique(cbind(fb_cand$Row[is_chk], fb_cand$Column[is_chk]))
      score  <- if (nrow(chk_rc) <= 1L) 0 else {
        d <- as.matrix(stats::dist(chk_rc)); diag(d) <- Inf
        mean(apply(d, 1, min))
      }
      if (score > best_score) { best_score <- score; best_blocks <- cand_blocks }
    }
    blocks <- best_blocks
  }

  final_data   <- do.call(rbind, blocks)
  n_assigned   <- nrow(final_data)
  rownames(final_data) <- NULL

  pos_mat            <- build_positions(n_rows, n_cols, order, serpentine)
  final_data$Plot    <- seq_len(n_assigned)
  final_data$Row     <- pos_mat$Row[seq_len(n_assigned)]
  final_data$Column  <- pos_mat$Column[seq_len(n_assigned)]

  # ============================================================
  # 8. LAYOUT MATRIX
  # ============================================================
  layout_matrix <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(n_assigned))
    layout_matrix[final_data$Row[i], final_data$Column[i]] <- final_data$Treatment[i]

  # ============================================================
  # 9. OPTIONAL GENETIC DISPERSION
  # ============================================================
  if (isTRUE(use_dispersion)) {
    Kdisp <- switch(dispersion_source, "K" = K, "A" = A, "GRM" = GRM)
    if (is.null(Kdisp))
      stop("use_dispersion = TRUE but selected dispersion matrix is NULL.")
    if (is.null(rownames(Kdisp)) || is.null(colnames(Kdisp)))
      stop("Dispersion matrix must have rownames and colnames.")

    seed_disp <- if (!is.null(dispersion_seed)) dispersion_seed else seed_used
    trt       <- as.character(final_data$Treatment)
    is_check  <- trt %in% check_treatments
    movable   <- which(!is_check)

    if (length(movable) >= 2L && dispersion_iters > 0L) {
      non_trt <- unique(trt[!is_check])

      if (is.null(line_id_map)) {
        line_ids <- setNames(non_trt, non_trt)
      } else {
        if (!is.data.frame(line_id_map) ||
            !all(c("Treatment", "LineID") %in% names(line_id_map)))
          stop("line_id_map must be a data.frame with columns: Treatment, LineID.")
        line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
      }

      ids  <- unname(line_ids[non_trt])
      miss <- setdiff(ids, rownames(Kdisp))
      if (length(miss) > 0)
        stop("Some non-check LineIDs are missing in the dispersion matrix.")

      Ksub           <- Kdisp[ids, ids, drop = FALSE]
      colnames(Ksub) <- non_trt
      rownames(Ksub) <- non_trt

      pairs    <- .build_neighbor_pairs(final_data$Row, final_data$Column,
                                         radius = dispersion_radius)
      best_trt <- trt
      best_sc  <- .score_dispersion(best_trt, is_check, Ksub, pairs)

      .with_local_seed(seed_disp, {
        for (it in seq_len(dispersion_iters)) {
          ij   <- sample(movable, 2L, replace = FALSE)
          cand <- best_trt; cand[ij] <- cand[rev(ij)]
          sc   <- .score_dispersion(cand, is_check, Ksub, pairs)
          if (sc < best_sc) { best_sc <- sc; best_trt <- cand }
        }
      })

      final_data$Treatment <- best_trt
      final_data$Family    <- unname(family_lookup[best_trt])
      final_data$Gcluster  <- unname(gcluster_lookup[best_trt])

      layout_matrix <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
      for (i in seq_len(n_assigned))
        layout_matrix[final_data$Row[i], final_data$Column[i]] <- final_data$Treatment[i]
    }
  }

  # ============================================================
  # 10. OUTPUT
  # ============================================================
  field_book <- final_data[, c("Treatment", "Family", "Gcluster",
                                "Block", "Plot", "Row", "Column")]

  list(
    layout_matrix = layout_matrix,
    field_book    = field_book,
    seed_used     = seed_used
  )
}
