#' Classify treatments into replication roles based on seed availability
#'
#' `assign_replication_by_seed()` evaluates seed feasibility for each
#' non-check treatment assigned to a single environment and partitions those
#' treatments into replication roles -- replicated, unreplicated, or excluded --
#' according to the requested design mode and the available seed quantity per
#' treatment. The output is intended as direct input to [met_prep_famoptg()],
#' which constructs the within-environment field layout from the resulting role
#' assignments.
#'
#' @description
#' Seed availability is the binding constraint that determines whether a target
#' replication level is deployable for any given treatment. This function makes
#' that constraint explicit by checking each treatment against a per-plot seed
#' requirement, an optional minimum buffer, and the requested replication level,
#' then assigning roles accordingly. Treatments that cannot support the target
#' replication may be downgraded to unreplicated or excluded from the design,
#' depending on `shortage_action`.
#'
#' Three design modes are supported. In `"augmented"` mode, all treatments are
#' unreplicated; the function identifies which ones have sufficient seed for a
#' single plot and flags the rest for exclusion. In `"p_rep"` mode, a subset of
#' treatments is selected for replication at `desired_replications` and the
#' remainder are unreplicated; selection is governed by `priority`, restricted
#' to `candidate_prep`, and optionally capped at `max_prep`. In `"rcbd_type"`
#' mode, all treatments are targeted for `desired_replications`; those with
#' insufficient seed are handled according to `shortage_action`.
#'
#' @details
#' Let \eqn{s_i} be the seed available for treatment \eqn{i}, \eqn{q} the seed
#' required per plot, \eqn{r} the requested replication count
#' (`desired_replications`), and \eqn{b} the minimum seed buffer
#' (`minimum_seed_buffer`). The feasibility conditions are:
#'
#' \deqn{s_i \geq q + b \quad \text{(feasible for one unreplicated plot)}}
#' \deqn{s_i \geq r \times q + b \quad \text{(feasible for } r \text{ replicated plots)}}
#'
#' A treatment that fails the first condition cannot be placed in any plot and
#' is excluded regardless of mode. A treatment that satisfies the first but not
#' the second condition can support one plot but not the target replication; its
#' fate under shortage depends on `replication_mode` and `shortage_action`.
#'
#' ## Replication modes
#'
#' **`"augmented"`**: All treatments are targeted for a single unreplicated
#' plot. Treatments satisfying \eqn{s_i \geq q + b} are assigned the role
#' `"unreplicated"`. Those that fail are assigned `"excluded"` if
#' `shortage_action` is `"exclude"` or `"downgrade"`, or trigger an error if
#' `shortage_action = "error"`.
#'
#' **`"p_rep"`**: A selection of treatments is designated for replication at
#' `desired_replications`; the remaining treatments are unreplicated.
#' Candidates for replication are drawn from `candidate_prep` (all treatments
#' if `NULL`) restricted to those satisfying the replicated feasibility
#' condition. Among those candidates, selection priority is determined by
#' `priority`: `"seed_available"` ranks by descending seed quantity,
#' `"input_order"` preserves the order of `treatments`, and `"random"` draws
#' uniformly. The number of replicated treatments is optionally capped by
#' `max_prep`. Treatments not selected for replication are evaluated for
#' single-plot feasibility and assigned `"unreplicated"` or `"excluded"`.
#'
#' **`"rcbd_type"`**: All treatments are targeted for `desired_replications`.
#' Treatments satisfying the replicated feasibility condition receive role
#' `"p_rep"` at the requested replication level. For treatments that fail,
#' behavior depends on `shortage_action`: `"error"` stops immediately;
#' `"exclude"` drops the treatment from the design; `"downgrade"` reduces the
#' treatment to a single unreplicated plot if the single-plot feasibility
#' condition is met, or excludes it otherwise.
#'
#' ## Role vocabulary
#'
#' The function assigns each treatment one of four roles, stored in the
#' `Role` column of `seed_summary`:
#'
#' - `"p_rep"`: assigned `desired_replications` plots.
#' - `"unreplicated"`: assigned exactly one plot.
#' - `"excluded"`: assigned no plots; insufficient seed.
#' - `"unused"`: treatment not matched in the seed data frame (internal flag).
#'
#' @param treatments Character vector of non-check treatment IDs assigned to
#'   the environment. Check treatments should not appear here; they are managed
#'   separately within the field design functions. Duplicate values are
#'   silently deduplicated.
#'
#' @param check_treatments Optional character vector of check IDs. If
#'   supplied, the function verifies that no treatment in `treatments` is also
#'   listed as a check, stopping with an error if overlap is found.
#'
#' @param seed_available Data frame with columns `Treatment` and
#'   `SeedAvailable`. `SeedAvailable` is the total seed quantity available for
#'   each treatment. Every treatment in `treatments` must have a matching row;
#'   the function stops if any are missing.
#'
#' @param seed_required_per_plot Positive numeric scalar. Seed quantity
#'   consumed by a single plot of one treatment in this environment. Used
#'   together with `desired_replications` and `minimum_seed_buffer` to compute
#'   feasibility thresholds.
#'
#' @param replication_mode Character. Design mode controlling how the function
#'   partitions treatments into roles. One of:
#' \describe{
#'   \item{`"augmented"`}{Targets one unreplicated plot per treatment.}
#'   \item{`"p_rep"`}{Targets replication for a feasible subset; unreplicated
#'     for the remainder.}
#'   \item{`"rcbd_type"`}{Targets `desired_replications` plots for all
#'     treatments.}
#' }
#'   The selected mode must be consistent with the field design intended for the
#'   environment. [met_prep_famoptg()] interprets the output roles directly.
#'
#' @param desired_replications Positive integer. Number of plots requested for
#'   replicated treatments. In `"augmented"` mode this value is not used for
#'   role assignment but is stored for reference. In `"p_rep"` mode it applies
#'   to the selected replicated subset. In `"rcbd_type"` mode it applies to all
#'   feasible treatments.
#'
#' @param candidate_prep Optional character vector. Treatments eligible for
#'   replication in `"p_rep"` mode. Must be a subset of `treatments`; values
#'   not present in `treatments` are silently dropped. If `NULL`, all
#'   treatments are candidates. Ignored in `"augmented"` and `"rcbd_type"`
#'   modes.
#'
#' @param priority Character. Selection rule for determining which candidate
#'   treatments are replicated in `"p_rep"` mode when the feasible candidate
#'   pool exceeds `max_prep`. One of:
#' \describe{
#'   \item{`"seed_available"`}{Ranks by descending seed quantity, breaking ties
#'     by position in `treatments`.}
#'   \item{`"input_order"`}{Follows the order of `treatments` as supplied.}
#'   \item{`"random"`}{Draws uniformly at random from feasible candidates.}
#' }
#'   Ignored in `"augmented"` and `"rcbd_type"` modes.
#'
#' @param minimum_seed_buffer Non-negative numeric scalar. Additional seed
#'   quantity subtracted from the available supply before feasibility is
#'   assessed. Use this when a safety reserve must be maintained, for example
#'   to allow for germination losses or re-planting. The buffer is applied
#'   uniformly across all treatments.
#'
#' @param max_prep Optional non-negative integer. Maximum number of treatments
#'   assigned the `"p_rep"` role in `"p_rep"` mode. If `NULL`, all treatments
#'   satisfying the replicated feasibility condition and present in
#'   `candidate_prep` are replicated. Use this when the field layout imposes an
#'   upper bound on the number of replicated plots independent of seed
#'   availability. Ignored in `"augmented"` and `"rcbd_type"` modes.
#'
#' @param shortage_action Character. Action taken when a treatment has
#'   insufficient seed for its target role. One of:
#' \describe{
#'   \item{`"error"`}{Stops immediately with an informative message listing the
#'     affected treatments.}
#'   \item{`"downgrade"`}{Reduces the role to `"unreplicated"` if the
#'     single-plot feasibility condition is met; otherwise assigns `"excluded"`.
#'     Applies in `"rcbd_type"` mode only.}
#'   \item{`"exclude"`}{Removes the treatment from the design.}
#' }
#'
#' @param seed Optional integer. Random seed for reproducibility. Affects
#'   random candidate selection in `"p_rep"` mode when `priority = "random"`.
#'   If `NULL`, no seed is set and results may differ across runs. The seed
#'   used internally is returned as `seed_used`.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`p_rep_treatments`}{Character vector of treatments assigned the
#'     `"p_rep"` role, in the order they appear in `treatments`. Empty if no
#'     treatments are replicated.}
#'   \item{`p_rep_reps`}{Integer vector of the same length as
#'     `p_rep_treatments` giving the replication count for each. All values
#'     equal `desired_replications`.}
#'   \item{`unreplicated_treatments`}{Character vector of treatments assigned
#'     exactly one plot, in the order they appear in `treatments`.}
#'   \item{`excluded_treatments`}{Character vector of treatments removed from
#'     the design due to insufficient seed, in the order they appear in
#'     `treatments`.}
#'   \item{`seed_summary`}{Data frame with one row per treatment in
#'     `treatments`. Contains columns `Treatment`, `SeedAvailable`,
#'     `Required_Unrep` (seed threshold for one plot), `Required_Rep` (seed
#'     threshold for `desired_replications` plots), `Can_Unrep` (logical),
#'     `Can_Rep` (logical), and `Role` (the assigned role string).}
#'   \item{`replication_mode`}{Character scalar. The mode used, returned for
#'     traceability.}
#'   \item{`seed_used`}{The integer seed passed to the function, or `NULL` if
#'     none was supplied.}
#' }
#'
#' @seealso [met_prep_famoptg()] to construct the within-environment field
#'   layout using the role assignments produced by this function.
#'   [plan_sparse_met_design()] for the end-to-end pipeline which calls this
#'   function internally when seed-aware replication is requested.
#'
#' @examples
#' treatments <- paste0("L", sprintf("%03d", 1:10))
#'
#' seed_df <- data.frame(
#'   Treatment     = treatments,
#'   SeedAvailable = c(50, 45, 40, 35, 30, 20, 18, 15, 12, 10),
#'   stringsAsFactors = FALSE
#' )
#'
#' ## Example 1: p_rep mode -- replicate the 4 best-seeded candidates
#' ## at 2 reps; remaining feasible treatments receive one plot.
#' out1 <- assign_replication_by_seed(
#'   treatments             = treatments,
#'   seed_available         = seed_df,
#'   seed_required_per_plot = 10,
#'   replication_mode       = "p_rep",
#'   desired_replications   = 2,
#'   priority               = "seed_available",
#'   max_prep               = 4
#' )
#'
#' out1$p_rep_treatments       # L001-L004: sufficient seed for 2 plots each
#' out1$unreplicated_treatments # remaining feasible treatments
#' out1$excluded_treatments     # treatments with seed < 10
#' out1$seed_summary
#'
#' ## Example 2: rcbd_type mode with downgrade -- all treatments targeted for
#' ## 3 reps; those below threshold are downgraded to 1 rep or excluded.
#' out2 <- assign_replication_by_seed(
#'   treatments             = treatments,
#'   seed_available         = seed_df,
#'   seed_required_per_plot = 10,
#'   replication_mode       = "rcbd_type",
#'   desired_replications   = 3,
#'   shortage_action        = "downgrade"
#' )
#'
#' out2$p_rep_treatments        # treatments with seed >= 30
#' out2$unreplicated_treatments # treatments with 10 <= seed < 30
#' out2$excluded_treatments     # treatments with seed < 10
#' out2$seed_summary
#'
#' ## Example 3: augmented mode with a seed buffer -- each treatment gets one
#' ## plot; a buffer of 5 seeds is reserved per treatment before feasibility
#' ## is assessed.
#' out3 <- assign_replication_by_seed(
#'   treatments             = treatments,
#'   seed_available         = seed_df,
#'   seed_required_per_plot = 10,
#'   replication_mode       = "augmented",
#'   minimum_seed_buffer    = 5
#' )
#'
#' out3$unreplicated_treatments # treatments with seed >= 15
#' out3$excluded_treatments     # treatments with seed < 15
#' out3$seed_summary
#'
#' @export
assign_replication_by_seed <- function(
    treatments,
    check_treatments = NULL,
    seed_available,
    seed_required_per_plot,
    replication_mode = c("augmented", "p_rep", "rcbd_type"),
    desired_replications = 2,
    candidate_prep = NULL,
    priority = c("seed_available", "input_order", "random"),
    minimum_seed_buffer = 0,
    max_prep = NULL,
    shortage_action = c("error", "downgrade", "exclude"),
    seed = NULL
) {

  replication_mode <- match.arg(replication_mode)
  priority         <- match.arg(priority)
  shortage_action  <- match.arg(shortage_action)

  seed_used <- seed
  if (!is.null(seed_used)) set.seed(seed_used)

  if (length(treatments) == 0)
    stop("`treatments` must contain at least one non-check treatment.")
  treatments <- unique(as.character(treatments))

  if (!is.null(check_treatments)) {
    check_treatments <- unique(as.character(check_treatments))
    overlap <- intersect(treatments, check_treatments)
    if (length(overlap) > 0)
      stop("Some treatments are also listed as checks: ", paste(overlap, collapse = ", "))
  }

  if (!is.data.frame(seed_available) ||
      !all(c("Treatment", "SeedAvailable") %in% names(seed_available)))
    stop("`seed_available` must be a data frame with columns `Treatment` and `SeedAvailable`.")

  if (!is.numeric(seed_required_per_plot) || length(seed_required_per_plot) != 1 ||
      is.na(seed_required_per_plot) || seed_required_per_plot <= 0)
    stop("`seed_required_per_plot` must be a single positive number.")

  if (!is.numeric(desired_replications) || length(desired_replications) != 1 ||
      is.na(desired_replications) || desired_replications < 1)
    stop("`desired_replications` must be a single integer >= 1.")
  desired_replications <- as.integer(desired_replications)

  if (!is.numeric(minimum_seed_buffer) || length(minimum_seed_buffer) != 1 ||
      is.na(minimum_seed_buffer) || minimum_seed_buffer < 0)
    stop("`minimum_seed_buffer` must be a single numeric value >= 0.")

  if (!is.null(max_prep)) {
    if (!is.numeric(max_prep) || length(max_prep) != 1 || is.na(max_prep) || max_prep < 0)
      stop("`max_prep` must be NULL or a single integer >= 0.")
    max_prep <- as.integer(max_prep)
  }

  seed_df <- seed_available[, c("Treatment", "SeedAvailable")]
  seed_df$Treatment <- as.character(seed_df$Treatment)
  seed_df <- seed_df[seed_df$Treatment %in% treatments, , drop = FALSE]

  missing_seed <- setdiff(treatments, seed_df$Treatment)
  if (length(missing_seed) > 0)
    stop("Missing seed availability for treatments: ", paste(missing_seed, collapse = ", "))

  seed_df <- seed_df[match(treatments, seed_df$Treatment), , drop = FALSE]

  req_unrep <- seed_required_per_plot + minimum_seed_buffer
  req_rep   <- desired_replications * seed_required_per_plot + minimum_seed_buffer

  seed_df$Required_Unrep <- req_unrep
  seed_df$Required_Rep   <- req_rep
  seed_df$Can_Unrep      <- seed_df$SeedAvailable >= req_unrep
  seed_df$Can_Rep        <- seed_df$SeedAvailable >= req_rep
  seed_df$Role           <- NA_character_

  if (is.null(candidate_prep)) {
    candidate_prep <- treatments
  } else {
    candidate_prep <- intersect(as.character(candidate_prep), treatments)
  }

  rank_candidates <- function(df, treatment_order, mode) {
    if (nrow(df) == 0) return(df)
    if (mode == "seed_available") {
      ord <- order(-df$SeedAvailable, match(df$Treatment, treatment_order))
      return(df[ord, , drop = FALSE])
    }
    if (mode == "input_order") {
      ord <- match(df$Treatment, treatment_order)
      return(df[order(ord), , drop = FALSE])
    }
    if (mode == "random")
      return(df[sample(seq_len(nrow(df))), , drop = FALSE])
    df
  }

  p_rep_treatments        <- character(0)
  p_rep_reps              <- integer(0)
  unreplicated_treatments <- character(0)
  excluded_treatments     <- character(0)

  if (replication_mode == "augmented") {
    feasible_unrep          <- seed_df$Treatment[seed_df$Can_Unrep]
    infeasible              <- seed_df$Treatment[!seed_df$Can_Unrep]
    unreplicated_treatments <- feasible_unrep
    excluded_treatments     <- infeasible
    seed_df$Role[seed_df$Treatment %in% unreplicated_treatments] <- "unreplicated"
    seed_df$Role[seed_df$Treatment %in% excluded_treatments]     <- "excluded"
  }

  if (replication_mode == "p_rep") {
    cand_df <- seed_df[
      seed_df$Treatment %in% candidate_prep & seed_df$Can_Rep, , drop = FALSE
    ]
    cand_df <- rank_candidates(cand_df, treatments, priority)
    if (!is.null(max_prep)) cand_df <- head(cand_df, max_prep)

    p_rep_treatments <- cand_df$Treatment
    p_rep_reps       <- rep(desired_replications, length(p_rep_treatments))

    remaining    <- setdiff(treatments, p_rep_treatments)
    remaining_df <- seed_df[seed_df$Treatment %in% remaining, , drop = FALSE]

    feasible_unrep <- remaining_df$Treatment[remaining_df$Can_Unrep]
    infeasible     <- remaining_df$Treatment[!remaining_df$Can_Unrep]

    if (length(infeasible) > 0 && shortage_action == "error")
      stop("Some remaining treatments do not have enough seed for one plot: ",
           paste(infeasible, collapse = ", "))
    if (length(infeasible) > 0) excluded_treatments <- infeasible
    unreplicated_treatments <- feasible_unrep

    seed_df$Role[seed_df$Treatment %in% p_rep_treatments]        <- "p_rep"
    seed_df$Role[seed_df$Treatment %in% unreplicated_treatments] <- "unreplicated"
    seed_df$Role[seed_df$Treatment %in% excluded_treatments]     <- "excluded"
  }

  if (replication_mode == "rcbd_type") {
    feasible_rep <- seed_df$Treatment[seed_df$Can_Rep]
    not_rep      <- seed_df$Treatment[!seed_df$Can_Rep]

    if (length(not_rep) > 0 && shortage_action == "error")
      stop("Some treatments do not have enough seed for the requested replication: ",
           paste(not_rep, collapse = ", "))

    if (length(not_rep) > 0 && shortage_action == "exclude") {
      excluded_treatments <- not_rep
      p_rep_treatments    <- feasible_rep
    }

    if (length(not_rep) > 0 && shortage_action == "downgrade") {
      downgradable        <- not_rep[seed_df$Can_Unrep[match(not_rep, seed_df$Treatment)]]
      excluded_treatments <- setdiff(not_rep, downgradable)
      unreplicated_treatments <- downgradable
      p_rep_treatments    <- feasible_rep
    }

    if (length(not_rep) == 0) p_rep_treatments <- feasible_rep
    p_rep_reps <- rep(desired_replications, length(p_rep_treatments))

    seed_df$Role[seed_df$Treatment %in% p_rep_treatments]        <- "p_rep"
    seed_df$Role[seed_df$Treatment %in% unreplicated_treatments] <- "unreplicated"
    seed_df$Role[seed_df$Treatment %in% excluded_treatments]     <- "excluded"
  }

  p_rep_treatments <- treatments[treatments %in% p_rep_treatments]
  p_rep_reps <- if (length(p_rep_treatments) > 0)
    rep(desired_replications, length(p_rep_treatments)) else integer(0)

  unreplicated_treatments <- treatments[treatments %in% unreplicated_treatments]
  excluded_treatments     <- treatments[treatments %in% excluded_treatments]
  seed_df$Role[is.na(seed_df$Role)] <- "unused"

  list(
    p_rep_treatments        = p_rep_treatments,
    p_rep_reps              = p_rep_reps,
    unreplicated_treatments = unreplicated_treatments,
    excluded_treatments     = excluded_treatments,
    seed_summary            = seed_df,
    replication_mode        = replication_mode,
    seed_used               = seed_used
  )
}
