#' Plot a cost-minimisation Pareto frontier
#'
#' @description
#' Plots one statistical endpoint against normalised cost minimisation. The
#' bottom axis runs from 0% (the most expensive observed design) to 100% (the
#' least expensive observed design). A secondary top axis reports the
#' corresponding real cost in decreasing order. Therefore, a movement to the
#' right always means stronger cost minimisation, whereas a movement upwards
#' always means stronger statistical performance.
#'
#' Only feasible non-dominated designs are joined, using straight segments
#' between the observed design points. When a `feasible` column is present,
#' infeasible designs remain visible as grey crosses but cannot define the
#' frontier. Dominated feasible designs remain visible as grey circles. Genetic
#' gain is expressed as a percentage of the maximum observed gain; accuracy
#' and reliability are expressed directly as percentages. Gain and accuracy
#' are never overlaid on different vertical scales.
#'
#' @param pf A [pareto_designs()] result, its `frontier` data frame, a
#'   [benchmark_met_designs()] result, or its `design_summary` data frame.
#' @param x Real resource column: `"cost"` (default) or `"plots"`.
#' @param highlight_pareto Logical. Highlight and join non-dominated designs.
#' @param col_gain,col_acc Colours retained for backwards compatibility.
#'   `col_gain` is used for gain and `col_acc` for accuracy or reliability.
#' @param metric One endpoint: `"accuracy"`, `"gain"`, or `"reliability"`.
#'   The default selects accuracy when available, then reliability, then gain.
#' @param reference_cost Optional two-element vector `c(minimum, maximum)` that
#'   defines the normalisation and the top axis. By default the observed range
#'   is used.
#' @param col_dominated Colour for dominated designs.
#' @param label_designs Logical. Label points when a `design` column is present.
#' @param ... Passed to [graphics::plot()].
#'
#' @return The frontier data frame ordered by real cost, invisibly. Added
#'   columns record `cost_minimisation_percent`,
#'   `performance_maximisation_percent`, `frontier_eligible`, and
#'   `pareto_plotted`.
#' @seealso [pareto_designs()], [benchmark_met_designs()].
#' @examples
#' fr <- data.frame(
#'   design = c("Full", "M4", "M3"),
#'   accuracy = c(0.82, 0.79, 0.72),
#'   cost = c(10, 6, 1)
#' )
#' fr$pareto_accuracy <- c(TRUE, TRUE, TRUE)
#' plot_pareto_frontier(fr, metric = "accuracy")
#' @export
plot_pareto_frontier <- function(
    pf, x = c("cost", "plots"), highlight_pareto = TRUE,
    col_gain = "#0C4C80", col_acc = "#087E8B",
    metric = NULL, reference_cost = NULL,
    col_dominated = "#B7BDC6", label_designs = FALSE, ...) {
  x <- match.arg(x)
  fr <- if (is.data.frame(pf)) {
    pf
  } else if (!is.null(pf$frontier)) {
    pf$frontier
  } else if (!is.null(pf$design_summary)) {
    pf$design_summary
  } else {
    stop("`pf` must be a frontier data frame or a supported result list.")
  }
  if (!x %in% names(fr))
    stop("`pf` does not contain the requested resource column `", x, "`.")
  if (is.null(metric)) {
    available <- intersect(c("accuracy", "reliability", "gain"), names(fr))
    if (!length(available))
      stop("`pf` must contain accuracy, reliability, or gain.")
    metric <- available[1L]
  } else {
    metric <- match.arg(metric, c("accuracy", "gain", "reliability"))
  }
  if (!metric %in% names(fr))
    stop("`pf` does not contain the requested endpoint `", metric, "`.")
  raw_cost <- as.numeric(fr[[x]])
  performance <- as.numeric(fr[[metric]])
  if (any(!is.finite(raw_cost)) || any(raw_cost < 0) ||
      any(!is.finite(performance)))
    stop("Resource and performance columns must contain finite values.")

  if (is.null(reference_cost)) {
    limits <- range(raw_cost)
  } else {
    if (!is.numeric(reference_cost) || length(reference_cost) != 2L ||
        any(!is.finite(reference_cost)) ||
        reference_cost[1L] < 0 ||
        reference_cost[2L] <= reference_cost[1L])
      stop("`reference_cost` must be `c(minimum, maximum)` with maximum > minimum.")
    limits <- as.numeric(reference_cost)
    if (any(raw_cost < limits[1L] - 1e-8 |
            raw_cost > limits[2L] + 1e-8))
      stop("All observed costs must lie within `reference_cost`.")
  }
  cost_span <- diff(limits)
  cost_min <- if (cost_span > 0)
    100 * (limits[2L] - raw_cost) / cost_span else rep(0, nrow(fr))

  performance_percent <- switch(
    metric,
    accuracy = 100 * performance,
    reliability = 100 * performance,
    gain = {
      mx <- max(performance)
      if (mx <= 0)
        stop("Gain must have a positive maximum for percentage scaling.")
      100 * performance / mx
    }
  )
  if (any(performance_percent < -1e-8 |
          performance_percent > 100 + 1e-8))
    stop("The selected endpoint cannot be represented on a 0--100% scale.")
  performance_percent <- pmin(100, pmax(0, performance_percent))

  pareto_col <- paste0("pareto_", metric)
  pareto <- if (pareto_col %in% names(fr)) {
    as.logical(fr[[pareto_col]])
  } else if ("pareto" %in% names(fr)) {
    as.logical(fr$pareto)
  } else {
    .pareto_flag(performance, raw_cost)
  }
  if (anyNA(pareto))
    stop("Pareto indicators must not contain missing values.")
  feasible <- if ("feasible" %in% names(fr)) {
    if (anyNA(fr$feasible))
      stop("Feasibility indicators must not contain missing values.")
    as.logical(fr$feasible)
  } else {
    rep(TRUE, nrow(fr))
  }
  pareto <- pareto & feasible

  fr$cost_minimisation_percent <- cost_min
  fr$performance_maximisation_percent <- performance_percent
  fr$frontier_eligible <- feasible
  fr$pareto_plotted <- pareto
  ordering <- order(fr[[x]])
  fr <- fr[ordering, , drop = FALSE]
  cost_min <- fr$cost_minimisation_percent
  performance_percent <- fr$performance_maximisation_percent
  pareto <- pareto[ordering]
  feasible <- feasible[ordering]
  point_colour <- ifelse(pareto, if (metric == "gain") col_gain else col_acc,
                         col_dominated)

  op <- graphics::par(mar = c(4.6, 4.7, 4.2, 1.2))
  on.exit(graphics::par(op))
  graphics::plot(
    cost_min, performance_percent, type = "n",
    xlim = c(0, 100), ylim = c(0, 100),
    xlab = "Cost minimisation (%)",
    ylab = switch(
      metric,
      accuracy = "Prediction accuracy maximisation (%)",
      reliability = "Prediction reliability maximisation (%)",
      gain = "Genetic gain maximisation (% of maximum)"
    ),
    xaxt = "n", ...
  )
  bottom_ticks <- seq(0, 100, by = 20)
  graphics::axis(1, at = bottom_ticks,
                 labels = paste0(bottom_ticks, "%"))
  if (highlight_pareto && any(pareto)) {
    ii <- which(pareto)
    ii <- ii[order(cost_min[ii])]
    if (length(ii) > 1L)
      graphics::lines(
        cost_min[ii], performance_percent[ii],
        col = if (metric == "gain") col_gain else col_acc, lwd = 2
      )
  }
  graphics::points(
    cost_min[feasible], performance_percent[feasible], pch = 21,
    bg = point_colour[feasible], col = point_colour[feasible], cex = 1.25
  )
  if (any(!feasible))
    graphics::points(
      cost_min[!feasible], performance_percent[!feasible], pch = 4,
      col = col_dominated, lwd = 1.5, cex = 1.25
    )
  if (highlight_pareto && any(pareto))
    graphics::points(
      cost_min[pareto], performance_percent[pareto],
      pch = 21, bg = point_colour[pareto], col = "#FFFFFF",
      lwd = 1.2, cex = 1.45
    )
  if (isTRUE(label_designs) && "design" %in% names(fr))
    graphics::text(
      cost_min, performance_percent, labels = fr$design,
      pos = 3, cex = 0.72
    )

  top_ticks <- seq(0, 100, by = 20)
  actual_cost <- limits[2L] - top_ticks / 100 * cost_span
  cost_digits <- if (all(abs(actual_cost - round(actual_cost)) < 1e-8)) 0 else 2
  graphics::axis(
    3, at = top_ticks,
    labels = formatC(actual_cost, format = "f", digits = cost_digits)
  )
  graphics::mtext(
    if (x == "cost") "Actual cost (decreasing)" else
      "Actual physical plots (decreasing)",
    side = 3, line = 2.5
  )
  graphics::title(main = "Feasible cost--performance Pareto frontier")
  if (highlight_pareto) {
    legend_text <- c("Feasible non-dominated", "Feasible dominated")
    legend_pch <- c(21, 21)
    legend_bg <- c(if (metric == "gain") col_gain else col_acc, col_dominated)
    legend_col <- legend_bg
    if (any(!feasible)) {
      legend_text <- c(legend_text, "Infeasible comparator")
      legend_pch <- c(legend_pch, 4)
      legend_bg <- c(legend_bg, NA)
      legend_col <- c(legend_col, col_dominated)
    }
    graphics::legend(
      "bottomleft", bty = "n",
      legend = legend_text, pch = legend_pch,
      pt.bg = legend_bg, col = legend_col
    )
  }
  invisible(fr)
}
