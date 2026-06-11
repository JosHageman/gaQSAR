#' Plot training metrics (R2, R2adj, Q2) versus model size for nested CV runs
#'
#' Summarizes training metrics across outer folds for one or multiple
#' `gaQSAR_dcv` objects and plots mean +/- SE versus the number of predictors in
#' a single panel. Metrics share axes and are distinguished by color/shape.
#'
#' SE is computed across outer folds with status "ok".
#'
#' @param dcvResults A `gaQSAR_dcv` object or a list of `gaQSAR_dcv` objects.
#'   Use a list when you ran nested CV for multiple `numberOfVariables`.
#' @param metrics Character vector of metrics to plot. Allowed values:
#'   `"R2"`, `"R2adj"`, `"Q2"` (inner). Use `includeOuterQ2 = TRUE` to add the
#'   outer-fold Q2 as a separate series.
#' @param title Optional character plot title.
#' @param includeOuterQ2 Logical; if `TRUE`, append the outer Q2 value from each
#'   `gaQSAR_dcv` as an additional metric in the plot. Default `FALSE`.
#'
#' @return A `ggplot2` object.
#'
#' @export
createDCVTrainingMetricsPlot <- function(dcvResults,
                                        metrics = c("R2", "R2adj", "Q2"),
                                        title = "",
                                        includeOuterQ2 = FALSE) {

  allowedMetrics <- c("R2", "R2adj", "Q2")
  if (length(metrics) == 0) {
    stop("`metrics` must contain at least one of: R2, R2adj, Q2.", call. = FALSE)
  }
  metrics <- intersect(allowedMetrics, metrics)
  if (length(metrics) == 0) {
    stop("No valid `metrics` selected. Allowed: R2, R2adj, Q2.", call. = FALSE)
  }

  # Allow a single gaQSAR_dcv or a list of them
  if (inherits(dcvResults, "gaQSAR_dcv")) {
    dcvList <- list(dcvResults)
  } else if (is.list(dcvResults)) {
    dcvList <- dcvResults
  } else {
    stop("`dcvResults` must be a gaQSAR_dcv object or a list of gaQSAR_dcv objects.", call. = FALSE)
  }

  if (!all(vapply(dcvList, inherits, logical(1), what = "gaQSAR_dcv"))) {
    stop("All elements of `dcvResults` must be gaQSAR_dcv objects.", call. = FALSE)
  }

  meanSe <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n == 0) return(c(mean = NA_real_, se = NA_real_, n = 0))
    se <- if (n >= 2) stats::sd(x) / sqrt(n) else NA_real_
    c(mean = mean(x), se = se, n = n)
  }

  inferModelSize <- function(fs) {
    if (!("nPredictors" %in% names(fs))) return(NA_real_)
    ok <- if ("status" %in% names(fs)) fs$status == "ok" else rep(TRUE, nrow(fs))
    vals <- fs$nPredictors[ok]
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(NA_real_)
    # If all folds use the same model size, take it; otherwise take the median
    if (length(unique(vals)) == 1) unique(vals) else stats::median(vals)
  }

  extractOne <- function(dcvObj) {
    fs <- dcvObj$foldSummaries
    if (is.null(fs) || nrow(fs) == 0) {
      return(NULL)
    }

    ok <- if ("status" %in% names(fs)) fs$status == "ok" else rep(TRUE, nrow(fs))
    nPredictors <- inferModelSize(fs)

    metricMap <- list(
      R2 = "trainingR2",
      R2adj = "trainingR2Adj",
      Q2 = "innerQ2"
    )

    out <- lapply(metrics, function(m) {
      colName <- metricMap[[m]]
      if (!(colName %in% names(fs))) {
        return(data.frame(
          nPredictors = nPredictors,
          metric = m,
          mean = NA_real_,
          se = NA_real_,
          nOk = 0,
          stringsAsFactors = FALSE
        ))
      }
      ms <- meanSe(fs[[colName]][ok])
      data.frame(
        nPredictors = nPredictors,
        metric = m,
        mean = as.numeric(ms["mean"]),
        se = as.numeric(ms["se"]),
        nOk = as.integer(ms["n"]),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, out)
  }

  plotDf <- do.call(rbind, lapply(dcvList, extractOne))
  if (isTRUE(includeOuterQ2)) {
    outerDf <- lapply(dcvList, function(dcvObj) {
      fs <- dcvObj$foldSummaries
      nPredictors <- inferModelSize(fs)
      if (!is.null(dcvObj$outerQ2) && is.finite(dcvObj$outerQ2) && !is.na(nPredictors)) {
        return(data.frame(
          nPredictors = nPredictors,
          metric = "OuterQ2",
          mean = as.numeric(dcvObj$outerQ2),
          se = NA_real_,
          nOk = NA_integer_,
          stringsAsFactors = FALSE
        ))
      }
      NULL
    })
    outerDf <- do.call(rbind, outerDf)
    if (!is.null(outerDf) && nrow(outerDf) > 0) {
      plotDf <- rbind(plotDf, outerDf)
    }
  }
  if (is.null(plotDf) || nrow(plotDf) == 0) {
    stop("No fold summary data available to plot.", call. = FALSE)
  }

  plotDf <- plotDf[order(plotDf$nPredictors, plotDf$metric), , drop = FALSE]

  # Only draw lines where there are 2+ points for that metric
  nPerMetric <- table(plotDf$metric)
  lineDf <- plotDf[plotDf$metric %in% names(nPerMetric)[nPerMetric >= 2], , drop = FALSE]

  metricLabels <- list(
    R2 = expression(R^2),
    R2adj = expression(R[adj]^2),
    Q2 = expression(Q[inner]^2),
    OuterQ2 = expression(Q[outer]^2)
  )

  p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = nPredictors, y = mean, color = metric)) +
    ggplot2::geom_line(data = lineDf, ggplot2::aes(group = metric)) +
    ggplot2::geom_point(size = 2, ggplot2::aes(shape = metric)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.15,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(plotDf$nPredictors))) +
    ggplot2::scale_color_discrete(labels = function(x) metricLabels[x]) +
    ggplot2::scale_shape_discrete(labels = function(x) metricLabels[x]) +
    ggplot2::labs(
      title = title,
      x = "Number of predictors",
      y = "Mean across outer folds (+/- SE)",
      color = "Metric",
      shape = "Metric"
    ) 

  p
}
