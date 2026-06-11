#' Williams plot for double cross-validation diagnostics
#'
#' Build a Williams plot for a `gaQSAR_dcv` run using diagnostics collected
#' across training folds. Optionally aggregate per object (mean/median) or
#' show raw points (no aggregation).
#'
#' @param dcvResult A `gaQSAR_dcv` object produced by [gaDoubleCrossValidation()].
#' @param residualThreshold Numeric threshold for standardized residual bands.
#'   Default is 2.5.
#' @param aggregation Character scalar: "mean", "median", or "none" (raw points).
#'   Default is "mean".
#' @param colorBy Character scalar controlling color mapping: "objectId",
#'   "fold", or "none". "object" is accepted as a backward-compatible alias
#'   for "objectId". Default is "objectId".
#' @param label Optional character string appended to the plot title.
#' @param labelOutliers Character scalar controlling outlier labeling: "rowNumber",
#'   "rowName", or "none". Outliers are defined as points with high leverage
#'   (> h*) and/or high residual (> residualThreshold). Default is "rowNumber".
#'
#' @return A `ggplot2` object.
#'
#' @export
createDCVWilliamsPlot <- function(dcvResult,
                                  residualThreshold = 2.5,
                                  aggregation = c("mean", "median", "none"),
                                  colorBy = c("objectId", "fold", "none", "object"),
                                  label = "",
                                  labelOutliers = c("rowNumber", "rowName", "none")) {
  if (!inherits(dcvResult, "gaQSAR_dcv")) {
    stop("`dcvResult` must be a gaQSAR_dcv object.", call. = FALSE)
  }

  aggregation <- match.arg(aggregation)
  colorBy <- match.arg(colorBy)
  labelOutliers <- match.arg(labelOutliers)
  if (identical(colorBy, "object")) {
    colorBy <- "objectId"
  }

  foldDiagnostics <- dcvResult$foldDiagnostics
  if (is.null(foldDiagnostics) || length(foldDiagnostics) == 0) {
    stop("`dcvResult` does not contain fold diagnostics.", call. = FALSE)
  }

  # Collect inner (training) diagnostics across folds
  innerRaw <- do.call(rbind, lapply(seq_along(foldDiagnostics), function(i) {
    diag <- foldDiagnostics[[i]]
    if (is.null(diag) || is.null(diag$williamsData)) return(NULL)
    candidate <- diag$williamsData
    if (is.null(candidate)) return(NULL)
    cbind(candidate, fold = i, stringsAsFactors = FALSE)
  }))

  rawAll <- innerRaw

  if (is.null(rawAll) || nrow(rawAll) == 0) {
    stop("No diagnostic data available to plot.", call. = FALSE)
  }

  rawAll <- rawAll[stats::complete.cases(rawAll[, c("objectIndex", "leverage", "standardizedResidual", "fold")]), , drop = FALSE]

  aggFun <- switch(aggregation,
                   mean = function(z) mean(z, na.rm = TRUE),
                   median = function(z) stats::median(z, na.rm = TRUE),
                   none = NULL)

  plottedDf <- if (aggregation == "none") {
    rawAll
  } else {
    stats::aggregate(
      cbind(leverage, standardizedResidual) ~ objectIndex,
      data = rawAll,
      FUN = aggFun
    )
  }

  if (is.null(plottedDf) || nrow(plottedDf) == 0) {
    stop("No diagnostic data available to plot after processing.", call. = FALSE)
  }

  # Add objectId labels
  obsLabels <- names(dcvResult$yObserved)
  maxIndex <- suppressWarnings(max(plottedDf$objectIndex, na.rm = TRUE))

  plottedDf$objectId <- if (!is.null(obsLabels) && is.finite(maxIndex) && length(obsLabels) >= maxIndex) {
    obsLabels[plottedDf$objectIndex]
  } else {
    as.character(plottedDf$objectIndex)
  }

  # Ensure fold is present for coloring/faceting by fold in raw mode
  if (!("fold" %in% names(plottedDf))) {
    plottedDf$fold <- NA_integer_
  }

  pointsDf <- data.frame(
    objectId = plottedDf$objectId,
    leverage = plottedDf$leverage,
    residual = plottedDf$standardizedResidual,
    fold = plottedDf$fold,
    stringsAsFactors = FALSE
  )
  pointsDf <- pointsDf[stats::complete.cases(pointsDf[, c("leverage", "residual")]), , drop = FALSE]

  # Leverage cutoff line (h*) from folds, if available
  hStarVals <- vapply(foldDiagnostics, function(diag) {
    if (is.null(diag) || is.null(diag$hStar)) return(NA_real_)
    diag$hStar
  }, numeric(1))
  hStarLine <- stats::median(hStarVals, na.rm = TRUE)

  plotTitle <- paste0(
    "Williams plot (double CV)",
    if (nzchar(label)) paste0(" - ", label)
  )

  # Base plot with requested coloring
  if (colorBy == "none") {
    p <- ggplot2::ggplot(pointsDf, ggplot2::aes(x = leverage, y = residual)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85)
  } else if (colorBy == "fold") {
    pointsDf$fold <- factor(pointsDf$fold)
    p <- ggplot2::ggplot(pointsDf, ggplot2::aes(x = leverage, y = residual, color = fold)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::labs(color = "Fold")
  } else {
    p <- ggplot2::ggplot(pointsDf, ggplot2::aes(x = leverage, y = residual, color = objectId)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::labs(color = "Object")
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = c(-residualThreshold, residualThreshold),
      linetype = "dashed",
      color = "red",
      linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      title = plotTitle,
      x = "Leverage (h)",
      y = "Standardized residuals"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  if (is.finite(hStarLine)) {
    p <- p + ggplot2::geom_vline(
      xintercept = hStarLine,
      linetype = "dashed",
      color = "blue",
      linewidth = 0.7
    )
  }

# Add labels for outliers if requested
  if (labelOutliers != "none") {
    # Identify outliers: high leverage and/or high residual
    outliersDf <- pointsDf[
      (is.finite(hStarLine) & pointsDf$leverage > hStarLine) |
      abs(pointsDf$residual) > residualThreshold,
      , drop = FALSE
    ]

    if (nrow(outliersDf) > 0) {
      # Determine label text
      outliersDf$labelText <- if (labelOutliers == "rowName") {
        outliersDf$objectId
      } else {
        # rowNumber: extract numeric index from objectId if possible, or use row position
        as.character(match(outliersDf$objectId, pointsDf$objectId))
      }

      p <- p + ggplot2::geom_text(
        data = outliersDf,
        ggplot2::aes(x = leverage, y = residual, label = labelText),
        hjust = -0.2, vjust = -0.2, size = 3, inherit.aes = FALSE
      )
    }
  }

  
  p
}
