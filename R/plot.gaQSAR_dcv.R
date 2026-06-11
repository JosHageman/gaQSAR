#' Plot method for gaQSAR_dcv objects
#'
#' Generate diagnostic plots for nested cross-validation results using ggplot2.
#'
#' @param x An object of class "gaQSAR_dcv" (returned from `gaDoubleCrossValidation()`).
#' @param type Character; plot type. Options:
#'   - "outerPredObs" (default): Observed vs predicted values with identity line.
#'   - "selectionFrequency": Bar plot of predictor selection frequencies.
#'   - "williamsFrequency": Per-object outlier and high-leverage frequencies.
#'   - "all": Generate all three plot types sequentially.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object (or invisible NULL for type = "all").
#'
#' @seealso [print.gaQSAR_dcv()], [summary.gaQSAR_dcv()], [gaDoubleCrossValidation()]
#'
#' @export
plot.gaQSAR_dcv <- function(x, type = "outerPredObs", ...) {

  if (!inherits(x, "gaQSAR_dcv")) {
    stop("`x` must be a gaQSAR_dcv object.", call. = FALSE)
  }

  type <- match.arg(type, c("outerPredObs", "selectionFrequency", "williamsFrequency", "all"))

  if (type == "all") {
    # Generate all plot types sequentially
    cat("\n=== Plot 1/3: Outer Predictions vs Observed ===")
    p1 <- plotOuterPredObs(x)
    print(p1)

    cat("\n=== Plot 2/3: Predictor Selection Frequency ===")
    p2 <- plotSelectionFrequency(x)
    print(p2)

    cat("\n=== Plot 3/3: Williams Diagnostic Frequencies ===")
    p3 <- plotWilliamsFrequency(x)
    print(p3)

    return(invisible(NULL))
  }

  p <- if (type == "outerPredObs") {
    plotOuterPredObs(x)
  } else if (type == "selectionFrequency") {
    plotSelectionFrequency(x)
  } else if (type == "williamsFrequency") {
    plotWilliamsFrequency(x)
  }

  print(p)
  invisible(p)
}


#' Plot observed vs predicted for outer predictions
#'
#' @keywords internal
#' @noRd 
plotOuterPredObs <- function(x) {
  # Use stored observed y when available; fallback to reconstruction if absent
  yObsFull <- if (!is.null(x$yObserved)) x$yObserved else x$outerResiduals + x$outerPredictions

  validIdx <- !is.na(x$outerPredictions) & !is.na(yObsFull)
  yObs <- yObsFull[validIdx]
  yPred <- x$outerPredictions[validIdx]

  # Check for sufficient data
  if (length(yObs) < 2) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "Insufficient valid predictions for plot\n(< 2 complete cases)",
                       size = 5, hjust = 0.5, vjust = 0.5)
    return(p)
  }

  # Create data frame for ggplot
  plotData <- data.frame(
    observed = yObs,
    predicted = yPred
  )

  # Create ggplot
  p <- ggplot2::ggplot(plotData, ggplot2::aes(x = observed, y = predicted)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "Nested CV: Observed vs Outer Predictions",
      x = "Observed",
      y = "Predicted"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  p
}


#' Plot predictor selection frequency
#'
#' @keywords internal
#' @noRd
plotSelectionFrequency <- function(x) {
  if (length(x$selectionFrequency) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "No predictor selection frequency data available\n(all folds may have failed)",
                       size = 5, hjust = 0.5, vjust = 0.5)
    return(p)
  }

  # Sort by frequency (descending)
  sortedFreq <- sort(x$selectionFrequency, decreasing = TRUE)

  # Use all predictors
  topFreq <- sortedFreq
  nTop <- length(topFreq)

  posShare <- x$signStability[names(topFreq)]
  posShare[is.na(posShare)] <- 1
  posShare <- pmin(pmax(posShare, 0), 1)

  positiveFreq <- as.numeric(topFreq * posShare)
  negativeFreq <- as.numeric(topFreq * (1 - posShare))

  stackedData <- data.frame(
    predictor = factor(rep(names(topFreq), times = 2), levels = names(topFreq)),
    part = rep(c("Positive", "Negative"), each = nTop),
    frequency = c(positiveFreq, negativeFreq)
  )
  stackedData <- stackedData[stackedData$frequency > 0, , drop = FALSE]

  p <- ggplot2::ggplot(stackedData, ggplot2::aes(x = predictor, y = frequency, fill = part)) +
    ggplot2::geom_col(color = "black", alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.5) +
    ggplot2::labs(
      title = "Predictor Selection Frequency (Nested CV)",
      x = "",
      y = "Selection Frequency",
      fill = "Coefficient sign"
    ) +
    ggplot2::scale_fill_manual(values = c(Positive = "forestgreen", Negative = "firebrick")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    ggplot2::ylim(0, 1)

  p
}


#' Plot per-object outlier and high-leverage frequencies
#'
#' @keywords internal
#' @noRd
plotWilliamsFrequency <- function(x) {
  # Check if diagnostic data is available
  if (is.null(x$perObjectOutlierFrequency) || is.null(x$perObjectHighLeverageFrequency)) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "No Williams diagnostic data available\n(all folds may have failed)",
                       size = 5, hjust = 0.5, vjust = 0.5)
    return(p)
  }

  # Create data frame with both frequency metrics
  n <- length(x$perObjectOutlierFrequency)
  plotData <- data.frame(
    objectIndex = seq_len(n),
    outlierFrequency = x$perObjectOutlierFrequency,
    leverageFrequency = x$perObjectHighLeverageFrequency
  )

  # Reshape to long format for faceting
  plotDataLong <- data.frame(
    objectIndex = c(plotData$objectIndex, plotData$objectIndex),
    frequency = c(plotData$outlierFrequency, plotData$leverageFrequency),
    type = c(rep("Outlier", n), rep("High-Leverage", n))
  )

  # Create ggplot with facets
  p <- ggplot2::ggplot(plotDataLong, ggplot2::aes(x = objectIndex, y = frequency)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    ggplot2::facet_wrap(~type, nrow = 2, scales = "free_y") +
    ggplot2::labs(
      title = "Per-Object Diagnostic Frequencies (Nested CV)",
      x = "Object Index",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::ylim(0, 1)

  p
}
