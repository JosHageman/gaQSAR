#' Summary method for gaQSAR_dcv objects
#'
#' Produce a detailed summary of nested cross-validation results, including
#' per-fold summaries and stability metrics.
#'
#' @param object An object of class "gaQSAR_dcv" (returned from
#'   `gaDoubleCrossValidation()`).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "summary.gaQSAR_dcv" (invisibly).
#'
#' @seealso [print.gaQSAR_dcv()], [gaDoubleCrossValidation()]
#'
#' @export
summary.gaQSAR_dcv <- function(object, ...) {

  if (!inherits(object, "gaQSAR_dcv")) {
    stop("`object` must be a gaQSAR_dcv object.", call. = FALSE)
  }

  cat("\n=== Detailed Nested Cross-Validation Summary ===\n\n")

  # Overall metrics
  n <- length(object$outerPredictions)

  cat("Dataset: n =", n, "objects\n")
  cat("Outer CV:", toupper(object$outer$method))
  if (object$outer$method == "kfold") cat(" (k =", object$outer$k, ")")
  cat("\nInner CV: LOO (singleCV fitness)\n\n")

  # Outer performance metrics
  cat("Outer Performance Metrics:\n")
  cat(sprintf("  Q2: %.4f\n", object$outerQ2))
  cat(sprintf("  RMSE: %.4f\n", sqrt(mean(object$outerResiduals^2, na.rm = TRUE))))
  cat(sprintf("  MAE:  %.4f\n", mean(abs(object$outerResiduals), na.rm = TRUE)))
  cat("\n")

  # Per-fold summary table
  cat("Per-Fold Summary:\n")
  summaryTable <- object$foldSummaries[, c("foldIndex", "nPredictors", "selectedPredictors",
                                           "innerQ2", "trainingR2", "maxVif")]
  print(summaryTable, row.names = FALSE)
  cat("\n")

  # Selection frequency and sign stability
  if (length(object$selectionFrequency) > 0) {
    cat("Predictor Selection Frequency and Sign Stability:\n")
    selFreqTable <- data.frame(
      Predictor = names(object$selectionFrequency),
      Frequency = object$selectionFrequency,
      PositiveSignFreq = object$signStability[names(object$selectionFrequency)],
      stringsAsFactors = FALSE
    )
    rownames(selFreqTable) <- NULL
    print(selFreqTable, digits = 3)
    cat("\n")
  }

  # Model size distribution
  if (length(object$modelSizeDistribution) > 0) {
    cat("Model Size Distribution:\n")
    modelSizeTable <- data.frame(
      ModelSize = as.integer(names(object$modelSizeDistribution)),
      Frequency = as.integer(object$modelSizeDistribution)
    )
    print(modelSizeTable, row.names = FALSE)
    cat("\n")
  }

  # Outlier and leverage frequencies
  perObjectOutlierMean <- mean(object$perObjectOutlierFrequency, na.rm = TRUE)
  perObjectLeverageMean <- mean(object$perObjectHighLeverageFrequency, na.rm = TRUE)
  cat("Per-Object Diagnostic Frequencies (mean across all objects):\n")
  cat(sprintf("  Mean outlier frequency: %.3f\n", perObjectOutlierMean))
  cat(sprintf("  Mean high-leverage frequency: %.3f\n", perObjectLeverageMean))
  cat("\n")

  invisible(structure(object, class = c("summary.gaQSAR_dcv", "list")))
}
