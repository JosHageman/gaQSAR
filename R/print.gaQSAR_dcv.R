#' Print method for gaQSAR_dcv objects
#'
#' Display a brief overview of nested cross-validation results.
#'
#' @param x An object of class "gaQSAR_dcv" (returned from `gaDoubleCrossValidation()`).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the object.
#'
#' @seealso [summary.gaQSAR_dcv()], [gaDoubleCrossValidation()]
#'
#' @export
print.gaQSAR_dcv <- function(x, ...) {

  if (!inherits(x, "gaQSAR_dcv")) {
    stop("`x` must be a gaQSAR_dcv object.", call. = FALSE)
  }

  cat("\n=== Nested Cross-Validation Results (gaQSAR_dcv) ===\n\n")

  # Dataset and fold information
  n <- length(x$outerPredictions)
  cat(sprintf("Dataset: n = %d objects\n", n))
  cat(sprintf("Outer CV: %s", toupper(x$outer$method)))
  if (x$outer$method == "kfold") cat(sprintf(" (k = %d)", x$outer$k))
  cat("\n")

  cat("Inner CV: LOO (singleCV fitness)\n\n")

  # Outer performance
  cat(sprintf("Outer Q2: %.4f\n", x$outerQ2))

  # Model sizes
  avgModelSize <- mean(x$foldSummaries$nPredictors)
  cat(sprintf("Average model size: %.1f predictors\n", avgModelSize))
  cat(sprintf("Model size range: [%d, %d]\n\n",
              min(x$foldSummaries$nPredictors),
              max(x$foldSummaries$nPredictors)))

  # Top selected predictors
  topN <- min(5, length(x$selectionFrequency))
  topPredictors <- utils::head(sort(x$selectionFrequency, decreasing = TRUE), topN)

  cat("Top selected predictors (frequency):\n")
  for (i in seq_len(length(topPredictors))) {
    cat(sprintf("  %d. %s (%.1f%%)\n",
                i, names(topPredictors)[i], 100 * topPredictors[i]))
  }

  cat("\n")
  invisible(x)
}
