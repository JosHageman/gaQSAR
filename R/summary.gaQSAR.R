#' Summary method for gaQSAR objects
#'
#' Produce a concise summary of a gaQSAR object returned from
#' `gaVariableSelection()`. Displays key model configuration and performance
#' metrics.
#'
#' @param object An object of class "gaQSAR" (returned from `gaVariableSelection()`).
#' @param ... Additional arguments (currently unused).
#'
#' @details Prints a summary including:
#' - Number of selected predictors
#' - Selected predictor indices and names (if available)
#' - Training R2, adjusted R2, and LOOCV Q2
#' - External validation Q2 (if available)
#' - Model coefficients
#' - Variance Inflation Factors (VIF) for each selected predictor
#' - Residual statistics for training and validation sets
#'
#' @return An object of class "summary.gaQSAR" (invisibly).
#'
#' @seealso [print.gaQSAR()], [gaVariableSelection()], [plot.gaQSAR()]
#'
#' @export
summary.gaQSAR <- function(object, ...) {

  if (!inherits(object, "gaQSAR")) {
    stop("`object` must be a gaQSAR object.", call. = FALSE)
  }

  models <- list(object)
  summary_list <- list()

  safe_get <- function(m, field, default = NA) {
    if (!is.null(m[[field]])) m[[field]] else default
  }

  for (i in seq_along(models)) {
    model <- models[[i]]

    if (!is.list(model)) {
      warning(sprintf("Skipping element %d: not a list/gaQSAR model", i), call. = FALSE)
      next
    }

    numVar <- safe_get(model, "numVar")
    impPred <- safe_get(model, "importantPredictors", numeric(0))
    R2Train <- safe_get(model, "R2Train")
    R2AdjTrain <- safe_get(model, "R2AdjTrain")
    Q2Loocv <- safe_get(model, "Q2Loocv")
    Q2Ext <- safe_get(model, "Q2Ext")
    coefs <- safe_get(model, "model", numeric(0))
    vif_vals <- safe_get(model, "VIF", numeric(0))
    pred_names <- if (length(coefs)) names(coefs)[-1] else character(0)

    summary_info <- list(
      numVar = numVar,
      importantPredictors = impPred,
      R2Train = R2Train,
      R2AdjTrain = R2AdjTrain,
      Q2Loocv = Q2Loocv,
      Q2Ext = Q2Ext,
      coefficients = coefs,
      VIF = vif_vals,
      pred_names = pred_names
    )

    if (!is.null(model$yTrain)) {
      summary_info$training_residuals <- list(
        mean = mean(model$yTrain$residual),
        sd = stats::sd(model$yTrain$residual),
        min = min(model$yTrain$residual),
        max = max(model$yTrain$residual)
      )
    }

    if (!is.null(model$yExt)) {
      summary_info$external_residuals <- list(
        mean = mean(model$yExt$residual),
        sd = stats::sd(model$yExt$residual),
        min = min(model$yExt$residual),
        max = max(model$yExt$residual)
      )
    }

    summary_list[[paste0("Model_", i)]] <- summary_info
  }

  class(summary_list) <- c("summary.gaQSAR", "list")
  print.summary.gaQSAR(summary_list)
  invisible(summary_list)
}

#' Print method for summary.gaQSAR objects
#'
#' @keywords internal
#' @noRd
#' @export
print.summary.gaQSAR <- function(x, ...) {

  cat("\n")
  cat("============= gaQSAR Model Summary =============\n")

  for (i in seq_along(x)) {
    model_summary <- x[[i]]

    cat("\n--- Model---\n")

    # Model configuration
    cat(sprintf("Number of predictors: %d\n", model_summary$numVar))
    cat(sprintf("Predictor indices: %s\n",
                paste(model_summary$importantPredictors, collapse = ", ")))

    if (length(model_summary$pred_names) > 0 &&
        model_summary$pred_names[1] != "") {
      cat(sprintf("Predictor names: %s\n",
                  paste(model_summary$pred_names, collapse = ", ")))
    }

    # Performance metrics
    cat("\nPerformance Metrics:\n")
    cat(sprintf("  R2 (Training): %.4f\n", model_summary$R2Train))
    if (!is.null(model_summary$R2AdjTrain)) {
      cat(sprintf("  Adjusted R2 (Training): %.4f\n", model_summary$R2AdjTrain))
    }
    cat(sprintf("  Q2 (LOOCV):    %.4f\n", model_summary$Q2Loocv))

    if (!is.null(model_summary$Q2Ext)) {
      cat(sprintf("  Q2 (External): %.4f\n", model_summary$Q2Ext))
    }

    # Model coefficients
    cat("\nModel Coefficients:\n")
    coef_names <- names(model_summary$coefficients)
    for (j in seq_along(model_summary$coefficients)) {
      cat(sprintf("  %s: %.6f\n", coef_names[j], model_summary$coefficients[j]))
    }

    # VIF values if available
    if (!is.null(model_summary$VIF) && length(model_summary$VIF) > 0) {
      cat("\nVariance Inflation Factors (VIF):\n")
      vif_names <- names(model_summary$VIF)
      for (j in seq_along(model_summary$VIF)) {
        cat(sprintf("  %s: %.4f\n", vif_names[j], model_summary$VIF[j]))
      }
    }

    # Residual statistics
    if (!is.null(model_summary$training_residuals)) {
      cat("\nTraining Residuals:\n")
      cat(sprintf("  Mean: %.4f, Std Dev: %.4f\n",
                  model_summary$training_residuals$mean,
                  model_summary$training_residuals$sd))
      cat(sprintf("  Range: [%.4f, %.4f]\n",
                  model_summary$training_residuals$min,
                  model_summary$training_residuals$max))
    }

    if (!is.null(model_summary$external_residuals)) {
      cat("\nExternal Validation Residuals:\n")
      cat(sprintf("  Mean: %.4f, Std Dev: %.4f\n",
                  model_summary$external_residuals$mean,
                  model_summary$external_residuals$sd))
      cat(sprintf("  Range: [%.4f, %.4f]\n",
                  model_summary$external_residuals$min,
                  model_summary$external_residuals$max))
    }
  }

  cat("\n==============================================\n\n")
  invisible(x)
}
