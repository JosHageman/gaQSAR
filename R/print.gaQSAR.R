#' Print method for gaQSAR objects
#'
#' Print comprehensive information about a gaQSAR object returned from
#' `gaVariableSelection()`. Displays model configuration, performance metrics,
#' selected predictors, model coefficients, and data summaries.
#'
#' @param x An object of class "gaQSAR" or a list of such objects (returned
#'   from `gaVariableSelection()`).
#' @param ... Additional arguments (currently unused).
#'
#' @details Prints detailed information including:
#' - Number of selected predictors and their indices/names
#' - Training R2, adjusted R2, and LOOCV Q2 metrics
#' - External validation Q2 (if available)
#' - OLS model coefficients (intercept and slopes)
#' - Variance Inflation Factors (VIF) for each selected predictor
#' - Summary statistics for training set residuals
#' - Summary statistics for external validation residuals (if available)
#' - Information about attached plots
#'
#' @return Invisibly returns the object.
#'
#' @seealso [summary.gaQSAR()], [gaVariableSelection()], [plot.gaQSAR()]
#' 
#' @export
print.gaQSAR <- function(x, ...) {

  # Normalize input: single gaQSAR object or list of gaQSAR objects
  if (inherits(x, "gaQSAR")) {
    models <- list(x)
  } else if (is.list(x) && all(vapply(x, inherits, logical(1), "gaQSAR"))) {
    models <- x
  } else {
    stop("`x` must be a gaQSAR object or a list of gaQSAR objects.", call. = FALSE)
  }

  cat("\n")
  cat("================== gaQSAR Model Results ==================\n")

  # Iterate through models
  for (i in seq_along(models)) {
    model <- models[[i]]

    cat("\n--- Model---\n")

    # Model configuration
    cat(sprintf("Number of predictors selected: %d\n", model$numVar))
    cat(sprintf("Selected predictor indices: %s\n",
                paste(model$importantPredictors, collapse = ", ")))

    # Try to get predictor names from model coefficients
    pred_names <- names(model$model)[-1]  # exclude intercept
    if (length(pred_names) > 0 && pred_names[1] != "") {
      cat(sprintf("Selected predictor names: %s\n",
                  paste(pred_names, collapse = ", ")))
    }

    # Performance metrics
    cat("\n--- Training Set Performance ---\n")
    cat(sprintf("R2 (Training): %.4f\n", model$R2Train))
    if (!is.null(model$R2AdjTrain)) {
      cat(sprintf("Adjusted R2 (Training): %.4f\n", model$R2AdjTrain))
    }
    cat(sprintf("Q2 (LOOCV): %.4f\n", model$Q2Loocv))

    # External validation metrics if available
    if (!is.null(model$Q2Ext)) {
      cat("\n--- External Validation Performance ---\n")
      cat(sprintf("Q2 (External): %.4f\n", model$Q2Ext))
    }

    # Model coefficients
    cat("\n--- Model Coefficients (OLS) ---\n")
    coef_df <- data.frame(
      Coefficient = names(model$model),
      Value = as.numeric(model$model),
      row.names = NULL
    )
    print(coef_df, row.names = FALSE)

    # VIF values if available
    if (!is.null(model$VIF) && length(model$VIF) > 0) {
      cat("\n--- Variance Inflation Factors (VIF) ---\n")
      vif_df <- data.frame(
        Predictor = names(model$VIF),
        VIF = as.numeric(model$VIF),
        row.names = NULL
      )
      print(vif_df, row.names = FALSE)
    }

    # Training residuals summary
    if (!is.null(model$yTrain)) {
      cat("\n--- Training Set Residuals ---\n")
      cat(sprintf("Mean: %.4f\n", mean(model$yTrain$residual)))
      cat(sprintf("Std Dev: %.4f\n", stats::sd(model$yTrain$residual)))
      cat(sprintf("Min: %.4f\n", min(model$yTrain$residual)))
      cat(sprintf("Max: %.4f\n", max(model$yTrain$residual)))
    }

    # External validation residuals summary if available
    if (!is.null(model$yExt)) {
      cat("\n--- External Validation Residuals ---\n")
      cat(sprintf("Mean: %.4f\n", mean(model$yExt$residual)))
      cat(sprintf("Std Dev: %.4f\n", stats::sd(model$yExt$residual)))
      cat(sprintf("Min: %.4f\n", min(model$yExt$residual)))
      cat(sprintf("Max: %.4f\n", max(model$yExt$residual)))
    }

    # Williams plot data if available
    if (!is.null(model$williamsData)) {
      cat("\n--- Williams Plot Data Available ---\n")
      cat(sprintf("Training samples: %d\n",
                  sum(model$williamsData$type == "Training")))
      cat(sprintf("Validation samples: %d\n",
                  sum(model$williamsData$type == "Validation")))
    }

    # Plots information
    cat("\n--- Attached Plots ---\n")
    plot_names <- names(model)[grepl("Plot$", names(model))]
    if (length(plot_names) > 0) {
      for (pname in plot_names) {
        cat(sprintf("  - %s\n", pname))
      }
    } else {
      cat("  No plots attached. Run createWilliamsPlot().\n")
    }
  }

  cat("\n=========================================================\n\n")
  invisible(x)
}
