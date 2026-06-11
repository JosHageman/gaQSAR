#' Y-scrambling permutation test for GA-based variable selection
#'
#' Performs a Y-scrambling permutation test for objects produced by gaQSAR
#' or gaQSAR_dcv. In each permutation, the response vector is randomly
#' permuted while the descriptor matrix is kept unchanged. The GA variable
#' selection procedure is then repeated using the same settings as in the
#' original analysis.
#'
#' @param object A `gaQSAR` or `gaQSAR_dcv` object representing the true (unpermuted)
#'   result.
#' @param x A data.frame or matrix of predictors (same as used in the original analysis).
#' @param nPermutations Integer; number of permutations to perform. Default is 500.
#' @param seed Integer; random seed for reproducibility of permutations. If `NULL`,
#'   permutations are non-deterministic.
#' @param validateSettings Logical; if `TRUE`, re-runs the true (unpermuted) analysis
#'   using the stored settings to verify they reproduce the original result. If the
#'   re-run metric differs from the stored metric by more than 0.001, a warning is issued.
#'   Default is `FALSE`.
#' @param verbose Logical; if `TRUE`, print progress information for each permutation.
#'   Default is `FALSE`.
#' @param workers Integer; number of parallel workers to use for permutations. Use `1`
#'   to run sequentially. When `workers > 1`, the function uses a temporary
#'   `future::multisession` plan and `future.apply::future_lapply()` with
#'   `future.seed = TRUE`, restoring the previous plan on exit.
#' @param ... Currently unused.
#'
#' @details
#' The permutation test works by:
#'   1. Extracting the true performance metric (Q2 for gaQSAR, outer Q2 for gaQSAR_dcv).
#'   2. For each permutation: scrambling y, re-running the analysis, and recording the 
#'      permuted performance metric.
#'   3. The empirical p-value is computed with a plus-one correction:
#'      `(1 + nGreaterEqual) / (1 + nValid)`, where `nGreaterEqual` is the number of 
#'      valid permutations with a performance metric greater than or equal to the 
#'      observed metric, and `nValid` is the number of valid permutation runs.
#'
#' @return A list of class "gaQSAR_permTest" containing:
#'   - trueMetric: Numeric; the true performance metric (Q2 or outer Q2).
#'   - permutedMetrics: Numeric vector of length `nPermutations` with performance
#'     metrics from permuted runs.
#'   - pValue: Numeric; empirical p-value with plus-one correction.
#'   - nPermutations: Integer; number of permutations performed.
#'   - objectType: Character; either "gaQSAR" or "gaQSAR_dcv".
#'   - seed: Integer or NULL; random seed used.
#'
#' @seealso [gaVariableSelection()], [gaDoubleCrossValidation()]
#'
#' @examples
#' \dontrun{
#'   # Simple example with gaQSAR object
#'   set.seed(42)
#'   n <- 50
#'   p <- 20
#'   x <- matrix(rnorm(n * p), nrow = n)
#'   colnames(x) <- paste0("X", seq_len(p))
#'   y <- 2 * x[, 2] - 1.5 * x[, 5] + rnorm(n, sd = 0.5)
#'
#'   fit <- gaVariableSelection(
#'     x = x, y = y,
#'     numberOfVariables = 3,
#'     popSize = 30,
#'     maxIter = 100,
#'     seeds = 1:3,
#'     verbose = TRUE
#'   )
#'
#'   permTest <- gaPermutationTest(
#'     object = fit,
#'     x = x,
#'     nPermutations = 100,
#'     verbose = TRUE
#'   )
#'   print(permTest)
#'
#'   # Simple example with gaQSAR_dcv object
#'   dcvFit <- gaDoubleCrossValidation(
#'     x = x, y = y,
#'     outerMethod = "kfold",
#'     outerK = 5,
#'     numberOfVariables = 3,
#'     popSize = 30,
#'     maxIter = 100,
#'     seed = 1,
#'     verbose = TRUE 
#'   )
#'
#'   permTestDCV <- gaPermutationTest(
#'     object = dcvFit,
#'     x = x,
#'     nPermutations = 50,
#'     verbose = TRUE
#'   )
#'   print(permTestDCV)
#' }
#'
#' @export
gaPermutationTest <- function(object, x, nPermutations = 100, seed = NULL, 
                              validateSettings = FALSE, verbose = FALSE, 
                              workers = 1, ...) {

  if (!inherits(object, c("gaQSAR", "gaQSAR_dcv"))) {
    stop("`object` must be a gaQSAR or gaQSAR_dcv object.", call. = FALSE)
  }

  if (length(workers) != 1 || is.na(workers) || !is.numeric(workers)) {
    stop("`workers` must be a single numeric value.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (workers < 1) {
    stop("`workers` must be >= 1.", call. = FALSE)
  }

  objectType <- if (inherits(object, "gaQSAR_dcv")) "gaQSAR_dcv" else "gaQSAR"

  # Extract y and true metric based on object type
  if (objectType == "gaQSAR") {
    y <- object$yTrain$y
    trueMetric <- object$Q2Loocv
    metricName <- "Q2"
  } else {
    y <- object$yObserved
    trueMetric <- object$outerQ2
    metricName <- "Outer Q2"
  }

  if (is.null(y) || length(y) == 0) {
    stop("Could not extract response variable from object.", call. = FALSE)
  }

  if (is.na(trueMetric) || !is.finite(trueMetric)) {
    stop(sprintf("True %s is NA or non-finite. Cannot perform permutation test.", metricName),
         call. = FALSE)
  }

  # Extract GA settings based on object type
  if (objectType == "gaQSAR") {
    # Use stored settings from gaQSAR object
    if (is.null(object$gaSettings)) {
      stop("gaQSAR object does not contain gaSettings. Please re-run gaVariableSelection to store settings.", call. = FALSE)
    }
    gaSettings <- object$gaSettings
  } else {
    # For gaQSAR_dcv: extract settings from stored call
    gaSettings <- as.list(object$call)[-1]
  }

  # Print GA settings
  if (verbose) {
    message(sprintf("=== Y-Scrambling Permutation Test ==="))
    message(sprintf("Object type: %s", objectType))
    message(sprintf("True %s: %.4f", metricName, trueMetric))
    message(sprintf("Number of permutations: %d", nPermutations))
    message("")
    message("GA settings to be used:")
    
    if (objectType == "gaQSAR") {
      crossoverFuncName <- if (is.function(gaSettings$crossoverFunc)) "<custom function>" else as.character(gaSettings$crossoverFunc)
      message(sprintf("  Number of variables: %s", gaSettings$numberOfVariables))
      message(sprintf("  Population size: %s", gaSettings$popSize))
      message(sprintf("  Max iterations: %s", gaSettings$maxIter))
      message(sprintf("  Seeds: %s", paste(gaSettings$seeds, collapse = ", ")))
      message(sprintf("  Mutation rate: %s", gaSettings$pMutation))
      message(sprintf("  Crossover rate: %s", gaSettings$pCrossover))
      message(sprintf("  Crossover function: %s", crossoverFuncName))
      message(sprintf("  Elitism: %s", gaSettings$elitism))
    } else {
      message(sprintf("  Number of variables: %s", 
                     if ("numberOfVariables" %in% names(gaSettings)) gaSettings$numberOfVariables else "4 (default)"))
      message(sprintf("  Outer method: %s", 
                     if ("outerMethod" %in% names(gaSettings)) gaSettings$outerMethod else "loo (default)"))
      if ("outerK" %in% names(gaSettings)) {
        message(sprintf("  Outer K: %s", gaSettings$outerK))
      }
      message(sprintf("  Population size: %s", 
                     if ("popSize" %in% names(gaSettings)) gaSettings$popSize else "100 (default)"))
      message(sprintf("  Max iterations: %s", 
                     if ("maxIter" %in% names(gaSettings)) gaSettings$maxIter else "1300 (default)"))
      message(sprintf("  Seeds: %s", 
                     if ("seeds" %in% names(gaSettings)) paste(gaSettings$seeds, collapse = ", ") else "1:5 (default)"))
      message(sprintf("  Mutation rate: %s", 
                     if ("pMutation" %in% names(gaSettings)) gaSettings$pMutation else "0.2 (default)"))
      message(sprintf("  Crossover rate: %s", 
                     if ("pCrossover" %in% names(gaSettings)) gaSettings$pCrossover else "0.7 (default)"))
    }
    message("")
  }

  # Validate settings by re-calculating Q2 using the `singleCV` fitness function
  if (validateSettings) {
    if (verbose) {
      if (objectType == "gaQSAR") {
        message("Validating stored settings by re-computing LOOCV Q2 with singleCV (no full re-run)...")
      } else {
        message("Validating stored settings by recomputing Outer Q2 from saved fold models (no full re-run)...")
      }
    }

    # determine vifThreshold from stored settings if present
    vifThreshold <- if ("vifThreshold" %in% names(gaSettings)) gaSettings$vifThreshold else 5

    rerunMetric <- NA_real_

    if (objectType == "gaQSAR") {
      # For gaQSAR: compute LOOCV Q2 using the selected predictors and singleCV
      if (!is.null(object$importantPredictors) && length(object$importantPredictors) > 0) {
        rerunMetric <- tryCatch(
          singleCV(predictors = object$importantPredictors, x = x, y = y, vifThreshold = vifThreshold),
          error = function(e) {
            warning(sprintf("singleCV validation failed for gaQSAR: %s", conditionMessage(e)), call. = FALSE)
            NA_real_
          }
        )
      } else {
        warning("gaQSAR object does not contain selected predictors to validate.", call. = FALSE)
      }
    } else {
      # For gaQSAR_dcv: recompute the stored Outer Q2 from the saved fold models
      # (this mirrors how `gaDoubleCrossValidation()` computes `outerQ2` and
      # avoids comparing the stored outer Q2 to mean(inner Q2), which are
      # different quantities).
      if (!is.null(object$foldModels) && !is.null(object$outer) && !is.null(object$outer$folds)) {
        nObs <- length(object$yObserved)
        recomputedPreds <- rep(NA_real_, nObs)
        foldIdxs <- seq_along(object$outer$folds)

        for (ii in foldIdxs) {
          fm <- object$foldModels[[ii]]
          folds_i <- object$outer$folds[[ii]]
          testIdx <- folds_i$test
          trainIdx <- folds_i$train

          if (is.null(fm) || is.null(fm$importantPredictors) || length(fm$importantPredictors) == 0) {
            next
          }

          preds <- tryCatch({
            sel <- fm$importantPredictors
            trainModel <- stats::lm(y[trainIdx] ~ ., data = as.data.frame(x[trainIdx, sel, drop = FALSE]))
            stats::predict(trainModel, newdata = as.data.frame(x[testIdx, sel, drop = FALSE]))
          }, error = function(e) {
            if (verbose) warning(sprintf("Fold %d prediction re-computation failed: %s", ii, conditionMessage(e)), call. = FALSE)
            rep(NA_real_, length(testIdx))
          })

          if (length(preds) == length(testIdx)) {
            recomputedPreds[testIdx] <- preds
          }
        }

        validIdx <- !is.na(recomputedPreds) & !is.na(object$yObserved)
        if (sum(validIdx) >= 2) {
          # Use the same Q2() helper as the original DCV implementation
          rerunMetric <- Q2(object$yObserved[validIdx], recomputedPreds[validIdx])
        } else {
          warning("Could not recompute outer predictions for validation; insufficient valid predictions.", call. = FALSE)
        }
      } else {
        warning("gaQSAR_dcv object lacks fold models/fold indices for validation.", call. = FALSE)
      }
    }

    if (is.finite(rerunMetric)) {
      metricDiff <- abs(rerunMetric - trueMetric)
      # label the recomputed metric depending on object type
      recomputedLabel <- if (objectType == "gaQSAR") "singleCV-based recalculated value" else "recomputed outer Q2"

      if (verbose) {
        message(sprintf("  Original %s: %.4f", metricName, trueMetric))
        message(sprintf("  %s: %.4f", recomputedLabel, rerunMetric))
        message(sprintf("  Difference: %.4f", metricDiff))
      }

      if (metricDiff > 0.001) {
        warning(sprintf("%s (%.4f) differs from stored %s (%.4f) by %.4f. Settings/metrics may not fully reproduce the original analysis.",
                       recomputedLabel, rerunMetric, metricName, trueMetric, metricDiff), call. = FALSE)
      } else if (verbose) {
        okLabel <- if (objectType == "gaQSAR") "singleCV" else "outer Q2 reconstruction"
        message(sprintf("  Settings validation (%s): PASSED\n", okLabel))
      }
    }
  }

  # init random seed for reproducibility of permutations
  if (!is.null(seed)) set.seed(seed)

  # Pre-generate unique seeds for each permutation
  permSeeds <- sample.int(.Machine$integer.max, nPermutations)

  runPermutation <- function(i) {
    set.seed(permSeeds[i])
    if (verbose && workers <= 1) {
      cat(sprintf("Permutation %d/%d...\n", i, nPermutations))
    }

    # Scramble y
    yPerm <- sample(y)

    # Create a copy of gaSettings and modify seeds for this permutation
    # to ensure different GA runs for each permutation
    permSettings <- gaSettings
    if (objectType == "gaQSAR") {
      # For gaQSAR: use only the best seed for all permutations (no variation)
      if (!is.null(object$bestSeed)) {
        permSettings$seeds <- c(object$bestSeed)  # Ensure vector, not scalar
      }
      # If bestSeed is not available, fall back to generating new seeds
      else if ("seeds" %in% names(permSettings)) {
        nSeeds <- length(permSettings$seeds)
        permSettings$seeds <- sample.int(.Machine$integer.max, nSeeds)
      }
    } else {
      # For gaQSAR_dcv: keep outer seed fixed (for reproducible fold splits)
      # but generate new GA seeds for robustness variation
      if (!is.null(object$outer$seed)) {
        permSettings$seed <- object$outer$seed
      }
      if ("seeds" %in% names(permSettings)) {
        nSeeds <- length(permSettings$seeds)
        permSettings$seeds <- sample.int(.Machine$integer.max, nSeeds)
      }
    }

    # Run analysis with scrambled y
    permResult <- tryCatch(
      {
        if (objectType == "gaQSAR") {
          # For gaQSAR: use stored settings
          do.call(gaVariableSelection, c(list(x = x, y = yPerm, verbose = FALSE), permSettings))
        } else {
          # For gaQSAR_dcv: use settings from stored call
          permSettings$x <- x
          permSettings$y <- yPerm
          permSettings$verbose <- FALSE
          do.call(gaDoubleCrossValidation, permSettings)
        }
      },
      error = function(e) {
        if (verbose) {
          warning(sprintf("Permutation %d failed: %s", i, conditionMessage(e)), call. = FALSE)
        }
        NULL
      }
    )

    # Extract metric from permuted result
    if (!is.null(permResult)) {
      if (objectType == "gaQSAR") {
        permResult$Q2Loocv
      } else {
        permResult$outerQ2
      }
    } else {
      NA_real_
    }
  }

  # Run permutations (sequential or parallel)
  if (workers <= 1) {
    permutedMetrics <- vapply(seq_len(nPermutations), runPermutation, numeric(1))
  } else {
    if (verbose) {
      message(sprintf("Running permutations in parallel using %d workers...", workers))
    }
    oldPlan <- future::plan()
    on.exit(future::plan(oldPlan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    permutedMetrics <- unlist(
      future.apply::future_lapply(
        seq_len(nPermutations),
        runPermutation,
        future.seed = TRUE
      ),
      use.names = FALSE
    )
  }

  # Remove failed permutations (NA values)
  validPermutedMetrics <- permutedMetrics[is.finite(permutedMetrics)]
  nValid <- length(validPermutedMetrics)

  if (nValid == 0) {
    stop("All permutations failed. Cannot compute p-value.", call. = FALSE)
  }

  if (nValid < nPermutations) {
    warning(sprintf("%d/%d permutations failed and were excluded.", 
                    nPermutations - nValid, nPermutations), call. = FALSE)
  }

  # Compute p-value with plus-one correction
  pValue <- (1 + sum(validPermutedMetrics >= trueMetric)) / (1 + nValid)

  if (verbose) {
    message(sprintf("\nPermutation test complete:"))
    message(sprintf("  Valid permutations: %d/%d", nValid, nPermutations))
    message(sprintf("  True %s: %.4f", metricName, trueMetric))
    message(sprintf("  Mean permuted %s: %.4f", metricName, mean(validPermutedMetrics)))
    message(sprintf("  p-value: %.4f\n", pValue))
  }

  result <- list(
    trueMetric = trueMetric,
    permutedMetrics = validPermutedMetrics,
    pValue = pValue,
    nPermutations = nValid,
    nRequested = nPermutations,
    objectType = objectType,
    metricName = metricName,
    seed = seed
  )

  class(result) <- c("gaQSAR_permTest", "list")
  return(result)
}


#' Print method for gaQSAR_permTest objects
#'
#' @param x A `gaQSAR_permTest` object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.gaQSAR_permTest <- function(x, ...) {
  if (!inherits(x, "gaQSAR_permTest")) {
    stop("`x` must be a gaQSAR_permTest object.", call. = FALSE)
  }

  cat("\n=== Y-Scrambling Permutation Test ===\n\n")
  cat(sprintf("Object type: %s\n", x$objectType))
  cat(sprintf("Metric: %s\n", x$metricName))
  cat(sprintf("Permutations: %d", x$nPermutations))
  if (x$nPermutations < x$nRequested) {
    cat(sprintf(" (%d requested, %d failed)", x$nRequested, x$nRequested - x$nPermutations))
  }
  cat("\n\n")

  cat(sprintf("True %s: %.4f\n", x$metricName, x$trueMetric))
  cat(sprintf("Mean permuted %s: %.4f\n", x$metricName, mean(x$permutedMetrics)))
  cat(sprintf("SD permuted %s: %.4f\n", x$metricName, stats::sd(x$permutedMetrics)))
  cat(sprintf("p-value: %.4f", x$pValue))
  if (x$pValue < 0.001) {
    cat(" ***")
  } else if (x$pValue < 0.01) {
    cat(" **")
  } else if (x$pValue < 0.05) {
    cat(" *")
  }
  cat("\n\n")

  invisible(x)
}


#' Summary method for gaQSAR_permTest objects
#'
#' @param object A `gaQSAR_permTest` object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
summary.gaQSAR_permTest <- function(object, ...) {
  if (!inherits(object, "gaQSAR_permTest")) {
    stop("`object` must be a gaQSAR_permTest object.", call. = FALSE)
  }

  print(object)

  cat("Permuted metric distribution:\n")
  print(summary(object$permutedMetrics))
  cat("\n")

  invisible(object)
}


#' Plot method for gaQSAR_permTest objects
#'
#' Creates a histogram of permuted metrics with a vertical line indicating
#' the true metric value.
#'
#' @param x A `gaQSAR_permTest` object.
#' @param bins Integer; number of histogram bins. Default is 30.
#' @param title Optional character plot title.
#' @param ... Additional arguments (unused).
#'
#' @return A `ggplot2` object.
#'
#' @export
plot.gaQSAR_permTest <- function(x, bins = 30, title = NULL, ...) {
  if (!inherits(x, "gaQSAR_permTest")) {
    stop("`x` must be a gaQSAR_permTest object.", call. = FALSE)
  }

  permDf <- data.frame(
    metric = x$permutedMetrics,
    stringsAsFactors = FALSE
  )

  # Create metric label with proper superscript for Q2
  metricLabel <- if (grepl("Q2|Q\\^2", x$metricName)) {
    if (grepl("Outer", x$metricName)) {
      expression(Q[outer]^2)
    } else {
      expression(Q2)
    }
  } else {
    x$metricName
  }

  plotTitle <- if (is.null(title)) {
    if (grepl("Q2|Q\\^2", x$metricName)) {
      if (grepl("Outer", x$metricName)) {
        bquote("Permutation test:" ~ Q[outer]^2 ~ "(p =" ~ .(sprintf("%.4f", x$pValue)) * ")")
      } else {
        bquote("Permutation test:" ~ Q2 ~ "(p =" ~ .(sprintf("%.4f", x$pValue)) * ")")
      }
    } else {
      sprintf("Permutation test: %s (p = %.4f)", x$metricName, x$pValue)
    }
  } else {
    title
  }

  # Ensure the x-axis includes the trueMetric
  metric_min_val <- suppressWarnings(min(permDf$metric, na.rm = TRUE))
  x_min <- min(0, metric_min_val, x$trueMetric)
  # Add a small padding so vertical line/label are visible at the edge
  x_pad <- 0.02 * (1.0 - x_min)
  x_min_plot <- x_min - x_pad
  # Position label inside plot when trueMetric equals the left limit
  label_hjust <- ifelse(x$trueMetric <= x_min_plot + 1e-12, 0, -0.1)

  p <- ggplot2::ggplot(permDf, ggplot2::aes(x = metric)) +
    ggplot2::geom_histogram(bins = bins, fill = "lightgray", color = "black", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = x$trueMetric, color = "red", linewidth = 1.2, linetype = "solid") +
    ggplot2::annotate("text", 
                     x = x$trueMetric, 
                     y = Inf, 
                     label = if (is.expression(metricLabel)) {
                       sprintf("True = %.3f", x$trueMetric)
                     } else {
                       sprintf("True %s = %.3f", x$metricName, x$trueMetric)
                     },
                     hjust = label_hjust, vjust = 1.5, color = "red", size = 4) +
    ggplot2::labs(
      title = plotTitle,
      x = metricLabel,
      y = "Frequency"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min_plot, 1.0)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4)
    )

  p
}
