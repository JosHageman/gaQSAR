#' Nested (double) cross-validation for GA-based variable selection
#'
#' Performs nested cross-validation combining an outer CV loop (for final model
#' assessment) with an inner CV loop (used by the GA for fitness evaluation during
#' variable selection). The inner loop always uses LOOCV via `singleCV` as the
#' GA fitness. The outer loop may be leave-one-out or k-fold (default), providing
#' an unbiased estimate of outer-loop predictive performance and stability analysis
#' of selected predictors.
#'
#' @param x A data.frame or matrix of predictors (rows = observations,
#'   columns = candidate variables).
#' @param y Numeric response vector aligned with rows of `x`.
#' @param outerMethod Character; outer CV method: "loo" for leave-one-out
#'   or "kfold" for k-fold cross-validation (default).
#' @param outerK Integer; number of folds for outer k-fold CV (only used if
#'   `outerMethod = "kfold"`). Defaults to 5.
#' @param seed Integer; random seed for reproducible fold splits (especially relevant
#'   for k-fold CV). If `NULL`, fold splits are non-deterministic.
#' @param residualThreshold Numeric; threshold for standardized residuals in Williams plot diagnostics.
#' @param verbose Logical; if `TRUE`, print brief progress information per outer fold.
#'   Defaults to `FALSE` (silent operation).
#' @param ... Additional arguments forwarded to `gaVariableSelection()`, such as
#'   `numberOfVariables`, `popSize`, `maxIter`, `pMutation`, `pCrossover`,
#'   and other GA settings.
#'
#' @details
#' For each outer fold:
#'   1. Create outerTrain (all data except outerTest) and outerTest (hold-out fold).
#'   2. Run gaVariableSelection() on outerTrain using inner CV for fitness evaluation.
#'   3. Fit the selected lm model on outerTrain using selected predictors.
#'   4. Predict on outerTest and store predictions and residuals.
#'   5. Compute and store diagnostic measures (Williams plot data, VIF, etc.).
#'
#' Two types of seeds are employed:
#'   - The `seed` parameter (outer fold seed) controls reproducible creation of outer CV folds.
#'   - The `seeds` parameter (passed via `...` to gaVariableSelection) enables multiple GA runs
#'     for robustness within each fold. The best seed is tracked in `foldSummaries$bestSeed`.
#'
#' Across all outer folds, summaries are computed:
#'   - selectionFrequency: How often each predictor is selected across successful
#'     outer folds.
#'   - selectionFrequencyAllFolds: Optional frequency over all outer folds,
#'     including failed folds.
#'   - signStability: For each selected predictor, the proportion of times its
#'     coefficient is positive (vs negative) across folds.
#'   - modelSizeDistribution: Distribution of numbers of selected predictors.
#'   - perObject outlier/high-leverage frequencies based on Williams plot thresholds.
#'
#' @return An object of class "gaQSAR_dcv" (list) containing:
#'   - call: The matched call.
#'   - outer: List with elements `method`, `k`, `seed` (the outer fold seed for reproducibility),
#'     and `folds` (list of outer fold indices).
#'   - outerPredictions: Numeric vector of length n (full dataset) with out-of-fold predictions.
#'   - outerResiduals: Numeric vector of length n with outer residuals.
#'   - yObserved: Original observed response vector (aligned with rows of `x`).
#'   - outerQ2: Numeric; Q2 computed from outer predictions (NA if < 2 valid predictions).
#'   - foldModels: List of gaQSAR objects (one per outer fold; NULL for failed folds).
#'   - foldSummaries: Data.frame with one row per outer fold, including columns:
#'     nPredictors, selectedPredictors, innerQ2 (LOOCV fitness), trainingR2,
#'     trainingR2Adj, maxVif, bestSeed (the GA seed that produced the best model for that fold),
#'     status ("ok" or "failed"), errorMessage.
#'   - foldDiagnostics: List with Williams plot diagnostic data per fold (NULL for failed folds);
#'     each element includes williamsData, hStar, residualThreshold.
#'   - adThresholds: List describing applicability domain thresholds used:
#'     leverageCutoffMethod, leverageCutoff (NA if varies), residualCutoff, notes.
#'   - selectionFrequency: Named numeric vector of predictor frequencies (0 to 1)
#'     computed over successful folds only.
#'   - selectionFrequencyAllFolds: Named numeric vector of predictor frequencies
#'     (0 to 1) computed over all folds, including failed folds.
#'   - signStability: Named numeric vector of positive sign frequency for selected predictors.
#'   - modelSizeDistribution: Named integer vector of model size frequencies.
#'   - perObjectOutlierFrequency: Numeric vector of length n with outlier frequencies.
#'   - perObjectHighLeverageFrequency: Numeric vector of length n with high-leverage frequencies.
#'
#' @seealso [gaVariableSelection()], [singleCV()], [Q2()], [createWilliamsPlot()]
#'
#' @examples
#' set.seed(42)
#' n <- 50
#' p <- 20
#' x <- matrix(rnorm(n * p), nrow = n)
#' colnames(x) <- paste0("X", seq_len(p))
#' y <- 2 * x[, 2] - 1.5 * x[, 5] + rnorm(n, sd = 0.5)
#'
#' # Nested CV with LOO outer and inner LOOCV fitness
#' # NOTE: Settings below are example-only (fast runtime, not optimal).
#' # For real QSAR work, use larger population sizes, more generations, and multiple seeds.
#' # Good starting point are the default GA settings, or run an experimental design to
#' # tune GA settings for your dataset.
#' dcvFit <- gaDoubleCrossValidation(
#'   x = x,
#'   y = y,
#'   outerMethod = "kfold",
#'   outerK = 2,
#'   numberOfVariables = 3,
#'   popSize = 20,
#'   maxIter = 20,
#'   seeds = 1,
#'   interval = 5,
#'   verbose = TRUE
#' )
#'
#' print(dcvFit)
#' summary(dcvFit)
#' plot(dcvFit)
#' plot(dcvFit, type = "selectionFrequency")
#' plot(dcvFit, type = "williamsFrequency")
#'
#' @export
gaDoubleCrossValidation <- function(x, y,
                                    outerMethod = "kfold",
                                    outerK = 5,
                                    seed = NULL,
                                    residualThreshold = 2.5,
                                    verbose = FALSE,
                                    ...) {

  # Internal helper for verbose output
  .vcat <- function(...) {
    if (verbose) cat(...)
  }

  # Argument validation
  outerMethod <- match.arg(outerMethod, c("loo", "kfold"))

  if (!is.null(seed)) set.seed(seed)

  call <- match.call()
  n <- nrow(x)

  # Ensure predictor names exist; fallback to X1...Xp if missing
  xLocal <- x
  if (is.null(colnames(xLocal)) || any(colnames(xLocal) == "")) {
    colnames(xLocal) <- paste0("X", seq_len(ncol(xLocal)))
  }
  x <- xLocal

  # Create outer folds
  outerFolds <- if (outerMethod == "loo") {
    lapply(seq_len(n), function(i) list(train = (1:n)[-i], test = i))
  } else {
    createFolds(n, outerK, seed)
  }

  nOuterFolds <- length(outerFolds)

  # Storage for results
  foldModels <- list()
  foldSummaries <- data.frame()
  foldDiagnostics <- list()
  outerPredictions <- rep(NA_real_, n)
  outerResiduals <- rep(NA_real_, n)
  allSelectedPredictors <- list()
  allCoefficients <- list()
  perObjectOutlierFlag <- matrix(FALSE, nrow = n, ncol = nOuterFolds)
  perObjectHighLeverageFlag <- matrix(FALSE, nrow = n, ncol = nOuterFolds)

  colnames(perObjectOutlierFlag) <- paste0("fold", seq_len(nOuterFolds))
  colnames(perObjectHighLeverageFlag) <- paste0("fold", seq_len(nOuterFolds))

  # Outer loop
  for (foldIdx in seq_len(nOuterFolds)) {
    outerTrainIdx <- outerFolds[[foldIdx]]$train
    outerTestIdx <- outerFolds[[foldIdx]]$test

    # Print fold start info
    if (outerMethod == "loo") {
      .vcat(sprintf("Outer fold %d/%d (method=loo) ...", foldIdx, nOuterFolds))
    } else {
      .vcat(sprintf("Outer fold %d/%d (kfold, test=%d) ...", 
                    foldIdx, nOuterFolds, length(outerTestIdx)))
    }

    xTrain <- x[outerTrainIdx, , drop = FALSE]
    yTrain <- y[outerTrainIdx]
    xTest <- x[outerTestIdx, , drop = FALSE]
    yTest <- y[outerTestIdx]

    # mark test objects as NA for diagnostics to avoid counting them as non-outliers
    perObjectOutlierFlag[outerTestIdx, foldIdx] <- NA
    perObjectHighLeverageFlag[outerTestIdx, foldIdx] <- NA

    # Wrap fold processing in tryCatch for robustness
    foldResult <- tryCatch(
      {
        # Run GA variable selection on training fold
        gaModel <- gaVariableSelection(
          x = xTrain,
          y = yTrain,
          ...
        )

        # Extract selected predictors (indices in original x)
        selectedIdx <- gaModel$importantPredictors
        if (length(selectedIdx) == 0) {
          stop("GA returned no selected predictors.")
        }

        selectedNames <- colnames(xTrain)[selectedIdx]

        # Fit model on outer training data
        xTrainSub <- xTrain[, selectedIdx, drop = FALSE]
        trainModel <- stats::lm(yTrain ~ ., data = data.frame(yTrain, xTrainSub))

        # Check for singular fit
        if (any(is.na(stats::coef(trainModel)))) {
          stop("Model fit is singular or has NA coefficients.")
        }

        # Make predictions on outer test fold
        xTestSub <- xTest[, selectedIdx, drop = FALSE]
        testPreds <- stats::predict(trainModel, newdata = data.frame(xTestSub))

        if (any(is.na(testPreds))) {
          stop("Predictions contain NA values.")
        }

        testResids <- yTest - testPreds

        # Store coefficients for sign stability
        coefs <- stats::coef(trainModel)
        names(coefs)[1] <- "Intercept"

        # Compute Williams plot data (leverages, residuals)
        leverages <- stats::hatvalues(trainModel)
        stdResids <- stats::rstandard(trainModel)
        rawResids <- stats::residuals(trainModel)

        # Compute thresholds: h* = 3p/n (p = number of coefficients including intercept)
        numberOfPredictors <- length(selectedIdx)
        hStar <- 3 * (numberOfPredictors + 1) / nrow(xTrain)

        williamsData <- data.frame(
          objectIndex = outerTrainIdx,
          leverage = leverages,
          standardizedResidual = stdResids,
          rawResidual = rawResids,
          highLeverage = leverages > hStar,
          highResidual = abs(stdResids) > residualThreshold,
          stringsAsFactors = FALSE
        )

        diagnostics <- list(
          williamsData = williamsData,
          hStar = hStar,
          residualThreshold = residualThreshold
        )

        # Flag outliers and high leverage in training fold
        for (i in seq_len(nrow(williamsData))) {
          objIdx <- williamsData$objectIndex[i]
          if (williamsData$highResidual[i]) perObjectOutlierFlag[objIdx, foldIdx] <- TRUE
          if (williamsData$highLeverage[i]) perObjectHighLeverageFlag[objIdx, foldIdx] <- TRUE
        }

        # Return successful fold result
        list(
          status = "ok",
          errorMessage = NA_character_,
          gaModel = gaModel,
          selectedNames = selectedNames,
          selectedIdx = selectedIdx,
          testPreds = testPreds,
          testResids = testResids,
          coefs = coefs,
          diagnostics = diagnostics
        )
      },
      error = function(e) {
        # Return failed fold result
        list(
          status = "failed",
          errorMessage = substr(conditionMessage(e), 1, 100),
          gaModel = NULL,
          selectedNames = character(0),
          selectedIdx = integer(0),
          testPreds = NA_real_,
          testResids = NA_real_,
          coefs = NULL,
          diagnostics = NULL
        )
      }
    )

    # Process fold result
    if (foldResult$status == "ok") {
      # Store successful fold results
      foldModels[[foldIdx]] <- foldResult$gaModel
      allSelectedPredictors[[foldIdx]] <- foldResult$selectedNames
      allCoefficients[[foldIdx]] <- foldResult$coefs
      foldDiagnostics[[foldIdx]] <- foldResult$diagnostics

      outerPredictions[outerTestIdx] <- foldResult$testPreds
      outerResiduals[outerTestIdx] <- foldResult$testResids

      # Store fold summary
      foldSummary <- data.frame(
        foldIndex = foldIdx,
        nPredictors = length(foldResult$selectedIdx),
        selectedPredictors = paste(foldResult$selectedNames, collapse = ","),
        innerQ2 = foldResult$gaModel$Q2Loocv,
        trainingR2 = foldResult$gaModel$R2Train,
        trainingR2Adj = foldResult$gaModel$R2AdjTrain,
        maxVif = if (!is.null(foldResult$gaModel$VIF) && length(foldResult$gaModel$VIF) > 0) max(foldResult$gaModel$VIF) else NA_real_,
        bestSeed = if (!is.null(foldResult$gaModel$bestSeed)) foldResult$gaModel$bestSeed else NA_integer_,
        status = "ok",
        errorMessage = NA_character_,
        stringsAsFactors = FALSE
      )

      # Print fold success info
      validSoFar <- !is.na(outerPredictions) & !is.na(y)
      outerQ2SoFar <- if (sum(validSoFar) >= 2) {
        Q2(y[validSoFar], outerPredictions[validSoFar])
      } else {
        NA_real_
      }
      .vcat(sprintf(" OK: p=%d, innerQ2=%.2f, R2adj=%.2f, outerQ2_sofar=%.2f\n",
                    length(foldResult$selectedIdx),
                    foldResult$gaModel$Q2Loocv,
                    foldResult$gaModel$R2AdjTrain,
                    outerQ2SoFar))
    } else {
      # Store failed fold results
      foldModels[[foldIdx]] <- NULL
      allSelectedPredictors[[foldIdx]] <- character(0)
      allCoefficients[[foldIdx]] <- NULL
      foldDiagnostics[[foldIdx]] <- NULL

      outerPredictions[outerTestIdx] <- NA_real_
      outerResiduals[outerTestIdx] <- NA_real_

      # Set entire fold column to NA to avoid bias in frequency calculations
      perObjectOutlierFlag[, foldIdx] <- NA
      perObjectHighLeverageFlag[, foldIdx] <- NA

      # Store fold summary with failure info
      foldSummary <- data.frame(
        foldIndex = foldIdx,
        nPredictors = 0L,
        selectedPredictors = "",
        innerQ2 = NA_real_,
        trainingR2 = NA_real_,
        trainingR2Adj = NA_real_,
        maxVif = NA_real_,
        bestSeed = NA_integer_,
        status = "failed",
        errorMessage = foldResult$errorMessage,
        stringsAsFactors = FALSE
      )

      # Print fold failure info
      errMsg <- substr(foldResult$errorMessage, 1, 80)
      .vcat(sprintf(" FAILED: %s\n", errMsg))
    }

    foldSummaries <- rbind(foldSummaries, foldSummary)
  }

  # Compute selection frequency from successful folds only
  allPredNames <- unique(unlist(allSelectedPredictors))
  if (length(allPredNames) == 0) {
    allPredNames <- colnames(x)
  }

  okFoldIdx <- foldSummaries$status == "ok"
  selectedPredictorsOk <- allSelectedPredictors[okFoldIdx]

  if (!any(okFoldIdx)) {
    selectionFrequency <- rep(NA_real_, length(allPredNames))
    names(selectionFrequency) <- allPredNames
  } else {
    selectionFrequency <- vapply(allPredNames, function(pname) {
      mean(vapply(selectedPredictorsOk, function(selected) {
        pname %in% selected
      }, logical(1)))
    }, numeric(1))
  }

  # Optional transparency metric: frequency over all folds (including failed)
  selectionFrequencyAllFolds <- vapply(allPredNames, function(pname) {
    mean(vapply(allSelectedPredictors, function(selected) pname %in% selected, logical(1)))
  }, numeric(1))

  names(selectionFrequency) <- allPredNames
  names(selectionFrequencyAllFolds) <- allPredNames

  # Compute sign stability for selected predictors
  signStabilityList <- lapply(allPredNames, function(pname) {
    coefVec <- vapply(allCoefficients, function(coef) {
      # Handle NULL or empty coefficient objects from failed folds
      if (is.null(coef) || length(coef) == 0) return(NA_real_)
      if (pname %in% names(coef)) coef[[pname]] else NA_real_
    }, numeric(1))
    validCoefs <- coefVec[!is.na(coefVec)]
    nAvailable <- length(validCoefs)
    if (nAvailable == 0) {
      return(list(positive = NA_real_, negative = NA_real_, nAvailable = 0L))
    }
    list(
      positive = mean(validCoefs > 0),
      negative = mean(validCoefs < 0),
      nAvailable = nAvailable
    )
  })
  names(signStabilityList) <- allPredNames
  
  # Extract sign stability as proportion positive (backward compatible)
  signStability <- vapply(signStabilityList, function(x) x$positive, numeric(1))
  names(signStability) <- allPredNames

  # Model size distribution
  modelSizes <- vapply(foldSummaries$nPredictors, as.integer, integer(1))
  modelSizeDistribution <- table(modelSizes)

  # Per-object outlier and high-leverage frequencies
  perObjectOutlierFrequency <- rowMeans(perObjectOutlierFlag, na.rm = TRUE)
  perObjectHighLeverageFrequency <- rowMeans(perObjectHighLeverageFlag, na.rm = TRUE)
  perObjectNAvailable <- rowSums(!is.na(perObjectOutlierFlag))

  # Compute outer Q2 (require at least 2 valid predictions)
  validIdx <- !is.na(outerPredictions) & !is.na(y)
  outerQ2 <- if (sum(validIdx) >= 2) {
    Q2(y[validIdx], outerPredictions[validIdx])
  } else {
    NA_real_
  }

  # Print final summary
  nFailed <- sum(foldSummaries$status == "failed")
  successfulFolds <- foldSummaries[foldSummaries$status == "ok", ]
  meanP <- if (nrow(successfulFolds) > 0) {
    mean(successfulFolds$nPredictors)
  } else {
    NA_real_
  }
  .vcat(sprintf("Done: outerQ2=%.2f, failedFolds=%d/%d, meanP=%.1f\n",
                outerQ2, nFailed, nOuterFolds, meanP))

  # Build applicability domain thresholds summary
  adThresholds <- list(
    leverageCutoffMethod = "3(p+1)/n",
    leverageCutoff = NA_real_,  # varies per fold depending on p
    residualCutoff = residualThreshold,
    notes = "Leverage cutoff computed per fold; varies with number of selected predictors."
  )

  # Build result object
  result <- list(
    call = call,
    outer = list(
      method = outerMethod,
      k = if (outerMethod == "kfold") outerK else NA,
      seed = seed,
      folds = outerFolds
    ),
    yObserved = y,
    outerPredictions = outerPredictions,
    outerResiduals = outerResiduals,
    outerQ2 = outerQ2,
    foldModels = foldModels,
    foldSummaries = foldSummaries,
    foldDiagnostics = foldDiagnostics,
    adThresholds = adThresholds,
    selectionFrequency = selectionFrequency,
    selectionFrequencyAllFolds = selectionFrequencyAllFolds,
    signStability = signStability,
    signStabilityList = signStabilityList,
    modelSizeDistribution = modelSizeDistribution,
    perObjectOutlierFrequency = perObjectOutlierFrequency,
    perObjectHighLeverageFrequency = perObjectHighLeverageFrequency,
    perObjectNAvailable = perObjectNAvailable
  )

  class(result) <- c("gaQSAR_dcv", "list")
  return(result)
}


#' Create fold indices for k-fold cross-validation
#'
#' Internal helper function to generate deterministic k-fold split indices.
#'
#' @param n Integer; number of observations.
#' @param k Integer; number of folds.
#' @param seed Integer; random seed for reproducibility.
#'
#' @return A list of length k, where each element is a list with:
#'   - train: Integer vector of training indices
#'   - test: Integer vector of test indices
#'
#' @keywords internal
#' @noRd
createFolds <- function(n, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Assign each observation to a fold (1 to k)
  foldAssignment <- sample(rep(seq_len(k), length.out = n))

  folds <- lapply(seq_len(k), function(i) {
    testIdx <- which(foldAssignment == i)
    trainIdx <- which(foldAssignment != i)
    list(train = trainIdx, test = testIdx)
  })

  return(folds)
}
