#' Genetic algorithm based variable selection for QSAR
#'
#' Runs a genetic algorithm (GA) to select a fixed number of predictors.
#' For the requested model size, the GA is executed for each seed in `seeds`,
#' and the best run (highest Q2) is returned.
#'
#' @param x A data.frame or matrix of predictors (rows = observations,
#'   columns = candidate variables).
#' @param y Numeric response vector aligned with rows of `x`.
#' @param evalFunc Character or function; fitness function used by the GA.
#'   The default is "singleCV", which computes LOOCV-based $Q^2$.
#' @param numberOfVariables Integer; number of predictors to select.
#' @param pMutation Numeric; mutation probability.
#' @param pCrossover Numeric; crossover probability.
#' @param crossoverFunc Character or function; crossover function used by the GA.
#' @param popSize Integer; GA population size.
#' @param maxIter Integer; maximum number of GA iterations.
#' @param elitism Integer; number of elite individuals carried over.
#' @param seeds Integer vector; RNG seeds to repeat GA runs for robustness.
#' @param interval Integer; iteration interval for progress monitoring when `verbose = TRUE`.
#' @param verbose Logical; if `TRUE`, prints progress information and GA settings.
#'
#' @details For each seed, the GA is run with the given parameters and the
#' resulting Q2 values are collected. The run with the highest Q2 is selected
#' and returned.
#'
#' @return An object of class "gaQSAR" containing the best GA run for the 
#'   specified model size. The structure typically includes: `numVar`, 
#'   `importantPredictors`, `yTrain`, `model`, `R2Train`, and `Q2Loocv`.
#'   Use `plot.gaQSAR()` to visualize any attached plots.
#'
#' @seealso [singleCV()], [createQ2Plot()], [createWilliamsPlot()]
#'
#' @export
gaVariableSelection <- function(x, y,
                                evalFunc = "singleCV",
                                numberOfVariables = 4,
                                pMutation = 0.2,
                                pCrossover = 0.7,
                                crossoverFunc = "gaintegerOnePointCrossover",
                                popSize = 100,
                                maxIter = 1300,
                                elitism = 3,
                                seeds = 1:5,
                                interval = 50,
                                verbose = FALSE) {

  resolveFunction <- function(fun, argName, envir = parent.frame()) {
    if (is.function(fun)) return(fun)

    if (!is.character(fun) || length(fun) != 1L || is.na(fun)) {
      stop(sprintf("%s must be a function or a single function name (character).", argName),
           call. = FALSE)
    }

    # First: look in calling environment / package namespace chain
    if (exists(fun, mode = "function", inherits = TRUE)) {
      return(get(fun, mode = "function", inherits = TRUE))
    }

    stop(sprintf("Unknown %s: '%s'. Function not found on search path.", argName, fun),
         call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(sprintf("Method of evaluation: %s", evalFunc))
    message(sprintf("Number of objects: %d", nrow(x)))
    message(sprintf("Number of predictors: %d", ncol(x)))
    message(sprintf("Number of predictors in model: %d\n", numberOfVariables))
    message("")
    message("Genetic algorithm settings:")
    message(sprintf("Seeds: %s", paste(seeds, collapse = ", ")))
    message(sprintf("Population size: %d", popSize))
    message(sprintf("Mutation rate: %.2f", pMutation))
    message(sprintf("Crossover rate: %.2f", pCrossover))
    message(sprintf("Crossover function: %s", crossoverFunc))
    message(sprintf("Elitism: %d", elitism))
    message(sprintf("Maximum iterations: %d", maxIter))
  }

  # usage inside gaVariableSelection()
  evalFunc <- resolveFunction(evalFunc, "evalFunc")
  crossoverFunc <- resolveFunction(crossoverFunc, "crossoverFunc")

  nPred <- ncol(x)

  computeR2 <- function(yObs, residuals) {
    sse <- sum(residuals^2)
    tss <- sum((yObs - mean(yObs))^2)
    1 - (sse / tss)
  }

  gaResults <- lapply(seq_along(seeds), function(i) {
    sol <- GA::ga(
      type = "real-valued",
      fitness = evalFunc,
      x = x,
      y = y,
      lower = rep(1, numberOfVariables),
      upper = rep(nPred, numberOfVariables),
      population = gaintegerPopulation,
      selection = "ga_lrSelection",
      crossover = crossoverFunc,
      mutation = "gaintegerMutation",
      parallel = FALSE,
      monitor = if (isTRUE(verbose)) QSARMonitorFactory(interval) else FALSE,
      seed = seeds[i],
      pmutation = pMutation,
      pcrossover = pCrossover,
      elitism = elitism,
      popSize = popSize,
      maxiter = maxIter)

    importantPredictors <- sol@solution[1, ]
    Q2 <- sol@fitnessValue

    xFit <- cbind(1, x[, importantPredictors, drop = FALSE])
    mdl <- stats::lm.fit(y = y, x = xFit)
    R2 <- computeR2(y, mdl$residuals)

    yTrainDf <- data.frame(y = y, yhat = mdl$fitted.values, residual = mdl$residuals)

    if (isTRUE(verbose)) {
      cat(sprintf("Seed %d: Q2: %.2f  -  ", seeds[i], Q2))
      cat(sprintf("Predictors: %s\n", paste(colnames(xFit)[-1], collapse= ", ")))
    }

    list(numVar = numberOfVariables, importantPredictors = importantPredictors,
         yTrain = yTrainDf, model = mdl$coefficients, R2Train = R2, Q2Loocv = Q2)
  })

  idxMaxQ2 <- which.max(vapply(seq_along(gaResults), function(jj) {
    gaResults[[jj]]$Q2Loocv
  }, numeric(1)))

  if (isTRUE(verbose)) {
    message(sprintf("Highest Q2 with %d predictors in model: %.3f - ",
                    numberOfVariables, gaResults[[idxMaxQ2]]$Q2Loocv))
    idx2 <- sort(gaResults[[idxMaxQ2]]$importantPredictors)
    message(sprintf("Predictors: %s", paste(colnames(x)[idx2], collapse= ", ")))
  }

  result <- gaResults[[idxMaxQ2]]
  class(result) <- c("gaQSAR", "list")
  return(result)
}
