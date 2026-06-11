#' Genetic algorithm based variable selection for QSAR
#'
#' Runs a genetic algorithm (GA) to select a fixed number of predictors using
#' leave-one-out cross-validation (LOOCV) for stable optimization.
#' The GA is executed for each seed in `seeds` for robustness, and the best run 
#' (highest LOOCV Q2) is returned.
#'
#' @param x A data.frame or matrix of predictors (rows = observations,
#'   columns = candidate variables).
#' @param y Numeric response vector aligned with rows of `x`.
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
#' @details 
#' **Fitness function:** Uses LOOCV-based Q2 (via `singleCV`) as the GA fitness metric.
#' LOOCV provides a deterministic, non-random evaluation landscape that enables stable
#' variable selection, preventing fold-splitting variability from interfering with 
#' optimization. The nested CV structure (outer: validation, inner: LOOCV fitness) 
#' separates variable selection from performance assessment.
#'
#' **Robustness:** For each seed, the GA runs with identical parameters. The run with 
#' the highest LOOCV Q2 is selected and returned. The best seed is tracked for 
#' reproducibility in permutation tests.
#'
#' @return An object of class "gaQSAR" containing the best GA run for the 
#'   specified model size. The structure typically includes: `numVar`, 
#'   `importantPredictors`, `yTrain`, `model`, `R2Train`, `R2AdjTrain`, `Q2Loocv` 
#'   (LOOCV-based Q2 used by GA for fitness), `VIF` (variance inflation factors for 
#'   each selected predictor), `bestFitnessPerGeneration` (best GA fitness value at 
#'   each generation), `bestSeed` (the seed that produced the optimal solution), and 
#'   `gaSettings` (list of GA parameters used, for reproducibility in permutation 
#'   tests). Use `plot.gaQSAR()` to visualize any attached plots.
#'
#' @seealso [singleCV()], [createQ2Plot()], [createWilliamsPlot()]
#'
#' @examples
#' # This is a toy example for the documentation: GA settings are set so it runs fast.
#' # For real QSAR work: start with the  default GA settings, or run an experimental design
#' # to tune GA settings for your dataset.
#'
#' set.seed(1)
#'
#' # Create toy descriptor matrix (n compounds x p predictors)
#' n <- 40
#' p <- 20
#' x <- matrix(rnorm(n * p), nrow = n)
#'
#' # Add names to mimic typical QSAR data structures (molecules + descriptor names)
#' rownames(x) <- paste0("mol_", seq_len(n))
#' colnames(x) <- paste0("X", seq_len(p))
#'
#' # Create a synthetic response with a known signal in predictors 2 and 7
#' y <- 1.5 * x[, 2] - 0.8 * x[, 7] + rnorm(n, sd = 0.5)
#'
#' # Run GA variable selection
#' # NOTE: Settings below are example-only (fast runtime, not optimal).
#' fit <- gaVariableSelection(
#'   x = x,
#'   y = y,
#'   numberOfVariables = 2,       # target model size (number of selected predictors)
#'   popSize = 20,                # small population for speed (increase in real runs)
#'   maxIter = 20,                # few generations for speed (increase in real runs)
#'   seeds = 1,                   # single seed for a minimal example (use multiple seeds)
#'   interval = 5,                # print progress every 5 generations when verbose=TRUE
#'   verbose = FALSE
#' )
#'
#' # Methods for gaQSAR objects
#' print(fit)
#' summary(fit)
#' plot(fit, type = "predObs")  # observed vs predicted plot
#'
#' @export
gaVariableSelection <- function(x, y,
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
    message(sprintf("Method of evaluation: LOOCV (singleCV)"))
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

  # Resolve crossover function
  crossoverFunc <- resolveFunction(crossoverFunc, "crossoverFunc")
  
  # Use singleCV (LOOCV) as the fitness function
  evalFunc <- resolveFunction("singleCV", "singleCV")

  nPred <- ncol(x)

  computeR2 <- function(yObs, residuals) {
    sse <- sum(residuals^2)
    tss <- sum((yObs - mean(yObs))^2)
    1 - (sse / tss)
  }

  computeAdjR2 <- function(r2, n, p) {
    if (n <= (p + 1)) return(NA_real_)
    1 - (1 - r2) * ((n - 1) / (n - p - 1))
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

    gaSummary <- sol@summary
    bestFitnessPerGen <- NULL
    if (!is.null(gaSummary) && nrow(gaSummary) > 0) {
      bestCol <- intersect(c("max", "Best", "best", "Max"), colnames(gaSummary))
      if (length(bestCol) > 0) {
        bestFitnessPerGen <- gaSummary[, bestCol[1], drop = TRUE]
      }
    }

    importantPredictors <- sol@solution[1, ]
    Q2 <- sol@fitnessValue

    xFit <- cbind(1, x[, importantPredictors, drop = FALSE])
    mdl <- stats::lm.fit(y = y, x = xFit)
    R2 <- computeR2(y, mdl$residuals)
    R2Adj <- computeAdjR2(R2, length(y), ncol(xFit) - 1)

    yTrainDf <- data.frame(y = y, yhat = mdl$fitted.values, residual = mdl$residuals)

    # Calculate VIF values for selected predictors
    vifValues <- vapply(2:ncol(xFit), function(j) {
      xOther <- xFit[, -j, drop = FALSE]
      xTarget <- xFit[, j]
      
      qrObj <- qr(xOther)
      if (qrObj$rank < ncol(xOther)) return(Inf)
      
      beta <- qr.coef(qrObj, xTarget)
      resid <- xTarget - xOther %*% beta
      
      sse <- sum(resid * resid)
      tss <- sum((xTarget - mean(xTarget))^2)
      
      if (tss == 0) Inf else 1 / (sse / tss)
    }, numeric(1))
    
    # Name the VIF values according to predictor names
    names(vifValues) <- colnames(xFit)[-1]

    if (isTRUE(verbose)) {
      cat(sprintf("Seed %d: Q2: %.2f  -  ", seeds[i], Q2))
      cat(sprintf("Predictors: %s\n", paste(colnames(xFit)[-1], collapse= ", ")))
    }

        list(numVar = numberOfVariables, importantPredictors = importantPredictors,
          yTrain = yTrainDf, model = mdl$coefficients,
          R2Train = R2, R2AdjTrain = R2Adj, Q2Loocv = Q2,
          bestFitnessPerGeneration = bestFitnessPerGen,
          VIF = vifValues)
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
    
  # Store GA settings for reproducibility (e.g., permutation tests)
  # Store function objects directly, not their deparsed versions
  result$gaSettings <- list(
    numberOfVariables = numberOfVariables,
    pMutation = pMutation,
    pCrossover = pCrossover,
    crossoverFunc = crossoverFunc,
    popSize = popSize,
    maxIter = maxIter,
    elitism = elitism,
    seeds = seeds,
    interval = interval
  )
  
  # Store the best seed that produced the optimal solution
  result$bestSeed <- seeds[idxMaxQ2]
  
  class(result) <- c("gaQSAR", "list")
  return(result)
}
