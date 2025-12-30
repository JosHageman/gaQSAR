# gaQSAR

![R-CMD-check](https://github.com/joshageman/gaQSAR/actions/workflows/R-CMD-check.yaml/badge.svg)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![CRAN status](https://www.r-pkg.org/badges/version/gaQSAR)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/gaQSAR)

**Genetic Algorithm-based Variable Selection for QSAR Modeling**

[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`gaQSAR` implements genetic algorithm-based variable selection for building quantitative structure-activity relationship (QSAR) models. The package provides a comprehensive workflow for selecting optimal predictor subsets from large descriptor spaces using leave-one-out cross-validation (LOOCV) with Q2 as the fitness criterion.

### Key Features

- **Genetic Algorithm Optimization**: Select optimal predictor subsets from large descriptor spaces
- **Automatic Multicollinearity Handling**: VIF thresholding to avoid correlated predictors
- **Cross-Validation**: LOOCV-based Q2 for robust model selection
- **External Validation**: Test model performance on independent test sets
- **Diagnostic Plots**: Q2 curves and Williams plots for model evaluation
- **Flexible Data Splitting**: Kennard-Stone or random split methods

### Typical Workflow

1. Prepare descriptors + response.
2. Optionally run a small experimental design to tune GA hyperparameters.
3. Run the GA for a range of subset sizes (number of predictors).
4. Validate on an external test set.
5. Inspect diagnostics: Q2 curves and applicability domain (Williams plot).

## Installation

From GitHub (example):

```r
# install.packages("remotes")
remotes::install_github("joshageman/gaQSAR")
```

Load the package:

```r
library(gaQSAR)
```

## Quick start

If you already have `x` (descriptor matrix) and `y` (numeric response):

```r
library(gaQSAR)

fit <- gaVariableSelection(
  x = x,
  y = y,
  numberOfVariables = 4,
  popSize = 100,
  pMutation = 0.2,
  pCrossover = 0.7,
  crossoverFunc = "gaintegerOnePointCrossover",
  elitism = 3,
  maxIter = 300,
  seeds = 1:5,
  interval = 50,
  verbose = TRUE
)

fit$Q2Loocv
```

## Full example: Aquatic Toxicity dataset (QSARdata)

This is a complete, runnable example showing:

- data preparation (AquaticTox)
- optional experimental design for GA settings
- final GA run for multiple subset sizes
- external validation
- diagnostics (Q2 curve, Williams plot)
- optional parallel execution with `future`

### Example script

```r
# load necessary libraries
library(gaQSAR)
library(future)
library(future.apply)
library(parallel)

timestampFileName <- function(fileName, time = Sys.time()) {
  paste0(format(time, "%Y%m%d_%H%M%S_"), fileName)
}

EXPDES <- FALSE
nWorkers <- 6

# start the clock
ptm <- proc.time()

################################
# 1. Prepare the data          #
################################

# use Aquatic Toxicity data set
library(QSARdata)

data(AquaticTox)
QSARLabel <- "AquaticTox"
myTrainData <- cbind(
  activity = AquaticTox_Outcome$Activity,
  AquaticTox_moe2D[, -1],
  AquaticTox_moe3D[, -1]
)

idx <- which(colSums(is.na(myTrainData)) != 0)
myTrainData <- myTrainData[, -idx]

################################
# 2. Do experimental design    #
################################

if (EXPDES) {

  # start the clock
  ptm <- proc.time()

  # settings for full GA-run
  numberOfVariables <- 4
  pMutation         <- c(0.1, 0.2)
  pCrossover        <- c(0.6, 0.7, 0.8)
  popSize           <- c(100, 300)
  elitism           <- c(3, 30)
  mySeeds           <- 1:3
  maxIter           <- 300
  crossoverType     <- c("gaintegerOnePointCrossover", "gaintegerTwoPointCrossover")
  KSpercentage      <- c(0.75, 0.85, 0.95)

  experimentalDesign <- expand.grid(
    numberOfVariables = numberOfVariables,
    pMutation         = pMutation,
    pCrossover        = pCrossover,
    popSize           = popSize,
    elitism           = elitism,
    theSeed           = mySeeds,
    maxIter           = maxIter,
    crossoverType     = crossoverType,
    KSpercentage      = KSpercentage,
    stringsAsFactors  = FALSE
  )

  cat(sprintf("Number of experiments: %d\n", nrow(experimentalDesign)))

  # number of cores
  plan(multisession, workers = nWorkers)
  on.exit(plan(sequential), add = TRUE)

  # Run experiments
  output <- future_lapply(
    X = seq_len(nrow(experimentalDesign)),
    FUN = function(i) {

      expRow <- experimentalDesign[i, , drop = FALSE]

      # divide in training and test set (KS split)
      splitupMolecules <- splitUp(myTrainData, method = "KS", pc = expRow$KSpercentage)

      xtrain <- as.matrix(myTrainData[splitupMolecules$model, -1, drop = FALSE])
      ytrain <- myTrainData[splitupMolecules$model,  1]
      xtest  <- as.matrix(myTrainData[splitupMolecules$test,  -1, drop = FALSE])
      ytest  <- myTrainData[splitupMolecules$test,   1]

      L <- list(
        x = xtrain,
        y = ytrain,
        numberOfVariables = expRow$numberOfVariables,
        popSize           = expRow$popSize,
        pMutation         = expRow$pMutation,
        pCrossover        = expRow$pCrossover,
        crossoverFunc     = expRow$crossoverType,
        elitism           = expRow$elitism,
        maxIter           = expRow$maxIter,
        seeds             = expRow$theSeed,
        verbose           = FALSE
      )

      # Run GA
      fit <- do.call(gaQSAR::gaVariableSelection, L)

      return(fit$Q2Loocv)
    },
    future.seed = TRUE,
    future.packages = c("gaQSAR")
  )

  resultsExpDes <- cbind(experimentalDesign, Q2Loocv = unlist(output))

  cat(sprintf("Combination with highest Q2:\n "))
  print(resultsExpDes[which.max(unlist(output)), ])

  save(resultsExpDes, file = timestampFileName(paste0(QSARLabel, "ExpDes.Rdata")))
  elapsed <- (proc.time() - ptm)[["elapsed"]]
  cat(sprintf("Done. Elapsed time: %.1f sec\n", elapsed))
}

################################
# 3. Do the final run          #
################################

# set up the GA
if (EXPDES) {
  gaSettings <- resultsExpDes[which.max(resultsExpDes$Q2Loocv), ]
  gaSettings$interval <- 50
} else {
  # fixed settings = best settings from ExpDes
  gaSettings <- list(
    pMutation = 0.2, pCrossover = 0.7, popSize = 100, maxIter = 300,
    interval = 50, elitism = 3, crossoverType = "gaintegerOnePointCrossover",
    KSpercentage = 0.95)
  }

# divide in training and test set
splitupMolecules <- splitUp(myTrainData, method = "KS", pc = gaSettings$KSpercentage)
xtrain <- as.matrix(myTrainData[splitupMolecules$model, -1])
ytrain <- myTrainData[splitupMolecules$model, 1]
xtest  <- as.matrix(myTrainData[splitupMolecules$test,  -1])
ytest  <- myTrainData[splitupMolecules$test,   1]

numVars  <- 2:9
theSeeds <- 1:5

useFuture <- nWorkers > 1L

# Defaults: sequential behavior
applyFun <- sapply
applyArgs <- list(simplify = FALSE)
resetPlan <- NULL

# Tip: when prototyping, set nWorkers <- 1 so you do not start multisession.
# This keeps console output responsive while the GA is running.
if (useFuture) {
  oldPlan <- plan()
  plan(multisession, workers = nWorkers)
  resetPlan <- function() plan(oldPlan)
  on.exit(resetPlan(), add = TRUE)

  applyFun <- future_sapply
  applyArgs <- c(applyArgs, list(future.seed = TRUE, future.packages = "gaQSAR"))
}

baseArgs <- list(
  x = xtrain,
  y = ytrain,
  popSize = gaSettings$popSize,
  pMutation = gaSettings$pMutation,
  pCrossover = gaSettings$pCrossover,
  crossoverFunc = gaSettings$crossoverType,
  elitism = gaSettings$elitism,
  maxIter = gaSettings$maxIter,
  interval = gaSettings$interval,
  seeds = theSeeds,
  verbose = TRUE
)

fun <- function(numberOfVariables) {
  args <- baseArgs
  args$numberOfVariables <- numberOfVariables
  do.call(gaQSAR::gaVariableSelection, args)
}

output <- do.call(applyFun, c(list(X = numVars, FUN = fun), applyArgs))

# add predictions for validation sets
output <- predictOOBObjects(output, xtest, ytest)

# add Q2 curves for LOOCV and validation sets
output <- createQ2Plot(output, label = QSARLabel)

# add Williams plots for applicability domain
output <- createWilliamsPlot(output, xtrain, ytrain, xtest, ytest, residualThreshold = 2.5)

# save all the results
save(output, file = timestampFileName(paste0(QSARLabel, "Results.Rdata")))

print(output[[1]]$q2Plot)
plot(output[[2]])
summary(output[[2]])

# what's our runtime?
cat(sprintf("Done in %.1f mins.\n", (proc.time() - ptm)[3] / 60))
```

## What you get back

After the final run, `output` is a list of GA results objects (one per `numberOfVariables` in `numVars`), enriched with:

- external test-set predictions (`predictOOBObjects`)
- a Q2 curve plot (`createQ2Plot`)
- applicability domain diagnostics (Williams plot via `createWilliamsPlot`)

## Parallel execution notes

- Set `nWorkers <- 1` during prototyping for more responsive console output.
- Set `nWorkers` to a higher value for batch runs. In that case, the code uses `future_sapply()` and a `multisession` plan.

## References

This method is demonstrated in:

- Araya-Cloutier, C., Vincken, JP., van de Schans, M.G.M. et al. QSAR-based molecular signatures of prenylated (iso)flavonoids underlying antimicrobial potency against and membrane-disruption in Gram positive and Gram negative bacteria. Sci Rep 8, 9267 (2018). [doi:10.1038/s41598-018-27545-4](https://doi.org/10.1038/s41598-018-27545-4)

- Kalli, S., Araya-Cloutier, C., Hageman, J. et al. Insights into the molecular properties underlying antibacterial activity of prenylated (iso)flavonoids against MRSA. Sci Rep 11, 14180 (2021). [doi:10.1038/s41598-021-92964-9](https://doi.org/10.1038/s41598-021-92964-9)

## License

GPL-3

## Author

Jos Hageman (jos.hageman@wur.nl)  
Wageningen University & Research
