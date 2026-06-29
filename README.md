# gaQSAR

![R-CMD-check](https://github.com/joshageman/gaQSAR/actions/workflows/r.yml/badge.svg)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![CRAN status](https://www.r-pkg.org/badges/version/gaQSAR)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/gaQSAR)

**Genetic algorithm-based variable selection for QSAR modelling**

## Overview

`gaQSAR` performs variable selection for quantitative structure-activity relationship (QSAR) models using a genetic algorithm. Predictor subsets are evaluated with leave-one-out cross-validation (LOOCV), using Q2 as the fitness criterion.

`gaQSAR` contains functions for:

- genetic algorithm-based descriptor selection;
- LOOCV-based model evaluation during variable selection;
- external test set prediction;
- double cross-validation;
- permutation testing by y-scrambling;
- Williams plots and Q2 plots;
- inspection of descriptor selection frequency.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("joshageman/gaQSAR")
```

Load the package:

```r
library(gaQSAR)
```

## Two ways to use gaQSAR

There are two main workflows.

### 1. Train/test workflow

Use this route when you want to split the data into a training set and a separate test set. The genetic algorithm is run on the training set. The selected models are then evaluated on the test set.

This workflow is useful for a first analysis, for comparing model sizes, or when a fixed external test set is available.

Main functions:

```r
splitUp()
gaVariableSelection()
predictOOBObjects()
createQ2Plot()
createWilliamsPlot()
gaPermutationTest()
```

A full example is available in:

```text
vignettes/train-test-workflow.Rmd
```

### 2. Double cross-validation workflow

Use this route when you want a stricter estimate of predictive performance. The outer cross-validation loop is used for model evaluation. Inside each training fold, the genetic algorithm performs descriptor selection.

This workflow is usually the better choice when the goal is to report model performance, because model selection and model evaluation are more clearly separated.

Main functions:

```r
gaDoubleCrossValidation()
createDCVTrainingMetricsPlot()
createDCVWilliamsPlot()
createBestFitnessPlot()
gaPermutationTest()
```

A full example is available in:

```text
vignettes/double-cross-validation-workflow.Rmd
```

## Quick start: one GA run

If you already have a descriptor matrix `x` and a numeric response vector `y`, a single GA run can be started as follows:

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
fit$importantPredictors
summary(fit)
```

This gives one `gaQSAR` object. In practice, it is often useful to repeat the run for several numbers of predictors and compare the resulting models.

## Example: train/test workflow

The train/test workflow starts by splitting the compounds into a training set and a test set. Then `gaVariableSelection()` is run for several model sizes.

```r
library(gaQSAR)
library(QSARdata)

data(AquaticTox)

qsarData <- cbind(
  activity = AquaticTox_Outcome$Activity,
  AquaticTox_moe2D[, -1],
  AquaticTox_moe3D[, -1]
)

missingColumns <- which(colSums(is.na(qsarData)) != 0)

if (length(missingColumns) > 0) {
  qsarData <- qsarData[, -missingColumns]
}

splitupMolecules <- splitUp(qsarData, method = "KS", pc = 0.95)

trainData <- splitupMolecules$trainData
xtest <- as.matrix(splitupMolecules$testData[, -1, drop = FALSE])
ytest <- splitupMolecules$testData[, 1]

xtrain <- as.matrix(trainData[, -1, drop = FALSE])
ytrain <- trainData[, 1]

numVars <- 2:7

fitModels <- lapply(numVars, function(numberOfVariables) {
  gaVariableSelection(
    x = xtrain,
    y = ytrain,
    numberOfVariables = numberOfVariables,
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
})

fitModels <- predictOOBObjects(fitModels, xtest = xtest, ytest = ytest)

createQ2Plot(fitModels, label = "AquaticTox")
createWilliamsPlot(fitModels, xtrain, ytrain, xtest, ytest)
```

`fitModels` is a list of `gaQSAR` objects. Each element corresponds to one model size. This is why the plotting functions for this workflow take a list as input.

A final model can be selected from the list, for example by choosing the model with the highest LOOCV Q2:

```r
Q2Values <- vapply(fitModels, function(object) object$Q2Loocv, numeric(1))
bestIdx <- which.max(Q2Values)
bestModel <- fitModels[[bestIdx]]

summary(bestModel)
plot(bestModel)
```

## Example: double cross-validation workflow

The double cross-validation workflow uses all compounds. The external validation is created inside the outer cross-validation loop.

```r
library(gaQSAR)
library(QSARdata)

data(AquaticTox)

qsarData <- cbind(
  activity = AquaticTox_Outcome$Activity,
  AquaticTox_moe2D[, -1],
  AquaticTox_moe3D[, -1]
)

missingColumns <- which(colSums(is.na(qsarData)) != 0)

if (length(missingColumns) > 0) {
  qsarData <- qsarData[, -missingColumns]
}

xAll <- as.matrix(qsarData[, -1, drop = FALSE])
yAll <- qsarData[, 1]

numVars <- 1:10

dcvModels <- lapply(numVars, function(numberOfVariables) {
  gaDoubleCrossValidation(
    x = xAll,
    y = yAll,
    outerMethod = "kfold",
    outerK = 5,
    seed = 1,
    numberOfVariables = numberOfVariables,
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
})

createDCVTrainingMetricsPlot(dcvModels, label = "AquaticTox")
```

The result is a list of `gaQSAR_dcv` objects. Each element corresponds to one number of predictors.

A model size can be selected using the outer cross-validated Q2:

```r
outerQ2Values <- vapply(dcvModels, function(object) object$outerQ2, numeric(1))
bestIdx <- which.max(outerQ2Values)
bestDcvModel <- dcvModels[[bestIdx]]

summary(bestDcvModel)
plot(bestDcvModel, type = "all")

createDCVWilliamsPlot(
  bestDcvModel,
  label = "AquaticTox",
  colorBy = "fold",
  aggregation = "none",
  labelOutliers = "rowName"
)
```

## Permutation testing

`gaPermutationTest()` can be used to compare the observed model performance with the performance obtained after randomly permuting the response variable.

For a selected train/test model:

```r
permutationResult <- gaPermutationTest(
  object = bestModel,
  x = xtrain,
  nPermutations = 100,
  seed = 1,
  validateSettings = TRUE,
  verbose = TRUE
)

summary(permutationResult)
plot(permutationResult)
```

For a selected double cross-validation result:

```r
permutationResult <- gaPermutationTest(
  object = bestDcvModel,
  x = xAll,
  nPermutations = 100,
  seed = 1,
  validateSettings = TRUE,
  verbose = TRUE
)

summary(permutationResult)
plot(permutationResult)
```

Permutation testing can be computationally expensive. For checking code, use a small value of `nPermutations`. For a real analysis, use a larger number.

## Notes on computation

Genetic algorithm runs can take time, especially when many model sizes, seeds or permutations are used. For quick checks, use smaller values for `popSize`, `maxIter`, `seeds` and `nPermutations`.

For parallel execution, use the `future` and `future.apply` packages in the analysis script. The vignettes contain examples of this.

## References

This approach is related to the QSAR analyses described in:

- Araya-Cloutier, C., Vincken, J.P., van de Schans, M.G.M. et al. QSAR-based molecular signatures of prenylated (iso)flavonoids underlying antimicrobial potency against and membrane-disruption in Gram positive and Gram negative bacteria. *Scientific Reports* 8, 9267 (2018). [doi:10.1038/s41598-018-27545-4](https://doi.org/10.1038/s41598-018-27545-4)

- Kalli, S., Araya-Cloutier, C., Hageman, J. et al. Insights into the molecular properties underlying antibacterial activity of prenylated (iso)flavonoids against MRSA. *Scientific Reports* 11, 14180 (2021). [doi:10.1038/s41598-021-92964-9](https://doi.org/10.1038/s41598-021-92964-9)

## License

GPL-3

## Author

Jos Hageman  
Wageningen University & Research
