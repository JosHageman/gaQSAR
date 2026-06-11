## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----setup--------------------------------------------------------------------
# library(gaQSAR)
# library(QSARdata)
# 
# timestampFileName <- function(fileName, time = Sys.time()) {
#   paste0(format(time, "%Y%m%d_%H%M%S_"), fileName)
# }

## ----prepare-data-------------------------------------------------------------
# data(AquaticTox)
# 
# qsarLabel <- "AquaticTox"
# 
# qsarData <- cbind(
#   activity = AquaticTox_Outcome$Activity,
#   AquaticTox_moe2D[, -1],
#   AquaticTox_moe3D[, -1]
# )
# 
# missingColumns <- which(colSums(is.na(qsarData)) != 0)
# 
# if (length(missingColumns) > 0) {
#   qsarData <- qsarData[, -missingColumns]
# }
# 
# xAll <- as.matrix(qsarData[, -1, drop = FALSE])
# yAll <- qsarData[, 1]

## ----settings-----------------------------------------------------------------
# smokeTest <- FALSE
# permutationTest <- TRUE
# nWorkers <- 1L
# 
# numVars <- 1:10
# theSeeds <- 1:5
# 
# gaSettings <- list(
#   pMutation = 0.2,
#   pCrossover = 0.7,
#   popSize = 300,
#   maxIter = 300,
#   interval = 50,
#   elitism = 30,
#   crossoverType = "gaintegerOnePointCrossover"
# )
# 
# outerK <- 5
# outerSeed <- 1
# 
# if (smokeTest) {
#   gaSettings$popSize <- 25
#   gaSettings$maxIter <- 10
#   gaSettings$elitism <- 2
#   gaSettings$interval <- 5
#   numVars <- c(3, 5, 7)
#   theSeeds <- 1:2
#   outerK <- 3
# }

## ----run-dcv------------------------------------------------------------------
# baseArgs <- list(
#   x = xAll,
#   y = yAll,
#   outerMethod = "kfold",
#   outerK = outerK,
#   seed = outerSeed,
#   popSize = gaSettings$popSize,
#   pMutation = gaSettings$pMutation,
#   pCrossover = gaSettings$pCrossover,
#   crossoverFunc = gaSettings$crossoverType,
#   elitism = gaSettings$elitism,
#   maxIter = gaSettings$maxIter,
#   interval = gaSettings$interval,
#   seeds = theSeeds,
#   verbose = TRUE
# )
# 
# fitOneModelSize <- function(numberOfVariables) {
#   args <- baseArgs
#   args$numberOfVariables <- numberOfVariables
#   do.call(gaQSAR::gaDoubleCrossValidation, args)
# }
# 
# output <- lapply(numVars, fitOneModelSize)
# names(output) <- paste0("p", numVars)

## ----parallel-run-------------------------------------------------------------
# library(future)
# library(future.apply)
# 
# useFuture <- nWorkers > 1L
# 
# if (useFuture) {
#   oldPlan <- future::plan()
#   on.exit(future::plan(oldPlan), add = TRUE)
#   future::plan(future::multisession, workers = nWorkers)
# 
#   output <- future.apply::future_lapply(
#     X = numVars,
#     FUN = fitOneModelSize,
#     future.seed = TRUE,
#     future.packages = "gaQSAR",
#     future.chunk.size = 1
#   )
#   names(output) <- paste0("p", numVars)
# }

## ----compare-models-----------------------------------------------------------
# trainingPlot <- createDCVTrainingMetricsPlot(
#   output,
#   metrics = c("R2", "R2adj", "Q2"),
#   includeOuterQ2 = TRUE
# )
# 
# print(trainingPlot)

## ----select-model-------------------------------------------------------------
# outerQ2Values <- vapply(output, function(object) object$outerQ2, numeric(1))
# bestIdx <- which.max(outerQ2Values)
# bestObj <- output[[bestIdx]]
# 
# cat(sprintf("Selected model size: %d predictors\n", numVars[bestIdx]))
# cat(sprintf("Outer Q2: %.4f\n", bestObj$outerQ2))

## ----inspect-model------------------------------------------------------------
# summary(bestObj)
# plot(bestObj, type = "all")

## ----dcv-williams-------------------------------------------------------------
# williamsPlot <- createDCVWilliamsPlot(
#   bestObj,
#   label = "AquaticTox data",
#   colorBy = "fold",
#   aggregation = "none",
#   labelOutliers = "rowName"
# )
# 
# print(williamsPlot + ggplot2::facet_wrap(~ fold))

## ----best-fitness-------------------------------------------------------------
# fitnessPlot <- createBestFitnessPlot(bestObj)
# print(fitnessPlot)

## ----permutation-test---------------------------------------------------------
# if (permutationTest) {
#   nPermutations <- if (smokeTest) 20 else 100
# 
#   permutationResult <- gaPermutationTest(
#     bestObj,
#     x = xAll,
#     nPermutations = nPermutations,
#     seed = 1,
#     validateSettings = TRUE,
#     verbose = TRUE,
#     workers = nWorkers
#   )
# 
#   print(plot(permutationResult))
#   summary(permutationResult)
# }

## ----save-results-------------------------------------------------------------
# save(output, file = timestampFileName(paste0(qsarLabel, "NestedCV_Results.Rdata")))
# 
# if (exists("permutationResult")) {
#   save(
#     permutationResult,
#     file = timestampFileName(paste0(qsarLabel, "NestedCV_PermutationResults.Rdata"))
#   )
# }

