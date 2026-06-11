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

## ----settings-----------------------------------------------------------------
# smokeTest <- FALSE
# permutationTest <- TRUE
# nWorkers <- 1L
# 
# numVars <- 2:7
# theSeeds <- 1:5
# 
# gaSettings <- list(
#   pMutation = 0.2,
#   pCrossover = 0.7,
#   popSize = 100,
#   maxIter = 300,
#   interval = 50,
#   elitism = 3,
#   crossoverType = "gaintegerOnePointCrossover",
#   KSpercentage = 0.95
# )
# 
# if (smokeTest) {
#   gaSettings$popSize <- 25
#   gaSettings$maxIter <- 10
#   gaSettings$elitism <- 2
#   gaSettings$interval <- 5
#   numVars <- c(3, 5, 7)
#   theSeeds <- 1:2
# }

## ----split-data---------------------------------------------------------------
# splitupMolecules <- splitUp(qsarData, method = "KS", pc = gaSettings$KSpercentage)
# 
# xTrain <- as.matrix(qsarData[splitupMolecules$model, -1, drop = FALSE])
# yTrain <- qsarData[splitupMolecules$model, 1]
# 
# xTest <- as.matrix(qsarData[splitupMolecules$test, -1, drop = FALSE])
# yTest <- qsarData[splitupMolecules$test, 1]

## ----run-ga-------------------------------------------------------------------
# baseArgs <- list(
#   x = xTrain,
#   y = yTrain,
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
#   do.call(gaQSAR::gaVariableSelection, args)
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
#     future.packages = "gaQSAR"
#   )
#   names(output) <- paste0("p", numVars)
# }

## ----predict-test-------------------------------------------------------------
# output <- predictOOBObjects(output, xTest, yTest)

## ----Q2-plot------------------------------------------------------------------
# myQ2Plot <- createQ2Plot(output, label = qsarLabel)
# print(myQ2Plot)
# 

## ----diverse-plot-------------------------------------------------------------
# output <- createWilliamsPlot(
#   output,
#   xTrain,
#   yTrain,
#   xTest,
#   yTest,
#   residualThreshold = 2.5
# )
# 
# print(plot(output[[2]]))

## ----select-model-------------------------------------------------------------
# Q2Values <- vapply(output, function(object) object$Q2Loocv, numeric(1))
# bestIdx <- which.max(Q2Values)
# bestObj <- output[[bestIdx]]
# 
# cat(sprintf("Selected model size: %d predictors\n", numVars[bestIdx]))
# cat(sprintf("LOOCV Q2: %.4f\n", bestObj$Q2Loocv))
# 
# summary(bestObj)

## ----permutation-test---------------------------------------------------------
# if (permutationTest) {
#   nPermutations <- if (smokeTest) 20 else 100
# 
#   permutationResult <- gaPermutationTest(
#     bestObj,
#     x = xTrain,
#     y = yTrain,
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
# save(output, file = timestampFileName(paste0(qsarLabel, "Results.Rdata")))
# 
# if (exists("permutationResult")) {
#   save(
#     permutationResult,
#     file = timestampFileName(paste0(qsarLabel, "PermutationResults.Rdata"))
#   )
# }

